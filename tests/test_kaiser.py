"""Independent scalar/model/derivative checks for the built-in signal."""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from fishhighz.covariance import gaussian_covariance
from fishhighz.derivatives import evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.grids import IntegrationGrid
from fishhighz.kernels.kaiser import _coordinates, _damping, _scales
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.models.kaiser import KaiserModel, Scaling
from fishhighz.models.templates import prepare_template
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult


def polys(k, derivative=0):
    """Evaluate smooth and wiggle polynomials in log wavenumber.

    Parameters
    ----------
    k : ndarray of shape (n_nodes,)
        Comoving wavenumbers in h/Mpc.
    derivative : int, optional
        Derivative selector: zero returns power; nonzero returns its first
        derivative with respect to k. Default is 0.

    Returns
    -------
    components : ndarray of shape (n_nodes, 2)
        Smooth and wiggle power in (Mpc/h)^3, or their k derivatives when
        requested.
    """
    log_k = np.log(k)
    if derivative:
        return np.column_stack(((2 + 0.6 * log_k) / k, (1 + 0.3 * log_k * log_k) / k))
    return np.column_stack(
        (12 + 2 * log_k + 0.3 * log_k * log_k, 2 + log_k + 0.1 * log_k**3)
    )


def template(*, lo=0.001, hi=3, h_fid=0.7, smooth_only=False, wiggle_only=False):
    """Prepare a polynomial synthetic matter-power template.

    Parameters
    ----------
    lo : float, optional
        Minimum template wavenumber in h/Mpc. Default is 0.001.
    hi : float, optional
        Maximum template wavenumber in h/Mpc. Default is 3.
    h_fid : float, optional
        Dimensionless fiducial Hubble parameter H0/(100 km/s/Mpc). Default is
        0.7.
    smooth_only : bool, optional
        Whether to set the oscillatory template component to zero. Default is
        False.
    wiggle_only : bool, optional
        Whether to set the smooth template component to zero. Default is False.

    Returns
    -------
    template : PowerTemplate
        Smooth and wiggle spline components on a geometric k grid.
    """
    k_grid = np.geomspace(lo, hi, 15)
    template_components = polys(k_grid)
    if smooth_only:
        template_components[:, 1] = 0
    if wiggle_only:
        template_components[:, 0] = 0
    return prepare_template(
        k_grid,
        template_components.sum(axis=1),
        template_components[:, 0],
        z_ref=2.4,
        h_template=0.7,
        h_fid=h_fid,
    )


def fields():
    """Construct the ordered synthetic observed-field inventory.

    Returns
    -------
    fields : list of ObservedField
        Field identities and tracer kinds in the order used by pair indices.
    """
    return [
        ObservedField("F", "forest", "same", background="qso"),
        ObservedField("g", "galaxy", "same"),
    ]


def model(t=None, **kwargs):
    """Construct the default two-field Kaiser model with explicit overrides.

    Parameters
    ----------
    t : PowerTemplate or None, optional
        Matter-power template; None prepares the standard synthetic template.
        Default is None.
    **kwargs : dict
        KaiserModel overrides, including biases, forest beta, galaxy growth rate, damping widths in Mpc/h, and component dilations.

    Returns
    -------
    model : KaiserModel
        Forest/galaxy power model using the supplied or synthetic template.
    """
    settings = dict(
        biases={"F": -0.3, "g": 1.7},
        betas={"F": 1.2},
        widths={"F": (4, 2), "g": (3, 1)},
        f=0.8,
    )
    settings.update(kwargs)
    return KaiserModel(template() if t is None else t, fields(), **settings)


def scalar_oracle(
    k, mu, pairs, field_list, biases, betas, f, widths, scales, *, growth=1, h_fid=0.7
):
    # Direct per-node/pair formula, independent of production mapping/assembly.
    """Evaluate the Kaiser and AP calculation independently in scalar loops.

    Parameters
    ----------
    k : ndarray of shape (n_nodes,)
        Comoving wavenumbers in h/Mpc.
    mu : ndarray of shape (n_nodes,)
        Dimensionless line-of-sight direction cosines.
    pairs : ndarray of int, shape (n_pairs, 2)
        Observed-field indices defining the requested spectra.
    field_list : sequence of ObservedField
        Observed fields in the order used by pair indices.
    biases : mapping of str to float
        Dimensionless linear bias by field identifier.
    betas : mapping of str to float
        Dimensionless forest RSD parameter by field identifier.
    f : float
        Dimensionless logarithmic growth rate shared by galaxy fields.
    widths : mapping of str to tuple of float
        Parallel and transverse BAO broadening widths in Mpc/h by field.
    scales : sequence of tuple of float
        Parallel/transverse dilation factors for smooth and wiggle components.
    growth : float, optional
        Dimensionless power-growth normalization multiplying the spectrum.
        Default is 1.
    h_fid : float, optional
        Dimensionless fiducial Hubble parameter H0/(100 km/s/Mpc). Default is
        0.7.

    Returns
    -------
    power : ndarray of shape (n_nodes, n_pairs)
        Intrinsic power in (Mpc/h)^3 including component-specific dilation and
        wiggle damping.
    """
    result = []
    for q, m in zip(k, mu):
        row = []
        for i, j in pairs:
            total = 0
            for component, (ap, at) in enumerate(scales):
                # Apply each component dilation before evaluating its RSD factors.
                k_parallel = q * m / ap
                k_transverse = q * math.sqrt(1 - m * m) / at
                mapped = math.sqrt(
                    k_parallel * k_parallel + k_transverse * k_transverse
                )
                angle = k_parallel / mapped
                factors = []
                for index in (i, j):
                    field = field_list[index]
                    bias = biases[field.id]
                    factors.append(
                        bias * (1 + betas[field.id] * angle * angle)
                        if field.kind == "forest"
                        else bias + f * angle * angle
                    )

                # Only the wiggle component receives BAO broadening.
                damp = 1
                if component:
                    sigma_parallel_squared = (
                        widths[field_list[i].id][0] ** 2
                        + widths[field_list[j].id][0] ** 2
                    ) / 2
                    sigma_transverse_squared = (
                        widths[field_list[i].id][1] ** 2
                        + widths[field_list[j].id][1] ** 2
                    ) / 2
                    damp = math.exp(
                        -(
                            k_parallel * k_parallel * sigma_parallel_squared
                            + k_transverse * k_transverse * sigma_transverse_squared
                        )
                        / 2
                    )

                power = (
                    polys(np.array([mapped / (0.7 / h_fid)]))[0, component]
                    * (h_fid / 0.7) ** 3
                )
                total += factors[0] * factors[1] * damp * power / (ap * at * at)
            row.append(growth * total)
        result.append(row)
    return np.array(result)


def test_five_fields_all_pairs_and_subsets():
    """Check five fields all pairs and subsets."""
    observed_fields = [
        ObservedField("F", "forest", "same", background="qso"),
        ObservedField("L", "forest", "same", background="lbg"),
    ] + [ObservedField(x, "galaxy", "same") for x in ("g", "q", "a")]
    biases = dict(zip([f.id for f in observed_fields], [-0.3, -0.4, 0, 1.7, 2.1]))
    betas = {"F": 1.2, "L": 0.5}
    widths = {f.id: (0, 0) for f in observed_fields}
    kaiser_model = KaiserModel(
        template(), observed_fields, biases=biases, betas=betas, widths=widths, f=0.8
    )
    pairs = PairSelection(observed_fields).required_pairs
    k_grid, mu_grid = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    expected = scalar_oracle(
        k_grid,
        mu_grid,
        pairs,
        observed_fields,
        biases,
        betas,
        0.8,
        widths,
        [(1, 1), (1, 1)],
    )
    value = kaiser_model([], 2.4, k_grid, mu_grid, pairs)
    assert value.shape == (3, 15)
    np.testing.assert_allclose(value, expected, rtol=3e-14, atol=2e-13)
    subset = pairs[[12, 1, 7]][:, ::-1]
    np.testing.assert_array_equal(
        kaiser_model([], 2.4, k_grid, mu_grid, subset), value[:, [12, 1, 7]]
    )
    np.testing.assert_array_equal(
        np.concatenate(
            [
                kaiser_model([], 2.4, k_grid[:1], mu_grid[:1], pairs),
                kaiser_model([], 2.4, k_grid[1:], mu_grid[1:], pairs),
            ]
        ),
        value,
    )
    assert value[0, 2] == 0  # F x zero-bias galaxy at mu=0
    assert value[1, 2] < 0


@pytest.mark.parametrize("ap,at", [(1.13, 0.87), (0.93, 1.07), (1.2, 1.2)])
def test_equivalent_bases_coordinates_and_Q(ap, at):
    """Check equivalent bases coordinates and Q.

    Parameters
    ----------
    ap : float
        Parallel dilation factor, supplied by pytest parametrization.
    at : float
        Transverse dilation factor, supplied by pytest parametrization.
    """
    alpha, phi = np.sqrt(ap * at), at / ap
    aiso, epsilon = np.cbrt(ap * at * at), np.cbrt(ap / at) - 1
    k_grid, mu_grid = np.array([0.12, 0.24, 0.3]), np.array([0, 0.4, 1])
    expected = None
    for code, (basis, coordinates) in enumerate(
        (
            ("ap_at", {"ap": ap, "at": at}),
            ("alpha_phi", {"alpha": alpha, "phi": phi}),
            ("alpha_iso_epsilon", {"alpha_iso": aiso, "epsilon": epsilon}),
            ("alpha_iso_phi", {"alpha_iso": aiso, "phi": phi}),
        )
    ):
        parallel_scale, transverse_scale, volume_normalization = _scales(
            np.array(list(coordinates.values())), code
        )
        np.testing.assert_allclose(
            [parallel_scale, transverse_scale, volume_normalization],
            [ap, at, 1 / (ap * at * at)],
            rtol=6e-16,
        )
        mapped, angle, par, per = _coordinates(
            k_grid, mu_grid, parallel_scale, transverse_scale
        )
        np.testing.assert_allclose(
            mapped,
            k_grid * np.sqrt(mu_grid**2 / ap**2 + (1 - mu_grid**2) / at**2),
            rtol=5e-16,
        )
        np.testing.assert_allclose(
            angle,
            (mu_grid / ap) / np.sqrt(mu_grid**2 / ap**2 + (1 - mu_grid**2) / at**2),
            rtol=1e-15,
        )
        scale = Scaling(basis, **coordinates)
        kaiser_model = model(smooth=scale, wiggle=scale)
        value = kaiser_model([], 2.4, k_grid, mu_grid, [[0, 0], [0, 1], [1, 1]])
        if expected is None:
            expected = value
        np.testing.assert_allclose(value, expected, rtol=2e-14, atol=2e-13)


def test_separate_bases_shapes_transformed_angles_units_and_G():
    """Check separate bases shapes transformed angles units and G."""
    power_template = template(h_fid=0.5)
    kaiser_model = model(
        power_template,
        smooth=Scaling("alpha_phi", alpha=1.05, phi=0.81),
        wiggle=Scaling("alpha_iso_epsilon", alpha_iso=0.96, epsilon=0.1),
        z=3,
        growth=0.64,
    )
    k_grid, mu_grid = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    pairs = [[1, 1], [0, 1], [0, 0]]
    expected = scalar_oracle(
        k_grid,
        mu_grid,
        pairs,
        fields(),
        {"F": -0.3, "g": 1.7},
        {"F": 1.2},
        0.8,
        {"F": (4, 2), "g": (3, 1)},
        [(1.05 / 0.9, 1.05 * 0.9), (0.96 * 1.1**2, 0.96 / 1.1)],
        growth=0.64,
        h_fid=0.5,
    )
    np.testing.assert_allclose(
        kaiser_model([], 3, k_grid, mu_grid, pairs), expected, rtol=3e-14, atol=2e-13
    )


def test_constant_power_isolates_separate_Q():
    """Check constant power isolates separate Q."""
    k_grid = np.geomspace(0.001, 3, 5)
    power_template = prepare_template(
        k_grid, np.full(5, 7.0), np.full(5, 3.0), z_ref=0, h_template=1, h_fid=1
    )
    kaiser_model = model(
        power_template,
        biases={"F": 1, "g": 1},
        betas={"F": 0},
        f=0,
        widths={"F": (0, 0), "g": (0, 0)},
        smooth=Scaling("ap_at", ap=1.2, at=0.8),
        wiggle=Scaling("ap_at", ap=0.9, at=1.1),
    )
    expected = 3 / (1.2 * 0.8**2) + 4 / (0.9 * 1.1**2)
    np.testing.assert_allclose(
        kaiser_model([], 0, [0.1, 0.2], [0.2, 0.8], [[0, 1]]), expected, rtol=5e-16
    )


def test_damping_axes_crosses_smooth_only_and_strong_limit():
    """Check damping axes crosses smooth only and strong limit."""
    pairs = np.array([[0, 0], [0, 1], [1, 1]])
    widths = np.array([[4.0, 2.0], [3.0, 1.0]])
    par, per = np.array([0, 0.2, 0.3]), np.array([0.1, 0.2, 0])
    damp, exponent = _damping(par, per, widths**2, pairs)
    np.testing.assert_allclose(damp[:, 1], np.sqrt(damp[:, 0] * damp[:, 2]), rtol=3e-16)
    np.testing.assert_allclose(
        damp[:, 0], np.exp(-0.5 * ((par * 4) ** 2 + (per * 2) ** 2)), rtol=3e-16
    )
    zero, _ = _damping(par, per, np.zeros((2, 2)), pairs)
    np.testing.assert_array_equal(zero, 1)
    huge_zero, _ = _damping(np.array([1e200]), np.array([0.0]), np.zeros((2, 2)), pairs)
    np.testing.assert_array_equal(huge_zero, 1)
    strong, exponent = _damping(par, per, np.full((2, 2), 1e100), pairs)
    assert np.isfinite(exponent).all()
    np.testing.assert_array_equal(strong, 0)
    args = ([], 2.4, [0.1, 0.2, 0.3], [0, 0.4, 1], pairs)
    smooth = template(smooth_only=True)
    np.testing.assert_array_equal(
        model(smooth)(*args), model(smooth, widths={"F": (0, 0), "g": (0, 0)})(*args)
    )
    wiggle = template(wiggle_only=True)
    damped = model(wiggle)(*args)
    plain = model(wiggle, widths={"F": (0, 0), "g": (0, 0)})(*args)
    par = np.array(args[2]) * args[3]
    per = np.array(args[2]) * np.sqrt(1 - np.array(args[3]) ** 2)
    damp, _ = _damping(par, per, widths**2, pairs)
    np.testing.assert_allclose(damped, plain * damp, rtol=3e-16, atol=1e-15)


def test_slots_partial_scaling_and_explicit_ties():
    """Check slots partial scaling and explicit ties."""
    kaiser_model = model(
        local_names=["phi", "b1", "b2"],
        biases={"F": "b1", "g": "b2"},
        smooth=Scaling("alpha_phi", alpha=1, phi="phi"),
    )
    np.testing.assert_allclose(
        kaiser_model([0.81, -0.3, 1.7], 2.4, [0.1], [0.6], [[0, 1]]),
        model(smooth=Scaling("alpha_phi", alpha=1, phi=0.81))(
            [], 2.4, [0.1], [0.6], [[0, 1]]
        ),
        rtol=1e-15,
    )
    observed_fields = [
        ObservedField("a", "forest", "same", background="q"),
        ObservedField("b", "forest", "same", background="q"),
    ]
    tied = KaiserModel(
        template(),
        observed_fields,
        biases={"a": "ba", "b": "bb"},
        betas={"a": 0, "b": 0},
        widths={"a": (0, 0), "b": (0, 0)},
        local_names=["ba", "bb"],
    )
    reg = ParameterRegistry([Parameter("shared", -0.3, "nuisance", step=0.001)])
    collection = PreparedP3D(
        reg,
        PairSelection(observed_fields),
        [
            P3DProvider(
                "tied",
                tied,
                BoundParameters(
                    reg, tied.local_names, {"ba": "shared", "bb": "shared"}
                ),
                [(0, 0), (0, 1), (1, 1)],
            )
        ],
    )
    derivatives = evaluate_derivatives(collection, reg.fiducials, 2.4, [0.1], [0.5])
    np.testing.assert_allclose(
        derivatives.jacobian, 2 * (-0.3) * polys(np.array([0.1])).sum(), atol=2e-13
    )
    assert derivatives.calls[0].model == 3


@pytest.mark.parametrize(
    "kwargs",
    [
        {"biases": {"F": -0.3}},
        {"betas": {"F": 1, "g": 0.5}},
        {"widths": {"F": (1, 2)}},
        {"widths": {"F": (-1, 2), "g": (1, 2)}},
        {"widths": {"F": ("width", 2), "g": (1, 2)}},
        {"widths": {"F": (1e200, 2), "g": (1, 2)}},
        {"f": None},
        {"f": np.inf},
        {"f": "x"},
        {"local_names": ["unused"]},
        {"local_names": ["x", "x"]},
        {"local_names": ["x"], "biases": {"F": "x", "g": "x"}},
        {"smooth": {}},
        {"z": 3},
        {"growth": 0.5},
        {"z": -1},
        {"z": 3, "growth": "G"},
    ],
)
def test_invalid_preparation(kwargs):
    """Check invalid preparation.

    Parameters
    ----------
    kwargs : dict
        Keyword arguments selecting the parametrized case, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError):
        model(**kwargs)


@pytest.mark.parametrize(
    "basis,kwargs",
    [
        ("ap_at", {"ap": 1}),
        ("ap_at", {"ap": 1, "at": 1, "phi": 1}),
        ("aiso_aap", {"aiso": 1, "aap": 1}),
        ("alpha_phi", {"alpha": 1, "phi": 0}),
        ("alpha_iso_epsilon", {"alpha_iso": 1, "epsilon": -1}),
        ("ap_at", {"ap": True, "at": 1}),
    ],
)
def test_invalid_scaling(basis, kwargs):
    """Check invalid scaling.

    Parameters
    ----------
    basis : str
        Dilation parameterization, supplied by pytest parametrization.
    kwargs : dict
        Keyword arguments selecting the parametrized case, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError):
        Scaling(basis, **kwargs)


@pytest.mark.parametrize(
    "theta,z,k,mu,pairs",
    [
        ([1], 2.4, [0.1], [0.5], [[0, 1]]),
        ([], 3, [0.1], [0.5], [[0, 1]]),
        ([], 2.4, [0], [0.5], [[0, 1]]),
        ([], 2.4, [0.1], [1.1], [[0, 1]]),
        ([], 2.4, [0.1], [0.5, 0.6], [[0, 1]]),
        ([], 2.4, [[0.1]], [0.5], [[0, 1]]),
        ([], 2.4, [0.1], [0.5], [[True, 1]]),
        ([], 2.4, [0.1], [0.5], [[0.0, 1.0]]),
        ([], 2.4, [0.1], [0.5], [[0, 2]]),
        ([], 2.4, [0.1], [0.5], [[-1, 1]]),
        ([], 2.4, [0.1], [0.5], []),
        ([], 2.4, [np.nan], [0.5], [[0, 1]]),
    ],
)
def test_invalid_calls(theta, z, k, mu, pairs):
    """Check invalid calls.

    Parameters
    ----------
    theta : list
        Model parameter values, supplied by pytest parametrization.
    z : int or float
        Dimensionless redshift test input, supplied by pytest parametrization.
    k : list
        Wavenumber test input in the convention stated by the tested function,
        supplied by pytest parametrization.
    mu : list
        Dimensionless direction-cosine input, supplied by pytest
        parametrization.
    pairs : list
        Pairs of tracer indices or identities, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError):
        model()(theta, z, k, mu, pairs)


def test_ownership_f_variation_A_B_A_and_no_preparation(monkeypatch):
    """Check ownership f variation A B A and no preparation.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    widths = {"F": [4, 2], "g": [3, 1]}
    kaiser_model = model(
        widths=widths, local_names=["rate"], f="rate", z=3, growth=0.64
    )
    widths["F"][0] = 999
    saved = kaiser_model.widths.copy()
    k_grid = np.array([0.1, 0.9, 0.2, 0.9, 0.3, 0.9])[::2]
    mu_grid = np.array([0, 0.9, 0.4, 0.9, 1, 0.9])[::2]
    k_grid.flags.writeable = mu_grid.flags.writeable = False

    def forbidden(*a, **kw):
        """Fail if a supposedly frozen or unused operation is invoked.

        Parameters
        ----------
        *a : tuple
            Positional arguments forwarded to the original callable or accepted by
            the test callback.
        **kw : dict
            Keyword options forwarded to the original callable or inspected by the
            test callback.

        Raises
        ------
        pytest.fail.Exception
            Raised if the forbidden operation is reached.
        """
        pytest.fail("preparation or FITS called")

    monkeypatch.setattr("scipy.interpolate.CubicSpline", forbidden)
    monkeypatch.setattr("astropy.io.fits.open", forbidden)
    args = (3, k_grid, mu_grid, [[0, 0], [0, 1], [1, 1]])
    fiducial_power = kaiser_model([0.8], *args)
    shifted_growth_power = kaiser_model([0.5], *args)
    np.testing.assert_array_equal(kaiser_model([0.8], *args), fiducial_power)
    np.testing.assert_array_equal(fiducial_power[:, 0], shifted_growth_power[:, 0])
    expected = scalar_oracle(
        k_grid,
        mu_grid,
        args[-1],
        fields(),
        {"F": -0.3, "g": 1.7},
        {"F": 1.2},
        0.5,
        {"F": (4, 2), "g": (3, 1)},
        [(1, 1), (1, 1)],
        growth=0.64,
    )
    np.testing.assert_allclose(shifted_growth_power, expected, rtol=4e-14, atol=1e-13)
    np.testing.assert_array_equal(kaiser_model.widths, saved)
    assert kaiser_model.growth == 0.64 and not kaiser_model.widths.flags.writeable
    assert fiducial_power.flags.owndata and fiducial_power.flags.c_contiguous
    with pytest.raises(ValueError):
        kaiser_model([np.nan], *args)
    with pytest.raises(ValueError, match="nonfinite"):
        kaiser_model([1e308], *args)


def test_domains_identity_and_both_components():
    """Check domains identity and both components."""
    power_template = template(lo=0.1, hi=0.3)
    k_grid = np.array([power_template.k[0], power_template.k[-1]])
    for mu in ([0, 1], [0.37, 0.62], [1, 0]):
        assert np.isfinite(model(power_template)([], 2.4, k_grid, mu, [[0, 1]])).all()
    for name in ("smooth", "wiggle"):
        with pytest.raises(ValueError, match=f"{name}.*mapped range.*template domain"):
            model(power_template, **{name: Scaling("ap_at", ap=1.1, at=1.1)})(
                [], 2.4, k_grid, [0.3, 0.7], [[0, 1]]
            )
    with pytest.raises(ValueError, match="smooth.*mapped range"):
        model(power_template)(
            [], 2.4, [np.nextafter(power_template.k[0], 0)], [0.3], [[0, 1]]
        )


def collection_for(m, registry, bindings):
    """Bind a Kaiser model to a global registry and its complete pair selection.

    Parameters
    ----------
    m : KaiserModel
        Built-in provider whose fields and local parameters are bound.
    registry : ParameterRegistry or None
        Global parameter definitions; None constructs the helper default.
    bindings : mapping of str to str
        Local-to-global parameter-name bindings.

    Returns
    -------
    prepared : PreparedP3D
        Bound built-in provider and covariance-required routing.
    """
    selection = PairSelection(m.fields)
    return PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "builtin",
                m,
                BoundParameters(registry, m.local_names, bindings),
                selection.required_pairs,
            )
        ],
    )


def test_stencil_domain_failure_padding_and_fixed_grid():
    """Check stencil domain failure padding and fixed grid."""
    grid = IntegrationGrid(
        [0.1, 0.3], [0.1, 0.1], [0, 1], [0.5, 0.5], k_min=0.1, k_max=0.3, h_fid=0.7
    )
    snapshots = {
        name: getattr(grid, name).copy()
        for name in ("k_flat", "mu_flat", "weights", "q_mode")
    }
    registry = ParameterRegistry([Parameter("a", 1, "target", step=0.001)])
    for pad in (False, True):
        kaiser_model = model(
            template(lo=0.01 if pad else 0.1, hi=1 if pad else 0.3),
            local_names=["a"],
            wiggle=Scaling("alpha_iso_epsilon", alpha_iso="a", epsilon=0),
        )
        collection = collection_for(kaiser_model, registry, {"a": "a"})
        args = (collection, registry.fiducials, 2.4, grid.k_flat, grid.mu_flat)
        if pad:
            assert np.isfinite(evaluate_derivatives(*args).jacobian).all()
        else:
            with pytest.raises(
                ValueError, match="builtin.*parameter 'a'.*wiggle.*mapped range"
            ):
                evaluate_derivatives(*args)
    for name, snapshot in snapshots.items():
        np.testing.assert_array_equal(getattr(grid, name), snapshot)
    assert grid.k_min == 0.1 and grid.k_max == 0.3


def test_bias_beta_f_product_rule_including_zero_power():
    """Check bias beta f product rule including zero power."""
    reg = ParameterRegistry(
        [
            Parameter(name, value, "nuisance", step=1e-4)
            for name, value in [("bf", -0.3), ("bg", 0), ("beta", 1.2), ("f", 0.8)]
        ]
    )
    kaiser_model = model(
        local_names=reg.ids,
        biases={"F": "bf", "g": "bg"},
        betas={"F": "beta"},
        f="f",
        widths={"F": (0, 0), "g": (0, 0)},
    )
    collection = collection_for(kaiser_model, reg, {p: p for p in reg.ids})
    k_grid, mu_grid = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    derivatives = evaluate_derivatives(collection, reg.fiducials, 2.4, k_grid, mu_grid)
    kaiser_factors = np.column_stack((-0.3 * (1 + 1.2 * mu_grid**2), 0.8 * mu_grid**2))
    local = np.zeros((3, 2, 4))
    local[:, 0, 0] = 1 + 1.2 * mu_grid**2
    local[:, 1, 1] = 1
    local[:, 0, 2] = -0.3 * mu_grid**2
    local[:, 1, 3] = mu_grid**2
    expected = np.empty((3, 3, 4))
    for p, (i, j) in enumerate(collection.selection.required_pairs):
        expected[:, p] = (
            local[:, i] * kaiser_factors[:, j, None]
            + kaiser_factors[:, i, None] * local[:, j]
        ) * polys(k_grid).sum(axis=1)[:, None]
    np.testing.assert_allclose(derivatives.jacobian, expected, rtol=3e-11, atol=3e-11)
    assert derivatives.power[0, 1] == 0 and derivatives.jacobian[0, 1, 1] != 0


@pytest.mark.parametrize(
    "component,damped",
    [("smooth", False), ("wiggle", False), ("wiggle", True), ("tied", True)],
)
def test_isotropic_dilation_analytic_oracle(component, damped):
    """Check isotropic dilation analytic oracle.

    Parameters
    ----------
    component : str
        Spectrum component under examination, supplied by pytest
        parametrization.
    damped : bool
        Whether BAO broadening is enabled, supplied by pytest parametrization.
    """
    reg = ParameterRegistry([Parameter("d", 1, "target", step=2e-5)])
    names = (
        ["s", "w"] if component == "tied" else ["s" if component == "smooth" else "w"]
    )
    kwargs = {}
    for name in names:
        kwargs["smooth" if name == "s" else "wiggle"] = Scaling(
            "alpha_iso_epsilon", alpha_iso=name, epsilon=0
        )
    widths = {"F": (4, 2), "g": (3, 1)} if damped else {"F": (0, 0), "g": (0, 0)}
    kaiser_model = model(local_names=names, z=3, growth=0.64, widths=widths, **kwargs)
    collection = collection_for(kaiser_model, reg, {name: "d" for name in names})
    k_grid, mu_grid = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    result = evaluate_derivatives(collection, reg.fiducials, 3, k_grid, mu_grid)
    powers, slopes = polys(k_grid), polys(k_grid, 1)
    kaiser_factors = np.column_stack(
        (-0.3 * (1 + 1.2 * mu_grid**2), 1.7 + 0.8 * mu_grid**2)
    )
    expected = np.zeros((3, 3))
    for column, (i, j) in enumerate(collection.selection.required_pairs):
        for name in names:
            component_index = 0 if name == "s" else 1
            sigma_parallel_squared = (
                widths[fields()[i].id][0] ** 2 + widths[fields()[j].id][0] ** 2
            ) / 2
            sigma_transverse_squared = (
                widths[fields()[i].id][1] ** 2 + widths[fields()[j].id][1] ** 2
            ) / 2
            sigma = (
                mu_grid**2 * sigma_parallel_squared
                + (1 - mu_grid**2) * sigma_transverse_squared
            )
            damping = np.exp(-0.5 * k_grid * k_grid * sigma) if component_index else 1
            expected[:, column] += (
                0.64
                * kaiser_factors[:, i]
                * kaiser_factors[:, j]
                * damping
                * (
                    -3 * powers[:, component_index]
                    - k_grid * slopes[:, component_index]
                    + (
                        k_grid * k_grid * sigma * powers[:, component_index]
                        if component_index
                        else 0
                    )
                )
            )
    np.testing.assert_allclose(result.jacobian[:, :, 0], expected, rtol=3e-8, atol=2e-8)
    assert result.calls[0].model == 3


def test_basis_jacobian_and_fixed_fisher_transform():
    """Check basis jacobian and fixed fisher transform."""
    reg = ParameterRegistry(
        [
            Parameter("one", 1, "target", step=1e-5),
            Parameter("two", 1, "target", step=1e-5),
        ]
    )
    k_grid, mu_grid = np.array([0.1, 0.15, 0.2, 0.25]), np.array([0, 0.3, 0.7, 1])
    results = []
    for basis, coord, point in [
        ("ap_at", {"ap": "a", "at": "b"}, [1, 1]),
        ("alpha_phi", {"alpha": "a", "phi": "b"}, [1, 1]),
        ("alpha_iso_epsilon", {"alpha_iso": "a", "epsilon": "b"}, [1, 0]),
    ]:
        kaiser_model = model(local_names=["a", "b"], wiggle=Scaling(basis, **coord))
        collection = collection_for(kaiser_model, reg, {"a": "one", "b": "two"})
        results.append(evaluate_derivatives(collection, point, 2.4, k_grid, mu_grid))
    total = results[0].power + np.array([1, 0, 2])
    factors = factor_covariance(
        gaussian_covariance(total, np.ones(4), collection.selection)
    )
    fisher = fisher_from_factors(results[0].jacobian, factors)
    for result, chain in zip(
        results[1:], [np.array([[1, -0.5], [1, 0.5]]), np.array([[1, 2], [1, -1]])]
    ):
        expected = results[0].jacobian @ chain
        np.testing.assert_allclose(result.jacobian, expected, rtol=2e-7, atol=2e-8)
        np.testing.assert_allclose(
            fisher_from_factors(result.jacobian, factors),
            chain.T @ fisher @ chain,
            rtol=3e-8,
            atol=2e-9,
        )


def test_forest_only_f_is_unconstrained():
    """Check forest only f is unconstrained."""
    observed_fields = [ObservedField("F", "forest", "same", background="qso")]
    kaiser_model = KaiserModel(
        template(),
        observed_fields,
        biases={"F": "b"},
        betas={"F": 1},
        widths={"F": (0, 0)},
        local_names=["b"],
    )
    reg = ParameterRegistry(
        [Parameter("b", -0.3, "nuisance", step=0.001), Parameter("f", 0.8, "target")]
    )
    collection = collection_for(kaiser_model, reg, {"b": "b"})
    derivatives = evaluate_derivatives(
        collection, reg.fiducials, 2.4, [0.1, 0.2], [0.3, 0.7]
    )
    np.testing.assert_array_equal(derivatives.jacobian[:, :, 1], 0)
    factors = factor_covariance(
        gaussian_covariance(derivatives.power + 1, np.ones(2), collection.selection)
    )
    forecast = FisherResult(reg, fisher_from_factors(derivatives.jacobian, factors))
    assert forecast.diagnostics.rank == 1 and np.isinf(
        forecast.conditional_errors(["f"])[0]
    )
    with pytest.raises(ValueError, match="rank"):
        forecast.marginalized_errors(["f"])
    with pytest.raises(ValueError, match="forbidden"):
        KaiserModel(
            template(),
            observed_fields,
            biases={"F": -0.3},
            betas={"F": 1},
            widths={"F": (0, 0)},
            f=0.8,
        )


def test_synthetic_example_fixed_factors_and_step_convergence(monkeypatch):
    """Check synthetic example fixed factors and step convergence.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    spec = importlib.util.spec_from_file_location(
        "builtin_example",
        Path(__file__).resolve().parents[1] / "examples/builtin_forecast.py",
    )
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    calls = []
    original = example.factor_covariance

    def count(covariance):
        """Record a covariance factorization and call the original implementation.

        Parameters
        ----------
        covariance : ndarray of shape (n_nodes, n_pairs, n_pairs)
            Selected-spectrum covariance blocks in power-squared units.

        Returns
        -------
        factors : ndarray
            Covariance Cholesky factors.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        calls.append(covariance.copy())
        return original(covariance)

    monkeypatch.setattr(example, "factor_covariance", count)
    report = example.run()
    assert len(calls) == 2
    for result in report.values():
        assert result["convergence_passed"]
        assert np.isfinite(result["target_errors_by_step_scale"]).all()
        assert result["fisher_max_changes"][1] < result["fisher_max_changes"][0] / 2
        assert result["jacobian_max_changes"][1] < result["jacobian_max_changes"][0] / 2


@pytest.mark.parametrize(
    "basis,coordinates,theta",
    [
        ("ap_at", {"ap": "a", "at": 1}, [0]),
        ("alpha_phi", {"alpha": 1, "phi": "a"}, [-1]),
        ("alpha_iso_epsilon", {"alpha_iso": 1, "epsilon": "a"}, [-1]),
        ("ap_at", {"ap": "a", "at": 1e200}, [1e200]),
    ],
)
def test_invalid_free_scales_and_volume(basis, coordinates, theta):
    """Check invalid free scales and volume.

    Parameters
    ----------
    basis : str
        Dilation parameterization, supplied by pytest parametrization.
    coordinates : dict
        Dilation coordinates, supplied by pytest parametrization.
    theta : list
        Model parameter values, supplied by pytest parametrization.
    """
    kaiser_model = model(local_names=["a"], wiggle=Scaling(basis, **coordinates))
    with pytest.raises(ValueError, match="wiggle.*scal"):
        kaiser_model(theta, 2.4, [0.1], [0.5], [[0, 1]])


def test_template_evaluated_once_per_component(monkeypatch):
    """Check template evaluated once per component.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    kaiser_model = model()
    original = type(kaiser_model.template).evaluate
    calls = []

    def counted(self, k, **kwargs):
        """Record template queries and call the original evaluator.

        Parameters
        ----------
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        **kwargs : dict
            Keyword options forwarded to the original callable or inspected by the
            test callback.

        Returns
        -------
        components : ndarray
            Template component values from the original evaluator.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        calls.append(np.array(k))
        return original(self, k, **kwargs)

    monkeypatch.setattr(type(kaiser_model.template), "evaluate", counted)
    kaiser_model([], 2.4, [0.1, 0.2], [0.3, 0.7], [[0, 0], [0, 1], [1, 1]])
    assert len(calls) == 2


def test_builtin_import_keeps_preparation_libraries_lazy(tmp_path):
    """Check builtin import keeps preparation libraries lazy.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    import subprocess
    import sys

    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-c",
            "import sys; import fishhighz.models.kaiser; assert not {'astropy','scipy','vega','camb'} & set(sys.modules)",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == result.stderr == ""
