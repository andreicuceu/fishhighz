"""Individual-spectrum Fisher matrices from the joint pass, and active columns."""

import numpy as np
import pytest

from fishhighz.covariance import gaussian_covariance, gaussian_variances
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.forecast import prepare_bin, run_bin, run_forecast
from fishhighz.geometry import prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.public import _pair_spec
from fishhighz.response import InstrumentResponse
from fishhighz.survey import BinSpec, ForestInput

BASE = np.array([[4, -1, 0.5], [-1, 3, -0.2], [0.5, -0.2, 2]])


def _geometry():
    """Prepare a synthetic common-volume redshift-bin geometry.

    Returns
    -------
    geometry : BinGeometry
        Distances in Mpc/h, volume in (Mpc/h)^3, and velocity conversion in km/s
        per Mpc/h.
    """
    return prepare_geometry(
        2,
        3,
        z_eval=2.5,
        area_deg2=10,
        h_fid=0.7,
        z_order=4,
        hubble=lambda z: np.full_like(z, 200),
        transverse_distance=lambda z: 1000 * (1 + z),
    )


def _registry():
    # "c" is bound to no provider, so its global column stays exactly zero.
    """Construct the shared target/nuisance registry for spectrum extraction.

    Returns
    -------
    registry : ParameterRegistry
        Two target parameters and one nuisance parameter in fixed order.
    """
    return ParameterRegistry(
        [
            Parameter("a", 1.5, "target", step=0.01),
            Parameter("c", 0.0, "target", step=0.01),
            Parameter("b", 0.2, "nuisance", step=0.01),
        ]
    )


def _model(t, z, k, mu, pairs):
    """Evaluate the synthetic spectrum used by the enclosing regression test.

    Parameters
    ----------
    t : ndarray of shape (n_parameters,)
        Local model parameters in the provider binding order.
    z : float or ndarray
        Dimensionless redshift.
    k : ndarray of shape (n_nodes,)
        Comoving wavenumbers in h/Mpc.
    mu : ndarray of shape (n_nodes,)
        Dimensionless line-of-sight direction cosines.
    pairs : ndarray of int, shape (n_pairs, 2)
        Observed-field indices defining the requested spectra.

    Returns
    -------
    power : ndarray of shape (n_nodes, n_pairs)
        Synthetic intrinsic power in (Mpc/h)^3 before response and noise.
    """
    shape = (1 + k[:, None]) * (1 + t[1] * mu[:, None] ** 2)
    return t[0] * shape * BASE[pairs[:, 0], pairs[:, 1]]


def _spec(fields, selected, *, noise="galaxy"):
    """Construct a forecast bin with the requested field, pair, and noise choices.

    Parameters
    ----------
    fields : sequence of ObservedField
        Observed fields in the order used by pair indices.
    selected : sequence of pair or None
        Selected field pairs; None selects every unique pair.
    noise : str, optional
        Noise construction: galaxy, forest, or explicitly supplied full noise.
        Default is 'galaxy'.

    Returns
    -------
    spec : BinSpec
        Synthetic bin with fixed response and either full, galaxy, or forest
        noise.
    """
    registry = _registry()
    selection = PairSelection(fields, selected)
    binding = BoundParameters(registry, ("a", "b"), {"a": "a", "b": "b"})
    p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("matrix", _model, binding, selection.required_pairs)],
    )
    grid = gauss_legendre_grid([0.02, 0.1, 0.2], k_order=3, mu_order=4, h_fid=0.7)
    responses = {f.id: InstrumentResponse(30, 10) for f in fields}
    if noise == "full":
        i, j = selection.required_pairs.T
        full = 0.5 + 0.1 * (i == j) + 0.01 * np.arange(len(grid.k_flat))[:, None]
        return BinSpec(
            "full", _geometry(), grid, p3d, responses, full_noise=full * (i == j)
        )
    forests = {}
    if noise == "forest":
        forests["f"] = ForestInput(
            dict(
                z_source=4,
                magnitudes=[20, 21],
                quadrature=[1, 1],
                rho=[0.01, 0.02],
                variance=[1, 2],
                length_velocity=10000,
                method="legacy",
                iterations=3,
            ),
            lambda t, z, k: np.ones_like(k) * 2,
            BoundParameters(registry, (), {}),
            registry.fiducials,
            auxiliary_coordinates=(24, 0.00035),
        )
    active = set(np.unique(selection.selected_pairs).tolist())
    galaxies = {
        f.id: 0.7 + index
        for index, f in enumerate(fields)
        if f.kind == "galaxy" and index in active
    }
    return BinSpec(
        noise,
        _geometry(),
        grid,
        p3d,
        responses,
        forests=forests,
        galaxies=galaxies,
        independent_sampling=True,
    )


GALAXIES = [ObservedField(str(i), "galaxy", str(i)) for i in range(3)]
MIXED = [
    ObservedField("f", "forest", "lya", background="qso"),
    ObservedField("g", "galaxy", "g"),
    ObservedField("h", "galaxy", "h"),
]


@pytest.mark.parametrize(
    "fields,selected,noise",
    [
        (GALAXIES, None, "galaxy"),
        (GALAXIES, [("2", "0"), ("1", "1"), ("0", "0")], "galaxy"),
        (GALAXIES, None, "full"),
        (MIXED, None, "forest"),
        (MIXED, [("f", "f"), ("f", "g"), ("g", "g")], "forest"),
    ],
)
@pytest.mark.parametrize("batch", [None, 7])
def test_individual_equals_independent_one_spectrum_bins(
    fields, selected, noise, batch
):
    """Check individual equals independent one spectrum bins.

    Parameters
    ----------
    fields : list
        Observed-field definitions, supplied by pytest parametrization.
    selected : list or None
        Selected spectrum definitions, supplied by pytest parametrization.
    noise : str
        Known-noise test input, supplied by pytest parametrization.
    batch : int or None
        Fourier-node batch size, supplied by pytest parametrization.
    """
    spec = _spec(fields, selected, noise=noise)
    joint = prepare_bin(spec)
    run = run_bin(joint, batch_size=batch, individual=True)
    pairs = joint.p3d.selection.selected_pairs.tolist()
    assert len(run.individual) == len(pairs)
    for pair, own in zip(pairs, run.individual):
        alone = prepare_bin(_pair_spec(spec, tuple(pair)))
        reference = run_forecast([alone], batch_size=batch).bins[0].result
        np.testing.assert_array_equal(own.data_fisher, reference.data_fisher)
        assert own.registry is reference.registry


def test_variances_are_the_covariance_diagonal():
    """Check variances are the covariance diagonal."""
    rng = np.random.default_rng(7)
    selection = PairSelection(GALAXIES, [("2", "0"), ("1", "1"), ("0", "0")])
    raw = rng.normal(size=(40, 3, 3))
    total = raw @ raw.swapaxes(1, 2) + 3 * np.eye(3)
    i, j = selection.required_pairs.T
    modes = rng.uniform(0.2, 4, 40)
    covariance = gaussian_covariance(total[:, i, j], modes, selection)
    variances = gaussian_variances(total[:, i, j], modes, selection)
    np.testing.assert_array_equal(variances, np.diagonal(covariance, axis1=1, axis2=2))
    with pytest.raises(ValueError, match="nonpositive variance"):
        gaussian_variances(np.zeros_like(total[:, i, j]), modes, selection)


def test_inactive_columns_are_exactly_zero_and_joint_matches_direct_solve():
    """Check inactive columns are exactly zero and joint matches direct solve."""
    spec = _spec(GALAXIES, None)
    prepared = prepare_bin(spec)
    data = run_bin(prepared, batch_size=5).result.data_fisher
    unused = spec.p3d.registry.ids.index("c")
    assert np.count_nonzero(data[unused]) == 0
    assert np.count_nonzero(data[:, unused]) == 0
    covariance = prepared.factors @ prepared.factors.swapaxes(1, 2)
    selected = prepared.p3d.selection.selected_to_required
    step = 1e-6
    theta = prepared.theta.copy()
    jac = []
    for index in (0, 2):
        up, down = theta.copy(), theta.copy()
        up[index] += step
        down[index] -= step
        nodes = (2.5, prepared.k, prepared.mu, spec.p3d.selection.required_pairs)
        local = (_model(up[[0, 2]], *nodes) - _model(down[[0, 2]], *nodes)) / (2 * step)
        jac.append((prepared.products * local)[:, selected])
    jac = np.stack(jac, axis=-1)
    direct = np.einsum("nsi,nsj->ij", jac, np.linalg.solve(covariance, jac))
    np.testing.assert_allclose(data[np.ix_([0, 2], [0, 2])], direct, rtol=1e-7)


@pytest.mark.parametrize("nodes", [1, 7, 300])
def test_one_spectrum_numpy_contraction_matches_reference_loop(nodes, monkeypatch):
    """Check one spectrum numpy contraction matches reference loop.

    Parameters
    ----------
    nodes : int
        Evaluation nodes, supplied by pytest parametrization.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    from fishhighz.fisher import factor_covariance, fisher_from_factors

    monkeypatch.setenv("FISHHIGHZ_FISHER_BACKEND", "numpy")
    rng = np.random.default_rng(nodes)
    factors = factor_covariance(rng.uniform(0.1, 10, (nodes, 1, 1)))
    jac = rng.normal(size=(nodes, 1, 3)) * np.geomspace(1e-3, 1e3, 3)
    fast = fisher_from_factors(jac, factors)
    with np.errstate(under="warn"):  # selects the reference node loop
        loop = fisher_from_factors(jac, factors)
    np.testing.assert_array_equal(fast, loop)
    bad = jac.copy()
    bad[-1, 0, 1] = np.inf
    messages = []
    for policy in ("ignore", "warn"):
        with np.errstate(under=policy), pytest.raises(ValueError) as error:
            fisher_from_factors(bad, factors)
        messages.append(str(error.value))
    assert messages[0] == messages[1]
