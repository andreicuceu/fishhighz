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
    x = np.log(k)
    if derivative:
        return np.column_stack(((2 + 0.6 * x) / k, (1 + 0.3 * x * x) / k))
    return np.column_stack((12 + 2 * x + 0.3 * x * x, 2 + x + 0.1 * x**3))


def template(*, lo=0.001, hi=3, h_fid=0.7, smooth_only=False, wiggle_only=False):
    k = np.geomspace(lo, hi, 15)
    p = polys(k)
    if smooth_only:
        p[:, 1] = 0
    if wiggle_only:
        p[:, 0] = 0
    return prepare_template(
        k, p.sum(axis=1), p[:, 0], z_ref=2.4, h_template=0.7, h_fid=h_fid
    )


def fields():
    return [
        ObservedField("F", "forest", "same", background="qso"),
        ObservedField("g", "galaxy", "same"),
    ]


def model(t=None, **kwargs):
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
    result = []
    for q, m in zip(k, mu):
        row = []
        for i, j in pairs:
            total = 0
            for component, (ap, at) in enumerate(scales):
                kp = q * m / ap
                kt = q * math.sqrt(1 - m * m) / at
                mapped = math.sqrt(kp * kp + kt * kt)
                angle = kp / mapped
                factors = []
                for index in (i, j):
                    field = field_list[index]
                    b = biases[field.id]
                    factors.append(
                        b * (1 + betas[field.id] * angle * angle)
                        if field.kind == "forest"
                        else b + f * angle * angle
                    )
                damp = 1
                if component:
                    sp = (
                        widths[field_list[i].id][0] ** 2
                        + widths[field_list[j].id][0] ** 2
                    ) / 2
                    st = (
                        widths[field_list[i].id][1] ** 2
                        + widths[field_list[j].id][1] ** 2
                    ) / 2
                    damp = math.exp(-(kp * kp * sp + kt * kt * st) / 2)
                power = (
                    polys(np.array([mapped / (0.7 / h_fid)]))[0, component]
                    * (h_fid / 0.7) ** 3
                )
                total += factors[0] * factors[1] * damp * power / (ap * at * at)
            row.append(growth * total)
        result.append(row)
    return np.array(result)


def test_five_fields_all_pairs_and_subsets():
    fs = [
        ObservedField("F", "forest", "same", background="qso"),
        ObservedField("L", "forest", "same", background="lbg"),
    ] + [ObservedField(x, "galaxy", "same") for x in ("g", "q", "a")]
    biases = dict(zip([f.id for f in fs], [-0.3, -0.4, 0, 1.7, 2.1]))
    betas = {"F": 1.2, "L": 0.5}
    widths = {f.id: (0, 0) for f in fs}
    m = KaiserModel(template(), fs, biases=biases, betas=betas, widths=widths, f=0.8)
    pairs = PairSelection(fs).required_pairs
    k, mu = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    expected = scalar_oracle(
        k, mu, pairs, fs, biases, betas, 0.8, widths, [(1, 1), (1, 1)]
    )
    value = m([], 2.4, k, mu, pairs)
    assert value.shape == (3, 15)
    np.testing.assert_allclose(value, expected, rtol=3e-14, atol=2e-13)
    subset = pairs[[12, 1, 7]][:, ::-1]
    np.testing.assert_array_equal(m([], 2.4, k, mu, subset), value[:, [12, 1, 7]])
    np.testing.assert_array_equal(
        np.concatenate(
            [m([], 2.4, k[:1], mu[:1], pairs), m([], 2.4, k[1:], mu[1:], pairs)]
        ),
        value,
    )
    assert value[0, 2] == 0  # F x zero-bias galaxy at mu=0
    assert value[1, 2] < 0


@pytest.mark.parametrize("ap,at", [(1.13, 0.87), (0.93, 1.07), (1.2, 1.2)])
def test_equivalent_bases_coordinates_and_Q(ap, at):
    alpha, phi = np.sqrt(ap * at), at / ap
    aiso, epsilon = np.cbrt(ap * at * at), np.cbrt(ap / at) - 1
    k, mu = np.array([0.12, 0.24, 0.3]), np.array([0, 0.4, 1])
    expected = None
    for code, (basis, coordinates) in enumerate(
        (
            ("ap_at", {"ap": ap, "at": at}),
            ("alpha_phi", {"alpha": alpha, "phi": phi}),
            ("alpha_iso_epsilon", {"alpha_iso": aiso, "epsilon": epsilon}),
        )
    ):
        p, t, q = _scales(np.array(list(coordinates.values())), code)
        np.testing.assert_allclose([p, t, q], [ap, at, 1 / (ap * at * at)], rtol=6e-16)
        mapped, angle, par, per = _coordinates(k, mu, p, t)
        np.testing.assert_allclose(
            mapped, k * np.sqrt(mu**2 / ap**2 + (1 - mu**2) / at**2), rtol=5e-16
        )
        np.testing.assert_allclose(
            angle, (mu / ap) / np.sqrt(mu**2 / ap**2 + (1 - mu**2) / at**2), rtol=1e-15
        )
        scale = Scaling(basis, **coordinates)
        m = model(smooth=scale, wiggle=scale)
        value = m([], 2.4, k, mu, [[0, 0], [0, 1], [1, 1]])
        if expected is None:
            expected = value
        np.testing.assert_allclose(value, expected, rtol=2e-14, atol=2e-13)


def test_separate_bases_shapes_transformed_angles_units_and_G():
    t = template(h_fid=0.5)
    m = model(
        t,
        smooth=Scaling("alpha_phi", alpha=1.05, phi=0.81),
        wiggle=Scaling("alpha_iso_epsilon", alpha_iso=0.96, epsilon=0.1),
        z=3,
        growth=0.64,
    )
    k, mu = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    pairs = [[1, 1], [0, 1], [0, 0]]
    expected = scalar_oracle(
        k,
        mu,
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
    np.testing.assert_allclose(m([], 3, k, mu, pairs), expected, rtol=3e-14, atol=2e-13)


def test_constant_power_isolates_separate_Q():
    k = np.geomspace(0.001, 3, 5)
    t = prepare_template(
        k, np.full(5, 7.0), np.full(5, 3.0), z_ref=0, h_template=1, h_fid=1
    )
    m = model(
        t,
        biases={"F": 1, "g": 1},
        betas={"F": 0},
        f=0,
        widths={"F": (0, 0), "g": (0, 0)},
        smooth=Scaling("ap_at", ap=1.2, at=0.8),
        wiggle=Scaling("ap_at", ap=0.9, at=1.1),
    )
    expected = 3 / (1.2 * 0.8**2) + 4 / (0.9 * 1.1**2)
    np.testing.assert_allclose(
        m([], 0, [0.1, 0.2], [0.2, 0.8], [[0, 1]]), expected, rtol=5e-16
    )


def test_damping_axes_crosses_smooth_only_and_strong_limit():
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
    m = model(
        local_names=["phi", "b1", "b2"],
        biases={"F": "b1", "g": "b2"},
        smooth=Scaling("alpha_phi", alpha=1, phi="phi"),
    )
    np.testing.assert_allclose(
        m([0.81, -0.3, 1.7], 2.4, [0.1], [0.6], [[0, 1]]),
        model(smooth=Scaling("alpha_phi", alpha=1, phi=0.81))(
            [], 2.4, [0.1], [0.6], [[0, 1]]
        ),
        rtol=1e-15,
    )
    fs = [
        ObservedField("a", "forest", "same", background="q"),
        ObservedField("b", "forest", "same", background="q"),
    ]
    tied = KaiserModel(
        template(),
        fs,
        biases={"a": "ba", "b": "bb"},
        betas={"a": 0, "b": 0},
        widths={"a": (0, 0), "b": (0, 0)},
        local_names=["ba", "bb"],
    )
    reg = ParameterRegistry([Parameter("shared", -0.3, "nuisance", step=0.001)])
    collection = PreparedP3D(
        reg,
        PairSelection(fs),
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
    d = evaluate_derivatives(collection, reg.fiducials, 2.4, [0.1], [0.5])
    np.testing.assert_allclose(
        d.jacobian, 2 * (-0.3) * polys(np.array([0.1])).sum(), atol=2e-13
    )
    assert d.calls[0].model == 3


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
    with pytest.raises(ValueError):
        model()(theta, z, k, mu, pairs)


def test_ownership_f_variation_A_B_A_and_no_preparation(monkeypatch):
    widths = {"F": [4, 2], "g": [3, 1]}
    m = model(widths=widths, local_names=["rate"], f="rate", z=3, growth=0.64)
    widths["F"][0] = 999
    saved = m.widths.copy()
    k = np.array([0.1, 0.9, 0.2, 0.9, 0.3, 0.9])[::2]
    mu = np.array([0, 0.9, 0.4, 0.9, 1, 0.9])[::2]
    k.flags.writeable = mu.flags.writeable = False

    def forbidden(*a, **kw):
        pytest.fail("preparation or FITS called")

    monkeypatch.setattr("scipy.interpolate.CubicSpline", forbidden)
    monkeypatch.setattr("astropy.io.fits.open", forbidden)
    args = (3, k, mu, [[0, 0], [0, 1], [1, 1]])
    a = m([0.8], *args)
    b = m([0.5], *args)
    np.testing.assert_array_equal(m([0.8], *args), a)
    np.testing.assert_array_equal(a[:, 0], b[:, 0])
    expected = scalar_oracle(
        k,
        mu,
        args[-1],
        fields(),
        {"F": -0.3, "g": 1.7},
        {"F": 1.2},
        0.5,
        {"F": (4, 2), "g": (3, 1)},
        [(1, 1), (1, 1)],
        growth=0.64,
    )
    np.testing.assert_allclose(b, expected, rtol=4e-14, atol=1e-13)
    np.testing.assert_array_equal(m.widths, saved)
    assert m.growth == 0.64 and not m.widths.flags.writeable
    assert a.flags.owndata and a.flags.c_contiguous
    with pytest.raises(ValueError):
        m([np.nan], *args)
    with pytest.raises(ValueError, match="nonfinite"):
        m([1e308], *args)


def test_domains_identity_and_both_components():
    t = template(lo=0.1, hi=0.3)
    k = np.array([t.k[0], t.k[-1]])
    for mu in ([0, 1], [0.37, 0.62], [1, 0]):
        assert np.isfinite(model(t)([], 2.4, k, mu, [[0, 1]])).all()
    for name in ("smooth", "wiggle"):
        with pytest.raises(ValueError, match=f"{name}.*mapped range.*template domain"):
            model(t, **{name: Scaling("ap_at", ap=1.1, at=1.1)})(
                [], 2.4, k, [0.3, 0.7], [[0, 1]]
            )
    with pytest.raises(ValueError, match="smooth.*mapped range"):
        model(t)([], 2.4, [np.nextafter(t.k[0], 0)], [0.3], [[0, 1]])


def collection_for(m, registry, bindings):
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
    grid = IntegrationGrid(
        [0.1, 0.3], [0.1, 0.1], [0, 1], [0.5, 0.5], k_min=0.1, k_max=0.3, h_fid=0.7
    )
    snapshots = {
        name: getattr(grid, name).copy()
        for name in ("k_flat", "mu_flat", "weights", "q_mode")
    }
    registry = ParameterRegistry([Parameter("a", 1, "target", step=0.001)])
    for pad in (False, True):
        m = model(
            template(lo=0.01 if pad else 0.1, hi=1 if pad else 0.3),
            local_names=["a"],
            wiggle=Scaling("alpha_iso_epsilon", alpha_iso="a", epsilon=0),
        )
        collection = collection_for(m, registry, {"a": "a"})
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
    reg = ParameterRegistry(
        [
            Parameter(name, value, "nuisance", step=1e-4)
            for name, value in [("bf", -0.3), ("bg", 0), ("beta", 1.2), ("f", 0.8)]
        ]
    )
    m = model(
        local_names=reg.ids,
        biases={"F": "bf", "g": "bg"},
        betas={"F": "beta"},
        f="f",
        widths={"F": (0, 0), "g": (0, 0)},
    )
    collection = collection_for(m, reg, {p: p for p in reg.ids})
    k, mu = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    d = evaluate_derivatives(collection, reg.fiducials, 2.4, k, mu)
    B = np.column_stack((-0.3 * (1 + 1.2 * mu**2), 0.8 * mu**2))
    local = np.zeros((3, 2, 4))
    local[:, 0, 0] = 1 + 1.2 * mu**2
    local[:, 1, 1] = 1
    local[:, 0, 2] = -0.3 * mu**2
    local[:, 1, 3] = mu**2
    expected = np.empty((3, 3, 4))
    for p, (i, j) in enumerate(collection.selection.required_pairs):
        expected[:, p] = (
            local[:, i] * B[:, j, None] + B[:, i, None] * local[:, j]
        ) * polys(k).sum(axis=1)[:, None]
    np.testing.assert_allclose(d.jacobian, expected, rtol=3e-11, atol=3e-11)
    assert d.power[0, 1] == 0 and d.jacobian[0, 1, 1] != 0


@pytest.mark.parametrize(
    "component,damped",
    [("smooth", False), ("wiggle", False), ("wiggle", True), ("tied", True)],
)
def test_isotropic_dilation_analytic_oracle(component, damped):
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
    m = model(local_names=names, z=3, growth=0.64, widths=widths, **kwargs)
    collection = collection_for(m, reg, {name: "d" for name in names})
    k, mu = np.array([0.1, 0.2, 0.3]), np.array([0, 0.4, 1])
    result = evaluate_derivatives(collection, reg.fiducials, 3, k, mu)
    powers, slopes = polys(k), polys(k, 1)
    B = np.column_stack((-0.3 * (1 + 1.2 * mu**2), 1.7 + 0.8 * mu**2))
    expected = np.zeros((3, 3))
    for column, (i, j) in enumerate(collection.selection.required_pairs):
        for name in names:
            c = 0 if name == "s" else 1
            sp = (widths[fields()[i].id][0] ** 2 + widths[fields()[j].id][0] ** 2) / 2
            st = (widths[fields()[i].id][1] ** 2 + widths[fields()[j].id][1] ** 2) / 2
            sigma = mu**2 * sp + (1 - mu**2) * st
            damping = np.exp(-0.5 * k * k * sigma) if c else 1
            expected[:, column] += (
                0.64
                * B[:, i]
                * B[:, j]
                * damping
                * (
                    -3 * powers[:, c]
                    - k * slopes[:, c]
                    + (k * k * sigma * powers[:, c] if c else 0)
                )
            )
    np.testing.assert_allclose(result.jacobian[:, :, 0], expected, rtol=3e-8, atol=2e-8)
    assert result.calls[0].model == 3


def test_basis_jacobian_and_fixed_fisher_transform():
    reg = ParameterRegistry(
        [
            Parameter("one", 1, "target", step=1e-5),
            Parameter("two", 1, "target", step=1e-5),
        ]
    )
    k, mu = np.array([0.1, 0.15, 0.2, 0.25]), np.array([0, 0.3, 0.7, 1])
    results = []
    for basis, coord, point in [
        ("ap_at", {"ap": "a", "at": "b"}, [1, 1]),
        ("alpha_phi", {"alpha": "a", "phi": "b"}, [1, 1]),
        ("alpha_iso_epsilon", {"alpha_iso": "a", "epsilon": "b"}, [1, 0]),
    ]:
        m = model(local_names=["a", "b"], wiggle=Scaling(basis, **coord))
        collection = collection_for(m, reg, {"a": "one", "b": "two"})
        results.append(evaluate_derivatives(collection, point, 2.4, k, mu))
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
    fs = [ObservedField("F", "forest", "same", background="qso")]
    m = KaiserModel(
        template(),
        fs,
        biases={"F": "b"},
        betas={"F": 1},
        widths={"F": (0, 0)},
        local_names=["b"],
    )
    reg = ParameterRegistry(
        [Parameter("b", -0.3, "nuisance", step=0.001), Parameter("f", 0.8, "target")]
    )
    collection = collection_for(m, reg, {"b": "b"})
    d = evaluate_derivatives(collection, reg.fiducials, 2.4, [0.1, 0.2], [0.3, 0.7])
    np.testing.assert_array_equal(d.jacobian[:, :, 1], 0)
    factors = factor_covariance(
        gaussian_covariance(d.power + 1, np.ones(2), collection.selection)
    )
    forecast = FisherResult(reg, fisher_from_factors(d.jacobian, factors))
    assert forecast.diagnostics.rank == 1 and np.isinf(
        forecast.conditional_errors(["f"])[0]
    )
    with pytest.raises(ValueError, match="rank"):
        forecast.marginalized_errors(["f"])
    with pytest.raises(ValueError, match="forbidden"):
        KaiserModel(
            template(),
            fs,
            biases={"F": -0.3},
            betas={"F": 1},
            widths={"F": (0, 0)},
            f=0.8,
        )


def test_synthetic_example_fixed_factors_and_step_convergence(monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "builtin_example",
        Path(__file__).resolve().parents[1] / "examples/builtin_forecast.py",
    )
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    calls = []
    original = example.factor_covariance

    def count(covariance):
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
    m = model(local_names=["a"], wiggle=Scaling(basis, **coordinates))
    with pytest.raises(ValueError, match="wiggle.*scal"):
        m(theta, 2.4, [0.1], [0.5], [[0, 1]])


def test_template_evaluated_once_per_component(monkeypatch):
    m = model()
    original = type(m.template).evaluate
    calls = []

    def counted(self, k, **kwargs):
        calls.append(np.array(k))
        return original(self, k, **kwargs)

    monkeypatch.setattr(type(m.template), "evaluate", counted)
    m([], 2.4, [0.1, 0.2], [0.3, 0.7], [[0, 0], [0, 1], [1, 1]])
    assert len(calls) == 2


def test_builtin_import_keeps_preparation_libraries_lazy(tmp_path):
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
