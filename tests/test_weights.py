"""Independent scalar checks of finite cumulative weighting and input contracts."""

import numpy as np
import pytest

from fishhighz.fields import ObservedField
from fishhighz.geometry import prepare_geometry
from fishhighz.response import InstrumentResponse
from fishhighz.weights import density_per_velocity, prepare_forest_weights


def geometry(h=0.7, area=100):
    return prepare_geometry(
        2,
        3,
        z_eval=2.4,
        area_deg2=area,
        h_fid=h,
        hubble=lambda z: np.full_like(z, 250),
        transverse_distance=lambda z: np.full_like(z, 5500),
        z_order=4,
    )


FIELD = ObservedField("forest", "forest", "lya", "qso")
RESPONSE = InstrumentResponse(0.5, 0)


def prepare(**kwargs):
    args = dict(
        z_source=3,
        magnitudes=[20, 21, 23],
        quadrature=[0.2, 0.5, 1.3],
        rho=np.array([1, 2, 4]) / [0.2, 0.5, 1.3],
        variance=[1, 4, 9],
        length_velocity=10,
        method="legacy",
        iterations=3,
        signal=3,
        alias=2,
    )
    args.update(kwargs)
    return prepare_forest_weights(FIELD, geometry(), RESPONSE, **args)


def oracle(r, v, count):
    w = [2 / (2 + 0.5 * x) if mass else 0 for mass, x in zip(r, v)]
    changes = []
    for _ in range(count):
        next_w = []
        for m in range(len(r)):
            prefix = sum(r[j] * w[j] for j in range(m + 1))
            next_w.append(3 / (3 + v[m] / (prefix * 10 / 0.5)) if r[m] else 0)
        changes.append(max(abs(a - b) for a, b in zip(w, next_w)))
        w = next_w
    integrals = [
        [
            sum(
                r[j] * w[j] ** power * (v[j] if power == 3 else 1) for j in range(m + 1)
            )
            for m in range(len(r))
        ]
        for power in (1, 2)
    ]
    integrals.append(
        [sum(r[j] * w[j] ** 2 * v[j] for j in range(m + 1)) for m in range(len(r))]
    )
    return w, integrals, changes


@pytest.mark.parametrize("count", [0, 1, 2, 3, 6])
def test_scalar_prefix_oracle(count):
    result = prepare(iterations=count)
    w, integrals, changes = oracle([1, 2, 4], [1, 4, 9], count)
    np.testing.assert_allclose(result.weights, w, rtol=5e-15)
    np.testing.assert_allclose([result.I1, result.I2, result.I3], integrals, rtol=5e-15)
    np.testing.assert_allclose(result.weight_changes, changes, rtol=5e-14, atol=2e-16)
    assert result.A == pytest.approx(
        integrals[1][-1] / integrals[0][-1] ** 2 / 10, rel=5e-15
    )
    assert result.P_pixel == pytest.approx(
        integrals[2][-1] * 0.5 / integrals[0][-1] ** 2 / 10, rel=5e-15
    )


def test_distinguishes_wrong_iterations():
    w = np.array(oracle([1, 2, 4], [1, 4, 9], 0)[0])
    total = 3 / (3 + np.array([1, 4, 9]) / (sum(w * [1, 2, 4]) * 20))
    sweep = w.copy()
    for m in range(3):
        sweep[m] = 3 / (
            3 + [1, 4, 9][m] / (sum(sweep[: m + 1] * np.array([1, 2, 4])[: m + 1]) * 20)
        )
    assert not np.allclose(prepare(iterations=1).weights, total)
    assert not np.allclose(prepare(iterations=1).weights, sweep)


@pytest.mark.parametrize("variance", [0, 7])
@pytest.mark.parametrize("weight", [0.01, 1, 19])
def test_single_population(variance, weight):
    p = prepare(
        magnitudes=[21],
        quadrature=[0.3],
        rho=[10],
        variance=[variance],
        method="supplied",
        weights=[weight],
        iterations=None,
        signal=None,
        alias=None,
    )
    assert p.A == pytest.approx(1 / 30, rel=5e-15)
    assert p.P_pixel == pytest.approx(variance * 0.5 / 30, rel=5e-15)


def test_scaling_and_immutable():
    kw = dict(
        method="supplied",
        weights=[0.3, 0.7, 2],
        iterations=None,
        signal=None,
        alias=None,
    )
    p = prepare(**kw)
    for scale in [0.001, 100]:
        other = prepare(**dict(kw, weights=np.array(kw["weights"]) * scale))
        np.testing.assert_allclose(
            [p.A, p.P_pixel], [other.A, other.P_pixel], rtol=5e-15
        )
    other = prepare(**kw, rho=np.array([1, 2, 4]) / [0.2, 0.5, 1.3] * 2)
    np.testing.assert_allclose(
        [p.A, p.P_pixel], 2 * np.array([other.A, other.P_pixel]), rtol=5e-15
    )
    other = prepare(**kw, length_velocity=20)
    np.testing.assert_allclose(
        [p.A, p.P_pixel], 2 * np.array([other.A, other.P_pixel]), rtol=5e-15
    )
    for name in (
        "weights",
        "I1",
        "I2",
        "I3",
        "rho",
        "variance",
        "magnitudes",
        "quadrature",
        "weight_changes",
    ):
        with pytest.raises(ValueError):
            getattr(p, name).flags.writeable = True
    uniform = prepare(**dict(kw, weights=[1, 1, 1]), variance=[2, 2, 2])
    assert uniform.A == pytest.approx(1 / 70, rel=5e-15)
    assert uniform.P_pixel / uniform.A == pytest.approx(1, rel=5e-15)


@pytest.mark.parametrize("r,v", [([0, 2, 4], [9, 0, 1]), ([1, 0, 4], [0, 9, 1])])
def test_zero_support_prefixes(r, v):
    p = prepare(rho=np.array(r) / [0.2, 0.5, 1.3], variance=v)
    w, ints, _ = oracle(r, v, 3)
    np.testing.assert_allclose(p.weights, w, rtol=5e-15)
    np.testing.assert_allclose([p.I1, p.I2, p.I3], ints, rtol=5e-15)
    assert np.all(p.weights[np.array(r) == 0] == 0)


@pytest.mark.parametrize(
    "change",
    [
        dict(z_source=2.4),
        dict(magnitudes=[21, 20, 23]),
        dict(magnitudes=[]),
        dict(quadrature=[1, 0, 1]),
        dict(rho=[0, 0, 0]),
        dict(rho=[-1, 1, 1]),
        dict(variance=[1, np.nan, 1]),
        dict(iterations=True),
        dict(iterations=-1),
        dict(iterations=1.2),
        dict(signal=0),
        dict(alias=-1),
        dict(weights=[1, 1, 1]),
        dict(method="guess"),
        dict(rho=[1e308] * 3),
        dict(
            method="supplied",
            iterations=None,
            signal=None,
            alias=None,
            weights=[0, 0, 0],
        ),
        dict(
            method="supplied",
            iterations=None,
            signal=None,
            alias=None,
            weights=[1e200] * 3,
        ),
    ],
)
def test_invalid(change):
    with pytest.raises(ValueError):
        prepare(**change)


def test_density_conversion_and_context():
    row = np.array([0, 10, 3.2])
    np.testing.assert_allclose(
        density_per_velocity(row, z_source=3.1), row * 4.1 / 299792.458, rtol=5e-15
    )
    p = prepare()
    p.validate_context(FIELD, geometry(area=200), RESPONSE)
    for f, g, response in [
        (FIELD, geometry(h=0.8), RESPONSE),
        (FIELD, geometry(), InstrumentResponse(1, 0)),
        (ObservedField("other", "forest", "lya", "qso"), geometry(), RESPONSE),
    ]:
        with pytest.raises(ValueError, match="mismatch"):
            p.validate_context(f, g, response)


@pytest.mark.parametrize(
    "row,z",
    [
        ([-1], 3),
        ([np.nan], 3),
        ([1], -1),
        ([1], np.inf),
        ([], 3),
        ([[1]], 3),
        ([1e-323], 0),
    ],
)
def test_density_invalid(row, z):
    with pytest.raises(ValueError):
        density_per_velocity(row, z_source=z)


@pytest.mark.parametrize(
    "change",
    [
        dict(variance=[-1, 1, 1]),
        dict(quadrature=[1, 1]),
        dict(quadrature=[1e-300] * 3, rho=[1e-300] * 3),
        dict(length_velocity=0),
        dict(magnitudes=[1, np.inf, 3]),
        dict(alias=np.inf),
        dict(signal=np.nan),
        dict(method="supplied", weights=[1, 1, 1]),
        dict(
            method="supplied",
            weights=[1e-300] * 3,
            iterations=None,
            signal=None,
            alias=None,
        ),
        dict(
            method="supplied",
            weights=[-1, 1, 1],
            iterations=None,
            signal=None,
            alias=None,
        ),
        dict(
            method="supplied", weights=[1, 1], iterations=None, signal=None, alias=None
        ),
    ],
)
def test_additional_invalid_boundaries(change):
    with pytest.raises(ValueError):
        prepare(**change)
