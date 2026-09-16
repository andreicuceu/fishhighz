"""Fixed-reference coefficients, retained arithmetic guards and provider ownership."""

from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
from math import exp, sin

import numpy as np
import pytest
from test_forecast import forest_spec
from test_weights import FIELD, geometry

from fishhighz.fields import ObservedField
from fishhighz.forecast import prepare_bin, run_bin
from fishhighz.models.external import P3DProvider, PreparedP3D
from fishhighz.noise import forest_noise
from fishhighz.response import InstrumentResponse
from fishhighz.weights import ForestWeights, prepare_forest_weights

RESPONSE = InstrumentResponse(1, 0)


def prepare(field=FIELD, response=RESPONSE, **changes):
    options = dict(
        z_source=3,
        magnitudes=[20, 21],
        quadrature=[0.5, 2],
        rho=[2, 0.5],
        variance=[1, 4],
        length_velocity=1,
        method="inverse_variance",
        alias=1,
    )
    options.update(changes)
    return prepare_forest_weights(field, geometry(), response, **options)


def agree(actual, expected):
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    np.testing.assert_allclose(actual, expected, rtol=5e-12, atol=0)


def test_rational_coefficients_and_legacy_seed(monkeypatch):
    legacy = prepare(method="legacy", iterations=0, signal=7)
    # Independent exact rational control, including the nonuniform measure once.
    r = [Fraction(2) * Fraction(1, 2), Fraction(1, 2) * 2]
    nu = [Fraction(1, 2), Fraction(1, 5)]
    i1 = sum(m * w for m, w in zip(r, nu))
    i2 = sum(m * w**2 for m, w in zip(r, nu))
    i3 = sum(m * w**2 * v for m, w, v in zip(r, nu, [1, 4]))
    assert (i1, i2, i3) == (Fraction(7, 10), Fraction(29, 100), Fraction(41, 100))
    assert i2 / i1**2 == Fraction(29, 49)
    assert i3 / i1**2 == Fraction(41, 49)
    assert (i2 + i3) / i1**2 == Fraction(10, 7)
    import fishhighz.weights as module

    def forbidden(*args, **kwargs):
        pytest.fail("inverse_variance invoked an iteration or provider")

    for name in ["_iterate", "sample_auxiliary", "evaluate_p1d", "evaluate_p3d"]:
        monkeypatch.setattr(module, name, forbidden)
    result = prepare()
    assert isinstance(result, ForestWeights)
    supplied = prepare(method="supplied", alias=None, weights=[float(w) for w in nu])
    agree(result.weights, [float(w) for w in nu])
    agree(
        [result.A, result.P_pixel, result.A + result.P_pixel],
        [29 / 49, 41 / 49, 10 / 7],
    )
    for other in [legacy, supplied]:
        for name in ["weights", "I1", "I2", "I3", "A", "P_pixel"]:
            agree(getattr(result, name), getattr(other, name))
    agree(
        [result.I1[-1], result.I2[-1], result.I3[-1]], [float(i1), float(i2), float(i3)]
    )


@pytest.mark.parametrize("background", ["qso", "lbg"])
def test_support_metadata_immutability_and_context(background):
    field = ObservedField(f"lya({background})", "forest", "lya", background)
    variance = np.array([0.0, 4.0])
    result = prepare(field, variance=variance, rho=[2, 0])
    np.testing.assert_array_equal(result.weights, [1, 0])
    np.testing.assert_array_equal(result.I3, [0, 0])
    assert result.P_pixel == 0
    assert result.method == "inverse_variance" and result.alias == 1
    assert result.signal is result.iterations is result.auxiliary is None
    assert result.weight_changes.shape == (0,)
    variance[:] = 9
    np.testing.assert_array_equal(result.variance, [0, 4])
    for name in [
        "weights",
        "rho",
        "variance",
        "magnitudes",
        "quadrature",
        "I1",
        "I2",
        "I3",
        "weight_changes",
    ]:
        with pytest.raises(ValueError):
            getattr(result, name).flags.writeable = True
    with pytest.raises(FrozenInstanceError):
        result.alias = 2
    result.validate_context(field, geometry(), RESPONSE)
    for f, g, response in [
        (ObservedField("other", "forest", "lya", background), geometry(), RESPONSE),
        (field, geometry(h=0.8), RESPONSE),
        (field, geometry(), InstrumentResponse(2, 0)),
    ]:
        with pytest.raises(ValueError, match="mismatch"):
            result.validate_context(f, g, response)


@pytest.mark.parametrize(
    "change",
    [
        dict(weights=[1, 1]),
        dict(iterations=0),
        dict(signal=1),
        dict(auxiliary=object()),
    ],
)
def test_conflicting_settings(change):
    with pytest.raises(ValueError, match="inverse_variance requires.*None"):
        prepare(**change)


@pytest.mark.parametrize("alias", [None, 0, -1, np.inf, -np.inf, np.nan, [1], True, 1j])
def test_invalid_alias(alias):
    with pytest.raises(ValueError):
        prepare(alias=alias)


@pytest.mark.parametrize(
    "change",
    [
        dict(alias=1e-300, variance=[0, 1e100]),  # zero computed weight
        dict(alias=1e308, variance=[0, 1e308]),  # denominator overflow
        dict(response=InstrumentResponse(2, 0), variance=[0, 1e308]),
        dict(variance=[0, 1e200]),  # mixed r*w*w underflow, positive total
        dict(rho=[2, 1e-200], variance=[0, 1e200]),  # mixed r*w underflow
        dict(quadrature=[1e-300, 2], rho=[1e-300, 1]),
    ],
)
def test_unrepresentable_weights_and_mixed_products(change):
    with pytest.raises(ValueError, match="not representable"):
        prepare(**change)


def test_positive_instrumental_power_underflow():
    tiny = np.nextafter(0.0, 1.0)
    # R1 accepted weights=(1,1), A=4, P_pixel=tiny despite losing tiny/2.
    with pytest.raises(ValueError, match="forest:.*not representable.*underflow"):
        prepare(
            response=InstrumentResponse(tiny, 0),
            alias=tiny,
            quadrature=[1, 1],
            rho=[1 / 8, 1 / 8],
            variance=[0, 1 / 2],
        )


def test_exact_zero_instrumental_power():
    tiny = np.nextafter(0.0, 1.0)
    result = prepare(
        response=InstrumentResponse(tiny, 0),
        alias=tiny,
        quadrature=[1, 1],
        rho=[1 / 8, 1 / 8],
        variance=[0, 0],
    )
    np.testing.assert_array_equal(result.weights, [1, 1])
    np.testing.assert_array_equal(result.I3, [0, 0])
    assert result.A == 4 and result.P_pixel == 0


def test_mode_dependent_noise_scalar_response():
    response = InstrumentResponse(30, 10)
    result = prepare(response=response, alias=2)
    supplied = prepare(
        response=response, method="supplied", alias=None, weights=[2 / 32, 2 / 122]
    )
    g = geometry()
    q = np.array([0.003, 0.02])
    p1d = np.array([4.0, 9.0])  # intrinsic values distinct from B_star
    args = (FIELD, g, response, q * g.a_v, np.ones(2), p1d)
    noise = forest_noise(result, *args)
    other = forest_noise(supplied, *args)
    expected_alias = []
    for qi, pi in zip(q, p1d):
        x = qi * 30 / 2
        w = sin(x) / x * exp(-((qi * 10) ** 2) / 2)
        expected_alias.append(result.A * pi * w**2 * g.d_deg**2 / g.a_v)
    agree(noise.aliasing, expected_alias)
    agree(noise.pixel, np.full(2, result.P_pixel * g.d_deg**2 / g.a_v))
    agree(noise.total, np.array(expected_alias) + noise.pixel)
    agree(noise.total, other.total)


def test_survey_equivalence_and_frozen_provider_counts(monkeypatch):
    import fishhighz.forecast as module

    spec, p1d_calls = forest_spec("supplied", auxiliary=False)
    source = spec.forests["f"]
    base = {
        k: v for k, v in source.weight_options.items() if k not in ("method", "weights")
    }
    p3d_calls = []
    owner = spec.p3d.routes[0].provider

    def model(*args):
        p3d_calls.append(1)
        return owner.model(*args)

    provider = P3DProvider(owner.label, model, owner.parameters, [(0, 0)])
    spec = replace(
        spec, p3d=PreparedP3D(spec.p3d.registry, spec.p3d.selection, [provider])
    )
    bins = []
    for options in [
        dict(method="inverse_variance", alias=2),
        dict(method="supplied", weights=[2 / 32, 2 / 62]),
    ]:
        source_i = replace(source, weight_options=dict(base, **options))
        b = prepare_bin(replace(spec, forests={"f": source_i}))
        assert b.diagnostics["p1d_calls"] == {"f": 1}
        assert b.diagnostics["p3d_calls"] == {"scalar": 1}
        bins.append(b)
    assert len(p1d_calls) == len(p3d_calls) == 2
    for name in ["weights", "I1", "I2", "I3", "A", "P_pixel"]:
        agree(getattr(bins[0].weights["f"], name), getattr(bins[1].weights["f"], name))
    agree(bins[0].noise, bins[1].noise)
    inverse_source = replace(
        source, weight_options=dict(base, method="inverse_variance", alias=2)
    )
    with pytest.raises(ValueError, match="auxiliary coordinates conflict"):
        replace(inverse_source, auxiliary_coordinates=(24, 0.00035))

    def forbidden(*args, **kwargs):
        pytest.fail("fixed weights, noise, response or P1D reevaluated")

    for name in [
        "sample_auxiliary",
        "prepare_forest_weights",
        "forest_noise",
        "evaluate_p1d",
        "prepare_noise",
        "prepare_response",
    ]:
        monkeypatch.setattr(module, name, forbidden)
    snapshots = [
        (
            b.noise.tobytes(),
            b.weights["f"].weights.tobytes(),
            b.weights["f"].A,
            b.weights["f"].P_pixel,
        )
        for b in bins
    ]
    results = [run_bin(b) for b in bins]
    agree(results[0].result.data_fisher, results[1].result.data_fisher)
    assert len(p1d_calls) == 2
    assert len(p3d_calls) == 12  # two preparations, then 5 mean calls per bin
    assert snapshots == [
        (
            b.noise.tobytes(),
            b.weights["f"].weights.tobytes(),
            b.weights["f"].A,
            b.weights["f"].P_pixel,
        )
        for b in bins
    ]
