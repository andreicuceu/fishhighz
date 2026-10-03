"""Integrated forest-source weights: exact central limits and dense oracles.

The integrated source is a synthetic stand-in exposing the attributes of
fishhighz.forest_integration.IntegratedForestSource, so these tests do not
depend on that module.
"""

from dataclasses import FrozenInstanceError, fields, replace
from types import SimpleNamespace

import numpy as np
import pytest
from test_forecast import forest_spec
from test_weights import FIELD, RESPONSE, geometry

from fishhighz.forecast import prepare_bin, run_bin
from fishhighz.noise import forest_noise
from fishhighz.public import _convergence
from fishhighz.survey import ForestInput
from fishhighz.weights import (
    ForestWeights,
    IntegratedForestWeights,
    prepare_forest_weights,
    prepare_integrated_forest_weights,
)

SIGNAL, ALIAS = 3.0, 2.0
LENGTH = 10.0
MAGNITUDES = np.array([20.0, 21.0, 23.0])
QUADRATURE = np.array([0.2, 0.5, 1.3])
RHO = np.array([0.05, 0.1, 0.2])
VARIANCE = np.array([1.0, 4.0, 9.0])


def synthetic_source(measure, variance, z_pix, pixel_width_velocity=0.5):
    """Build a stand-in integrated forest source.

    Parameters
    ----------
    measure : array_like of shape (n_pixel, n_magnitude)
        Integration measure in deg^-2.
    variance : array_like of shape (n_pixel, n_magnitude)
        Dimensionless pixel-noise variance.
    z_pix : array_like of shape (n_pixel,)
        Pixel redshifts.
    pixel_width_velocity : float, default=0.5
        Pixel width in km/s.

    Returns
    -------
    source : types.SimpleNamespace
        Object with the IntegratedForestSource attribute names; the density
        carries the full measure and the quadrature is unity.
    """
    measure = np.asarray(measure, dtype=float)
    n_pixel, n_magnitude = measure.shape
    nodes = SimpleNamespace(
        geom=np.ones(n_pixel),
        y_index=np.zeros(n_pixel, dtype=int),
        z_q=np.full(n_pixel, 3.0),
        lam_obs=np.full(n_pixel, 4000.0),
        z_pix=np.asarray(z_pix, dtype=float),
        zq_nodes=np.array([3.0]),
        info={},
    )
    return SimpleNamespace(
        nodes=nodes,
        density=measure[:1].copy(),
        magnitudes=MAGNITUDES[:n_magnitude].copy(),
        quadrature=np.ones(n_magnitude),
        variance=np.asarray(variance, dtype=float),
        pixel_width_velocity=pixel_width_velocity,
        info={},
        measure=measure,
    )


def central(**kwargs):
    """Prepare the central early-lyaforecast reference weights.

    Parameters
    ----------
    **kwargs : dict
        Overrides to prepare_forest_weights.

    Returns
    -------
    weights : ForestWeights
        Central weights for (RHO, QUADRATURE, LENGTH).
    """
    args = dict(
        z_source=3,
        magnitudes=MAGNITUDES,
        quadrature=QUADRATURE,
        rho=RHO,
        variance=VARIANCE,
        length_velocity=LENGTH,
        method="early_lyaforecast",
        signal=SIGNAL,
        alias=ALIAS,
    )
    args.update(kwargs)
    return prepare_forest_weights(FIELD, geometry(), RESPONSE, **args)


def integrated(source, **kwargs):
    """Prepare integrated early-lyaforecast weights with fixed S and B.

    Parameters
    ----------
    source : object
        Integrated forest source.
    **kwargs : dict
        Overrides to prepare_integrated_forest_weights.

    Returns
    -------
    weights : IntegratedForestWeights
        Integrated weights for the standard test field and response.
    """
    args = dict(method="early_lyaforecast", signal=SIGNAL, alias=ALIAS)
    args.update(kwargs)
    return prepare_integrated_forest_weights(
        FIELD, geometry(), RESPONSE, source, **args
    )


def assert_matches_central(prepared, reference, n_pixel):
    """Compare integrated weights with central weights at rtol 1e-12.

    Parameters
    ----------
    prepared : IntegratedForestWeights
        Integrated result.
    reference : ForestWeights
        Central result.
    n_pixel : int
        Number of identical pixels in the integrated source.
    """
    tolerance = dict(rtol=1e-12, atol=0)
    np.testing.assert_allclose(prepared.A, reference.A, **tolerance)
    np.testing.assert_allclose(prepared.P_pixel, reference.P_pixel, **tolerance)
    np.testing.assert_allclose(prepared.N1, LENGTH * reference.I1[-1], **tolerance)
    np.testing.assert_allclose(prepared.N2, LENGTH * reference.I2[-1], **tolerance)
    np.testing.assert_allclose(prepared.N3, LENGTH * reference.I3[-1], **tolerance)
    assert prepared.weights.shape == (n_pixel, len(MAGNITUDES))
    for row in prepared.weights:
        np.testing.assert_allclose(row, reference.weights, **tolerance)

    integrated_record, central_record = prepared.convergence, reference.convergence
    assert set(integrated_record) == set(central_record)
    for key in ("status", "updates", "state_updates", "candidate"):
        if key in central_record:
            assert integrated_record[key] == central_record[key]
    if central_record["status"] == "converged":
        for key in ("last_step", "confirmation"):
            for metric, value in central_record[key].items():
                # The metrics are differences of states that agree to rounding;
                # the dedicated recurrence differs from the central kernel by a
                # few units in the last place per weight (about 1e-15 absolute).
                np.testing.assert_allclose(
                    integrated_record[key][metric], value, rtol=1e-9, atol=1e-13
                )
        np.testing.assert_allclose(
            integrated_record["coefficients"][3:],
            central_record["coefficients"][3:],
            **tolerance,
        )
        np.testing.assert_allclose(
            integrated_record["weights"][0], central_record["weights"], **tolerance
        )


@pytest.mark.parametrize("iterations", [None, 4])
def test_single_pixel_equals_central(iterations):
    """A single pixel with measure L rho q reproduces the central weights."""
    measure = (LENGTH * RHO * QUADRATURE)[None, :]
    source = synthetic_source(measure, VARIANCE[None, :], [2.4])
    prepared = integrated(source, iterations=iterations)
    reference = central(iterations=iterations)
    assert_matches_central(prepared, reference, 1)
    np.testing.assert_allclose(prepared.z_eff, 2.4, rtol=1e-12)
    assert prepared.iterations == iterations


@pytest.mark.parametrize("iterations", [None, 3])
def test_identical_pixels_equal_central(iterations):
    """Pixel copies whose measures sum to L rho q reproduce the central weights."""
    fractions = np.array([0.2, 0.3, 0.5])
    z_pix = np.array([2.2, 2.4, 2.6])
    measure = fractions[:, None] * (LENGTH * RHO * QUADRATURE)[None, :]
    variance = np.tile(VARIANCE, (3, 1))
    prepared = integrated(
        synthetic_source(measure, variance, z_pix), iterations=iterations
    )
    assert_matches_central(prepared, central(iterations=iterations), 3)
    np.testing.assert_allclose(prepared.z_eff, np.sum(fractions * z_pix), rtol=1e-12)


def dense_reference(measure, variance, pixel, updates):
    """Evaluate seed and updates of the integrated recurrence with explicit loops.

    Parameters
    ----------
    measure, variance : ndarray of shape (n_pixel, n_magnitude)
        Measure in deg^-2 and dimensionless pixel variance.
    pixel : float
        Pixel width in km/s.
    updates : int
        Number of updates after the seed.

    Returns
    -------
    weights : ndarray of shape (n_pixel, n_magnitude)
        Dimensionless weights.
    moments : tuple of float
        N1, N2, N3.
    """
    n_pixel, n_magnitude = measure.shape
    weights = np.empty(measure.shape)
    for i in range(n_pixel):
        for j in range(n_magnitude):
            weights[i, j] = (ALIAS / pixel) / (ALIAS / pixel + variance[i, j])

    def totals(w):
        n1 = n2 = n3 = 0.0
        for i in range(n_pixel):
            for j in range(n_magnitude):
                n1 += measure[i, j] * w[i, j]
                n2 += measure[i, j] * w[i, j] ** 2
                n3 += measure[i, j] * w[i, j] ** 2 * variance[i, j]
        return n1, n2, n3

    for _ in range(updates):
        n1, _, _ = totals(weights)
        signal = SIGNAL + ALIAS / n1
        new = np.empty_like(weights)
        for i in range(n_pixel):
            for j in range(n_magnitude):
                new[i, j] = signal / (signal + pixel * variance[i, j] / n1)
        weights = new
    return weights, totals(weights)


def random_source(seed=3, n_pixel=4, n_magnitude=3):
    """Draw a small random integrated source with one empty cell.

    Parameters
    ----------
    seed : int, default=3
        Seed of the random generator.
    n_pixel, n_magnitude : int
        Array dimensions.

    Returns
    -------
    source : types.SimpleNamespace
        Synthetic integrated source.
    """
    rng = np.random.default_rng(seed)
    measure = rng.uniform(0.5, 5.0, (n_pixel, n_magnitude))
    measure[1, 2] = 0.0
    variance = rng.uniform(0.5, 6.0, (n_pixel, n_magnitude))
    return synthetic_source(measure, variance, np.linspace(2.1, 2.7, n_pixel))


def test_dense_loop_oracle_fixed_and_adaptive():
    """The kernels match an explicit double loop on a random (n_p, n_m) source."""
    source = random_source()
    pixel = RESPONSE.pixel_width_velocity

    fixed = integrated(source, iterations=5)
    weights, (n1, n2, n3) = dense_reference(source.measure, source.variance, pixel, 5)
    np.testing.assert_allclose(fixed.weights, weights, rtol=1e-12)
    np.testing.assert_allclose([fixed.N1, fixed.N2, fixed.N3], [n1, n2, n3], rtol=1e-12)
    np.testing.assert_allclose(fixed.A, n2 / n1**2, rtol=1e-12)
    np.testing.assert_allclose(fixed.P_pixel, pixel * n3 / n1**2, rtol=1e-12)
    assert fixed.convergence == dict(status="fixed_count", updates=5)

    adaptive = integrated(source)
    assert adaptive.convergence["status"] == "converged"
    weights, (n1, n2, n3) = dense_reference(
        source.measure, source.variance, pixel, adaptive.convergence["updates"]
    )
    np.testing.assert_allclose(adaptive.weights, weights, rtol=1e-12)
    np.testing.assert_allclose(adaptive.A, n2 / n1**2, rtol=1e-12)
    np.testing.assert_allclose(
        adaptive.z_eff,
        np.sum(source.nodes.z_pix[:, None] * source.measure * weights) / n1,
        rtol=1e-12,
    )
    assert adaptive.weights.shape == source.measure.shape
    assert adaptive.weights[1, 2] >= 0


@pytest.mark.parametrize(
    "method", ["legacy", "mcdonald", "inverse_variance", "supplied", "other"]
)
def test_only_early_lyaforecast(method):
    """Any method other than early_lyaforecast raises."""
    with pytest.raises(ValueError, match="early_lyaforecast"):
        integrated(random_source(), method=method, iterations=3)


def test_invalid_inputs():
    """Invalid sources, pixel widths and S/B inputs raise ValueError."""
    source = random_source()
    with pytest.raises(ValueError, match="pixel width"):
        integrated(replace_source(source, pixel_width_velocity=0.6), iterations=2)
    with pytest.raises(ValueError, match="lacks attributes"):
        integrated(SimpleNamespace(measure=source.measure), iterations=2)
    with pytest.raises(ValueError, match="measure support"):
        integrated(replace_source(source, measure=np.zeros_like(source.measure)))
    with pytest.raises(ValueError):
        integrated(source, signal=None, iterations=2)
    with pytest.raises(ValueError):
        integrated(replace_source(source, variance=-source.variance), iterations=2)


def replace_source(source, **changes):
    """Copy a synthetic source with some attributes replaced.

    Parameters
    ----------
    source : types.SimpleNamespace
        Synthetic integrated source.
    **changes : dict
        Replacement attributes.

    Returns
    -------
    source : types.SimpleNamespace
        Modified shallow copy.
    """
    return SimpleNamespace(**{**vars(source), **changes})


def test_noise_context_and_immutability():
    """forest_noise accepts both classes, rejects context mismatch; arrays frozen."""
    prepared = integrated(random_source(), iterations=4)
    assert isinstance(prepared, IntegratedForestWeights)
    k, mu, p1d = [0.1, 0.2], [0.3, 0.6], [1.0, 2.0]

    noise = forest_noise(prepared, FIELD, geometry(), RESPONSE, k, mu, p1d)
    reference = forest_noise(
        _central_like(prepared), FIELD, geometry(), RESPONSE, k, mu, p1d
    )
    np.testing.assert_array_equal(noise.total, reference.total)
    assert np.all(np.isfinite(noise.total)) and np.all(noise.total > 0)

    with pytest.raises(ValueError, match="mismatch"):
        forest_noise(prepared, FIELD, geometry(h=0.8), RESPONSE, k, mu, p1d)
    other_field = replace(FIELD, id="other")
    with pytest.raises(ValueError, match="mismatch"):
        forest_noise(prepared, other_field, geometry(), RESPONSE, k, mu, p1d)
    with pytest.raises(ValueError, match="require ForestWeights"):
        forest_noise(object(), FIELD, geometry(), RESPONSE, k, mu, p1d)
    prepared.validate_context(FIELD, geometry(), RESPONSE)

    assert not prepared.weights.flags.writeable
    with pytest.raises(ValueError):
        prepared.weights[0, 0] = 1.0
    with pytest.raises(FrozenInstanceError):
        prepared.A = 1.0
    assert isinstance(prepared.convergence, dict)
    assert prepared.convergence["status"] == "fixed_count"
    assert isinstance(integrated(random_source()).convergence, dict)


def _central_like(prepared):
    """Build a ForestWeights carrying an integrated result's A, P_pixel, context.

    Parameters
    ----------
    prepared : IntegratedForestWeights
        Integrated result.

    Returns
    -------
    weights : ForestWeights
        Object whose noise coefficients equal those of ``prepared``.
    """
    result = object.__new__(ForestWeights)
    for name in ("context", "A", "P_pixel"):
        object.__setattr__(result, name, getattr(prepared, name))
    return result


def integrated_forest_spec(integrated_source):
    """Build a forest bin whose forest carries an integrated source.

    Parameters
    ----------
    integrated_source : object or None
        Integrated source, or None for the central ForestInput.

    Returns
    -------
    spec : BinSpec
        Synthetic one-forest bin; its pixel width is 30 km/s.
    """
    spec, _ = forest_spec("legacy", auxiliary=True)
    old = spec.forests["f"]
    options = dict(method="early_lyaforecast", rtol=1e-6)
    if integrated_source is None:
        options.update(
            z_source=4,
            magnitudes=[20, 21],
            quadrature=[1, 1],
            rho=[0.01, 0.02],
            variance=[1, 2],
            length_velocity=10000,
        )
    source = ForestInput(
        options,
        old.p1d_model,
        old.p1d_parameters,
        old.theta_p1d,
        auxiliary_coordinates=old.auxiliary_coordinates,
        integrated=integrated_source,
    )
    return replace(spec, forests={"f": source})


def test_prepare_bin_integrated_and_central_unchanged():
    """prepare_bin dispatches on ForestInput.integrated; central path is unchanged."""
    measure = np.array([[40.0, 60.0], [50.0, 50.0], [30.0, 70.0]])
    variance = np.array([[1.0, 2.0], [1.5, 2.5], [1.0, 3.0]])
    source = synthetic_source(measure, variance, [2.3, 2.5, 2.7], 30.0)
    spec = integrated_forest_spec(source)
    prepared = prepare_bin(spec)

    weights = prepared.weights["f"]
    assert isinstance(weights, IntegratedForestWeights)
    assert weights.convergence["status"] == "converged"
    assert np.all(np.isfinite(prepared.noise)) and np.all(prepared.noise > 0)
    assert 2.3 < weights.z_eff < 2.7
    run_bin(prepared, batch_size=1)

    # The public convergence summary omits the large weight array.
    record = _convergence([prepared])[prepared.id]["f"]
    assert record["status"] == "converged"
    assert "weights" not in record and record["weights_shape"] == (3, 2)

    # Central path: reconstruct the noise directly from the public kernels.
    central_spec = integrated_forest_spec(None)
    central_bin = prepare_bin(central_spec)
    assert isinstance(central_bin.weights["f"], ForestWeights)
    options = dict(central_spec.forests["f"].weight_options)
    central_weights = central_bin.weights["f"]
    expected = prepare_forest_weights(
        central_spec.p3d.selection.fields[0],
        central_spec.geometry,
        central_spec.responses["f"],
        auxiliary=central_weights.auxiliary,
        **options,
    )
    np.testing.assert_array_equal(expected.weights, central_weights.weights)
    assert expected.A == central_weights.A
    assert expected.P_pixel == central_weights.P_pixel
    assert central_spec.forests["f"].integrated is None


def test_forest_input_integrated_validation():
    """ForestInput validates the integrated attributes and weight options."""
    spec = integrated_forest_spec(None)
    old = spec.forests["f"]
    source = random_source()
    with pytest.raises(ValueError, match="lacks attributes"):
        replace(old, integrated=SimpleNamespace(measure=1))
    with pytest.raises(ValueError, match="early_lyaforecast"):
        replace(
            old,
            weight_options=dict(method="mcdonald"),
            integrated=source,
        )
    with pytest.raises(ValueError, match="unexpected keys"):
        replace(
            old,
            weight_options=dict(method="early_lyaforecast", rho=[1.0]),
            integrated=source,
        )
    assert (
        replace(old, weight_options=dict(method="early_lyaforecast")).integrated is None
    )
    names = [f.name for f in fields(ForestInput)]
    assert names[:6] == [
        "weight_options",
        "p1d_model",
        "p1d_parameters",
        "theta_p1d",
        "auxiliary_coordinates",
        "provenance",
    ]
    assert names[6:] == ["integrated"]
