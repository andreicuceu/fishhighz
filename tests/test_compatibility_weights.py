"""Analytic controls for validation-only signed compatibility weights."""

import numpy as np
import pytest

from fishhighz.validation.compatibility_weights import (
    VARIANTS,
    WeightInputs,
    coefficients,
    fixed_weights,
    seed,
    update,
)


def inputs(density=(2.0, 3.0), variance=(1.0, 4.0), signal=2.0, p1d=3.0):
    """Construct a normalized two-sample compatibility-weight fixture.

    Parameters
    ----------
    density : array_like of shape (n_magnitudes,), optional
        Source density per deg^2, velocity interval in km/s, and magnitude.
        Default is (2.0, 3.0).
    variance : array_like of shape (n_magnitudes,), optional
        Dimensionless pixel-noise variance for each magnitude sample. Default is
        (1.0, 4.0).
    signal : float, optional
        Auxiliary three-dimensional signal power in deg^2 km/s. Default is 2.0.
    p1d : float, optional
        Auxiliary one-dimensional forest power in km/s. Default is 3.0.

    Returns
    -------
    inputs : WeightInputs
        Magnitude-ordered source densities and noise with unit quadrature,
        forest length 2 km/s, and pixel width 1 km/s.
    """
    return WeightInputs(
        np.arange(len(density)),
        np.array(density),
        np.ones(len(density)),
        np.array(variance),
        2.0,
        1.0,
        signal,
        p1d,
    )


@pytest.mark.parametrize("variant", VARIANTS)
def test_two_cell_simultaneous(variant):
    """Check two cell simultaneous.

    Parameters
    ----------
    variant : str
        Forest-weight recurrence variant, supplied by pytest parametrization.
    """
    weight_inputs = inputs(density=(-1.0, 3.0))
    source_weights = np.array([0.5, 0.25])
    # Independent explicit two-cell moments, including a negative first prefix.
    j1 = np.array([-0.5, 0.25])
    j2 = np.array([-0.25, -0.0625])
    if variant.startswith("sum"):
        j1, j2 = j1[-1], j2[-1]
    effective_signal = 2.0
    if variant.endswith("aliasing"):
        effective_signal += 3.0 * j2 / (2.0 * j1**2)
    elif variant.endswith("historical"):
        effective_signal += 3.0 / (2.0 * j1)
    expected = effective_signal / (effective_signal + np.array([1.0, 4.0]) / (2.0 * j1))
    np.testing.assert_allclose(
        update(weight_inputs, source_weights, variant), expected, rtol=2e-15
    )


@pytest.mark.parametrize("variant", VARIANTS)
def test_single_cell_fixed_point(variant):
    """Check single cell fixed point.

    Parameters
    ----------
    variant : str
        Forest-weight recurrence variant, supplied by pytest parametrization.
    """
    weight_inputs = inputs(density=(2.0,), variance=(1.0,), signal=2.0, p1d=3.0)
    # s=P*L*n=8; moment-aliasing s+B=11; historical positive quadratic root.
    if variant.endswith("intrinsic"):
        source_weights = 1 - 1 / 8
    elif variant.endswith("aliasing"):
        source_weights = 1 - 1 / 11
    else:
        source_weights = (4 + np.sqrt(112)) / 16
    np.testing.assert_allclose(
        update(weight_inputs, np.array([source_weights]), variant),
        [source_weights],
        rtol=2e-15,
    )


@pytest.mark.parametrize("variant", VARIANTS[1::2] + ("sum_historical",))
def test_equal_noise_and_permutation(variant):
    """Check equal noise and permutation.

    Parameters
    ----------
    variant : str
        Forest-weight recurrence variant, supplied by pytest parametrization.
    """
    weight_inputs = inputs(variance=(2.0, 2.0))
    source_weights = fixed_weights(weight_inputs, variant)
    np.testing.assert_allclose(source_weights[0], source_weights[1], rtol=2e-15)
    # Constant weights minimize both coefficients for this positive population.
    np.testing.assert_allclose(
        coefficients(weight_inputs, source_weights)[-2:], [0.1, 0.2]
    )
    permuted_inputs = inputs(density=(3.0, 2.0), variance=(4.0, 1.0))
    original_inputs = inputs()
    np.testing.assert_allclose(
        fixed_weights(permuted_inputs, variant)[::-1],
        fixed_weights(original_inputs, variant),
        rtol=2e-15,
    )


def test_seed_count_and_prefix_difference():
    """Check seed count and prefix difference."""
    weight_inputs = inputs(variance=(2.0, 2.0))
    np.testing.assert_array_equal(
        fixed_weights(weight_inputs, updates=0), seed(weight_inputs)
    )
    source_weights = fixed_weights(weight_inputs)
    assert source_weights[0] != source_weights[1]
    expected = seed(weight_inputs)
    for _ in range(3):
        expected = update(weight_inputs, expected, "prefix_intrinsic")
    np.testing.assert_array_equal(source_weights, expected)


def test_invalid_moment_is_not_regularized():
    """Check invalid moment is not regularized."""
    with pytest.raises(FloatingPointError):
        update(inputs(density=(-1.0, 1.0)), np.ones(2), "sum_intrinsic")


def test_noise_response_units_and_unfiltered_pixel_term():
    """Check noise response units and unfiltered pixel term."""
    from fishhighz.validation.compatibility_weights import forest_noise

    row = dict(
        _distance_to_velocity=100.0,
        _angle_to_distance=50.0,
        _pix_kms=20.0,
        _res_kms=15.0,
        _z_mean=3.0,
    )
    k_grid, mu_grid = np.array([0.1, 0.2]), np.array([0.5, 0.8])
    velocity_k = k_grid * mu_grid / 100.0
    floor = 0.009 * np.exp((-0.5 * -2.55 - 1) / -0.1)
    scaled_wavenumber = np.maximum(velocity_k, floor) / 0.009
    p1d = (
        np.pi
        * 0.064
        / 0.009
        * scaled_wavenumber ** (-0.55 - 0.1 * np.log(scaled_wavenumber))
    )
    window_squared = (np.sin(velocity_k * 10.0) / (velocity_k * 10.0)) ** 2 * np.exp(
        -((velocity_k * 15.0) ** 2)
    )
    expected = (2.0 * p1d * window_squared + 3.0) * 25.0
    np.testing.assert_allclose(
        forest_noise(row, k_grid, mu_grid, 2.0, 3.0), expected, rtol=2e-15
    )


@pytest.mark.parametrize("variant", ["sum_intrinsic", "sum_aliasing"])
def test_homogeneous_subthreshold_decay(variant):
    """Check homogeneous subthreshold decay.

    Parameters
    ----------
    variant : str
        Forest-weight recurrence variant, supplied by pytest parametrization.
    """
    weight_inputs = inputs(density=(1.0,), variance=(10.0,), signal=1.0, p1d=1.0)
    source_weights = seed(weight_inputs)
    for _ in range(12):
        next_w = update(weight_inputs, source_weights, variant)
        assert 0 < next_w[0] < source_weights[0]
        np.testing.assert_allclose(coefficients(weight_inputs, next_w)[-2:], [0.5, 5.0])
        source_weights = next_w
