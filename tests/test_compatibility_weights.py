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
    x = inputs(density=(-1.0, 3.0))
    w = np.array([0.5, 0.25])
    # Independent explicit two-cell moments, including a negative first prefix.
    j1 = np.array([-0.5, 0.25])
    j2 = np.array([-0.25, -0.0625])
    if variant.startswith("sum"):
        j1, j2 = j1[-1], j2[-1]
    s = 2.0
    if variant.endswith("aliasing"):
        s += 3.0 * j2 / (2.0 * j1**2)
    elif variant.endswith("historical"):
        s += 3.0 / (2.0 * j1)
    expected = s / (s + np.array([1.0, 4.0]) / (2.0 * j1))
    np.testing.assert_allclose(update(x, w, variant), expected, rtol=2e-15)


@pytest.mark.parametrize("variant", VARIANTS)
def test_single_cell_fixed_point(variant):
    x = inputs(density=(2.0,), variance=(1.0,), signal=2.0, p1d=3.0)
    # s=P*L*n=8; moment-aliasing s+B=11; historical positive quadratic root.
    if variant.endswith("intrinsic"):
        w = 1 - 1 / 8
    elif variant.endswith("aliasing"):
        w = 1 - 1 / 11
    else:
        w = (4 + np.sqrt(112)) / 16
    np.testing.assert_allclose(update(x, np.array([w]), variant), [w], rtol=2e-15)


@pytest.mark.parametrize("variant", VARIANTS[1::2] + ("sum_historical",))
def test_equal_noise_and_permutation(variant):
    x = inputs(variance=(2.0, 2.0))
    w = fixed_weights(x, variant)
    np.testing.assert_allclose(w[0], w[1], rtol=2e-15)
    # Constant weights minimize both coefficients for this positive population.
    np.testing.assert_allclose(coefficients(x, w)[-2:], [0.1, 0.2])
    y = inputs(density=(3.0, 2.0), variance=(4.0, 1.0))
    z = inputs()
    np.testing.assert_allclose(
        fixed_weights(y, variant)[::-1], fixed_weights(z, variant), rtol=2e-15
    )


def test_seed_count_and_prefix_difference():
    x = inputs(variance=(2.0, 2.0))
    np.testing.assert_array_equal(fixed_weights(x, updates=0), seed(x))
    w = fixed_weights(x)
    assert w[0] != w[1]
    expected = seed(x)
    for _ in range(3):
        expected = update(x, expected, "prefix_intrinsic")
    np.testing.assert_array_equal(w, expected)


def test_invalid_moment_is_not_regularized():
    with pytest.raises(FloatingPointError):
        update(inputs(density=(-1.0, 1.0)), np.ones(2), "sum_intrinsic")


def test_noise_response_units_and_unfiltered_pixel_term():
    from fishhighz.validation.compatibility_weights import forest_noise

    row = dict(
        _distance_to_velocity=100.0,
        _angle_to_distance=50.0,
        _pix_kms=20.0,
        _res_kms=15.0,
        _z_mean=3.0,
    )
    k, mu = np.array([0.1, 0.2]), np.array([0.5, 0.8])
    q = k * mu / 100.0
    floor = 0.009 * np.exp((-0.5 * -2.55 - 1) / -0.1)
    x = np.maximum(q, floor) / 0.009
    p1d = np.pi * 0.064 / 0.009 * x ** (-0.55 - 0.1 * np.log(x))
    window_squared = (np.sin(q * 10.0) / (q * 10.0)) ** 2 * np.exp(-((q * 15.0) ** 2))
    expected = (2.0 * p1d * window_squared + 3.0) * 25.0
    np.testing.assert_allclose(forest_noise(row, k, mu, 2.0, 3.0), expected, rtol=2e-15)


@pytest.mark.parametrize("variant", ["sum_intrinsic", "sum_aliasing"])
def test_homogeneous_subthreshold_decay(variant):
    x = inputs(density=(1.0,), variance=(10.0,), signal=1.0, p1d=1.0)
    w = seed(x)
    for _ in range(12):
        next_w = update(x, w, variant)
        assert 0 < next_w[0] < w[0]
        np.testing.assert_allclose(coefficients(x, next_w)[-2:], [0.5, 5.0])
        w = next_w
