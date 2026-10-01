"""Decimal oracles for mixed contributions concealed by positive totals (R1)."""

from decimal import Decimal, localcontext

import numpy as np
import pytest
from test_weights import FIELD, geometry

from fishhighz.response import InstrumentResponse
from fishhighz.weights import prepare_forest_weights


def prepare(rho, weights, variance):
    """Prepare supplied weights over extreme numerical ranges.

    Parameters
    ----------
    rho : array_like of shape (n_magnitudes,)
        Source density per deg^2, velocity interval in km/s, and magnitude.
    weights : array_like of shape (n_magnitudes,)
        Dimensionless source weights in magnitude order.
    variance : array_like of shape (n_magnitudes,)
        Dimensionless pixel-noise variance for each magnitude sample.

    Returns
    -------
    weights : ForestWeights
        Fixed source weights, integrals, and noise coefficients.
    """
    return prepare_forest_weights(
        FIELD,
        geometry(),
        InstrumentResponse(1, 0),
        z_source=3,
        magnitudes=[20, 21],
        quadrature=[1, 1],
        rho=rho,
        variance=variance,
        length_velocity=1,
        method="supplied",
        weights=weights,
    )


def oracle(rho, weights, variance):
    """Compute two-sample noise normalization with 500-digit arithmetic.

    Parameters
    ----------
    rho : array_like of shape (n_magnitudes,)
        Source density per deg^2, velocity interval in km/s, and magnitude.
    weights : array_like of shape (n_magnitudes,)
        Dimensionless source weights in magnitude order.
    variance : array_like of shape (n_magnitudes,)
        Dimensionless pixel-noise variance for each magnitude sample.

    Returns
    -------
    reference : tuple of ndarray
        Three cumulative density/weight integrals and the final alias/pixel
        coefficients.
    """
    with localcontext() as context:
        context.prec = 500
        decimal_density, decimal_weights, decimal_variance = (
            [Decimal(float(x)) for x in a] for a in (rho, weights, variance)
        )

        # Decimal products retain contributions below float64 intermediate range.
        prefixes = [
            [
                sum(
                    decimal_density[j]
                    * decimal_weights[j] ** power
                    * (decimal_variance[j] if noise else 1)
                    for j in range(i + 1)
                )
                for i in range(2)
            ]
            for power, noise in [(1, False), (2, False), (2, True)]
        ]
        alias_coefficient = prefixes[1][-1] / prefixes[0][-1] ** 2
        pixel_coefficient = prefixes[2][-1] / prefixes[0][-1] ** 2
        return np.array(prefixes, dtype=float), np.array(
            [alias_coefficient, pixel_coefficient], dtype=float
        )


@pytest.mark.parametrize("reverse", [False, True])
def test_r1_mixed_underflow_and_normalization(reverse):
    """Check r1 mixed underflow and normalization.

    Parameters
    ----------
    reverse : bool
        Whether to reverse input ordering, supplied by pytest parametrization.
    """
    accepted = []
    for scale in (1.0, 1e100):
        source_weights = np.array([1.0, 1e-200]) * scale
        pixel_variance = np.array([1e-150, 1e300])
        if reverse:
            source_weights, pixel_variance = source_weights[::-1], pixel_variance[::-1]
        prefixes, coefficients = oracle([1, 1], source_weights, pixel_variance)
        assert coefficients[1] == pytest.approx(1e-100, rel=5e-13, abs=0)
        if scale == 1:
            with pytest.raises(
                ValueError, match="forest.*not representable.*underflow"
            ):
                prepare([1, 1], source_weights, pixel_variance)
        else:
            result = prepare([1, 1], source_weights, pixel_variance)
            np.testing.assert_allclose(
                [result.I1, result.I2, result.I3], prefixes, rtol=5e-13, atol=0
            )
            np.testing.assert_allclose(
                [result.A, result.P_pixel], coefficients, rtol=5e-13, atol=0
            )
            accepted.append(result.P_pixel)
    # A second representable normalization must preserve the same coefficient.
    result2 = prepare([1, 1], source_weights * 2, pixel_variance)
    np.testing.assert_allclose(result2.P_pixel, accepted[0], rtol=5e-13, atol=0)


@pytest.mark.parametrize(
    "rho,w,v",
    [
        ([0, 1], [1e-200, 1], [1e300, 2]),
        ([1, 1], [0, 1], [1e300, 2]),
        ([1, 1], [1, 1], [0, 2]),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_exact_zero_contributions_with_positive_total(rho, w, v, reverse):
    """Check exact zero contributions with positive total.

    Parameters
    ----------
    rho : list
        Source density per velocity and magnitude, supplied by pytest
        parametrization.
    w : list
        Source-weight test input, supplied by pytest parametrization.
    v : list
        Parametrized numeric input, including invalid values where specified,
        supplied by pytest parametrization.
    reverse : bool
        Whether to reverse input ordering, supplied by pytest parametrization.
    """
    if reverse:
        rho, w, v = rho[::-1], w[::-1], v[::-1]
    prefixes, coefficients = oracle(rho, w, v)
    result = prepare(rho, w, v)
    np.testing.assert_allclose(
        [result.I1, result.I2, result.I3], prefixes, rtol=5e-13, atol=0
    )
    np.testing.assert_allclose(
        [result.A, result.P_pixel], coefficients, rtol=5e-13, atol=0
    )


@pytest.mark.parametrize(
    "rho,w,v",
    [
        ([1, 1e-200], [1, 1e-200], [1, 1e300]),  # r*w boundary
        ([1, 1], [1, 1e-150], [1, 1e-100]),  # r*w*w*variance boundary
    ],
)
def test_individual_loss_at_other_product_boundaries(rho, w, v):
    """Check individual loss at other product boundaries.

    Parameters
    ----------
    rho : list
        Source density per velocity and magnitude, supplied by pytest
        parametrization.
    w : list
        Source-weight test input, supplied by pytest parametrization.
    v : list
        Parametrized numeric input, including invalid values where specified,
        supplied by pytest parametrization.
    """
    with pytest.raises(ValueError, match="not representable.*underflow"):
        prepare(rho, w, v)
