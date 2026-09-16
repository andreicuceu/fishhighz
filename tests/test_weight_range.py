"""Decimal oracles for mixed contributions concealed by positive totals (R1)."""

from decimal import Decimal, localcontext

import numpy as np
import pytest
from test_weights import FIELD, geometry

from fishhighz.response import InstrumentResponse
from fishhighz.weights import prepare_forest_weights


def prepare(rho, weights, variance):
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
    with localcontext() as context:
        context.prec = 500
        r, w, v = ([Decimal(float(x)) for x in a] for a in (rho, weights, variance))
        prefixes = [
            [
                sum(r[j] * w[j] ** power * (v[j] if noise else 1) for j in range(i + 1))
                for i in range(2)
            ]
            for power, noise in [(1, False), (2, False), (2, True)]
        ]
        a = prefixes[1][-1] / prefixes[0][-1] ** 2
        p = prefixes[2][-1] / prefixes[0][-1] ** 2
        return np.array(prefixes, dtype=float), np.array([a, p], dtype=float)


@pytest.mark.parametrize("reverse", [False, True])
def test_r1_mixed_underflow_and_normalization(reverse):
    accepted = []
    for scale in (1.0, 1e100):
        w = np.array([1.0, 1e-200]) * scale
        v = np.array([1e-150, 1e300])
        if reverse:
            w, v = w[::-1], v[::-1]
        prefixes, coefficients = oracle([1, 1], w, v)
        assert coefficients[1] == pytest.approx(1e-100, rel=5e-13, abs=0)
        if scale == 1:
            with pytest.raises(
                ValueError, match="forest.*not representable.*underflow"
            ):
                prepare([1, 1], w, v)
        else:
            result = prepare([1, 1], w, v)
            np.testing.assert_allclose(
                [result.I1, result.I2, result.I3], prefixes, rtol=5e-13, atol=0
            )
            np.testing.assert_allclose(
                [result.A, result.P_pixel], coefficients, rtol=5e-13, atol=0
            )
            accepted.append(result.P_pixel)
    # A second representable normalization must preserve the same coefficient.
    result2 = prepare([1, 1], w * 2, v)
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
    with pytest.raises(ValueError, match="not representable.*underflow"):
        prepare(rho, w, v)
