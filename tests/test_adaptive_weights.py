"""Focused controls for conditional stopping and unchanged fixed counts."""

from dataclasses import replace

import numpy as np
import pytest

from fishhighz.validation.adaptive_weights import AUTO_CONTEXTS, adaptive_weights
from fishhighz.validation.compatibility_weights import (
    VARIANTS,
    WeightInputs,
    fixed_weights,
    seed,
    update,
)


def homogeneous(signal=2.0):
    return WeightInputs(*(np.ones(1) for _ in range(4)), 1.0, 1.0, signal, 1.0)


def test_exact_fixed_point_stops_after_actual_double_count():
    result = adaptive_weights(homogeneous(), "sum_intrinsic", context=AUTO_CONTEXTS[0])
    assert result["status"] == "converged"
    assert result["candidate"] == 3
    assert result["updates"] == result["state_updates"] == 6
    np.testing.assert_array_equal(result["weights"], [0.5])
    assert set(result["last_step"].values()) == {0.0}
    assert result["forward_residual"] is None


def test_decay_and_late_candidate_report_cap():
    for inputs, cap in ((homogeneous(0.99), 96), (homogeneous(), 5)):
        result = adaptive_weights(
            inputs, "sum_intrinsic", context=AUTO_CONTEXTS[0], max_updates=cap
        )
        assert result["status"] == "capped"
        assert result["updates"] == cap
        np.testing.assert_array_equal(
            result["weights"], fixed_weights(inputs, "sum_intrinsic", cap)
        )


@pytest.mark.parametrize("variant", VARIANTS)
def test_exactly_three_updates_unchanged(variant):
    inputs = homogeneous()
    expected = seed(inputs)
    for _ in range(3):
        expected = update(inputs, expected, variant)
    np.testing.assert_array_equal(fixed_weights(inputs, variant), expected)


def test_ineligible_prefix_cross_and_nonpositive_signal():
    for variant, context, signal in (
        ("prefix_intrinsic", AUTO_CONTEXTS[0], 2),
        ("prefix_aliasing", AUTO_CONTEXTS[0], 2),
        ("sum_historical", "lya(qso)_lya(lbg)", 2),
        ("sum_intrinsic", AUTO_CONTEXTS[0], 0),
        ("sum_intrinsic", AUTO_CONTEXTS[0], -1),
    ):
        result = adaptive_weights(homogeneous(signal), variant, context=context)
        assert result["status"] == "ineligible"
        assert result["updates"] == 0 and result["weights"] is None


def test_arithmetic_failure_retains_last_finite_state():
    # Signed density gives J1=-1/2: P + noise = 2 - 2 = 0.
    inputs = replace(homogeneous(), density=np.array([-1.0]))
    result = adaptive_weights(inputs, "sum_intrinsic", context=AUTO_CONTEXTS[0])
    assert result["status"] == "arithmetic_failure"
    assert result["updates"] == result["state_updates"] == 0
    np.testing.assert_array_equal(result["weights"], seed(inputs))
    assert result["reason"]
