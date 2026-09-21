"""Analytic amplitude decay and transient plateau controls."""

import numpy as np

from fishhighz.validation.compatibility_weights import WeightInputs
from fishhighz.validation.weight_convergence import (
    classify,
    relative_change,
    trajectory,
)


def test_relative_metric_does_not_accept_tiny_weights():
    assert relative_change(np.array([1e-200]), np.array([2e-200])) == 0.5
    assert relative_change(np.zeros(1), np.zeros(1)) == 0
    assert np.isinf(relative_change(np.ones(1), np.zeros(1)))


def test_subthreshold_homogeneous_decay_has_constant_coefficients():
    inputs = WeightInputs(
        np.array([0.0]),
        np.array([1.0]),
        np.array([1.0]),
        np.array([1.0]),
        1.0,
        1.0,
        0.99,
        0.01,
    )
    data = trajectory(inputs, "sum_intrinsic")
    np.testing.assert_allclose(data["coefficients"][:, -2:], 1.0, rtol=1e-14)
    assert classify(data, 1e-3)["status"] == "amplitude_decay_stable_normalized"
    assert data["amplitude"][-1] < data["amplitude"][0]


def test_transient_plateau_is_not_convergence():
    # A false fixed plateau until t=20 then order-unity departure.
    w = np.ones((97, 1))
    w[21:, 0] = np.linspace(1.0, 2.0, 76)
    c = np.ones((97, 5))
    from fishhighz.validation.weight_convergence import changes

    data = dict(
        weights=w,
        coefficients=c,
        amplitude=w[:, 0],
        failure=None,
        changes=np.array(
            [changes(w[t], w[t - 1], c[t], c[t - 1]) for t in range(1, 97)]
        ),
        residual=np.r_[np.abs(np.diff(w[:, 0])) / w[:-1, 0], np.nan],
    )
    assert classify(data, 1e-4)["status"] == "transient_plateau_then_drift"
