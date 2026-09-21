"""Finite-count diagnostics for signed compatibility recurrences."""

import numpy as np

from .compatibility_weights import coefficients, seed, update


def relative_change(new, old):
    """Infinity-norm relative change without an absolute denominator floor."""
    scale = np.max(np.abs(old))
    if scale == 0:
        return 0.0 if np.max(np.abs(new)) == 0 else np.inf
    return float(np.max(np.abs(new - old)) / scale)


def changes(new, old, new_c, old_c):
    """Amplitude, signed normalized shape, A, and pixel-power changes."""
    a, b = np.max(np.abs(new)), np.max(np.abs(old))
    if a == 0 or b == 0:
        return np.full(4, np.inf)
    return np.array(
        [
            abs(a / b - 1),
            relative_change(new / a, old / b),
            relative_change(new_c[-2:-1], old_c[-2:-1]),
            relative_change(new_c[-1:], old_c[-1:]),
        ]
    )


def trajectory(inputs, variant, cap=96):
    """Save all finite states; stop arithmetic failure without substitution."""
    weights, coeffs, failure = [], [], None
    w = seed(inputs)
    for t in range(cap + 1):
        try:
            c = coefficients(inputs, w)
        except FloatingPointError as error:
            failure = dict(iteration=t, reason=str(error), operation="coefficients")
            break
        weights.append(w.copy())
        coeffs.append(c)
        if t == cap:
            break
        try:
            w = update(inputs, w, variant)
        except FloatingPointError as error:
            failure = dict(iteration=t + 1, reason=str(error), operation="update")
            break
    w, c = np.asarray(weights), np.asarray(coeffs)
    delta = np.array(
        [changes(w[t], w[t - 1], c[t], c[t - 1]) for t in range(1, len(w))]
    )
    residual = np.full(len(w), np.nan)
    for t in range(len(w) - 1):
        residual[t] = relative_change(w[t + 1], w[t])
    # The final forward residual is unavailable: evaluating it would require
    # an update beyond the cap (or retrying the failed transition).
    return dict(
        weights=w,
        coefficients=c,
        changes=delta,
        residual=residual,
        amplitude=np.max(np.abs(w), axis=1),
        failure=failure,
    )


def classify(data, tolerance):
    """Three stable updates, then candidate-to-double-count confirmation."""
    w, c, d = data["weights"], data["coefficients"], data["changes"]
    candidates = []
    confirmed = []
    for t in range(3, len(w)):
        if np.all(d[t - 3 : t] <= tolerance) and np.all(
            data["residual"][t - 3 : t] <= tolerance
        ):
            candidates.append(t)
            if 2 * t < len(w):
                difference = changes(w[2 * t], w[t], c[2 * t], c[t])
                if (
                    np.all(difference <= tolerance)
                    and relative_change(w[2 * t], w[t]) <= tolerance
                ):
                    confirmed.append((t, 2 * t, difference))
    # A transient plateau can pass doubling and subsequently leave it.
    # The full saved trajectory must continue to support the candidate.
    doubled_candidates = confirmed.copy()
    confirmed = [
        (t, double, diff)
        for t, double, diff in confirmed
        if all(
            np.all(changes(w[j], w[t], c[j], c[t]) <= tolerance)
            and relative_change(w[j], w[t]) <= tolerance
            for j in range(t, len(w))
        )
    ]
    if data["failure"] is not None:
        status = "invalid_arithmetic"
    elif np.any(data["amplitude"] == 0):
        status = "exact_zero"
    elif confirmed:
        status = "finite_nonzero_convergence"
    elif doubled_candidates:
        status = "transient_plateau_then_drift"
    elif candidates:
        status = "unconfirmed_candidate"
    elif (
        len(w) >= 25
        and np.all(np.diff(data["amplitude"][-25:]) < 0)
        and np.max(d[-1, 1:]) <= tolerance
    ):
        status = "amplitude_decay_stable_normalized"
    elif len(w) >= 25 and np.all(np.diff(data["amplitude"][-25:]) < 0):
        status = "continued_decay"
    else:
        status = "continued_drift_or_oscillation"
    item = dict(
        status=status,
        first_candidate=candidates[0] if candidates else None,
        candidate=None,
        confirmed_at=None,
        confirmation_changes=None,
        first_doubled_candidate=doubled_candidates[0][0]
        if doubled_candidates
        else None,
    )
    if status == "finite_nonzero_convergence":
        t, double, diff = confirmed[0]
        item.update(
            candidate=t, confirmed_at=double, confirmation_changes=diff.tolist()
        )
    return item
