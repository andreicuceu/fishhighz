"""Finite-count diagnostics for signed compatibility recurrences."""

import numpy as np

from .compatibility_weights import coefficients, seed, update


def relative_change(new, old):
    """Infinity-norm relative change without an absolute denominator floor.

    Parameters
    ----------
    new : array_like, shape (n_magnitude,)
        Updated dimensionless source weights.
    old : array_like, shape (n_magnitude,)
        Reference dimensionless source weights.

    Returns
    -------
    relative_difference : float
        Dimensionless infinity-norm change; infinity if only the reference is
        exactly zero.
    """
    scale = np.max(np.abs(old))
    if scale == 0:
        return 0.0 if np.max(np.abs(new)) == 0 else np.inf
    return float(np.max(np.abs(new - old)) / scale)


def changes(new, old, new_c, old_c):
    """Amplitude, signed normalized shape, A, and pixel-power changes.

    Parameters
    ----------
    new : array_like, shape (n_magnitude,)
        Updated dimensionless source weights.
    old : array_like, shape (n_magnitude,)
        Reference dimensionless source weights.
    new_c : array_like, shape (5,)
        Updated I1, I2, I3, aliasing coefficient A and pixel-noise power.
    old_c : array_like, shape (5,)
        Reference I1, I2, I3, aliasing coefficient A and pixel-noise power.

    Returns
    -------
    changes : ndarray, shape (4,)
        Relative amplitude, normalized shape, A and pixel-power changes.
    """
    new_amplitude, old_amplitude = np.max(np.abs(new)), np.max(np.abs(old))
    if new_amplitude == 0 or old_amplitude == 0:
        return np.full(4, np.inf)
    return np.array(
        [
            abs(new_amplitude / old_amplitude - 1),
            relative_change(new / new_amplitude, old / old_amplitude),
            relative_change(new_c[-2:-1], old_c[-2:-1]),
            relative_change(new_c[-1:], old_c[-1:]),
        ]
    )


def trajectory(inputs, variant, cap=96):
    """Save all finite states; stop arithmetic failure without substitution.

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.
    variant : str
        Named compatibility recurrence, selecting prefix or full-sample moments
        and the auxiliary signal prescription.
    cap : int
        Maximum number of transitions following the initial state. Default is
        ``96``.

    Returns
    -------
    trajectory : dict
        Finite weight states with shape (n_state, n_magnitude), five noise
        coefficients per state, four changes per transition, residuals,
        amplitudes and any failure.
    """
    weights, coefficient_history, failure = [], [], None
    weight_states = seed(inputs)
    for t in range(cap + 1):
        try:
            coefficient_states = coefficients(inputs, weight_states)
        except FloatingPointError as error:
            failure = dict(iteration=t, reason=str(error), operation="coefficients")
            break
        weights.append(weight_states.copy())
        coefficient_history.append(coefficient_states)
        if t == cap:
            break
        try:
            weight_states = update(inputs, weight_states, variant)
        except FloatingPointError as error:
            failure = dict(iteration=t + 1, reason=str(error), operation="update")
            break

    weight_states, coefficient_states = (
        np.asarray(weights),
        np.asarray(coefficient_history),
    )
    delta = np.array(
        [
            changes(
                weight_states[t],
                weight_states[t - 1],
                coefficient_states[t],
                coefficient_states[t - 1],
            )
            for t in range(1, len(weight_states))
        ]
    )
    residual = np.full(len(weight_states), np.nan)
    for t in range(len(weight_states) - 1):
        residual[t] = relative_change(weight_states[t + 1], weight_states[t])
    # The final forward residual is unavailable: evaluating it would require
    # an update beyond the cap (or retrying the failed transition).
    return dict(
        weights=weight_states,
        coefficients=coefficient_states,
        changes=delta,
        residual=residual,
        amplitude=np.max(np.abs(weight_states), axis=1),
        failure=failure,
    )


def classify(data, tolerance):
    """Three stable updates, then candidate-to-double-count confirmation.

    Parameters
    ----------
    data : dict
        Saved finite trajectory returned by trajectory.
    tolerance : float
        Dimensionless threshold for all relative convergence metrics.

    Returns
    -------
    classification : dict
        Convergence status and the first supported candidate, doubling count and
        confirmation changes.
    """
    weight_states, coefficient_states, step_changes = (
        data["weights"],
        data["coefficients"],
        data["changes"],
    )
    candidates = []
    confirmed = []
    for t in range(3, len(weight_states)):
        if np.all(step_changes[t - 3 : t] <= tolerance) and np.all(
            data["residual"][t - 3 : t] <= tolerance
        ):
            candidates.append(t)
            if 2 * t < len(weight_states):
                difference = changes(
                    weight_states[2 * t],
                    weight_states[t],
                    coefficient_states[2 * t],
                    coefficient_states[t],
                )
                if (
                    np.all(difference <= tolerance)
                    and relative_change(weight_states[2 * t], weight_states[t])
                    <= tolerance
                ):
                    confirmed.append((t, 2 * t, difference))

    # A transient plateau can pass doubling and subsequently leave it.
    # The full saved trajectory must continue to support the candidate.
    doubled_candidates = confirmed.copy()
    confirmed = [
        (t, double, diff)
        for t, double, diff in confirmed
        if all(
            np.all(
                changes(
                    weight_states[j],
                    weight_states[t],
                    coefficient_states[j],
                    coefficient_states[t],
                )
                <= tolerance
            )
            and relative_change(weight_states[j], weight_states[t]) <= tolerance
            for j in range(t, len(weight_states))
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
        len(weight_states) >= 25
        and np.all(np.diff(data["amplitude"][-25:]) < 0)
        and np.max(step_changes[-1, 1:]) <= tolerance
    ):
        status = "amplitude_decay_stable_normalized"
    elif len(weight_states) >= 25 and np.all(np.diff(data["amplitude"][-25:]) < 0):
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
