"""Full-sample forest recurrences and confirmed stopping, array-only.

Arithmetic adapted from lyaforecast weights.py (GPLv3). Signed compatibility
preparation is separate from strict public preparation. No normalization occurs
inside the nonlinear update.
"""

import numpy as np

VARIANTS = ("sum_historical", "sum_aliasing")
METHODS = {"early_lyaforecast": "sum_historical", "mcdonald": "sum_aliasing"}


def seed(inputs):
    """Construct the common legacy forest-weight seed.

    Parameters
    ----------
    inputs : WeightInputs
        Prepared one-dimensional magnitude arrays: density in deg^-2 mag^-1
        (km/s)^-1, quadrature in mag, and dimensionless pixel variance; length
        and pixel in km/s, signal in deg^2 km/s, and p1d in km/s.

    Returns
    -------
    weights : ndarray of shape (n_magnitude,)
        Dimensionless P1D/pixel divided by P1D/pixel plus pixel variance.

    Raises
    ------
    FloatingPointError
        If the controlled division, invalid, or overflow arithmetic fails, or
        checked output is nonfinite.

    Notes
    -----
    Keep the saved division order; do not renormalize the seed.
    """
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        power = inputs.p1d / inputs.pixel
        return power / (power + inputs.variance)


def moments(inputs, weights):
    """Integrate the cumulative density-weight and noise moments.

    Parameters
    ----------
    inputs : WeightInputs
        Prepared one-dimensional magnitude arrays: density in deg^-2 mag^-1
        (km/s)^-1, quadrature in mag, and dimensionless pixel variance; length
        and pixel in km/s, signal in deg^2 km/s, and p1d in km/s.
    weights : ndarray of shape (n_magnitude,)
        Dimensionless source weights; signed arithmetic is retained.

    Returns
    -------
    first_moment, second_moment, noise_moment : ndarray of shape (n_magnitude,)
        Cumulative I1, I2, and I3 in deg^-2 (km/s)^-1.

    Notes
    -----
    The multiplication and cumulative reduction order are the literal legacy
    order. I2 and I3 use squared weights; I3 also includes pixel variance.
    """
    density, magnitude_weights, pixel_variance = (
        inputs.density,
        inputs.quadrature,
        inputs.variance,
    )
    return (
        np.cumsum(density * weights * magnitude_weights),
        np.cumsum(density * weights**2 * magnitude_weights),
        np.cumsum(density * weights**2 * pixel_variance * magnitude_weights),
    )


def update(inputs, weights, variant):
    """Apply one simultaneous full-sample forest-weight update.

    Parameters
    ----------
    inputs : WeightInputs
        Prepared one-dimensional magnitude arrays: density in deg^-2 mag^-1
        (km/s)^-1, quadrature in mag, and dimensionless pixel variance; length
        and pixel in km/s, signal in deg^2 km/s, and p1d in km/s.
    weights : ndarray of shape (n_magnitude,)
        Dimensionless source weights; signed arithmetic is retained.
    variant : {'sum_historical', 'sum_aliasing'}
        Full-sample recurrence prescription.

    Returns
    -------
    weights : ndarray of shape (n_magnitude,)
        Updated dimensionless weights, retaining their absolute amplitude.

    Raises
    ------
    ValueError
        If the recurrence variant is unknown.
    FloatingPointError
        If the controlled division, invalid, or overflow arithmetic fails, or
        checked output is nonfinite.

    Notes
    -----
    Full-sample sums use the same cumulative reduction as prefix controls.
    Floating-point failures propagate so the caller can retain the last finite
    state; there is no positivity replacement or weight normalization.
    """
    if variant not in VARIANTS:
        raise ValueError(f"unknown compatibility weighting variant: {variant}")
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        first_moment, second_moment, _ = moments(inputs, weights)
        if variant.startswith("sum_"):
            first_moment, second_moment = first_moment[-1], second_moment[-1]
        signal = inputs.signal
        if variant.endswith("aliasing"):
            signal = signal + inputs.p1d * second_moment / (
                first_moment**2 * inputs.length
            )
        elif variant == "sum_historical":
            signal = signal + inputs.p1d / (first_moment * inputs.length)
        noise = inputs.variance / (first_moment * (inputs.length / inputs.pixel))
        result = signal / (signal + noise)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite updated weights")
    return result


def coefficients(inputs, weights):
    """Normalize the full-sample moments into forest-noise coefficients.

    Parameters
    ----------
    inputs : WeightInputs
        Prepared one-dimensional magnitude arrays: density in deg^-2 mag^-1
        (km/s)^-1, quadrature in mag, and dimensionless pixel variance; length
        and pixel in km/s, signal in deg^2 km/s, and p1d in km/s.
    weights : ndarray of shape (n_magnitude,)
        Dimensionless source weights; signed arithmetic is retained.

    Returns
    -------
    coefficients : ndarray of shape (5,)
        I1, I2, I3 in deg^-2 (km/s)^-1, A in deg^2, and pixel power in deg^2
        km/s, in that order.

    Raises
    ------
    FloatingPointError
        If the controlled division, invalid, or overflow arithmetic fails, or
        checked output is nonfinite.

    Notes
    -----
    No positivity substitution is applied to signed compatibility inputs.
    """
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        first_moment, second_moment, noise_moment = (
            a[-1] for a in moments(inputs, weights)
        )
        denominator = first_moment**2 * inputs.length
        result = np.array(
            [
                first_moment,
                second_moment,
                noise_moment,
                second_moment / denominator,
                noise_moment * inputs.pixel / denominator,
            ]
        )
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite final moments or coefficients")
    return result


def fixed_weights(inputs, variant="prefix_intrinsic", updates=3):
    """Return weights after a prescribed number of updates from the seed.

    Parameters
    ----------
    inputs : WeightInputs
        Prepared one-dimensional magnitude arrays: density in deg^-2 mag^-1
        (km/s)^-1, quadrature in mag, and dimensionless pixel variance; length
        and pixel in km/s, signal in deg^2 km/s, and p1d in km/s.
    variant : str, default='prefix_intrinsic'
        Recurrence name. The retained historical default is not in VARIANTS and
        raises ValueError; callers must select a supported full-sample variant.
    updates : int, default=3
        Nonnegative number of simultaneous updates after the seed.

    Returns
    -------
    weights : ndarray of shape (n_magnitude,)
        Dimensionless weights after exactly the requested number of updates.

    Raises
    ------
    ValueError
        If the update count or variant is invalid.
    FloatingPointError
        If the controlled division, invalid, or overflow arithmetic fails, or
        checked output is nonfinite.
    """
    if isinstance(updates, bool) or not isinstance(updates, int) or updates < 0:
        raise ValueError("updates must be a nonnegative integer")
    if variant not in VARIANTS:
        raise ValueError("unknown compatibility weighting variant")
    weights = seed(inputs)
    for _ in range(updates):
        weights = update(inputs, weights, variant)
    return weights


def relative_change(new, old):
    """Measure an infinity-norm relative change without a denominator floor.

    Parameters
    ----------
    new, old : ndarray
        New and previous values with identical shape and units.

    Returns
    -------
    change : float
        Dimensionless max(abs(new-old))/max(abs(old)); zero for two zero arrays
        and infinity for a nonzero new array with zero old amplitude.
    """
    scale = np.max(np.abs(old))
    if scale == 0:
        return 0.0 if np.max(np.abs(new)) == 0 else np.inf
    return float(np.max(np.abs(new - old)) / scale)


def changes(new, old, new_c, old_c):
    """Compare weight amplitude, signed shape, and both noise coefficients.

    Parameters
    ----------
    new, old : ndarray
        New and previous values with identical shape and units.
    new_c, old_c : ndarray of shape (5,)
        New and previous integrated I1, I2, I3, aliasing coefficient A, and
        pixel power; comparisons are dimensionless.

    Returns
    -------
    changes : ndarray of shape (4,)
        Dimensionless relative changes in amplitude, normalized shape, A, and
        pixel power. All entries are infinite if either weight amplitude is
        zero.
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


METRICS = ("amplitude", "shape", "A", "P_pixel", "vector")


def residuals(new, old, new_c, old_c):
    """Append the full weight-vector change to the four convergence metrics.

    Parameters
    ----------
    new, old : ndarray
        New and previous values with identical shape and units.
    new_c, old_c : ndarray of shape (5,)
        New and previous integrated I1, I2, I3, aliasing coefficient A, and
        pixel power; comparisons are dimensionless.

    Returns
    -------
    residuals : ndarray of shape (5,)
        Dimensionless amplitude, shape, A, pixel-power, and vector changes.

    Notes
    -----
    These compare the supplied states; they are not forward residuals and
    do not evaluate another recurrence update.
    """
    return np.r_[changes(new, old, new_c, old_c), relative_change(new, old)]


def solve(
    inputs,
    variant,
    *,
    eligible=True,
    rtol=1e-4,
    min_updates=3,
    stable_steps=3,
    max_updates=96,
):
    """Iterate until full-sample forest weights satisfy confirmed convergence.

    Parameters
    ----------
    inputs : WeightInputs
        Prepared one-dimensional magnitude arrays: density in deg^-2 mag^-1
        (km/s)^-1, quadrature in mag, and dimensionless pixel variance; length
        and pixel in km/s, signal in deg^2 km/s, and p1d in km/s.
    variant : {'sum_historical', 'sum_aliasing'}
        Full-sample recurrence prescription.
    eligible : bool, default=True
        Whether the physical spectrum permits this stopping prescription.
    rtol : float, default=1e-4
        Positive dimensionless tolerance for all five relative-change metrics.
    min_updates : int, default=3
        Earliest update count eligible to nominate a candidate.
    stable_steps : int, default=3
        Required number of consecutive stable transitions.
    max_updates : int, default=96
        Maximum number of completed updates.

    Returns
    -------
    result : dict
        Status, reason, last finite weights and coefficients, update counts,
        candidate count, and preceding/confirmation metrics. Unavailable values
        remain None.

    Raises
    ------
    ValueError
        If the variant, integer controls, or relative tolerance is invalid.

    Notes
    -----
    Only the named forest autos with positive intrinsic signal and the two
    full-sum variants are eligible. Eligibility does not imply convergence.
    Three stable transitions nominate a candidate by default. Every subsequent
    transition and every candidate-to-current comparison must remain within rtol
    through twice the candidate count. Return that actually computed state.

    ``weights`` and ``coefficients`` are the last finite state on cap/failure,
    never a converged substitute. ``updates`` counts completed transitions;
    ``state_updates`` identifies the returned state if its coefficients failed.
    The forward residual at that state is unavailable: no hidden extra update
    is evaluated. ``last_step`` measures its preceding transition instead.
    """
    if variant not in VARIANTS:
        raise ValueError("unknown compatibility weighting variant")
    for value in (min_updates, stable_steps, max_updates):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("update controls must be positive integers")
    if not np.isfinite(rtol) or rtol <= 0:
        raise ValueError("rtol must be finite and positive")
    result = dict(
        status="ineligible",
        reason=None,
        weights=None,
        coefficients=None,
        updates=0,
        state_updates=None,
        candidate=None,
        last_step=None,
        confirmation=None,
        forward_residual=None,
    )
    if not eligible or variant not in ("sum_historical", "sum_aliasing"):
        result["reason"] = "adaptive stopping requires a full-sum forest-auto context"
        return result
    if (
        not np.isfinite(inputs.signal)
        or inputs.signal <= 0
        or not np.isfinite(inputs.p1d)
        or inputs.p1d <= 0
    ):
        result["reason"] = (
            "adaptive stopping requires finite positive intrinsic signal and P1D reference"
        )
        return result
    candidates = []
    stable = 0
    try:
        weights = seed(inputs)
        noise_coefficients = coefficients(inputs, weights)
        result.update(weights=weights, coefficients=noise_coefficients, state_updates=0)
        for t in range(1, max_updates + 1):
            new = update(inputs, weights, variant)
            result["updates"] = t
            new_coefficients = coefficients(inputs, new)
            step = residuals(new, weights, new_coefficients, noise_coefficients)
            result.update(
                weights=new,
                coefficients=new_coefficients,
                state_updates=t,
                last_step=dict(zip(METRICS, step.tolist(), strict=True)),
            )
            stable = stable + 1 if np.all(step <= rtol) else 0
            candidates = [
                (n, candidate_weights, candidate_coefficients)
                for n, candidate_weights, candidate_coefficients in candidates
                if stable
                and np.all(
                    residuals(
                        new, candidate_weights, new_coefficients, candidate_coefficients
                    )
                    <= rtol
                )
            ]
            for n, candidate_weights, candidate_coefficients in candidates:
                if t == 2 * n:
                    result.update(
                        status="converged",
                        candidate=n,
                        confirmation=dict(
                            zip(
                                METRICS,
                                residuals(
                                    new,
                                    candidate_weights,
                                    new_coefficients,
                                    candidate_coefficients,
                                ).tolist(),
                                strict=True,
                            )
                        ),
                    )
                    return result
            if t >= min_updates and stable >= stable_steps and 2 * t <= max_updates:
                candidates.append((t, new.copy(), new_coefficients.copy()))
            weights, noise_coefficients = new, new_coefficients
    except FloatingPointError as error:
        result.update(status="arithmetic_failure", reason=str(error))
        return result
    result.update(
        status="capped", reason="no confirmed finite nonzero convergence within cap"
    )
    return result
