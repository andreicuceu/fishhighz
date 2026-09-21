"""Full-sample forest recurrences and confirmed stopping, array-only.

Arithmetic adapted from lyaforecast weights.py (GPLv3). Signed compatibility
preparation is separate from strict public preparation. No normalization occurs
inside the nonlinear update.
"""

import numpy as np

VARIANTS = ("sum_historical", "sum_aliasing")
METHODS = {"early_lyaforecast": "sum_historical", "mcdonald": "sum_aliasing"}


def seed(inputs):
    """Common legacy seed, preserving the saved division order."""
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        power = inputs.p1d / inputs.pixel
        return power / (power + inputs.variance)


def moments(inputs, weights):
    """Cumulative I1/I2/I3 with literal legacy multiplication/reduction order."""
    d, q, v = inputs.density, inputs.quadrature, inputs.variance
    return (
        np.cumsum(d * weights * q),
        np.cumsum(d * weights**2 * q),
        np.cumsum(d * weights**2 * v * q),
    )


def update(inputs, weights, variant):
    """One simultaneous update, retaining amplitude and signed arithmetic.

    Floating-point failures propagate; callers must report the last finite state.
    Full-sample sums use the same cumulative reduction as the prefix controls.
    """
    if variant not in VARIANTS:
        raise ValueError(f"unknown compatibility weighting variant: {variant}")
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        j1, j2, _ = moments(inputs, weights)
        if variant.startswith("sum_"):
            j1, j2 = j1[-1], j2[-1]
        signal = inputs.signal
        if variant.endswith("aliasing"):
            signal = signal + inputs.p1d * j2 / (j1**2 * inputs.length)
        elif variant == "sum_historical":
            signal = signal + inputs.p1d / (j1 * inputs.length)
        noise = inputs.variance / (j1 * (inputs.length / inputs.pixel))
        result = signal / (signal + noise)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite updated weights")
    return result


def coefficients(inputs, weights):
    """Full-sample I1, I2, I3, A and pixel power; no positivity substitution."""
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        i1, i2, i3 = (a[-1] for a in moments(inputs, weights))
        denominator = i1**2 * inputs.length
        result = np.array(
            [i1, i2, i3, i2 / denominator, i3 * inputs.pixel / denominator]
        )
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite final moments or coefficients")
    return result


def fixed_weights(inputs, variant="prefix_intrinsic", updates=3):
    """Return weights after exactly updates transitions following the seed."""
    if isinstance(updates, bool) or not isinstance(updates, int) or updates < 0:
        raise ValueError("updates must be a nonnegative integer")
    if variant not in VARIANTS:
        raise ValueError("unknown compatibility weighting variant")
    weights = seed(inputs)
    for _ in range(updates):
        weights = update(inputs, weights, variant)
    return weights


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


METRICS = ("amplitude", "shape", "A", "P_pixel", "vector")


def residuals(new, old, new_c, old_c):
    """Relative changes from old to new; vector change is not a forward residual."""
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
    """Return an explicit status, last finite state, count and relative metrics.

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
        w = seed(inputs)
        c = coefficients(inputs, w)
        result.update(weights=w, coefficients=c, state_updates=0)
        for t in range(1, max_updates + 1):
            new = update(inputs, w, variant)
            result["updates"] = t
            new_c = coefficients(inputs, new)
            step = residuals(new, w, new_c, c)
            result.update(
                weights=new,
                coefficients=new_c,
                state_updates=t,
                last_step=dict(zip(METRICS, step.tolist(), strict=True)),
            )
            stable = stable + 1 if np.all(step <= rtol) else 0
            candidates = [
                (n, cw, cc)
                for n, cw, cc in candidates
                if stable and np.all(residuals(new, cw, new_c, cc) <= rtol)
            ]
            for n, cw, cc in candidates:
                if t == 2 * n:
                    result.update(
                        status="converged",
                        candidate=n,
                        confirmation=dict(
                            zip(
                                METRICS,
                                residuals(new, cw, new_c, cc).tolist(),
                                strict=True,
                            )
                        ),
                    )
                    return result
            if t >= min_updates and stable >= stable_steps and 2 * t <= max_updates:
                candidates.append((t, new.copy(), new_c.copy()))
            w, c = new, new_c
    except FloatingPointError as error:
        result.update(status="arithmetic_failure", reason=str(error))
        return result
    result.update(
        status="capped", reason="no confirmed finite nonzero convergence within cap"
    )
    return result
