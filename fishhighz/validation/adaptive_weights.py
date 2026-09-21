"""Opt-in finite-grid stopping for compatibility full-sum forest autos.

Empirically checked on the three stage-2 grids, not a convergence theorem for
arbitrary signed measures. Fixed-count preparation remains a separate API.
"""

import numpy as np

from .compatibility_weights import VARIANTS, coefficients, seed, update
from .weight_convergence import changes, relative_change

AUTO_CONTEXTS = ("lya(qso)_lya(qso)", "lya(lbg)_lya(lbg)")
ADAPTIVE_VARIANTS = ("sum_intrinsic", "sum_aliasing", "sum_historical")
METRICS = ("amplitude", "shape", "A", "P_pixel", "vector")


def residuals(new, old, new_c, old_c):
    """Relative changes from old to new; vector change is not a forward residual."""
    return np.r_[changes(new, old, new_c, old_c), relative_change(new, old)]


def adaptive_weights(
    inputs,
    variant,
    *,
    context,
    rtol=1e-4,
    min_updates=3,
    stable_steps=3,
    max_updates=96,
):
    """Return an explicit status, last finite state, count and relative metrics.

    Only the named forest autos with positive intrinsic signal and the three
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
    if variant in ("sum_historical", "sum_aliasing"):
        from ..kernels.full_sum_weights import solve

        return solve(
            inputs,
            variant,
            eligible=context in AUTO_CONTEXTS,
            rtol=rtol,
            min_updates=min_updates,
            stable_steps=stable_steps,
            max_updates=max_updates,
        )
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
    if context not in AUTO_CONTEXTS or variant not in ADAPTIVE_VARIANTS:
        result["reason"] = "adaptive stopping requires a full-sum forest-auto context"
        return result
    if not np.isfinite(inputs.signal) or inputs.signal <= 0:
        result["reason"] = "adaptive stopping requires finite positive intrinsic signal"
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
