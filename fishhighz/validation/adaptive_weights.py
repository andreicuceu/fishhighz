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
    """Relative changes from old to new; vector change is not a forward residual.

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
    residuals : ndarray, shape (5,)
        Relative amplitude, shape, A, pixel-power and full-vector changes.
    """
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

    Parameters
    ----------
    inputs : WeightInputs
        Magnitude-dependent source density, integration measure, pixel variance
        and auxiliary signal in angular/velocity units.
    variant : str
        Named compatibility recurrence, selecting prefix or full-sample moments
        and the auxiliary signal prescription.
    context : str
        Forest-auto pair identifier used to determine adaptive-stopping
        eligibility.
    rtol : float
        Dimensionless relative numerical tolerance. Default is ``0.0001``.
    min_updates : int
        Minimum update count eligible to nominate a candidate. Default is ``3``.
    stable_steps : int
        Required consecutive stable transitions. Default is ``3``.
    max_updates : int
        Maximum number of weight transitions. Default is ``96``.

    Returns
    -------
    result : dict
        Eligibility/convergence status, last finite weights and coefficients,
        completed update count and five relative convergence metrics.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.

    Notes
    -----
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
        weights = seed(inputs)
        noise_coefficients = coefficients(inputs, weights)
        result.update(weights=weights, coefficients=noise_coefficients, state_updates=0)
        for t in range(1, max_updates + 1):
            new = update(inputs, weights, variant)
            result["updates"] = t
            new_c = coefficients(inputs, new)
            step = residuals(new, weights, new_c, noise_coefficients)
            result.update(
                weights=new,
                coefficients=new_c,
                state_updates=t,
                last_step=dict(zip(METRICS, step.tolist(), strict=True)),
            )

            stable = stable + 1 if np.all(step <= rtol) else 0
            candidates = [
                (n, candidate_weights, candidate_coefficients)
                for n, candidate_weights, candidate_coefficients in candidates
                if stable
                and np.all(
                    residuals(new, candidate_weights, new_c, candidate_coefficients)
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
                                    new_c,
                                    candidate_coefficients,
                                ).tolist(),
                                strict=True,
                            )
                        ),
                    )
                    return result
            if t >= min_updates and stable >= stable_steps and 2 * t <= max_updates:
                candidates.append((t, new.copy(), new_c.copy()))
            weights, noise_coefficients = new, new_c
    except FloatingPointError as error:
        result.update(status="arithmetic_failure", reason=str(error))
        return result
    result.update(
        status="capped", reason="no confirmed finite nonzero convergence within cap"
    )
    return result
