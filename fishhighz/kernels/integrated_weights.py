"""Full-sum 'sum_historical' recurrence for integrated forest sources, array-only.

The recurrence, stopping metrics and confirmation logic reproduce
``kernels.full_sum_weights.solve`` for ``variant='sum_historical'`` with
``length = 1`` and unit quadrature (the flattened-measure mapping used by
``weights.prepare_integrated_forest_weights``), which stays the reference
implementation: same seed, update, stopping metrics and order of evaluation,
doubled-count confirmation, cap and convergence record. Stopping decisions are
identical except when ``rtol`` coincides with a metric to within its rounding
(~1e-16 absolute), where a different state still within ``rtol`` may be
accepted. The differences are

* every state of the recurrence is the one-parameter family
  w = c / (c + v) of the pixel variances v, with c = S' N1 / pixel and
  S' = S + B / N1 (the reference update S' / (S' + v pixel / N1) is
  algebraically the same, to a few units in the last place per element), so
  states, candidates and the previous state are carried as the scalar c and
  evaluated on the fly in streaming passes; no weight array is stored until
  the final state is written out, and the extreme weight is c / (min(v) + c);
* the full sums N1 = sum(mu w), N2 = sum(mu w^2) and N3 = sum(mu w^2 v) are
  single-pass block-wise (pairwise-like) reductions instead of cumulative
  prefix sums over the flattened array, which are accurate to about
  n * eps (about 4e-13 relative for the 1e7 elements of a production source);
* the normalized-shape metric multiplies by the reciprocal amplitude instead
  of dividing elementwise (a difference of one unit in the last place);
* a candidate is kept without a pass over the arrays when the triangle
  inequality in the maximum norm, applied to the intermediate steps, proves
  that all its metrics are within the tolerance; the others, and the one that
  reaches twice its count (which supplies the confirmation record), are
  evaluated exactly.

The streaming passes over the arrays are the Numba kernels of
``_compiled_integrated_weights`` (no fastmath, no parallel loops), imported
lazily. Compiled arithmetic does not trap floating-point underflow; a
nonfinite full sum raises FloatingPointError, which the recurrence records as
an arithmetic failure. The final moments of :func:`integrated_moments` remain a
chunked NumPy reduction that traps underflow.
"""

import numpy as np

from .full_sum_weights import METRICS, relative_change

# Elements per streaming chunk of the final NumPy moments (256 kB per float64
# buffer, resident in the L2 cache).
CHUNK = 1 << 15


def _moments_chunk(measure, weights, variance, work):
    """Return the three weighted partial sums of one chunk.

    Parameters
    ----------
    measure, weights, variance : ndarray of shape (n,)
        Measure in deg^-2, dimensionless weights and pixel variances.
    work : ndarray of shape (n,)
        Scratch array, overwritten.

    Returns
    -------
    n1, n2, n3 : float
        Partial sums of mu w, mu w^2 and mu w^2 v, formed by the sequential
        products of ``kernels.full_sum_weights._integrals``.
    """
    np.multiply(measure, weights, out=work)
    n1 = np.add.reduce(work)
    np.multiply(work, weights, out=work)
    n2 = np.add.reduce(work)
    np.multiply(work, variance, out=work)
    n3 = np.add.reduce(work)
    return n1, n2, n3


class _CompiledPasses:
    """Streaming passes of the recurrence, compiled with Numba.

    The compiled kernels cannot raise NumPy floating-point errors; they flag
    nonfinite full sums with a nonzero status, which is raised here as
    FloatingPointError so that the recurrence records an arithmetic failure.
    Compiled arithmetic does not trap underflow.
    """

    def __init__(self):
        """Import the compiled kernels (lazily, keeping Numba out of import)."""
        from . import _compiled_integrated_weights as module

        self._module = module

    def moments(self, measure, variance, coefficient):
        """Return N1, N2, N3 of the weights w = c / (c + v).

        Parameters
        ----------
        measure, variance : ndarray of shape (n,)
            Flattened float64 measure in deg^-2 and dimensionless pixel
            variance.
        coefficient : float
            Scalar c of the state.

        Returns
        -------
        n1, n2, n3 : float
            Full sums of mu w, mu w^2 and mu w^2 v.

        Raises
        ------
        FloatingPointError
            If a full sum is nonfinite.
        """
        *values, status = self._module.moments_pass(measure, variance, coefficient)
        if status:
            raise FloatingPointError(
                "moments pass: nonfinite full sums N1, N2, N3 of the weights"
            )
        return tuple(values)

    def step(self, measure, variance, *scalars):
        """Return the new state's moments and its change from the previous state.

        Parameters
        ----------
        measure, variance : ndarray of shape (n,)
            Flattened float64 measure in deg^-2 and dimensionless pixel
            variance.
        *scalars : float
            coefficient_new, coefficient_old, reciprocal_new, reciprocal_old:
            scalars c of the two states and their 1 / max(w).

        Returns
        -------
        n1, n2, n3 : float
            Full sums of mu w, mu w^2 and mu w^2 v of the new weights.
        max_difference, max_shape_difference : float
            max|w_new - w_old| and the maximum change of the amplitude-
            normalized weights.

        Raises
        ------
        FloatingPointError
            If a full sum of the new state is nonfinite.
        """
        *values, status = self._module.step_pass(measure, variance, *scalars)
        if status:
            raise FloatingPointError(
                "step pass: nonfinite full sums N1, N2, N3 of the updated weights"
            )
        return tuple(values)

    def difference(self, variance, *scalars):
        """Return max|w_new - w_old| and the maximum normalized-shape change.

        Parameters
        ----------
        variance : ndarray of shape (n,)
            Flattened float64 dimensionless pixel variance.
        *scalars : float
            coefficient_new, coefficient_old, reciprocal_new, reciprocal_old.

        Returns
        -------
        max_difference, max_shape_difference : float
            Changes between the two states; the normalization multiplies by
            the reciprocal amplitude instead of dividing elementwise (a
            difference of one unit in the last place).
        """
        return self._module.difference_pass(variance, *scalars)

    def fill(self, variance, coefficient, out):
        """Write the weights w = c / (c + v) of a state into ``out``.

        Parameters
        ----------
        variance : ndarray of shape (n,)
            Flattened float64 dimensionless pixel variance.
        coefficient : float
            Scalar c of the state.
        out : ndarray of shape (n,)
            Receives the dimensionless weights.
        """
        self._module.fill_pass(variance, coefficient, out)


def coefficients_from_moments(n1, n2, n3, pixel):
    """Normalize full-sum moments as ``full_sum_weights.coefficients`` does.

    Parameters
    ----------
    n1, n2, n3 : float
        N1, N2, N3 in deg^-2.
    pixel : float
        Pixel width in km/s.

    Returns
    -------
    coefficients : ndarray of shape (5,)
        N1, N2, N3, A = N2/N1^2 in deg^2 and P_pixel = pixel N3/N1^2.

    Raises
    ------
    FloatingPointError
        If the controlled arithmetic fails or an output is nonfinite.
    """
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        first, second, noise = np.float64(n1), np.float64(n2), np.float64(n3)
        denominator = first**2 * 1.0
        result = np.array(
            [first, second, noise, second / denominator, noise * pixel / denominator]
        )
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("nonfinite final moments or coefficients")
    return result


def _amplitude(variance_min, coefficient):
    """Return max(w) = c / (c + min(v)) of the state with parameter c.

    Parameters
    ----------
    variance_min : float
        Smallest pixel variance.
    coefficient : float
        Positive scalar c of the state.

    Returns
    -------
    amplitude : float
        w is a nonincreasing function of v in floating-point arithmetic as
        well (the sum and the quotient are monotone), so the maximum over the
        elements is attained at the smallest variance and equals the array
        maximum exactly.

    Raises
    ------
    FloatingPointError
        Under the caller's NumPy error policy.
    """
    return np.float64(coefficient) / (np.float64(variance_min) + coefficient)


def _metrics(max_difference, max_shape, amplitude_new, amplitude_old, new_c, old_c):
    """Assemble the five metrics of ``full_sum_weights.residuals``.

    Parameters
    ----------
    max_difference, max_shape : float
        max|w_new - w_old| and the maximum difference of the amplitude-
        normalized weights.
    amplitude_new, amplitude_old : float
        Positive max(w) of the two states.
    new_c, old_c : ndarray of shape (5,)
        Their coefficients (N1, N2, N3, A, P_pixel).

    Returns
    -------
    residuals : ndarray of shape (5,)
        Amplitude, shape, A, P_pixel and vector changes (``METRICS`` order).
    """
    return np.array(
        [
            abs(amplitude_new / amplitude_old - 1),
            # max|w_old / amplitude_old| is exactly 1 (x / x = 1).
            max_shape,
            relative_change(new_c[-2:-1], old_c[-2:-1]),
            relative_change(new_c[-1:], old_c[-1:]),
            max_difference / amplitude_old,
        ]
    )


def _proved_within(entry, new_c, amplitude_new, rtol, slack):
    """Decide from the intermediate steps that a candidate stays within rtol.

    Parameters
    ----------
    entry : list
        Candidate record [count, coefficient, coefficients, amplitude,
        summed_difference, summed_shape], the last two being the sums over the
        steps since the candidate of max|w_s - w_{s-1}| and of the maximum
        differences of the amplitude-normalized weights.
    new_c : ndarray of shape (5,)
        Coefficients (N1, N2, N3, A, P_pixel) of the current state.
    amplitude_new : float
        max(w) of the current state.
    rtol : float
        Tolerance of every metric.
    slack : float
        Absolute allowance for floating-point rounding of the sums.

    Returns
    -------
    proved : bool
        True only if all five metrics of the candidate-to-current comparison
        are guaranteed to be at most rtol: the amplitude, A and P_pixel
        metrics are evaluated exactly, and the weight-vector and shape metrics
        are bounded by the triangle inequality in the maximum norm. False means
        undecided; the caller then evaluates the metrics exactly.
    """
    _, _, candidate_c, amplitude, summed_difference, summed_shape = entry
    return bool(
        abs(amplitude_new / amplitude - 1) <= rtol
        and relative_change(new_c[-2:-1], candidate_c[-2:-1]) <= rtol
        and relative_change(new_c[-1:], candidate_c[-1:]) <= rtol
        and summed_difference / amplitude + slack <= rtol
        and summed_shape + slack <= rtol
    )


def _update_coefficient(first_moment, pixel, signal, alias):
    """Return the scalar c of the 'sum_historical' update w = c / (c + v).

    Parameters
    ----------
    first_moment : float
        N1 of the previous state in deg^-2.
    pixel, signal, alias : float
        Pixel width in km/s, S in deg^2 km/s and B in km/s.

    Returns
    -------
    coefficient : float
        c = S' k with S' = S + B / N1 and k = N1 / pixel. The reference update
        S' / (S' + v / k) is algebraically c / (c + v); the two evaluations
        agree to a few units in the last place per element, and the second
        needs one division instead of two.

    Raises
    ------
    FloatingPointError
        If N1 vanishes or the arithmetic overflows or is invalid.
    """
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        first = np.float64(first_moment)
        signal_total = signal + alias / (first * 1.0)
        scale = first * (1.0 / pixel)
        coefficient = signal_total * scale
    return coefficient


def integrated_moments(measure, weights, variance, pixel):
    """Normalize the full-sum moments of final weights, trapping underflow.

    Parameters
    ----------
    measure, weights, variance : ndarray of shape (n,)
        Flattened measure in deg^-2, dimensionless weights and pixel variances.
    pixel : float
        Pixel width in km/s.

    Returns
    -------
    n1, n2, n3 : float
        N1, N2, N3 in deg^-2.
    aliasing_coefficient : float
        A = N2 / N1^2 in deg^2.
    pixel_power : float
        P_pixel = pixel N3 / N1^2 in deg^2 km/s.

    Raises
    ------
    FloatingPointError
        If an individual weighted product underflows inexactly (as
        ``kernels.full_sum_weights._integrals``), or the arithmetic is invalid.
    """
    n_total = len(measure)
    work = np.empty(min(CHUNK, n_total))
    n1 = n2 = n3 = 0.0
    # Trap individual inexact underflows before a later positive sample can
    # conceal the lost contribution, as in kernels.full_sum_weights._integrals.
    with np.errstate(under="raise"):
        for start in range(0, n_total, CHUNK):
            stop = min(start + CHUNK, n_total)
            a1, a2, a3 = _moments_chunk(
                measure[start:stop],
                weights[start:stop],
                variance[start:stop],
                work[: stop - start],
            )
            n1, n2, n3 = n1 + a1, n2 + a2, n3 + a3
    first, second, noise = np.float64(n1), np.float64(n2), np.float64(n3)
    aliasing_coefficient = second / first / first / 1.0
    pixel_power = noise / first / first * pixel / 1.0
    return n1, n2, n3, aliasing_coefficient, pixel_power


def _run(measure, variance, pixel, signal, alias, controls, passes, result):
    """Iterate the recurrence, filling ``result``; see :func:`solve_integrated`.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened measure and pixel variance.
    pixel, signal, alias : float
        Pixel width, S and B.
    controls : tuple
        (rtol, min_updates, stable_steps, max_updates).
    passes : _CompiledPasses
        Streaming passes over the arrays.
    result : dict
        Convergence record, updated in place; its 'weights' entry is left
        unset.

    Returns
    -------
    coefficient : float or None
        Scalar c of the state that the record describes (None if the seed
        failed).
    """
    rtol, min_updates, stable_steps, max_updates = controls
    candidates = []
    stable = 0
    state = None
    margin = 16 * np.finfo(float).eps
    try:
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            variance_min = variance.min()
            coefficient = alias / pixel
            amplitude = _amplitude(variance_min, coefficient)
            coefficients = coefficients_from_moments(
                *passes.moments(measure, variance, coefficient), pixel
            )
            result.update(coefficients=coefficients, state_updates=0)
            state = coefficient
            for t in range(1, max_updates + 1):
                new_state = _update_coefficient(coefficients[0], pixel, signal, alias)
                new_amplitude = _amplitude(variance_min, new_state)
                n1, n2, n3, max_difference, max_shape = passes.step(
                    measure,
                    variance,
                    new_state,
                    state,
                    1.0 / new_amplitude,
                    1.0 / amplitude,
                )
                result["updates"] = t
                new_coefficients = coefficients_from_moments(n1, n2, n3, pixel)
                step = _metrics(
                    max_difference,
                    max_shape,
                    new_amplitude,
                    amplitude,
                    new_coefficients,
                    coefficients,
                )
                result.update(
                    coefficients=new_coefficients,
                    state_updates=t,
                    last_step=dict(zip(METRICS, step.tolist(), strict=True)),
                )
                state = new_state
                stable = stable + 1 if np.all(step <= rtol) else 0

                # Candidates must stay within rtol of the current state. Their
                # distance to it is bounded by the sum of the intermediate
                # steps (triangle inequality in the maximum norm, for the
                # weights and for the amplitude-normalized weights); a
                # candidate whose bound is within rtol passes without a pass
                # over the arrays, and only the others, and the one reaching
                # twice its count, are evaluated exactly.
                kept = []
                confirmation = None
                for entry in candidates:
                    if not stable:
                        break
                    count, candidate_state, candidate_c, candidate_amplitude = entry[:4]
                    entry[4] += max_difference
                    entry[5] += max_shape
                    if t != 2 * count and _proved_within(
                        entry, new_coefficients, new_amplitude, rtol, margin * (t + 2)
                    ):
                        kept.append(entry)
                        continue
                    exact_difference, exact_shape = passes.difference(
                        variance,
                        new_state,
                        candidate_state,
                        1.0 / new_amplitude,
                        1.0 / candidate_amplitude,
                    )
                    residual = _metrics(
                        exact_difference,
                        exact_shape,
                        new_amplitude,
                        candidate_amplitude,
                        new_coefficients,
                        candidate_c,
                    )
                    if np.all(residual <= rtol):
                        kept.append(entry)
                        if t == 2 * count and confirmation is None:
                            confirmation = (count, residual)
                candidates = kept
                if confirmation is not None:
                    result.update(
                        status="converged",
                        candidate=confirmation[0],
                        confirmation=dict(
                            zip(METRICS, confirmation[1].tolist(), strict=True)
                        ),
                    )
                    return state
                if t >= min_updates and stable >= stable_steps and 2 * t <= max_updates:
                    candidates.append(
                        [t, new_state, new_coefficients, new_amplitude, 0.0, 0.0]
                    )
                coefficients, amplitude = new_coefficients, new_amplitude
    except FloatingPointError as error:
        result.update(status="arithmetic_failure", reason=str(error))
        return state
    result.update(
        status="capped", reason="no confirmed finite nonzero convergence within cap"
    )
    return state


def solve_integrated(
    measure,
    variance,
    *,
    pixel,
    signal,
    alias,
    rtol=1e-4,
    min_updates=3,
    stable_steps=3,
    max_updates=96,
):
    """Iterate the 'sum_historical' recurrence to confirmed convergence.

    Reproduces ``full_sum_weights.solve(inputs, 'sum_historical', ...)`` for
    ``inputs.density = measure``, ``inputs.quadrature = 1``, ``inputs.length =
    1`` and ``inputs.p1d = alias`` with the same seed, update, stopping
    metrics, order of evaluation, doubled-count confirmation, cap and record;
    see the module documentation for the numerical differences.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened C-order float64 measure in deg^-2 and dimensionless pixel
        variance (nonnegative, finite).
    pixel : float
        Pixel width in km/s.
    signal : float
        Reference S in deg^2 km/s.
    alias : float
        Reference B in km/s.
    rtol : float, default=1e-4
        Positive tolerance for all five relative-change metrics.
    min_updates, stable_steps, max_updates : int
        Update controls as in the reference solver.

    Returns
    -------
    result : dict
        Keys as in the reference solver: status, reason, weights (n,),
        coefficients (5,), updates, state_updates, candidate, last_step,
        confirmation, forward_residual. ``weights`` is the (writeable)
        weight vector of the last finite state, written once at the end.

    Raises
    ------
    ValueError
        If the integer controls or the relative tolerance are invalid.
    """
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
    if not np.isfinite(signal) or signal <= 0 or not np.isfinite(alias) or alias <= 0:
        result["reason"] = (
            "adaptive stopping requires finite positive intrinsic signal and P1D reference"
        )
        return result

    passes = _CompiledPasses()
    controls = (rtol, min_updates, stable_steps, max_updates)
    state = _run(measure, variance, pixel, signal, alias, controls, passes, result)
    if state is not None:
        weights = np.empty_like(measure)
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            passes.fill(variance, state, weights)
        result["weights"] = weights
    return result
