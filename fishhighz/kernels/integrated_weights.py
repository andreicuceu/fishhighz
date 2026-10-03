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

Two interchangeable backends supply the streaming passes: a chunked NumPy
implementation (this module, the reference for the compiled one) and a Numba
implementation in ``_compiled_integrated_weights`` (no fastmath, no parallel
loops). Select with ``FISHHIGHZ_INTEGRATED_BACKEND`` = ``numba`` (default,
falling back to NumPy when Numba is unavailable) or ``numpy``.
"""

import os
from types import SimpleNamespace

import numpy as np

from .full_sum_weights import METRICS, relative_change

# Elements per streaming chunk of the NumPy backend (256 kB per float64 buffer,
# resident in the L2 cache).
CHUNK = 1 << 15

BACKENDS = ("numpy", "numba")


def _max_abs(chunk, buffer):
    """Return the maximum absolute value of a 1-D chunk without a temporary.

    Parameters
    ----------
    chunk : ndarray of shape (n,)
        Values; n > 0.
    buffer : ndarray of shape (n,)
        Scratch array of the same length, overwritten (may be ``chunk``).

    Returns
    -------
    value : float
        max(|chunk|).
    """
    np.abs(chunk, out=buffer)
    return buffer.max()


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
        products of ``kernels.weights._integrals``.
    """
    np.multiply(measure, weights, out=work)
    n1 = np.add.reduce(work)
    np.multiply(work, weights, out=work)
    n2 = np.add.reduce(work)
    np.multiply(work, variance, out=work)
    n3 = np.add.reduce(work)
    return n1, n2, n3


def _state_weights(variance, coefficient, out):
    """Evaluate w = c / (c + v) of one chunk.

    Parameters
    ----------
    variance : ndarray of shape (n,)
        Dimensionless pixel variances.
    coefficient : float
        Scalar c of the state.
    out : ndarray of shape (n,)
        Receives the weights.

    Raises
    ------
    FloatingPointError
        Under the caller's NumPy error policy (divide, invalid, over).
    """
    np.add(variance, coefficient, out=out)
    np.divide(coefficient, out, out=out)


def numpy_moments_pass(measure, variance, coefficient):
    """Evaluate the moments of the weights w = c / (c + variance), chunked.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened measure in deg^-2 and dimensionless pixel variance.
    coefficient : float
        Scalar c of the state (P1D/pixel for the seed).

    Returns
    -------
    n1, n2, n3 : float
        Full sums of mu w, mu w^2 and mu w^2 v.

    Raises
    ------
    FloatingPointError
        Under the caller's NumPy error policy.
    """
    n_total = len(measure)
    weights = np.empty(min(CHUNK, n_total))
    work = np.empty_like(weights)
    n1 = n2 = n3 = 0.0
    for start in range(0, n_total, CHUNK):
        stop = min(start + CHUNK, n_total)
        size = stop - start
        _state_weights(variance[start:stop], coefficient, weights[:size])
        a1, a2, a3 = _moments_chunk(
            measure[start:stop], weights[:size], variance[start:stop], work[:size]
        )
        n1, n2, n3 = n1 + a1, n2 + a2, n3 + a3
    return n1, n2, n3


def numpy_step_pass(
    measure, variance, coefficient_new, coefficient_old, reciprocal_new, reciprocal_old
):
    """Evaluate the new state's moments and its change from the previous state.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened measure in deg^-2 and dimensionless pixel variance.
    coefficient_new, coefficient_old : float
        Scalars c of the new and the previous state.
    reciprocal_new, reciprocal_old : float
        1 / max(w) of the two states.

    Returns
    -------
    n1, n2, n3 : float
        Full sums of mu w, mu w^2 and mu w^2 v of the new weights.
    max_difference : float
        max|w_new - w_old|.
    max_shape_difference : float
        max|w_new * reciprocal_new - w_old * reciprocal_old|; see
        :func:`numpy_difference_pass`.

    Raises
    ------
    FloatingPointError
        Under the caller's NumPy error policy.
    """
    n_total = len(measure)
    new = np.empty(min(CHUNK, n_total))
    old, work, other = np.empty_like(new), np.empty_like(new), np.empty_like(new)
    n1 = n2 = n3 = 0.0
    max_difference = max_shape = 0.0
    for start in range(0, n_total, CHUNK):
        stop = min(start + CHUNK, n_total)
        size = stop - start
        chunk_new, chunk_old = new[:size], old[:size]
        _state_weights(variance[start:stop], coefficient_new, chunk_new)
        _state_weights(variance[start:stop], coefficient_old, chunk_old)
        a1, a2, a3 = _moments_chunk(
            measure[start:stop], chunk_new, variance[start:stop], work[:size]
        )
        n1, n2, n3 = n1 + a1, n2 + a2, n3 + a3
        np.subtract(chunk_new, chunk_old, out=work[:size])
        max_difference = max(max_difference, _max_abs(work[:size], work[:size]))
        np.multiply(chunk_new, reciprocal_new, out=work[:size])
        np.multiply(chunk_old, reciprocal_old, out=other[:size])
        np.subtract(work[:size], other[:size], out=work[:size])
        max_shape = max(max_shape, _max_abs(work[:size], work[:size]))
    return n1, n2, n3, max_difference, max_shape


def numpy_difference_pass(
    variance, coefficient_new, coefficient_old, reciprocal_new, reciprocal_old
):
    """Evaluate the change between two states of the one-parameter family.

    Parameters
    ----------
    variance : ndarray of shape (n,)
        Flattened dimensionless pixel variance.
    coefficient_new, coefficient_old : float
        Scalars c of the two states.
    reciprocal_new, reciprocal_old : float
        1 / max(w) of the two states.

    Returns
    -------
    max_difference : float
        max|w_new - w_old|.
    max_shape_difference : float
        max|w_new * reciprocal_new - w_old * reciprocal_old|, the change of the
        amplitude-normalized weights.

    Notes
    -----
    The normalization multiplies by the reciprocal amplitude instead of
    dividing elementwise (the reference ``full_sum_weights.changes`` divides):
    the two differ by at most one unit in the last place of each normalized
    weight, i.e. by about 2e-16 in absolute terms against a tolerance of 1e-4,
    and the multiplication is several times cheaper than a division on the
    1e7-element arrays of an integrated source.
    """
    n_total = len(variance)
    new = np.empty(min(CHUNK, n_total))
    old, work = np.empty_like(new), np.empty_like(new)
    max_difference = max_shape = 0.0
    for start in range(0, n_total, CHUNK):
        stop = min(start + CHUNK, n_total)
        size = stop - start
        chunk_new, chunk_old = new[:size], old[:size]
        _state_weights(variance[start:stop], coefficient_new, chunk_new)
        _state_weights(variance[start:stop], coefficient_old, chunk_old)
        np.subtract(chunk_new, chunk_old, out=work[:size])
        max_difference = max(max_difference, _max_abs(work[:size], work[:size]))
        np.multiply(chunk_new, reciprocal_new, out=chunk_new)
        np.multiply(chunk_old, reciprocal_old, out=chunk_old)
        np.subtract(chunk_new, chunk_old, out=work[:size])
        max_shape = max(max_shape, _max_abs(work[:size], work[:size]))
    return max_difference, max_shape


def numpy_fill_pass(variance, coefficient, out):
    """Write the weights w = c / (c + variance) of a state, chunked.

    Parameters
    ----------
    variance : ndarray of shape (n,)
        Flattened dimensionless pixel variance.
    coefficient : float
        Scalar c of the state.
    out : ndarray of shape (n,)
        Receives the weights.
    """
    for start in range(0, len(variance), CHUNK):
        stop = min(start + CHUNK, len(variance))
        _state_weights(variance[start:stop], coefficient, out[start:stop])


NUMPY_BACKEND = SimpleNamespace(
    name="numpy",
    moments=numpy_moments_pass,
    step=numpy_step_pass,
    difference=numpy_difference_pass,
    fill=numpy_fill_pass,
)


class _CompiledBackend:
    """Compiled passes that replay a failed pass through the NumPy reference.

    The compiled kernels cannot raise NumPy floating-point errors; they flag
    nonfinite sums, and the host then repeats the pass with the NumPy
    implementation to reproduce the exact exception.

    Parameters
    ----------
    module : module
        ``_compiled_integrated_weights``.
    """

    name = "numba"

    def __init__(self, module):
        """Store the compiled passes."""
        self._module = module

    def moments(self, measure, variance, coefficient):
        """Compiled :func:`numpy_moments_pass` with NumPy replay on failure."""
        *values, status = self._module.moments_pass(measure, variance, coefficient)
        if status:
            return numpy_moments_pass(measure, variance, coefficient)
        return tuple(values)

    def step(self, measure, variance, *scalars):
        """Compiled :func:`numpy_step_pass` with NumPy replay on failure."""
        *values, status = self._module.step_pass(measure, variance, *scalars)
        if status:
            return numpy_step_pass(measure, variance, *scalars)
        return tuple(values)

    def difference(self, variance, *scalars):
        """Compiled :func:`numpy_difference_pass`."""
        return self._module.difference_pass(variance, *scalars)

    def fill(self, variance, coefficient, out):
        """Compiled :func:`numpy_fill_pass`."""
        self._module.fill_pass(variance, coefficient, out)


def select_backend(name=None):
    """Resolve the streaming-pass backend.

    Parameters
    ----------
    name : {'numpy', 'numba'} or None, optional
        Requested backend; None reads ``FISHHIGHZ_INTEGRATED_BACKEND`` and
        defaults to ``'numba'``.

    Returns
    -------
    backend : object
        Object with ``moments``, ``step``, ``difference`` and ``fill`` passes.
        Numba is used only when it imports and the NumPy underflow policy is
        'ignore' (compiled arithmetic cannot trap underflow); otherwise the
        NumPy backend is returned.

    Raises
    ------
    ValueError
        If the requested backend name is unknown.
    """
    if name is None:
        name = os.environ.get("FISHHIGHZ_INTEGRATED_BACKEND", "numba")
    if name not in BACKENDS:
        raise ValueError("FISHHIGHZ_INTEGRATED_BACKEND must be numpy or numba")
    if name == "numba" and np.geterr()["under"] == "ignore":
        try:
            from . import _compiled_integrated_weights as module
        except ImportError:
            # The NumPy-only installation stays executable.
            return NUMPY_BACKEND
        return _CompiledBackend(module)
    return NUMPY_BACKEND


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


def fixed_weights_integrated(
    measure, variance, *, pixel, signal, alias, updates, backend=None
):
    """Return the weights after a prescribed number of updates from the seed.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened C-order float64 measure in deg^-2 and dimensionless pixel
        variance.
    pixel : float
        Pixel width in km/s.
    signal : float
        Reference S in deg^2 km/s.
    alias : float
        Reference B in km/s.
    updates : int
        Nonnegative number of updates after the seed.
    backend : {'numpy', 'numba'} or None, optional
        See :func:`select_backend`.

    Returns
    -------
    weights : ndarray of shape (n,)
        Dimensionless weights, equal to ``full_sum_weights.fixed_weights`` for
        'sum_historical' on the flattened inputs up to the rounding of the
        update and the summation order of N1.

    Raises
    ------
    ValueError
        If the update count is invalid.
    FloatingPointError
        If the controlled arithmetic fails or the weights are nonfinite.
    """
    if isinstance(updates, bool) or not isinstance(updates, int) or updates < 0:
        raise ValueError("updates must be a nonnegative integer")
    passes = select_backend(backend)
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        coefficient = alias / pixel
        for _ in range(updates):
            n1, *_ = passes.moments(measure, variance, coefficient)
            coefficient = _update_coefficient(n1, pixel, signal, alias)
        weights = np.empty_like(measure)
        passes.fill(variance, coefficient, weights)
    return weights


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
        ``kernels.weights._integrals``), or the arithmetic is invalid.
    """
    n_total = len(measure)
    work = np.empty(min(CHUNK, n_total))
    n1 = n2 = n3 = 0.0
    # Trap individual inexact underflows before a later positive sample can
    # conceal the lost contribution, as in kernels.weights._integrals.
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
    passes : object
        Streaming-pass backend.
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
    backend=None,
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
    backend : {'numpy', 'numba'} or None, optional
        See :func:`select_backend`.

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

    passes = select_backend(backend)
    controls = (rtol, min_updates, stable_steps, max_updates)
    state = _run(measure, variance, pixel, signal, alias, controls, passes, result)
    if state is not None:
        weights = np.empty_like(measure)
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            passes.fill(variance, state, weights)
        result["weights"] = weights
    return result
