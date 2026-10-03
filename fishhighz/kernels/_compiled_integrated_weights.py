"""Optional compiled streaming passes of the integrated-source weight recurrence.

Every state of the 'sum_historical' recurrence is the one-parameter family
w = c / (c + v) of the pixel variances v, so the passes below take the scalar c
of a state instead of a weight array: no (pixel, magnitude)-sized array is
stored or re-read besides the measure and the variance. Fused loops process
blocks of ``BLOCK`` elements; the division of a block is a separate loop
without reductions, which LLVM vectorizes (``error_model='numpy'`` removes the
Python division-by-zero branch), and the sums are accumulated block-wise
(partial sums of one block added to the running total), mirroring the
pairwise-like NumPy reductions. No fastmath and no parallel loops; every
elementwise operation is the IEEE operation of the NumPy reference in
``integrated_weights``. The kernels cannot raise NumPy floating-point errors;
they return a nonzero status for nonfinite sums and the host replays the pass
through the NumPy reference.
"""

import numpy as np
from numba import njit

BLOCK = 2048


@njit(cache=True, fastmath=False, error_model="numpy")
def moments_pass(measure, variance, coefficient):
    """Evaluate the moments of the weights w = c / (c + variance).

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened float64 measure in deg^-2 and dimensionless pixel variance.
    coefficient : float
        Scalar c of the state (P1D/pixel for the seed).

    Returns
    -------
    n1, n2, n3 : float
        Full sums of mu w, mu w^2 and mu w^2 v.
    status : int
        Zero on success; one if a sum is nonfinite.
    """
    n_total = measure.shape[0]
    weights = np.empty(BLOCK)
    n1 = n2 = n3 = 0.0
    for start in range(0, n_total, BLOCK):
        size = min(BLOCK, n_total - start)
        for j in range(size):
            weights[j] = coefficient / (variance[start + j] + coefficient)
        b1 = b2 = b3 = 0.0
        for j in range(size):
            first = measure[start + j] * weights[j]
            second = first * weights[j]
            b1 += first
            b2 += second
            b3 += second * variance[start + j]
        n1 += b1
        n2 += b2
        n3 += b3
    status = 0 if np.isfinite(n1 + n2 + n3) else 1
    return n1, n2, n3, status


@njit(cache=True, fastmath=False, error_model="numpy")
def step_pass(
    measure, variance, coefficient_new, coefficient_old, reciprocal_new, reciprocal_old
):
    """Evaluate the new state's moments and its change from the previous state.

    Parameters
    ----------
    measure, variance : ndarray of shape (n,)
        Flattened float64 measure in deg^-2 and dimensionless pixel variance.
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
        max|w_new * reciprocal_new - w_old * reciprocal_old|.
    status : int
        Zero on success; one if a sum is nonfinite.
    """
    n_total = measure.shape[0]
    new = np.empty(BLOCK)
    old = np.empty(BLOCK)
    n1 = n2 = n3 = 0.0
    max_difference = 0.0
    max_shape = 0.0
    for start in range(0, n_total, BLOCK):
        size = min(BLOCK, n_total - start)
        for j in range(size):
            new[j] = coefficient_new / (variance[start + j] + coefficient_new)
            old[j] = coefficient_old / (variance[start + j] + coefficient_old)
        b1 = b2 = b3 = 0.0
        for j in range(size):
            first = measure[start + j] * new[j]
            second = first * new[j]
            b1 += first
            b2 += second
            b3 += second * variance[start + j]
        n1 += b1
        n2 += b2
        n3 += b3
        for j in range(size):
            difference = abs(new[j] - old[j])
            if difference > max_difference:
                max_difference = difference
            shape = abs(new[j] * reciprocal_new - old[j] * reciprocal_old)
            if shape > max_shape:
                max_shape = shape
    status = 0 if np.isfinite(n1 + n2 + n3) else 1
    return n1, n2, n3, max_difference, max_shape, status


@njit(cache=True, fastmath=False, error_model="numpy")
def difference_pass(
    variance, coefficient_new, coefficient_old, reciprocal_new, reciprocal_old
):
    """Evaluate the change between two states of the one-parameter family.

    Parameters
    ----------
    variance : ndarray of shape (n,)
        Flattened float64 dimensionless pixel variance.
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
        amplitude-normalized weights (the reference divides by the amplitude;
        see ``integrated_weights.numpy_difference_pass``).
    """
    n_total = variance.shape[0]
    new = np.empty(BLOCK)
    old = np.empty(BLOCK)
    max_difference = 0.0
    max_shape = 0.0
    for start in range(0, n_total, BLOCK):
        size = min(BLOCK, n_total - start)
        for j in range(size):
            new[j] = coefficient_new / (variance[start + j] + coefficient_new)
            old[j] = coefficient_old / (variance[start + j] + coefficient_old)
        for j in range(size):
            difference = abs(new[j] - old[j])
            if difference > max_difference:
                max_difference = difference
            shape = abs(new[j] * reciprocal_new - old[j] * reciprocal_old)
            if shape > max_shape:
                max_shape = shape
    return max_difference, max_shape


@njit(cache=True, fastmath=False, error_model="numpy")
def fill_pass(variance, coefficient, out):
    """Write the weights w = c / (c + variance) of a state.

    Parameters
    ----------
    variance : ndarray of shape (n,)
        Flattened float64 dimensionless pixel variance.
    coefficient : float
        Scalar c of the state.
    out : ndarray of shape (n,)
        Receives the weights.
    """
    for i in range(variance.shape[0]):
        out[i] = coefficient / (variance[i] + coefficient)
