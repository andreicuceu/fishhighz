"""Optional compiled contraction (default when Numba is installed); lazy import.

One cell of solve workspace, no fastmath, no parallel loops. Matrix products
retain the NumPy reference ordering. Status failures are replayed by the host
through the reference implementation to retain its exact diagnostic context.
"""

import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def contract(lower, jacobian):
    """Contract cell Jacobians after whitening with their Cholesky factors.

    Parameters
    ----------
    lower : ndarray of shape (n_node, n_pair, n_pair)
        Lower covariance factors in power units, with positive diagonals.
    jacobian : ndarray of shape (n_node, n_pair, n_parameter)
        Mean-power derivatives in power units per parameter unit.

    Returns
    -------
    result : ndarray of shape (n_parameter, n_parameter)
        Accumulated information in inverse products of parameter units.
        On failure this is the partial accumulation, not a valid forecast.
    status : int
        Zero on success; one for invalid factors or nonfinite arithmetic.

    Notes
    -----
    The host replays failures through the NumPy path for contextual errors.
    Numba compilation uses neither fastmath nor parallel loops.
    """
    n_node, n_pair, n_parameter = jacobian.shape
    result = np.zeros((n_parameter, n_parameter))
    solved = np.empty((n_pair, n_parameter))
    for node in range(n_node):
        for i in range(n_pair):
            for j in range(n_pair):
                value = lower[node, i, j]
                if not np.isfinite(value) or (j > i and value != 0):
                    return result, 1
            if lower[node, i, i] <= 0:
                return result, 1
            for col in range(n_parameter):
                if not np.isfinite(jacobian[node, i, col]):
                    return result, 1
        for row in range(n_pair):
            solved[row] = (
                jacobian[node, row] - lower[node, row, :row] @ solved[:row]
            ) / lower[node, row, row]
            for col in range(n_parameter):
                if not np.isfinite(solved[row, col]):
                    return result, 1
        result += solved.T @ solved
        for i in range(n_parameter):
            for j in range(n_parameter):
                if not np.isfinite(result[i, j]):
                    return result, 1
    return result, 0
