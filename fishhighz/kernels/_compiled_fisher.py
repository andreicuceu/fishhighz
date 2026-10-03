"""Compiled Fisher contraction (Numba); imported lazily by ``fisher``.

One cell of solve workspace, no fastmath, no parallel loops. Matrix products
retain the NumPy ordering of ``(L^-1 J).T @ (L^-1 J)``. Failures are reported
as a status code with the failing Fourier node and parameter column, from which
the host raises a contextual error.
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
        Zero on success; otherwise the first failed check of the first failing
        node, tested in this order: 1 nonfinite factor entry, 2 nonfinite
        Jacobian entry, 3 factor not exactly lower triangular or with a
        nonpositive diagonal, 4 nonfinite whitened Jacobian (solve),
        5 nonfinite accumulated information.
    node : int
        Index of the failing Fourier node; -1 on success.
    column : int
        Global parameter column of the first nonfinite solve entry, in row-major
        (pair, parameter) order, for status 4; -1 otherwise.

    Notes
    -----
    Numba compilation uses neither fastmath nor parallel loops. Compiled
    arithmetic does not trap floating-point underflow. Validation loops do not
    modify any array entering the solve or the accumulation.
    """
    n_node, n_pair, n_parameter = jacobian.shape
    result = np.zeros((n_parameter, n_parameter))
    solved = np.empty((n_pair, n_parameter))
    for node in range(n_node):
        # Validation order of the per-node checks: finite factor, finite
        # Jacobian, then factor structure.
        for i in range(n_pair):
            for j in range(n_pair):
                if not np.isfinite(lower[node, i, j]):
                    return result, 1, node, -1
        for i in range(n_pair):
            for col in range(n_parameter):
                if not np.isfinite(jacobian[node, i, col]):
                    return result, 2, node, -1
        for i in range(n_pair):
            for j in range(i + 1, n_pair):
                if lower[node, i, j] != 0:
                    return result, 3, node, -1
            if lower[node, i, i] <= 0:
                return result, 3, node, -1

        for row in range(n_pair):
            solved[row] = (
                jacobian[node, row] - lower[node, row, :row] @ solved[:row]
            ) / lower[node, row, row]
            for col in range(n_parameter):
                if not np.isfinite(solved[row, col]):
                    return result, 4, node, col
        result += solved.T @ solved
        for i in range(n_parameter):
            for j in range(n_parameter):
                if not np.isfinite(result[i, j]):
                    return result, 5, node, -1
    return result, 0, -1, -1
