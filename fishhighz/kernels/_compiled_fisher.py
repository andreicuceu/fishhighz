"""Explicitly requested optional contraction; imported lazily by the host.

One cell of solve workspace, no fastmath, no parallel loops. Matrix products
retain the NumPy reference ordering. Status failures are replayed by the host
through the reference implementation to retain its exact diagnostic context.
"""

import numpy as np
from numba import njit


@njit(cache=False, fastmath=False)
def contract(lower, jacobian):
    nodes, rows, parameters = jacobian.shape
    result = np.zeros((parameters, parameters))
    solved = np.empty((rows, parameters))
    for node in range(nodes):
        for i in range(rows):
            for j in range(rows):
                value = lower[node, i, j]
                if not np.isfinite(value) or (j > i and value != 0):
                    return result, 1
            if lower[node, i, i] <= 0:
                return result, 1
            for col in range(parameters):
                if not np.isfinite(jacobian[node, i, col]):
                    return result, 1
        for row in range(rows):
            solved[row] = (
                jacobian[node, row] - lower[node, row, :row] @ solved[:row]
            ) / lower[node, row, row]
            for col in range(parameters):
                if not np.isfinite(solved[row, col]):
                    return result, 1
        result += solved.T @ solved
        for i in range(parameters):
            for j in range(parameters):
                if not np.isfinite(result[i, j]):
                    return result, 1
    return result, 0
