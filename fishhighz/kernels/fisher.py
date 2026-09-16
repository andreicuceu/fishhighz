"""Numeric-only Cholesky, forward substitution, and Fisher contraction.

Callers validate shapes and attach numerical-failure context. No models or names
are accessed; no compilation or automatic parallelism is enabled.
"""

import numpy as np


def _cholesky(matrix):
    """Return lower Cholesky factor of one validated SPD matrix."""
    return np.linalg.cholesky(matrix)


def _forward_substitute(lower, rhs, out):
    """Overwrite out with L^-1 rhs for all RHS columns; return None.

    lower is (n,n), rhs/out are (n,p), all float64, with positive L diagonal.
    out must not alias inputs. Work beyond out is at most one p-vector.
    """
    for row in range(len(lower)):
        out[row] = (rhs[row] - lower[row, :row] @ out[:row]) / lower[row, row]


def _accumulate_fisher(solved, out):
    """Add solved.T @ solved to out; return None, using one p-by-p buffer."""
    out += solved.T @ solved
