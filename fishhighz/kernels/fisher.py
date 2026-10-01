"""Numeric-only Cholesky, forward substitution, and Fisher contraction.

Callers validate shapes and attach numerical-failure context. No models or names
are accessed; no compilation or automatic parallelism is enabled.
"""

import numpy as np


def _cholesky(matrix):
    """Factor one validated positive-definite covariance matrix.

    Parameters
    ----------
    matrix : ndarray of shape (n_pair, n_pair)
        Symmetric positive-definite covariance, in power units squared.

    Returns
    -------
    lower : ndarray of shape (n_pair, n_pair)
        Lower Cholesky factor, in power units.

    Raises
    ------
    numpy.linalg.LinAlgError
        If the covariance is not positive definite.
    """
    return np.linalg.cholesky(matrix)


def _forward_substitute(lower, rhs, out):
    """Solve the lower-triangular system for every parameter column.

    Parameters
    ----------
    lower : ndarray of shape (n_pair, n_pair)
        Float64 lower Cholesky factor with positive diagonal, in power units.
    rhs : ndarray of shape (n_pair, n_parameter)
        Float64 Jacobian in power units per parameter unit.
    out : ndarray of shape (n_pair, n_parameter)
        Float64 destination, which must not alias either input.

    Returns
    -------
    None
        Overwrite out with the whitened Jacobian, in inverse parameter units.

    Notes
    -----
    The caller validates inputs. At most one parameter-vector temporary is
    needed beyond the output array.
    """
    for row in range(len(lower)):
        out[row] = (rhs[row] - lower[row, :row] @ out[:row]) / lower[row, row]


def _accumulate_fisher(solved, out):
    """Add one cell's whitened-Jacobian contraction to the Fisher matrix.

    Parameters
    ----------
    solved : ndarray of shape (n_pair, n_parameter)
        Whitened Jacobian, in inverse parameter units.
    out : ndarray of shape (n_parameter, n_parameter)
        Accumulated Fisher information, in inverse products of parameter units.

    Returns
    -------
    None
        Add solved.T @ solved to out in place, using one matrix temporary.
    """
    out += solved.T @ solved
