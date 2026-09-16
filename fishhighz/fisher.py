"""Fisher information from supplied means and fixed fiducial covariance.

No additional volume, mode count, quadrature, or factor of two is applied.
Known subtracted noise affects covariance; its parameter response, if any, must
be supplied in the mean Jacobian. No covariance derivative information is added.
"""

import os

import numpy as np

from ._arrays import real_array
from ._information import inspect_information
from .kernels.fisher import _accumulate_fisher, _cholesky, _forward_substitute


def _numeric_input(value, name):
    # Preserve existing array views. Dtype conversion and finite checks happen
    # one node at a time, avoiding an all-node float64 Jacobian copy.
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name}: expected real numeric data")
    return array


def _blocks(value, name):
    array = _numeric_input(value, name)
    if array.ndim != 3 or 0 in array.shape or array.shape[1] != array.shape[2]:
        raise ValueError(f"{name}: expected nonempty (n_node,n_selected,n_selected)")
    return array


def _factor_covariance_scalar(covariance, offset=0):
    """Prepare reusable lower Cholesky factors for independent covariance blocks.

    Parameters
    ----------
    covariance : array_like, shape (n_node,n_selected,n_selected)
        Real finite fixed fiducial covariance in the selected mean-spectrum order.

    Returns
    -------
    factors : ndarray
        Owned C-contiguous float64 factors L, with C=L@L.T. Inputs are unchanged.
        Roundoff-level asymmetry is averaged only in the local factorization copy.

    Notes
    -----
    Strictly positive variances are required. In R=C/sqrt(diag(C))/sqrt(diag(C)),
    symmetry tolerance is elementwise 64*eps64*n*max(1,abs(Rij),abs(Rji)).
    All eigenvalues must exceed 64*eps64*n*max(1,max(abs(eigenvalues(R)))).
    Singular/numerically unresolved blocks fail with node and rank context.
    These tests are invariant under positive observable-unit rescaling. No
    jitter, clipping, node removal, or pseudoinverse is used.
    """
    covariance = _blocks(covariance, "covariance")
    factors = np.empty(covariance.shape, dtype=np.float64)
    for node, block in enumerate(covariance):
        context = f"covariance node {node + offset}"
        _, normalized, scales, values, _, tolerance = inspect_information(
            block, context
        )
        rank = np.count_nonzero(values > tolerance)
        if np.any(block.diagonal() <= 0) or rank < len(block):
            raise ValueError(
                f"{context}: singular or numerically unresolved, rank {rank}/{len(block)}, "
                f"normalized minimum eigenvalue {values[0]:.17g}, threshold {tolerance:.17g}; "
                "change the selection or supply an appropriate physical noise model"
            )
        try:
            with np.errstate(over="ignore", invalid="ignore"):
                factors[node] = scales[:, None] * _cholesky(normalized)
        except np.linalg.LinAlgError as error:
            raise ValueError(
                f"{context}: Cholesky failed after normalized rank validation"
            ) from error
        if not np.all(np.isfinite(factors[node])) or np.any(
            factors[node].diagonal() <= 0
        ):
            raise ValueError(f"{context}: nonfinite or unrepresentable Cholesky factor")
    return factors


def factor_covariance(covariance):
    """Return reusable factors with the scalar normalized-rank/error contract.

    Validation and factorization use batches of at most 256 cells. Exceptional
    batches are replayed in cell order through the scalar diagnostic, including
    cells near the rank threshold where eigvalsh/eigh rounding can differ.
    Inputs are unchanged; no eigenvalues, powers or thresholds are modified.
    """
    covariance = _blocks(covariance, "covariance")
    factors = np.empty(covariance.shape, dtype=np.float64)
    for start in range(0, len(covariance), 256):
        block = covariance[start : start + 256]
        try:
            factors[start : start + len(block)] = _factor_batch(block)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            factors[start : start + len(block)] = _factor_covariance_scalar(
                block, offset=start
            )
    return factors


def _factor_batch(block):
    from ._information import normalize_positive_batch

    normalized, scales = normalize_positive_batch(block)
    values = np.linalg.eigvalsh(normalized)
    tolerance = (
        64
        * np.finfo(np.float64).eps
        * block.shape[-1]
        * np.maximum(1.0, np.max(np.abs(values), axis=1))
    )
    # A conservative fast-path eligibility test, not a changed acceptance rule.
    if not np.all(np.isfinite(values)) or np.any(values[:, 0] <= 2 * tolerance):
        raise ValueError("scalar rank diagnostic required")
    factors = scales[:, :, None] * np.linalg.cholesky(normalized)
    if not np.all(np.isfinite(factors)) or np.any(
        np.diagonal(factors, axis1=1, axis2=2) <= 0
    ):
        raise ValueError("scalar factor diagnostic required")
    return factors


def fisher_from_factors(jacobian, factors):
    """Accumulate sum_q (L_q^-1 J_q).T @ (L_q^-1 J_q), reusing factors.

    Parameters
    ----------
    jacobian : array_like, shape (n_node,n_selected,n_global)
        Supplied selected mean derivatives in global parameter order. Zero
        columns are legal; all axes must be nonempty.
    factors : array_like, shape (n_node,n_selected,n_selected)
        Finite exactly lower-triangular factors with positive diagonals, normally
        from factor_covariance. Directly supplied factors receive structural
        validation only; their fiducial provenance cannot be inferred here.

    Returns
    -------
    matrix : ndarray, shape (n_global,n_global)
        Owned C-contiguous float64 data information. No inputs are changed.
        Solve/contraction workspace is bounded by one node, not all nodes.

    Raises
    ------
    ValueError
        For malformed inputs or nonfinite solve/accumulation with node context.
    """
    factors = _blocks(factors, "factors")
    jacobian = _numeric_input(jacobian, "jacobian")
    if (
        jacobian.ndim != 3
        or 0 in jacobian.shape
        or jacobian.shape[:2] != factors.shape[:2]
    ):
        raise ValueError(
            "jacobian: expected matching nonempty (n_node,n_selected,n_global)"
        )
    backend = os.environ.get("FISHHIGHZ_FISHER_BACKEND", "numpy")
    if backend not in ("numpy", "numba"):
        raise ValueError("FISHHIGHZ_FISHER_BACKEND must be numpy or numba")
    if (
        backend == "numba"
        and factors.dtype == jacobian.dtype == np.float64
        and np.geterr()["under"] == "ignore"
    ):
        try:
            from .kernels._compiled_fisher import contract
        except ImportError:
            # The NumPy-only installation remains executable even when an
            # optional backend was requested but is not installed.
            pass
        else:
            matrix, status = contract(factors, jacobian)
            if status == 0:
                return np.array(matrix, dtype=np.float64, order="C", copy=True)
            # Reproduce reference validation order and first-cell diagnostics.
    upper = np.triu_indices(factors.shape[1], 1)
    result = np.zeros((jacobian.shape[2], jacobian.shape[2]), dtype=np.float64)
    solved = np.empty(jacobian.shape[1:], dtype=np.float64)
    for node in range(len(jacobian)):
        lower = real_array(factors[node], f"factors node {node}")
        derivative = real_array(jacobian[node], f"jacobian node {node}")
        if np.any(lower[upper] != 0) or np.any(lower.diagonal() <= 0):
            raise ValueError(
                f"factors node {node}: require exact lower triangle and positive diagonal"
            )
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            _forward_substitute(lower, derivative, solved)
        if not np.all(np.isfinite(solved)):
            column = np.argwhere(~np.isfinite(solved))[0, 1]
            raise ValueError(
                f"nonfinite Fisher solve at node {node}, global parameter column {column}"
            )
        with np.errstate(over="ignore", invalid="ignore"):
            _accumulate_fisher(solved, result)
        if not np.all(np.isfinite(result)):
            raise ValueError(
                f"nonfinite Fisher accumulation at node {node}; check derivative magnitudes/units"
            )
    return result


def fisher_matrix(jacobian, covariance):
    """Factor fixed covariance once and assemble data Fisher from supplied J.

    Parameters
    ----------
    jacobian : array_like
        (node,selected,global) derivatives of the predicted mean.
    covariance : array_like
        (node,selected,selected) fixed fiducial covariance.

    Returns
    -------
    matrix : ndarray
        Owned C-contiguous float64 (global,global) Fisher matrix. See
        factor_covariance and fisher_from_factors for validation conventions.
    """
    return fisher_from_factors(jacobian, factor_covariance(covariance))
