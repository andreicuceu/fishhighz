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
    """Validate real numeric dtype without copying every Fourier node.

    Parameters
    ----------
    value : array_like
        Real numeric input of arbitrary shape and units.
    name : str
        Quantity name used in errors.

    Returns
    -------
    array : ndarray
        Input view when possible, preserving dtype, shape, and units.

    Raises
    ------
    ValueError
        If the input is not real integer or floating-point data.

    Notes
    -----
    Float64 conversion and finite checks are deferred to the node evaluation.
    """
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError(f"{name}: expected real numeric data")
    return array


def _blocks(value, name):
    """Validate the shape of a batch of square numerical matrices.

    Parameters
    ----------
    value : array_like of shape (n_node, n_selected, n_selected)
        Real numeric blocks; units are retained.
    name : str
        Quantity name used in errors.

    Returns
    -------
    array : ndarray
        Input array without an unconditional float64 copy.

    Raises
    ------
    ValueError
        If dtype or nonempty square-block shape is invalid.
    """
    array = _numeric_input(value, name)
    if array.ndim != 3 or 0 in array.shape or array.shape[1] != array.shape[2]:
        raise ValueError(f"{name}: expected nonempty (n_node,n_selected,n_selected)")
    return array


def _factor_covariance_scalar(covariance, offset=0):
    """Factor covariance blocks with scalar normalized-rank diagnostics.

    Parameters
    ----------
    covariance : array_like of shape (n_node, n_selected, n_selected)
        Fixed covariance blocks in (Mpc/h_fid)^6, in selected-spectrum order.
    offset : int, default=0
        Global Fourier-node index of the first block, used in errors.

    Returns
    -------
    factors : ndarray of shape (n_node, n_selected, n_selected)
        Owned float64 lower Cholesky factors in (Mpc/h_fid)^3.

    Raises
    ------
    ValueError
        If a block is invalid, rank deficient, or cannot be factored without
        regularization.

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

    Parameters
    ----------
    covariance : array_like of shape (n_node, n_selected, n_selected)
        Fixed covariance blocks in (Mpc/h_fid)^6, in selected-spectrum order.

    Returns
    -------
    factors : ndarray of shape (n_node, n_selected, n_selected)
        Owned float64 lower Cholesky factors in (Mpc/h_fid)^3.

    Raises
    ------
    ValueError
        If a covariance block is invalid, singular, or numerically unresolved.

    Notes
    -----
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
    """Factor a covariance batch sufficiently far from rank boundaries.

    Parameters
    ----------
    block : ndarray of shape (n_node, n_selected, n_selected)
        Covariance blocks in (Mpc/h_fid)^6.

    Returns
    -------
    factors : ndarray of shape (n_node, n_selected, n_selected)
        Owned float64 lower Cholesky factors in (Mpc/h_fid)^3.

    Raises
    ------
    ValueError
        If scalar diagonal, symmetry, rank, or factor diagnostics are required.
    numpy.linalg.LinAlgError
        If the eigensolve or Cholesky factorization fails.

    Notes
    -----
    The doubled rank tolerance selects the fast calculation only. A rejected
    batch is replayed through scalar diagnostics with the original acceptance
    threshold; no eigenvalue is modified.
    """
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


def _fisher_one_spectrum(jacobian, factors):
    """Vectorized one-spectrum (1x1 factor) form of the reference node loop.

    Parameters
    ----------
    jacobian : array_like of shape (n_node, 1, n_global)
        Mean derivatives in (Mpc/h_fid)^3 per global parameter unit.
    factors : ndarray of shape (n_node, 1, 1)
        Owned float64 lower Cholesky factors in (Mpc/h_fid)^3.

    Returns
    -------
    matrix : ndarray of shape (n_global, n_global) or None
        Float64 Fisher information in inverse products of global parameter
        units.
        Return None if any input or arithmetic check fails, requesting replay
        through the reference loop.

    Notes
    -----
    The one-row forward substitution is a single division and np.add.accumulate
    sums node contributions in the loop's order, so the result is identical to
    the loop. Workspace is O(n_node*n_global**2) for the supplied node batch.
    Returns None whenever any check or value would fail, so the caller replays
    the reference loop and its first-node diagnostics.
    """
    lower = factors[:, 0, 0]
    if (
        lower.dtype.kind != "f"
        or jacobian.dtype.kind != "f"
        or not np.all(np.isfinite(lower))
        or np.any(lower <= 0)
        or not np.all(np.isfinite(jacobian))
    ):
        return None
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        solved = (
            jacobian[:, 0, :].astype(np.float64) / lower.astype(np.float64)[:, None]
        )
        partial = np.add.accumulate(solved[:, :, None] * solved[:, None, :], axis=0)
    if not np.all(np.isfinite(solved)) or not np.all(np.isfinite(partial)):
        return None
    return np.array(partial[-1], dtype=np.float64, order="C", copy=True)


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
    # Unset selects the compiled contraction when Numba is installed; "numpy"
    # forces the reference loop below.
    backend = os.environ.get("FISHHIGHZ_FISHER_BACKEND", "numba")
    if backend not in ("numpy", "numba"):
        raise ValueError("FISHHIGHZ_FISHER_BACKEND must be numpy or numba")
    if (
        backend == "numba"
        and factors.dtype == jacobian.dtype == np.float64
        and np.geterr()["under"] == "ignore"
    ):
        try:
            # Numba matrix products call SciPy's BLAS; without it the compiled
            # kernel aborts the interpreter instead of raising.
            import scipy.linalg.cython_blas  # noqa: F401

            from .kernels._compiled_fisher import contract
        except ImportError:
            # The NumPy-only installation remains executable even when an
            # optional backend was requested but is not installed.
            pass
        else:
            # One writable C-contiguous signature avoids recompiling for sliced
            # or read-only (prepared, immutable) inputs.
            matrix, status = contract(
                np.require(factors, requirements=("C", "W")),
                np.require(jacobian, requirements=("C", "W")),
            )
            if status == 0:
                return np.array(matrix, dtype=np.float64, order="C", copy=True)
            # Reproduce reference validation order and first-cell diagnostics.
    if factors.shape[1] == 1 and np.geterr()["under"] == "ignore":
        matrix = _fisher_one_spectrum(jacobian, factors)
        if matrix is not None:
            return matrix
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
    """Factor the fixed covariance and contract the mean-power derivatives.

    Parameters
    ----------
    jacobian : array_like of shape (n_node, n_selected, n_global)
        Mean derivatives in (Mpc/h_fid)^3 per global parameter unit.
    covariance : array_like of shape (n_node, n_selected, n_selected)
        Fixed covariance blocks in (Mpc/h_fid)^6, in selected-spectrum order.

    Returns
    -------
    matrix : ndarray of shape (n_global, n_global)
        Float64 Fisher information in inverse products of global parameter
        units.

    Raises
    ------
    ValueError
        If covariance or Jacobian validation or numerical contraction fails.

    Notes
    -----
    No covariance derivatives, additional mode counts, or quadrature factors
    are added. See factor_covariance and fisher_from_factors for validation.
    """
    return fisher_from_factors(jacobian, factor_covariance(covariance))
