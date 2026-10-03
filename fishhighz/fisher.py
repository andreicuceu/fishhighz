"""Fisher information from supplied means and fixed fiducial covariance.

No additional volume, mode count, quadrature, or factor of two is applied.
Known subtracted noise affects covariance; its parameter response, if any, must
be supplied in the mean Jacobian. No covariance derivative information is added.
"""

import numpy as np

from ._information import inspect_information
from .kernels.fisher import _cholesky

# Diagnostics for the nonzero status codes of kernels._compiled_fisher.contract,
# formatted with the first failing Fourier node and global parameter column.
_CONTRACTION_ERRORS = {
    1: "factors node {node} must be finite",
    2: "jacobian node {node} must be finite",
    3: "factors node {node}: require exact lower triangle and positive diagonal",
    4: "nonfinite Fisher solve at node {node}, global parameter column {column}",
    5: (
        "nonfinite Fisher accumulation at node {node}; "
        "check derivative magnitudes/units"
    ),
}


def _numeric_input(value, name):
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
    Existing array views are preserved. Float64 conversion (only for
    non-float64 data) and finite checks are deferred to the Fisher contraction.
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
        Per node, the order of checks is: finite factors, finite Jacobian,
        exact lower triangle with positive diagonal, finite solve (reporting
        the first nonfinite global parameter column), finite accumulation.

    Notes
    -----
    The contraction is the compiled Numba kernel
    ``kernels._compiled_fisher.contract`` (no fastmath, no parallel loops).
    Integer inputs are converted to float64; float64 inputs are not copied.
    Compiled arithmetic does not trap floating-point underflow, whatever the
    NumPy ``errstate``: gradual underflow in the solve or accumulation is
    accepted silently, while nonfinite results raise as described above.
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
    # The compiled contraction has one float64 signature; integer inputs are
    # converted, while float64 inputs are passed on without a copy.
    if factors.dtype != np.float64:
        factors = np.asarray(factors, dtype=np.float64)
    if jacobian.dtype != np.float64:
        jacobian = np.asarray(jacobian, dtype=np.float64)

    # Lazy import keeps Numba compilation out of package import. One writable
    # C-contiguous signature avoids recompiling for sliced or read-only
    # (prepared, immutable) inputs.
    from .kernels._compiled_fisher import contract

    matrix, status, node, column = contract(
        np.require(factors, requirements=("C", "W")),
        np.require(jacobian, requirements=("C", "W")),
    )
    if status != 0:
        raise ValueError(_CONTRACTION_ERRORS[status].format(node=node, column=column))
    return np.array(matrix, dtype=np.float64, order="C", copy=True)


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
