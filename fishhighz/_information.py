"""Shared matrix-boundary checks in diagonally normalized coordinates.

Private to Fisher assembly/results; does not alter the Step 04 PSD contract.
"""

import numpy as np

from ._arrays import real_array


def inspect_information(value, context):
    """Return symmetric copy, normalized matrix, scales, eigensystem, tolerance.

    Negative diagonals fail exactly, and zero diagonals require exact zero rows
    and columns. scales=sqrt(diag), with scale=1 for zero rows. R=M/scales/scales.
    Symmetry uses elementwise 64*eps*n*max(1,abs(Rij),abs(Rji)). Accepted asymmetry
    is averaged only in returned copies. PSD/rank tolerance is
    64*eps*n*max(1,max(abs(eigenvalues))). No eigenvalues are modified.
    """
    matrix = real_array(value, context)
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{context}: expected nonempty square matrix")
    diagonal = matrix.diagonal()
    if np.any(diagonal < 0):
        index = np.flatnonzero(diagonal < 0)[0]
        raise ValueError(
            f"{context}: negative diagonal at index {index}: {diagonal[index]:.17g}"
        )
    zero = diagonal == 0
    if np.any(matrix[zero] != 0) or np.any(matrix[:, zero] != 0):
        raise ValueError(f"{context}: zero diagonal requires exact zero row and column")
    scales = np.sqrt(diagonal)
    scales[zero] = 1.0
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        normalized = matrix / scales[:, None] / scales[None, :]
    if not np.all(np.isfinite(normalized)):
        raise ValueError(f"{context}: nonfinite diagonally normalized matrix")
    eps_n = 64 * np.finfo(np.float64).eps * len(matrix)
    symmetry_tolerance = eps_n * np.maximum(
        1.0, np.maximum(np.abs(normalized), np.abs(normalized.T))
    )
    with np.errstate(over="ignore", invalid="ignore"):
        asymmetric = np.abs(normalized - normalized.T) > symmetry_tolerance
    if np.any(asymmetric):
        i, j = np.argwhere(asymmetric)[0]
        raise ValueError(
            f"{context}: asymmetric at ({i},{j}) in normalized basis; "
            f"entries {normalized[i, j]:.17g}, {normalized[j, i]:.17g}"
        )
    # Average off-diagonals only; preserve tiny diagonal values without halving.
    i, j = np.triu_indices(len(matrix), 1)
    different = matrix[i, j] != matrix[j, i]
    a, b = i[different], j[different]
    # Leave exactly symmetric entries untouched, including subnormal powers.
    matrix[a, b] = matrix[a, b] + 0.5 * (matrix[b, a] - matrix[a, b])
    matrix[b, a] = matrix[a, b]
    normalized[i, j] = 0.5 * normalized[i, j] + 0.5 * normalized[j, i]
    normalized[j, i] = normalized[i, j]
    try:
        values, vectors = np.linalg.eigh(normalized)
    except np.linalg.LinAlgError as error:
        raise ValueError(f"{context}: normalized eigensolve failed") from error
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{context}: nonfinite normalized eigenvalues")
    tolerance = eps_n * max(1.0, np.max(np.abs(values)))
    if values[0] < -tolerance:
        raise ValueError(
            f"{context}: indefinite; normalized minimum eigenvalue "
            f"{values[0]:.17g}, tolerance {tolerance:.17g}"
        )
    return matrix, normalized, scales, values, vectors, tolerance


def normalize_positive_batch(value):
    """Normalize a positive-diagonal batch; caller replays scalar failures.

    This private fast-path helper has no acceptance authority near boundaries.
    It preserves sequential divisions and off-diagonal averaging from
    inspect_information, but does not compute unused eigenvectors.
    """
    matrix = np.asarray(value, dtype=np.float64)
    diagonal = np.diagonal(matrix, axis1=1, axis2=2)
    if not np.all(np.isfinite(matrix)) or np.any(diagonal <= 0):
        raise ValueError("scalar diagonal diagnostic required")
    scales = np.sqrt(diagonal)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        normalized = matrix / scales[:, :, None] / scales[:, None, :]
        transpose = normalized.swapaxes(1, 2)
        tolerance = (
            64
            * np.finfo(np.float64).eps
            * matrix.shape[-1]
            * np.maximum(1.0, np.maximum(np.abs(normalized), np.abs(transpose)))
        )
        asymmetric = np.abs(normalized - transpose) > tolerance
    if not np.all(np.isfinite(normalized)) or np.any(asymmetric):
        raise ValueError("scalar symmetry diagnostic required")
    i, j = np.triu_indices(matrix.shape[-1], 1)
    normalized[:, i, j] = 0.5 * normalized[:, i, j] + 0.5 * normalized[:, j, i]
    normalized[:, j, i] = normalized[:, i, j]
    return normalized, scales
