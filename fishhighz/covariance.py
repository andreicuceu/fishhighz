"""Fixed-fiducial covariance from supplied observed field powers.

Work on one independent redshift block. Known noise belongs in total power even
when subtracted from the mean. No response, noise model, theta-dependent weight
update, or covariance-derivative information is computed here. Mode counts
include both conjugate hemispheres on mu in [0,1].
"""

import numpy as np

from ._arrays import real_array
from .fields import PairSelection
from .kernels.covariance import _gaussian_covariance_kernel


def _power_array(value, name):
    array = real_array(value, name)
    if array.ndim != 2 or 0 in array.shape:
        raise ValueError(f"{name} must have nonempty shape (n_node, n_required)")
    return array


def combine_observed_power(signal, noise):
    """Add supplied observed signal and noise without broadcasting.

    Parameters
    ----------
    signal, noise : array_like, shape (n_node, n_required)
        Real finite observed powers in (Mpc/h_fid)^3, in identical required-pair
        order. Signal already includes field responses W_i W_j. Noise includes
        any supplied off-diagonal terms. Components need exact matching shapes;
        component-specific physical validity remains the providers' responsibility.

    Returns
    -------
    total_power : ndarray
        Owned C-contiguous float64 sum. Inputs are unchanged. Physical validation
        of the total occurs in gaussian_covariance, not in this addition helper.

    Raises
    ------
    ValueError
        If types/shapes are invalid or addition gives nonfinite results.
    """
    signal = _power_array(signal, "signal")
    noise = _power_array(noise, "noise")
    if signal.shape != noise.shape:
        raise ValueError("signal and noise must have exactly matching shapes")
    with np.errstate(over="ignore", invalid="ignore"):
        total = signal + noise
    if not np.all(np.isfinite(total)):
        node, column = np.argwhere(~np.isfinite(total))[0]
        raise ValueError(
            f"nonfinite observed-power sum at node {node}, required column {column}; "
            "check component magnitudes and units"
        )
    return total


def _validate_field_power_scalar(power, selection, offset=0):
    active = np.unique(selection.selected_pairs)
    expected = np.array([(i, j) for i in active for j in active if i <= j])
    if not np.array_equal(selection.required_pairs, expected):
        raise ValueError(
            "required_pairs must contain every canonical pair of active fields"
        )
    row, column = np.searchsorted(active, selection.required_pairs).T
    ids = tuple(selection.fields[i].id for i in active)
    matrix = np.empty((len(active), len(active)), dtype=np.float64)
    eps = np.finfo(np.float64).eps
    for node, packed in enumerate(power, start=offset):
        matrix[row, column] = packed
        matrix[column, row] = packed
        auto = matrix.diagonal()
        negative = np.flatnonzero(auto < 0)
        if negative.size:
            i = negative[0]
            raise ValueError(
                f"node {node}, field {ids[i]!r}: negative auto power {auto[i]:.17g}"
            )
        zero = auto == 0
        if np.any(matrix[zero] != 0):
            i = np.flatnonzero(zero & np.any(matrix != 0, axis=1))[0]
            raise ValueError(
                f"node {node}, field {ids[i]!r}: zero auto power requires an exactly "
                f"zero row; max absolute cross power {np.max(np.abs(matrix[i])):.17g}"
            )
        positive = np.flatnonzero(~zero)
        if not positive.size:
            continue
        corr = matrix[np.ix_(positive, positive)]
        scale = np.sqrt(auto[positive])
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            corr /= scale[:, None]
            corr /= scale[None, :]
        # Sequential divisions avoid overflow/underflow of products of auto powers.
        n = len(positive)
        if not np.all(np.isfinite(corr)) or np.max(np.abs(corr)) > 1 + 64 * eps * n:
            raise ValueError(
                f"node {node}, active fields {ids!r}: invalid normalized correlation; "
                f"max absolute value {np.max(np.abs(corr)):.17g} (expected <= 1 "
                "within roundoff); check cross/auto powers"
            )
        try:
            eigenvalues = np.linalg.eigvalsh(corr)
        except np.linalg.LinAlgError as error:
            raise ValueError(
                f"node {node}, active fields {ids!r}: correlation eigensolve failed"
            ) from error
        tolerance = 64 * eps * n * max(1.0, np.max(np.abs(eigenvalues)))
        if eigenvalues[0] < -tolerance:
            raise ValueError(
                f"node {node}, active fields {ids!r}: total field power is not PSD; "
                f"normalized minimum eigenvalue {eigenvalues[0]:.17g}, "
                f"allowed negative tolerance {tolerance:.17g}"
            )


def _validate_field_power(power, selection):
    active = np.unique(selection.selected_pairs)
    expected = np.array([(i, j) for i in active for j in active if i <= j])
    if not np.array_equal(selection.required_pairs, expected):
        raise ValueError(
            "required_pairs must contain every canonical pair of active fields"
        )
    row, column = np.searchsorted(active, selection.required_pairs).T
    n = len(active)
    eps_n = 64 * np.finfo(np.float64).eps * n
    for start in range(0, len(power), 256):
        packed = power[start : start + 256]
        matrix = np.empty((len(packed), n, n), dtype=np.float64)
        matrix[:, row, column] = packed
        matrix[:, column, row] = packed
        auto = np.diagonal(matrix, axis1=1, axis2=2)
        try:
            if np.any(auto <= 0):
                raise ValueError("scalar zero/negative diagonal diagnostic required")
            scale = np.sqrt(auto)
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                corr = matrix / scale[:, :, None] / scale[:, None, :]
            if not np.all(np.isfinite(corr)) or np.any(np.abs(corr) > 1 + eps_n):
                raise ValueError("scalar correlation diagnostic required")
            values = np.linalg.eigvalsh(corr)
            tolerance = eps_n * np.maximum(1.0, np.max(np.abs(values), axis=1))
            if not np.all(np.isfinite(values)) or np.any(values[:, 0] <= 2 * tolerance):
                raise ValueError("scalar PSD diagnostic required")
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            _validate_field_power_scalar(packed, selection, offset=start)


def gaussian_covariance(total_power, mode_counts, selection):
    """Compute selected-spectrum Gaussian covariance at each independent node.

    Parameters
    ----------
    total_power : array_like, shape (n_node, n_required)
        Supplied fiducial observed clustering plus noise, in
        selection.required_pairs order, in (Mpc/h_fid)^3.
    mode_counts : array_like, shape (n_node,)
        Positive finite counts, e.g. V_fid * grid.q_mode. Fractional counts are
        allowed. Counts include both conjugate hemispheres; divide exactly once.
    selection : PairSelection
        Prepared field identities, dependencies, and selected-pair order.

    Returns
    -------
    covariance : ndarray, shape (n_node, n_selected, n_selected)
        Owned C-contiguous float64 blocks in selected-pair order, with squared
        power units. Inputs remain unchanged. Valid singular blocks are allowed;
        no invertibility claim, regularization, or solve is made.

    Notes
    -----
    C_AB = (T_im*T_jn + T_in*T_jm)/mode_counts for A=(i,j), B=(m,n).
    Nonnegative autos are required exactly; a zero auto requires its entire row
    to be exactly zero. Each active field matrix is checked separately. Positive
    autos normalize the matrix to correlations R. The dimensionless eigenvalue
    tolerance is 64*eps64*n_positive*max(1, max(abs(eigvalsh(R)))); only negative
    eigenvalues within this roundoff tolerance pass. Correlations also satisfy
    abs(R_ij) <= 1 + 64*eps64*n_positive. No clipping or power rescaling is applied
    to the supplied matrix or output. Unused fields are excluded, and the full
    active pair closure must be complete. Covariance eigensolves are not repeated.

    Raises
    ------
    ValueError
        For invalid shape/type, nonfinite data/arithmetic, nonpositive counts,
        incomplete dependencies, or unphysical total power (with node/field
        diagnostics). Covariance arithmetic overflow is an error even if a
        differently ordered calculation could avoid an intermediate overflow.
    """
    if not isinstance(selection, PairSelection):
        raise ValueError("selection must be a prepared PairSelection")
    power = _power_array(total_power, "total_power")
    if power.shape[1] != len(selection.required_pairs):
        raise ValueError(
            f"total_power needs exactly {len(selection.required_pairs)} required columns"
        )
    counts = real_array(mode_counts, "mode_counts")
    if counts.shape != (power.shape[0],) or np.any(counts <= 0):
        raise ValueError(
            "mode_counts must have shape (n_node,) and be strictly positive"
        )
    _validate_field_power(power, selection)
    n_selected = len(selection.selected_pairs)
    out = np.empty((len(power), n_selected, n_selected), dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        _gaussian_covariance_kernel(
            power, counts, selection.im, selection.jn, selection.in_, selection.jm, out
        )
    finite = np.all(np.isfinite(out), axis=(1, 2))
    if not np.all(finite):
        node = np.flatnonzero(~finite)[0]
        raise ValueError(
            f"nonfinite covariance arithmetic at node {node}; "
            "check power magnitudes, mode counts, and units"
        )
    return out
