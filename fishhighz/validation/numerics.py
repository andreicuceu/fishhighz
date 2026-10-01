"""Validation-only Wick, information summaries and literal legacy derivative.

Literal legacy algebra is labeled separately from physical field-PSD checks.
The independent checker uses direct solves rather than FishHighz factor kernels.
Peak extraction conventions are adapted from lyaforecast/fisher.py (GPLv3).
"""

import numpy as np

from .._arrays import real_array
from .._information import inspect_information
from ..fisher import factor_covariance, fisher_from_factors


def field_matrix(total, pairs, n_fields):
    """Unpack original field indices, preserving signed pair powers.

    Parameters
    ----------
    total : array_like, shape (n_cell, n_required_pair)
        Signed total field-pair powers, including noise, in (Mpc/h)^3.
    pairs : array_like, shape (n_pair, 2)
        Ordered pairs of integer field indices.
    n_fields : int
        Number of physical observed fields in the original field ordering.

    Returns
    -------
    matrix : ndarray, shape (n_cell, n_fields, n_fields)
        Symmetric total-power matrices in (Mpc/h)^3.
    """
    total = real_array(total, "total")
    matrix = np.zeros((len(total), n_fields, n_fields))
    i, j = np.asarray(pairs).T
    matrix[:, i, j] = total
    matrix[:, j, i] = total
    return matrix


def wick(total, modes, required, selected, n_fields):
    """Literal selected Wick covariance; caller separately reports field PSD.

    Parameters
    ----------
    total : array_like, shape (n_cell, n_required_pair)
        Signed total field-pair powers, including noise, in (Mpc/h)^3.
    modes : array_like, shape (n_cell,)
        Dimensionless independent-mode counts.
    required : array_like, shape (n_pair, 2)
        Ordered pairs of integer field indices.
    selected : array_like, shape (n_pair, 2)
        Ordered pairs of integer field indices.
    n_fields : int
        Number of physical observed fields in the original field ordering.

    Returns
    -------
    covariance : ndarray, shape (n_cell, n_selected_pair, n_selected_pair)
        Wick covariance in (Mpc/h)^6.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    field_power = field_matrix(total, required, n_fields)
    i, j = np.asarray(selected).T
    result = (
        field_power[:, i[:, None], i[None, :]] * field_power[:, j[:, None], j[None, :]]
        + field_power[:, i[:, None], j[None, :]]
        * field_power[:, j[:, None], i[None, :]]
    ) / np.asarray(modes)[:, None, None]
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite selected Wick covariance")
    return result


def information(fisher):
    """Rank-aware covariance/errors with explicit availability, no finite null bars.

    Parameters
    ----------
    fisher : array_like, shape (n_parameter, n_parameter)
        Fisher information in inverse products of parameter units.

    Returns
    -------
    information : dict
        Fisher matrix, covariance, marginal errors, correlations, constrained-
        coordinate mask and rank. Unavailable entries remain zero placeholders.

    Notes
    -----
    Uses the established diagonally normalized 64*eps*n symmetry/PSD/rank rule.
    Unavailable covariance/error/correlation entries are zero placeholders with
    a separate integer constrained mask. No pseudoinverse error is assigned to
    a coordinate containing a null component.
    """
    matrix, _, scale, eigenvalues, eigenvectors, rank_tolerance = inspect_information(
        fisher, "evidence Fisher"
    )
    constrained_mask = eigenvalues > rank_tolerance

    # A coordinate is estimable only when it has no component in the null space.
    null_basis = eigenvectors[:, ~constrained_mask]
    available = np.sum(null_basis**2, axis=1) <= 64 * np.finfo(float).eps * len(
        eigenvalues
    )
    covariance = np.zeros_like(matrix)
    if np.any(constrained_mask):
        basis = eigenvectors[:, constrained_mask] / scale[:, None]
        constrained_covariance = (basis / eigenvalues[constrained_mask]) @ basis.T
        mask = available[:, None] & available[None, :]
        covariance[mask] = constrained_covariance[mask]

    # Normalize only the covariance entries with available, nonzero errors.
    errors = np.sqrt(np.diag(covariance))
    correlation = np.zeros_like(matrix)
    usable = available & (errors > 0)
    usable_indices = np.ix_(usable, usable)
    correlation[usable_indices] = covariance[usable_indices] / np.outer(
        errors[usable], errors[usable]
    )
    return dict(
        fisher=matrix,
        covariance=covariance,
        errors=errors,
        correlation=correlation,
        constrained=available.astype(np.int64),
        rank=int(np.count_nonzero(constrained_mask)),
    )


def contract(covariance, observed_j, *, independent=False, batch_size=2048):
    """Combined correlated F and each independently selected spectrum's F.

    Parameters
    ----------
    covariance : ndarray, shape (n_cell, n_selected_pair, n_selected_pair)
        Wick covariance in squared power units.
    observed_j : ndarray, shape (n_cell, n_selected_pair, n_parameter)
        Observed mean-spectrum derivatives, including field responses, in power
        units per parameter unit.
    independent : bool
        Use direct Cholesky validation and solves as an independent contraction
        check. Default is ``False``.
    batch_size : int
        Maximum Fourier-cell count per contraction batch. Default is ``2048``.

    Returns
    -------
    fisher : ndarray, shape (n_parameter, n_parameter)
        Joint information from all selected spectra.
    single : ndarray, shape (n_selected_pair, n_parameter, n_parameter)
        Information from each spectrum considered independently.
    """
    # Accumulate the joint contraction and each independent spectrum together.
    n_parameters = observed_j.shape[-1]
    result = np.zeros((n_parameters, n_parameters))
    single = np.zeros((observed_j.shape[1], n_parameters, n_parameters))
    for start in range(0, len(covariance), batch_size):
        covariance_batch = covariance[start : start + batch_size]
        jacobian_batch = observed_j[start : start + batch_size]
        if independent:
            # Cholesky explicitly diagnoses SPD even if a generic solve exists.
            np.linalg.cholesky(covariance_batch)
            result += np.einsum(
                "nsi,nsj->ij",
                jacobian_batch,
                np.linalg.solve(covariance_batch, jacobian_batch),
            )
        else:
            result += fisher_from_factors(
                jacobian_batch, factor_covariance(covariance_batch)
            )
        single += np.einsum(
            "nsi,nsj,ns->sij",
            jacobian_batch,
            jacobian_batch,
            1 / np.diagonal(covariance_batch, axis1=1, axis2=2),
        )
    return result, single


def summaries(fisher, single):
    """Pack named combined and individual-spectrum derived quantities.

    Parameters
    ----------
    fisher : array_like, shape (n_parameter, n_parameter)
        Fisher information in inverse products of parameter units.
    single : array_like, shape (n_selected_pair, n_parameter, n_parameter)
        Independent-spectrum Fisher matrices.

    Returns
    -------
    summaries : dict of str to ndarray
        Joint and individual Fisher matrices, covariances, errors, correlations,
        availability masks and ranks.
    """
    main = information(fisher)
    rows = [information(f) for f in single]
    result = {k: np.asarray(v) for k, v in main.items() if k != "rank"}
    result["rank"] = np.array([main["rank"]], dtype=np.int64)
    for key in ("fisher", "covariance", "errors", "correlation", "constrained"):
        result["pair_" + key] = np.stack([r[key] for r in rows])
    result["pair_rank"] = np.array([r["rank"] for r in rows], dtype=np.int64)
    return result


def relative(a, b):
    """Zero-safe Frobenius discrepancy; exact-zero reference requires exact zero.

    Parameters
    ----------
    a : array_like
        Comparison quantity, with the same shape and units as b.
    b : array_like
        Reference quantity.

    Returns
    -------
    relative_difference : float
        Dimensionless Frobenius discrepancy relative to b; infinity denotes
        nonfinite inputs or a nonzero difference from an exactly zero reference.
    """
    a, b = np.asarray(a), np.asarray(b)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return float("inf")
    scale = max(float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    if scale == 0:
        return 0.0
    # Scale before subtraction and squaring, including subnormal operands.
    scaled_comparison, scaled_reference = a / scale, b / scale
    norm = np.linalg.norm(scaled_reference)
    if norm == 0:
        return float("inf")
    return float(np.linalg.norm(scaled_comparison - scaled_reference) / norm)


def change(a, b, pair_a, pair_b, volume_a, volume_b):
    """Information/error/individual-spectrum convergence without signed division.

    Parameters
    ----------
    a : array_like, shape (n_parameter, n_parameter)
        Fisher information in inverse products of parameter units.
    b : array_like, shape (n_parameter, n_parameter)
        Fisher information in inverse products of parameter units.
    pair_a : array_like, shape (n_pair, n_parameter, n_parameter)
        Comparison individual-spectrum Fisher matrices.
    pair_b : array_like, shape (n_pair, n_parameter, n_parameter)
        Reference individual-spectrum Fisher matrices.
    volume_a : float
        Comparison volume in (Mpc/h)^3.
    volume_b : float
        Reference volume in (Mpc/h)^3.

    Returns
    -------
    changes : dict of str to float
        Relative joint information, constrained errors, individual errors and
        volume changes.
    """

    def error_change(x, y):
        """Compare marginal errors on the same constrained parameter subspace.

        Parameters
        ----------
        x : array_like, shape (n_parameter, n_parameter)
            Fisher information in inverse products of parameter units.
        y : array_like, shape (n_parameter, n_parameter)
            Fisher information in inverse products of parameter units.

        Returns
        -------
        relative_difference : float
            Largest fractional constrained-error change, or infinity if constrained
            coordinates differ.
        """
        comparison_information, reference_information = information(x), information(y)
        if not np.array_equal(
            comparison_information["constrained"], reference_information["constrained"]
        ):
            return float("inf")
        mask = reference_information["constrained"].astype(bool)
        if not np.any(mask):
            return 0.0
        return float(
            np.max(
                abs(
                    comparison_information["errors"][mask]
                    / reference_information["errors"][mask]
                    - 1
                )
            )
        )

    return dict(
        fisher_relative=relative(a, b),
        error_relative=error_change(a, b),
        pair_error_relative=max(error_change(x, y) for x, y in zip(pair_a, pair_b)),
        volume_relative=abs(float(volume_a) / float(volume_b) - 1),
    )


def legacy_peak(model, k):
    """Literal degree-8 log-polynomial peak, 1e-12 regularizer, first-three weights.

    Parameters
    ----------
    model : array_like, shape (n_pair, n_k)
        Mean spectra including instrument response, in (Mpc/h)^3.
    k : array_like
        Comoving Fourier wavenumbers in h/Mpc; array shape follows the model or
        paired grid.

    Returns
    -------
    peak : ndarray, shape (n_pair, n_k)
        Signed residual after the literal legacy smooth fit, in the input power
        units.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    k = real_array(k, "legacy k")
    model = real_array(model, "legacy mean")
    if k.ndim != 1 or len(k) < 9 or model.shape[-1] != len(k) or np.any(k <= 0):
        raise ValueError("legacy peak requires >=9 positive k nodes and aligned mean")
    scaled_log_k = np.log(k)
    scaled_log_k = (scaled_log_k - scaled_log_k.mean()) / (
        scaled_log_k.max() - scaled_log_k.min()
    )
    weights = np.ones(len(k))
    weights[:3] *= 1e8
    return np.stack(
        [
            row
            - np.sign(row)
            * np.exp(
                np.polyval(
                    np.polyfit(scaled_log_k, np.log(abs(row) + 1e-12), 8, w=weights),
                    scaled_log_k,
                )
            )
            for row in model
        ]
    )


def legacy_jacobian(model, k, mu, *, widths):
    """Backward derivative including damping, zero first node, positive mu factors.

    Parameters
    ----------
    model : array_like, shape (n_pair, n_k)
        Mean spectra including instrument response, in (Mpc/h)^3.
    k : array_like
        Comoving Fourier wavenumbers in h/Mpc; array shape follows the model or
        paired grid.
    mu : float
        Direction cosine for this k slice.
    widths : array_like, shape (n_pair, 2)
        Fixed parallel and transverse BAO damping lengths in Mpc/h.

    Returns
    -------
    jacobian : ndarray, shape (n_k, n_pair, 2)
        Derivatives with respect to dimensionless parallel and transverse BAO
        dilations, in (Mpc/h)^3.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.

    Notes
    -----
    model is (pair,k), widths is explicit (pair,parallel/transverse). Mean already
    includes instrument response; this function never applies that response.
    """
    peak = legacy_peak(model, k)
    widths = real_array(widths, "legacy pair widths")
    if widths.shape != (len(model), 2):
        raise ValueError("legacy widths shape")
    damping = np.exp(
        -0.5
        * (
            (widths[:, 0, None] * mu * k) ** 2
            + (widths[:, 1, None] * np.sqrt(1 - mu**2) * k) ** 2
        )
    )
    peak *= damping

    # Retain the upstream backward difference and its explicitly zero first node.
    peak_derivative = np.zeros_like(peak)
    log_k_spacing = (k[1] - k[0]) / k
    peak_derivative[:, 1:] = np.diff(peak, axis=1) / log_k_spacing[1:]
    return peak_derivative.T[:, :, None] * np.array([mu**2, 1 - mu**2])
