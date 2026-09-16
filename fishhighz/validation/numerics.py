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
    """Unpack original field indices, preserving signed pair powers."""
    total = real_array(total, "total")
    matrix = np.zeros((len(total), n_fields, n_fields))
    i, j = np.asarray(pairs).T
    matrix[:, i, j] = total
    matrix[:, j, i] = total
    return matrix


def wick(total, modes, required, selected, n_fields):
    """Literal selected Wick covariance; caller separately reports field PSD."""
    t = field_matrix(total, required, n_fields)
    i, j = np.asarray(selected).T
    result = (
        t[:, i[:, None], i[None, :]] * t[:, j[:, None], j[None, :]]
        + t[:, i[:, None], j[None, :]] * t[:, j[:, None], i[None, :]]
    ) / np.asarray(modes)[:, None, None]
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite selected Wick covariance")
    return result


def information(fisher):
    """Rank-aware covariance/errors with explicit availability, no finite null bars.

    Uses the established diagonally normalized 64*eps*n symmetry/PSD/rank rule.
    Unavailable covariance/error/correlation entries are zero placeholders with
    a separate integer constrained mask. No pseudoinverse error is assigned to
    a coordinate containing a null component.
    """
    matrix, _, scale, eig, vec, tol = inspect_information(fisher, "evidence Fisher")
    good = eig > tol
    null = vec[:, ~good]
    available = np.sum(null**2, axis=1) <= 64 * np.finfo(float).eps * len(eig)
    covariance = np.zeros_like(matrix)
    if np.any(good):
        basis = vec[:, good] / scale[:, None]
        sub = (basis / eig[good]) @ basis.T
        mask = available[:, None] & available[None, :]
        covariance[mask] = sub[mask]
    errors = np.sqrt(np.diag(covariance))
    correlation = np.zeros_like(matrix)
    usable = available & (errors > 0)
    ix = np.ix_(usable, usable)
    correlation[ix] = covariance[ix] / np.outer(errors[usable], errors[usable])
    return dict(
        fisher=matrix,
        covariance=covariance,
        errors=errors,
        correlation=correlation,
        constrained=available.astype(np.int64),
        rank=int(np.count_nonzero(good)),
    )


def contract(covariance, observed_j, *, independent=False, batch_size=2048):
    """Combined correlated F and each independently selected spectrum's F."""
    npar = observed_j.shape[-1]
    result = np.zeros((npar, npar))
    single = np.zeros((observed_j.shape[1], npar, npar))
    for start in range(0, len(covariance), batch_size):
        c = covariance[start : start + batch_size]
        j = observed_j[start : start + batch_size]
        if independent:
            # Cholesky explicitly diagnoses SPD even if a generic solve exists.
            np.linalg.cholesky(c)
            result += np.einsum("nsi,nsj->ij", j, np.linalg.solve(c, j))
        else:
            result += fisher_from_factors(j, factor_covariance(c))
        single += np.einsum(
            "nsi,nsj,ns->sij", j, j, 1 / np.diagonal(c, axis1=1, axis2=2)
        )
    return result, single


def summaries(fisher, single):
    """Pack named combined and individual-spectrum derived quantities."""
    main = information(fisher)
    rows = [information(f) for f in single]
    result = {k: np.asarray(v) for k, v in main.items() if k != "rank"}
    result["rank"] = np.array([main["rank"]], dtype=np.int64)
    for key in ("fisher", "covariance", "errors", "correlation", "constrained"):
        result["pair_" + key] = np.stack([r[key] for r in rows])
    result["pair_rank"] = np.array([r["rank"] for r in rows], dtype=np.int64)
    return result


def relative(a, b):
    """Zero-safe Frobenius discrepancy; exact-zero reference requires exact zero."""
    a, b = np.asarray(a), np.asarray(b)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return float("inf")
    scale = max(float(np.max(np.abs(a))), float(np.max(np.abs(b))))
    if scale == 0:
        return 0.0
    # Scale before subtraction and squaring, including subnormal operands.
    x, y = a / scale, b / scale
    norm = np.linalg.norm(y)
    if norm == 0:
        return float("inf")
    return float(np.linalg.norm(x - y) / norm)


def change(a, b, pair_a, pair_b, volume_a, volume_b):
    """Information/error/individual-spectrum convergence without signed division."""

    def error_change(x, y):
        sx, sy = information(x), information(y)
        if not np.array_equal(sx["constrained"], sy["constrained"]):
            return float("inf")
        mask = sy["constrained"].astype(bool)
        if not np.any(mask):
            return 0.0
        return float(np.max(abs(sx["errors"][mask] / sy["errors"][mask] - 1)))

    return dict(
        fisher_relative=relative(a, b),
        error_relative=error_change(a, b),
        pair_error_relative=max(error_change(x, y) for x, y in zip(pair_a, pair_b)),
        volume_relative=abs(float(volume_a) / float(volume_b) - 1),
    )


def legacy_peak(model, k):
    """Literal degree-8 log-polynomial peak, 1e-12 regularizer, first-three weights."""
    k = real_array(k, "legacy k")
    model = real_array(model, "legacy mean")
    if k.ndim != 1 or len(k) < 9 or model.shape[-1] != len(k) or np.any(k <= 0):
        raise ValueError("legacy peak requires >=9 positive k nodes and aligned mean")
    x = np.log(k)
    x = (x - x.mean()) / (x.max() - x.min())
    weights = np.ones(len(k))
    weights[:3] *= 1e8
    return np.stack(
        [
            row
            - np.sign(row)
            * np.exp(
                np.polyval(np.polyfit(x, np.log(abs(row) + 1e-12), 8, w=weights), x)
            )
            for row in model
        ]
    )


def legacy_jacobian(model, k, mu, *, widths):
    """Backward derivative including damping, zero first node, positive mu factors.

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
    deriv = np.zeros_like(peak)
    dlogk = (k[1] - k[0]) / k
    deriv[:, 1:] = np.diff(peak, axis=1) / dlogk[1:]
    return deriv.T[:, :, None] * np.array([mu**2, 1 - mu**2])
