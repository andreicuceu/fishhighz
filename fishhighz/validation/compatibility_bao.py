"""Saved-input BAO contractions with explicit forest-auto noise dependencies."""

import numpy as np

from .compatibility_weights import forest_noise
from .numerics import contract, information, wick


def dependencies(task, pair):
    """Auto contexts needed by the individual spectrum's Wick variance.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    pair : sequence of int
        One pair of field indices.

    Returns
    -------
    contexts : list of str
        Sorted unique forest-auto context names required by the spectrum.
    """
    return sorted(
        {
            task["fields"][i]["id"] + "_" + task["fields"][i]["id"]
            for i in pair
            if task["fields"][i]["kind"] == "forest"
        }
    )


def plot_errors(errors):
    """Mask both components together; NaNs preserve gaps along redshift.

    Parameters
    ----------
    errors : array_like, shape (..., 2)
        Dimensionless parallel and transverse BAO uncertainties.

    Returns
    -------
    masked_errors : ndarray, shape (..., 2)
        Copy with both BAO-error components replaced by NaN if either is
        nonfinite or exceeds 0.2.
    """
    result = np.asarray(errors, dtype=float).copy()
    usable = np.all(np.isfinite(result) & (result <= 0.2), axis=-1)
    result[~usable] = np.nan
    return result


def plot_ratio(numerator, denominator):
    """Fractional changes require both operands to survive the joint cut.

    Parameters
    ----------
    numerator : array_like, shape (..., 2)
        Comparison parallel/transverse BAO uncertainties.
    denominator : array_like, shape (..., 2)
        Reference parallel/transverse BAO uncertainties.

    Returns
    -------
    fractional_change : ndarray, shape (..., 2)
        Dimensionless numerator/denominator minus one, retaining paired-
        component masking.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        return plot_errors(numerator) / plot_errors(denominator) - 1


def forecast(arrays, task, rows, auto_coefficients):
    """Reassemble only available auto noise and contract valid observables.

    Parameters
    ----------
    arrays : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    rows : dict
        Captured source and auxiliary-power inputs keyed by field-pair name.
    auto_coefficients : dict of str to array_like
        Available forest-auto coefficient vectors; final two entries are A and
        pixel power.

    Returns
    -------
    singles : list of dict or None
        Rank-aware information for each valid individual spectrum.
    joint : dict or None
        Joint information, available only when every required auto is present.
    total : ndarray, shape (n_cell, n_required_pair)
        Reconstructed total powers in (Mpc/h)^3.
    covariance : ndarray, shape (n_cell, n_selected_pair, n_selected_pair)
        Reconstructed Wick covariance in (Mpc/h)^6.
    valid : ndarray of bool, shape (n_selected_pair,)
        Availability of each spectrum under the supplied auto coefficients.

    Notes
    -----
    Missing autos use zero scratch entries solely to form unaffected variances.
    No spectrum depending on a missing auto, nor a partial joint forecast, is
    returned. Intrinsic covariance-redshift powers and observed mean derivatives
    are retained exactly. No priors are added to the saved two-parameter recipe.
    """
    total = arrays["total"].copy()
    missing = set()
    for column, (i, jacobian_batch) in enumerate(task["required_pairs"]):
        if i != jacobian_batch or task["fields"][i]["kind"] != "forest":
            continue
        name = dependencies(task, [i])[0]
        if name not in auto_coefficients:
            missing.add(name)
            total[:, column] = 0
            continue
        row = rows[name]
        baseline = forest_noise(
            row,
            arrays["k"],
            arrays["mu"],
            row["_aliasing_weights"][-1],
            row["_effective_noise_power"][-1],
        )
        changed = forest_noise(
            row, arrays["k"], arrays["mu"], *auto_coefficients[name][-2:]
        )
        total[:, column] = arrays["total"][:, column] - baseline + changed

    # Missing forest autos invalidate every selected spectrum that needs them.
    valid_mask = np.array(
        [
            not missing.intersection(dependencies(task, p))
            for p in task["selected_pairs"]
        ]
    )
    covariance = wick(
        total,
        arrays["modes"],
        task["required_pairs"],
        task["selected_pairs"],
        len(task["fields"]),
    )

    singles = [None] * len(valid_mask)
    # Same batching and reduction order as the immutable compatibility baseline.
    single = np.zeros((len(valid_mask), 2, 2))
    for start in range(0, len(covariance), 2048):
        jacobian_batch = arrays["observed_j"][start : start + 2048, valid_mask]
        diagonal = np.diagonal(covariance[start : start + 2048], axis1=1, axis2=2)[
            :, valid_mask
        ]
        single[valid_mask] += np.einsum(
            "nsi,nsj,ns->sij", jacobian_batch, jacobian_batch, 1 / diagonal
        )
    for index in np.flatnonzero(valid_mask):
        singles[index] = information(single[index])

    joint = (
        information(contract(covariance, arrays["observed_j"])[0])
        if np.all(valid_mask)
        else None
    )
    return singles, joint, total, covariance, valid_mask
