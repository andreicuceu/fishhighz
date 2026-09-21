"""Saved-input BAO contractions with explicit forest-auto noise dependencies."""

import numpy as np

from .compatibility_weights import forest_noise
from .numerics import contract, information, wick


def dependencies(task, pair):
    """Auto contexts needed by the individual spectrum's Wick variance."""
    return sorted(
        {
            task["fields"][i]["id"] + "_" + task["fields"][i]["id"]
            for i in pair
            if task["fields"][i]["kind"] == "forest"
        }
    )


def plot_errors(errors):
    """Mask both components together; NaNs preserve gaps along redshift."""
    result = np.asarray(errors, dtype=float).copy()
    usable = np.all(np.isfinite(result) & (result <= 0.2), axis=-1)
    result[~usable] = np.nan
    return result


def plot_ratio(numerator, denominator):
    """Fractional changes require both operands to survive the joint cut."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return plot_errors(numerator) / plot_errors(denominator) - 1


def forecast(arrays, task, rows, auto_coefficients):
    """Reassemble only available auto noise and contract valid observables.

    Missing autos use zero scratch entries solely to form unaffected variances.
    No spectrum depending on a missing auto, nor a partial joint forecast, is
    returned. Intrinsic covariance-redshift powers and observed mean derivatives
    are retained exactly. No priors are added to the saved two-parameter recipe.
    """
    total = arrays["total"].copy()
    missing = set()
    for column, (i, j) in enumerate(task["required_pairs"]):
        if i != j or task["fields"][i]["kind"] != "forest":
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
    valid = np.array(
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
    singles = [None] * len(valid)
    # Same batching and reduction order as the immutable compatibility baseline.
    single = np.zeros((len(valid), 2, 2))
    for start in range(0, len(covariance), 2048):
        j = arrays["observed_j"][start : start + 2048, valid]
        diagonal = np.diagonal(covariance[start : start + 2048], axis1=1, axis2=2)[
            :, valid
        ]
        single[valid] += np.einsum("nsi,nsj,ns->sij", j, j, 1 / diagonal)
    for index in np.flatnonzero(valid):
        singles[index] = information(single[index])
    joint = (
        information(contract(covariance, arrays["observed_j"])[0])
        if np.all(valid)
        else None
    )
    return singles, joint, total, covariance, valid
