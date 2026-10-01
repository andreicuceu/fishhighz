"""Literal upstream legacy estimator, independent FishHighz Wick/Fisher assembly."""

import json
from pathlib import Path

import numpy as np

from .evidence import digest
from .numerics import change, legacy_jacobian, relative
from .schema import assemble


def load_reference(root, case, index):
    """Validate captured case/bin/field identities and hash before reading upstream.

    Parameters
    ----------
    root : str or pathlib.Path
        Root directory of the input checkout or saved evidence bundle.
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.
    index : int
        Zero-based redshift-bin index.

    Returns
    -------
    arrays : dict of str to ndarray
        Exact captured upstream arrays for this bin.
    record : dict
        Bound bin identity, physical inputs and derivative conventions.
    metadata : dict
        Case-wide reference provenance and capture inventory.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    from .cases import bins, selection

    root = Path(root).resolve()
    case_root = root / case
    manifest = json.loads((root / "manifest.json").read_text())
    if not manifest["complete"] or manifest["requested"] != [
        r["case"] for r in manifest["records"]
    ]:
        raise ValueError("incomplete reference inventory")
    entry = next(r for r in manifest["records"] if r["case"] == case)
    metadata_path = root / entry["metadata"]
    if (
        metadata_path.resolve() != case_root / "metadata.json"
        or digest(metadata_path) != entry["sha256"]
    ):
        raise ValueError("reference metadata identity")
    meta = json.loads(metadata_path.read_text())
    row = meta["records"][index]
    if (
        row["case"] != case
        or row["bin"] != index
        or row["bounds"] != list(bins(case)[index])
        or row["fields"] != [f.id for f in selection(case).fields]
    ):
        raise ValueError("swapped reference case/bin/fields")
    path = case_root / row["array"]
    if path.resolve().parent != case_root or digest(path) != row["sha256"]:
        raise ValueError("reference array identity")
    with np.load(path, allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    return arrays, row, meta


def run(task, reference, provenance):
    """Never use upstream F as FishHighz information; use it solely as comparator.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    reference : str or pathlib.Path
        Root of the reference checkout or captured bundle.
    provenance : dict
        Input paths, hashes and software identities retained with the
        calculation.

    Returns
    -------
    arrays : dict of str to ndarray
        Independently assembled FishHighz information and retained upstream
        comparison arrays.
    report : dict
        Literal-estimator agreement, provenance and field-PSD diagnostics.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    reference_arrays, row, meta = load_reference(reference, task["case"], task["bin"])
    k_grid = reference_arrays["k_axis"]
    k_spacing = k_grid[1] - k_grid[0]
    volume = float(
        reference_arrays["modes"][0] * 2 * np.pi**2 / (k_grid[0] ** 2 * k_spacing * 0.1)
    )
    settings = dict(
        profile="compatibility",
        parameters=task["parameters"],
        bounds=task["bounds"],
        fields=row["fields"],
        grid=dict(kind="legacy", volume=volume),
        mean_z=row["mean_z"],
        covariance_z=row["covariance_z"],
        response_ownership=task["response_ownership"],
        pair_inputs=row["pair_inputs"],
        estimator="literal legacy degree-8 peak, backward dlogk including damping, first node zero",
        numerical_limits="algebraic reproduction; independent field PSD diagnostics, selected C must be SPD",
    )
    arrays, report = assemble(
        task,
        reference_arrays["total"],
        reference_arrays["observed_j"],
        settings,
        extra=dict(
            reference_fisher=reference_arrays["legacy_fisher"],
            reference_pair_fisher=reference_arrays["legacy_pair_fisher"],
            mean=reference_arrays["mean"],
        ),
    )
    # Independent literal derivative implementation against the captured upstream derivative.
    from .cases import recipe

    reconstruction = float(recipe(task["case"])["survey"]["reconstruction factor"])
    widths = []
    for name in row["selected_labels"]:
        sigma_transverse = (
            3.26 * row["growth_ratio"] / np.sqrt(1 if "lya" in name else reconstruction)
        )
        widths.append([(1 + row["growth_rate"]) * sigma_transverse, sigma_transverse])
    literal = np.concatenate(
        [
            legacy_jacobian(
                reference_arrays["mean"]
                .reshape(len(reference_arrays["mu_axis"]), len(k_grid), -1)[i]
                .T,
                k_grid,
                mu,
                widths=widths,
            )
            for i, mu in enumerate(reference_arrays["mu_axis"])
        ]
    )
    error = relative(literal, reference_arrays["observed_j"])
    if error > 5e-12:
        raise ValueError(f"legacy derivative mismatch {error}")
    report.update(
        provenance=provenance,
        reference_metadata_sha256=digest(
            Path(reference) / task["case"] / "metadata.json"
        ),
        comparison=change(
            arrays["fisher"],
            reference_arrays["legacy_fisher"],
            arrays["pair_fisher"],
            reference_arrays["legacy_pair_fisher"],
            1,
            1,
        ),
        upstream_j_relative=error,
        field_psd_violations=int(np.count_nonzero(arrays["field_min_eigenvalue"] < 0)),
    )
    comparison = report["comparison"]
    report["passed"] = (
        comparison["fisher_relative"] <= 5e-12
        and max(comparison["error_relative"], comparison["pair_error_relative"]) <= 1e-6
    )
    return arrays, report
