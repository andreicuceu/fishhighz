"""Exclusive validation JSON/NPZ bundles; no executable-object deserialization."""

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from ..adapters.legacy_compat import plain
from .cases import CASE_IDS, bins, recipe, selection


def digest(path):
    """Hash a saved input or evidence file.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the file to inspect.

    Returns
    -------
    digest : str
        SHA-256 hexadecimal digest of the file bytes.
    """
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record_report(root, record):
    """Read one bounded, hash-bound report; inline schema-2 fixtures stay valid.

    Parameters
    ----------
    root : str or pathlib.Path
        Root directory of the input checkout or saved evidence bundle.
    record : dict
        Manifest entry containing the inline report or its relative path and
        hashes.

    Returns
    -------
    report : dict
        Verified external report or the retained inline historical report.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if "report_file" not in record:
        return record["report"]
    if "report" in record:
        raise ValueError("ambiguous inline/external report")
    root = Path(root).resolve()
    path = (root / record["report_file"]).resolve()
    if path.parent != root or path.suffix != ".json":
        raise ValueError("invalid report path")
    if digest(path) != record["report_sha256"]:
        raise ValueError("stale report file hash")
    report = json.loads(path.read_text())
    if canonical(report) != record["effective_hash"]:
        raise ValueError("stale effective report hash")
    return report


def canonical(value):
    """Hash a deterministic JSON representation of scientific metadata.

    Parameters
    ----------
    value : object
        JSON-compatible metadata, allowing NumPy values through plain
        conversion.

    Returns
    -------
    digest : str
        SHA-256 hexadecimal digest of the sorted finite JSON representation.
    """
    return hashlib.sha256(
        json.dumps(plain(value), sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def requests(suite, cases=None, bin_indices=None):
    """Explicit quick/full work inventory; full execution requires caller opt-in.

    Parameters
    ----------
    suite : str
        Requested quick or explicitly selected full validation inventory.
    cases : sequence of str or None
        Case identifiers to include; None uses the suite-defined inventory.
        Default is ``None``.
    bin_indices : sequence of int or None
        Zero-based bins to include; None uses the suite-defined bin selection.
        Default is ``None``.

    Returns
    -------
    requests : list of dict
        Ordered case/bin requests with bounds and original selected pairs.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    if suite not in ("quick", "full"):
        raise ValueError("suite must be quick or full")
    cases = (
        tuple(cases)
        if cases is not None
        else (("lya_qso_lbg_lae_15x2pt",) if suite == "quick" else CASE_IDS)
    )
    if suite == "quick":
        if cases != ("lya_qso_lbg_lae_15x2pt",) or bin_indices not in (None, [2], (2,)):
            raise ValueError("quick coverage is exactly 15x2pt bin 2")
        bin_indices = (2,)
    if not cases or len(set(cases)) != len(cases):
        raise ValueError("duplicate or empty cases")
    work = []
    for case in cases:
        bounds = bins(case)
        indices = (
            tuple(range(len(bounds))) if bin_indices is None else tuple(bin_indices)
        )
        if not indices or len(set(indices)) != len(indices):
            raise ValueError("empty or duplicate bin selection")
        for index in indices:
            if (
                isinstance(index, bool)
                or not isinstance(index, int)
                or not 0 <= index < len(bounds)
            ):
                raise ValueError("invalid bin index")
            work.append(
                dict(
                    case=case,
                    bin=index,
                    bounds=bounds[index],
                    pairs=selection(case).selected_pairs.tolist(),
                )
            )
    return work


def modern_requests(
    suite,
    cases=None,
    bin_indices=None,
    *,
    profiles=("accuracy",),
    kind="real_bao",
    recipe_revision=None,
):
    """Exact primary inventory; full two-profile coverage is 78 records.

    Parameters
    ----------
    suite : str
        Requested quick or explicitly selected full validation inventory.
    cases : sequence of str or None
        Case identifiers to include; None uses the suite-defined inventory.
        Default is ``None``.
    bin_indices : sequence of int or None
        Zero-based bins to include; None uses the suite-defined bin selection.
        Default is ``None``.
    profiles : sequence of str
        Ordered scientific profiles to include. Default is ``('accuracy',)``.
    kind : str
        Scientific record type used to distinguish primary and diagnostic
        evidence. Default is ``'real_bao'``.
    recipe_revision : str or None
        Declared revision of the physical recipe; None retains historical
        request semantics. Default is ``None``.

    Returns
    -------
    requests : list of dict
        Ordered case/bin/profile requests with fixed physical and numerical
        identities.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    from .schema import PROFILES, request

    profiles = tuple(profiles)
    if (
        not profiles
        or len(set(profiles)) != len(profiles)
        or any(p not in PROFILES for p in profiles)
    ):
        raise ValueError("invalid profile inventory")
    return [
        request(
            t["case"], t["bin"], profile, kind=kind, recipe_revision=recipe_revision
        )
        for t in requests(suite, cases, bin_indices)
        for profile in profiles
    ]


def execute(
    output,
    *,
    suite,
    worker,
    cases=None,
    bin_indices=None,
    inputs=None,
    profiles=("accuracy",),
    kind="real_bao",
    diagnostic_requests=(),
    recipe_revision=None,
):
    """Write schema-2 semantic evidence, retaining execution/science failures.

    Parameters
    ----------
    output : str or pathlib.Path
        Destination directory for the generated evidence or figures.
    suite : str
        Requested quick or explicitly selected full validation inventory.
    worker : callable
        Function taking one request and returning its numerical arrays and
        report.
    cases : sequence of str or None
        Case identifiers to include; None uses the suite-defined inventory.
        Default is ``None``.
    bin_indices : sequence of int or None
        Zero-based bins to include; None uses the suite-defined bin selection.
        Default is ``None``.
    inputs : sequence of str or pathlib.Path or None
        Source paths whose hashes are recorded in the manifest. Default is
        ``None``.
    profiles : sequence of str
        Ordered scientific profiles to include. Default is ``('accuracy',)``.
    kind : str
        Scientific record type used to distinguish primary and diagnostic
        evidence. Default is ``'real_bao'``.
    diagnostic_requests : sequence of dict
        Separately declared diagnostic requests. Default is ``()``.
    recipe_revision : str or None
        Declared revision of the physical recipe; None retains historical
        request semantics. Default is ``None``.

    Returns
    -------
    manifest : dict
        Written inventory with every attempted outcome and separate
        execution/scientific status.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.

    Notes
    -----
    Each worker returns a payload built for its exact request. A worker that
    merely returns Fisher/errors is invalid. Diagnostics have a separately
    predeclared exact inventory; failed tasks never disappear from the manifest.

    Creates a new output directory and writes compressed arrays, JSON reports and a progressively updated manifest. Existing output directories are rejected.
    """
    from .schema import validate_payload, validate_request

    work = modern_requests(
        suite,
        cases,
        bin_indices,
        profiles=profiles,
        kind=kind,
        recipe_revision=recipe_revision,
    )
    diagnostics = list(diagnostic_requests)
    for task in diagnostics:
        validate_request(task)
        if task["kind"] != "diagnostic":
            raise ValueError("invalid diagnostic request")
    if len({canonical(t) for t in diagnostics}) != len(diagnostics):
        raise ValueError("duplicate diagnostics")
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    manifest = dict(
        schema=3,
        kind="fishhighz-validation",
        suite=suite,
        cases=list(cases) if cases is not None else None,
        bin_indices=list(bin_indices) if bin_indices is not None else None,
        profiles=list(profiles),
        payload_kind=kind,
        recipe_revision=recipe_revision,
        requested=work,
        diagnostics_requested=diagnostics,
        records=[],
        diagnostics=[],
        execution_finished=False,
        complete=False,
        not_run_cases=[c for c in CASE_IDS if c not in {w["case"] for w in work}],
        inputs={str(Path(p).resolve()): digest(p) for p in (inputs or [])},
    )

    def save():
        """Atomically replace the manifest with the current execution inventory.

        Notes
        -----
        Writes manifest.next.json and replaces manifest.json in the enclosing output directory.
        """
        pending = out / "manifest.next.json"
        pending.write_text(
            json.dumps(plain(manifest), indent=2, allow_nan=False) + "\n"
        )
        pending.replace(out / "manifest.json")

    save()
    for group, tasks in (("records", work), ("diagnostics", diagnostics)):
        for index, task in enumerate(tasks):
            start = time.monotonic()
            record = dict(
                task=task, original_recipe_hash=canonical(recipe(task["case"]))
            )
            try:
                arrays, report = worker(json.loads(json.dumps(task)))
                arrays = {k: np.asarray(v) for k, v in arrays.items()}
                name = f"{group}-{index:03d}.npz"
                # Preserve failed scientific evidence as well; never load object arrays.
                if any(a.dtype.hasobject for a in arrays.values()):
                    raise ValueError("object evidence forbidden")
                np.savez_compressed(out / name, **arrays)
                record.update(
                    arrays=name,
                    sha256=digest(out / name),
                    inventory={k: list(a.shape) for k, a in arrays.items()},
                    effective_hash=canonical(report),
                )
                if task["kind"] in ("real_bao", "diagnostic"):
                    report_name = f"{group}-{index:03d}.report.json"
                    (out / report_name).write_text(
                        json.dumps(plain(report), indent=2, allow_nan=False) + "\n"
                    )
                    record.update(
                        report_file=report_name, report_sha256=digest(out / report_name)
                    )
                else:
                    record["report"] = plain(report)
                validate_payload(task, arrays, report, require_pass=False)
                record["status"] = "completed"
                try:
                    validate_payload(task, arrays, report, require_pass=True)
                    record["scientific_passed"] = True
                except Exception as error:
                    record.update(scientific_passed=False, scientific_error=str(error))
            except Exception as error:
                record.update(
                    status="failed",
                    scientific_passed=False,
                    error=f"{type(error).__name__}: {error}",
                )
            if task.get("recipe_revision"):
                record["numerical_status"] = (
                    (
                        "qualified_finite_refinement"
                        if task["profile"] == "accuracy"
                        else "reproducibility_control"
                    )
                    if record.get("scientific_passed")
                    else "unresolved"
                )
            record["seconds"] = time.monotonic() - start
            manifest[group].append(record)
            save()
    manifest["execution_finished"] = True
    manifest["complete"] = all(
        r.get("scientific_passed") is True
        for group in ("records", "diagnostics")
        for r in manifest[group]
    )
    save()
    return manifest


def inspect_legacy(output):
    """Explicit limited inspection; schema 1 can never pass schema-2 acceptance.

    Parameters
    ----------
    output : str or pathlib.Path
        Destination directory for the generated evidence or figures.

    Returns
    -------
    inspection : dict
        Explicitly limited historical status; scientific completeness remains
        false.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    out = Path(output)
    manifest = json.loads((out / "manifest.json").read_text())
    if manifest.get("schema") != 1:
        raise ValueError("not historical schema 1")
    return dict(
        schema=1,
        limited=True,
        complete=False,
        execution_finished=manifest.get("complete") is True,
        limitation="Historical hash/shape checks do not establish scientific payload consistency",
    )


def check(output, *, verify_sources=True):
    """Independent offline semantic validation with exact primary/diagnostic inventories.

    Parameters
    ----------
    output : str or pathlib.Path
        Destination directory for the generated evidence or figures.
    verify_sources : bool
        Whether to rehash the original source files as well as saved evidence.
        Default is ``True``.

    Returns
    -------
    manifest : dict
        Checked inventory; schema-2 evidence remains explicitly limited.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
    from .schema import validate_payload, validate_request

    out = Path(output).resolve()
    manifest = json.loads((out / "manifest.json").read_text())
    if (
        manifest.get("schema") not in (2, 3)
        or manifest.get("kind") != "fishhighz-validation"
    ):
        raise ValueError(
            "schema 2 required; inspect_legacy exposes historical limitations"
        )
    expected = plain(
        modern_requests(
            manifest["suite"],
            manifest["cases"],
            manifest["bin_indices"],
            profiles=manifest["profiles"],
            kind=manifest["payload_kind"],
            recipe_revision=manifest.get("recipe_revision"),
        )
    )
    if manifest["requested"] != expected:
        raise ValueError("wrong requested case/profile/bin inventory")
    diagnostics = manifest["diagnostics_requested"]
    if len({canonical(t) for t in diagnostics}) != len(diagnostics):
        raise ValueError("duplicate diagnostic inventory")
    for task in diagnostics:
        validate_request(task)
        if task["kind"] != "diagnostic":
            raise ValueError("incorrect diagnostic kind")
    for group, tasks in (("records", expected), ("diagnostics", diagnostics)):
        if [r["task"] for r in manifest[group]] != tasks:
            raise ValueError("missing/extra/swapped case/profile/bin/diagnostic")
    if manifest["not_run_cases"] != [
        c for c in CASE_IDS if c not in {w["case"] for w in expected}
    ]:
        raise ValueError("incorrect not-run inventory")
    if (
        manifest.get("execution_finished") is not True
        or manifest.get("complete") is not True
    ):
        raise ValueError("partial or scientifically failed validation")
    if verify_sources:
        for path, sha in manifest["inputs"].items():
            if digest(path) != sha:
                raise ValueError("stale source hash")
    for group in ("records", "diagnostics"):
        for record in manifest[group]:
            if (
                record.get("status") != "completed"
                or record.get("scientific_passed") is not True
            ):
                raise ValueError("partial/failed validation record")
            report = record_report(out, record)
            if record["original_recipe_hash"] != canonical(
                recipe(record["task"]["case"])
            ) or record["effective_hash"] != canonical(report):
                raise ValueError("stale original/effective recipe/report hash")
            path = (out / record["arrays"]).resolve()
            if path.parent != out or not path.name.endswith(".npz"):
                raise ValueError("invalid evidence path")
            if digest(path) != record["sha256"]:
                raise ValueError("stale array hash")
            with np.load(path, allow_pickle=False) as data:
                if set(data.files) != set(record["inventory"]):
                    raise ValueError("wrong array inventory")
                arrays = {name: data[name] for name in data.files}
            for name, a in arrays.items():
                if list(a.shape) != record["inventory"][name]:
                    raise ValueError("wrong recorded shape")
            validate_payload(record["task"], arrays, report, schema=manifest["schema"])
    if manifest["schema"] == 2:
        manifest = {
            **manifest,
            "limited": True,
            "complete": False,
            "limitation": "Historical schema 2 does not bind convergence operands to actual trials",
        }
    return manifest
