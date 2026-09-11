#!/usr/bin/env python3
"""Capture and validate reproducible lyaforecast DESI-2 baselines."""

from __future__ import annotations

import argparse
import configparser
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKER = Path(__file__).with_name("lyaforecast_baseline_worker.py")
SCHEMA = "fishhighz.lyaforecast-baseline"
SCHEMA_VERSION = 2
LEGACY_SCHEMA_VERSION = 1
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
PRIMARY_CASES = {
    "lbg_lae_3x2pt.ini": 3,
    "lya_lbg_lae_3x2pt.ini": 3,
    "lya_lbg_lae_6x2pt.ini": 6,
    "lya_qso_2x2pt.ini": 2,
    "lya_qso_lbg_lae_4x2pt.ini": 4,
    "lya_qso_lbg_lae_8x2pt.ini": 8,
    "lya_qso_lbg_lae_15x2pt.ini": 15,
}
REPEAT_CASE = "lya_qso_lbg_lae_15x2pt.ini"
SUITE_CASES = {"quick": [REPEAT_CASE], "full": list(PRIMARY_CASES)}
PATH_OPTIONS = {("cosmo", "filename")}
BASE_CASE_ARTIFACTS = {
    "original_config",
    "effective_config",
    "configuration_comparison",
    "worker_request",
    "stdout",
    "stderr",
    "status",
    "worker_response",
    "result_pickle",
    "result_json",
    "metadata",
    "summary",
    "forecast_log",
}
LEGACY_SHARED_ARTIFACTS = {
    "environment",
    "probe_stdout",
    "probe_stderr",
    "reference_git",
    "reference_dirty_diff",
    "source_manifest",
    "input_manifest",
    "tool_state",
    "tool_git_diff",
}
CURRENT_SHARED_ARTIFACTS = LEGACY_SHARED_ARTIFACTS | {"configuration_manifest"}
COORDINATE_RTOL = 5e-14
COORDINATE_ATOL = 5e-14


class BaselineError(RuntimeError):
    """Raised when capture evidence is incomplete or inconsistent."""


def utc_now() -> str:
    """Return a UTC timestamp suitable for evidence metadata."""
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds")


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one regular file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    """Write deterministic, strict JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def read_json(path: Path) -> object:
    """Read JSON with duplicate-object-key rejection."""

    def reject_duplicates(items: list[tuple[str, object]]) -> dict[str, object]:
        result = {}
        for key, value in items:
            if key in result:
                raise BaselineError(f"duplicate JSON key {key!r} in {path}")
            result[key] = value
        return result

    try:
        with path.open(encoding="utf-8") as stream:
            return json.load(stream, object_pairs_hook=reject_duplicates)
    except (OSError, json.JSONDecodeError) as error:
        raise BaselineError(f"cannot read JSON artifact {path}: {error}") from error


def relative(bundle: Path, path: Path) -> str:
    """Return a portable path relative to the bundle root."""
    try:
        return path.resolve().relative_to(bundle.resolve()).as_posix()
    except ValueError as error:
        raise BaselineError(f"artifact is outside bundle: {path}") from error


def bundle_path(bundle: Path, path_value: object, owner: str) -> Path:
    """Resolve a nested metadata path and require bundle containment."""
    if not isinstance(path_value, str):
        raise BaselineError(f"{owner}: non-string artifact path")
    path = (bundle / path_value).resolve()
    try:
        path.relative_to(bundle.resolve())
    except ValueError as error:
        raise BaselineError(
            f"{owner}: artifact escapes bundle: {path_value}"
        ) from error
    return path


def artifact(bundle: Path, path: Path) -> dict[str, str]:
    """Describe one bundle artifact."""
    if not path.is_file():
        raise BaselineError(f"missing artifact while building manifest: {path}")
    return {"path": relative(bundle, path), "sha256": sha256(path)}


def parse_config(path: Path) -> configparser.ConfigParser:
    """Read an INI while preserving section and option case/order."""
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    loaded = parser.read(path, encoding="utf-8")
    if loaded != [str(path)]:
        raise BaselineError(f"could not parse configuration {path}")
    return parser


def config_items(
    parser: configparser.ConfigParser,
) -> list[tuple[str, list[tuple[str, str]]]]:
    """Return ordered explicit INI contents, excluding ConfigParser defaults."""
    return [
        (section, list(parser._sections[section].items()))  # noqa: SLF001
        for section in parser.sections()
    ]


def allowed_path_option(section: str, option: str) -> bool:
    """Return whether Step 02 permits changing this effective field."""
    return (section, option) in PATH_OPTIONS or (
        section.startswith("tracer ") and option in {"dn dz", "snr-file-dir"}
    )


def compare_configs(original: Path, effective: Path) -> dict[str, object]:
    """Compare logical INI contents against the Step 02 path-only allowlist."""
    original_parser = parse_config(original)
    effective_parser = parse_config(effective)
    original_items = config_items(original_parser)
    effective_items = config_items(effective_parser)
    original_sections = [section for section, _ in original_items]
    effective_sections = [section for section, _ in effective_items]
    option_order_matches = original_sections == effective_sections
    changes = []
    unauthorized = []

    if original_sections != effective_sections:
        unauthorized.append(
            {
                "kind": "section-order-or-membership",
                "original": original_sections,
                "effective": effective_sections,
            }
        )

    for section in set(original_sections) | set(effective_sections):
        before = dict(original_parser._sections.get(section, {}))  # noqa: SLF001
        after = dict(effective_parser._sections.get(section, {}))  # noqa: SLF001
        before_order = list(before)
        after_order = list(after)
        if before_order != after_order:
            option_order_matches = False
            unauthorized.append(
                {
                    "kind": "option-order-or-membership",
                    "section": section,
                    "original": before_order,
                    "effective": after_order,
                }
            )
        for option in set(before) | set(after):
            if before.get(option) == after.get(option):
                continue
            change = {
                "section": section,
                "option": option,
                "original": before.get(option),
                "effective": after.get(option),
            }
            if allowed_path_option(section, option) or (
                section == "output" and option == "filename"
            ):
                changes.append(change)
            else:
                unauthorized.append({"kind": "value", **change})

    return {
        "section_order_matches": original_sections == effective_sections,
        "option_order_matches": option_order_matches,
        "allowed_changes": sorted(
            changes, key=lambda item: (item["section"], item["option"])
        ),
        "unauthorized_changes": unauthorized,
    }


def worker_environment() -> dict[str, str]:
    """Return a minimal worker environment with deterministic thread settings."""
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    for variable in THREAD_VARIABLES:
        environment[variable] = "1"
    return environment


def worker_command(
    python: Path,
    worker: Path,
    action: str,
    reference: Path,
    response: Path,
    *extra: str,
) -> list[str]:
    """Construct an isolated worker command."""
    return [
        str(python),
        "-B",
        "-I",
        str(worker),
        action,
        "--reference-checkout",
        str(reference),
        "--response",
        str(response),
        *extra,
    ]


def run_worker(
    command: list[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> subprocess.CompletedProcess[str]:
    """Run one worker while retaining its complete output streams."""
    with (
        stdout_path.open("w", encoding="utf-8") as stdout,
        stderr_path.open("w", encoding="utf-8") as stderr,
    ):
        return subprocess.run(
            command,
            cwd=cwd,
            env=worker_environment(),
            stdout=stdout,
            stderr=stderr,
            text=True,
            check=False,
        )


def git_output(checkout: Path, *arguments: str) -> str:
    """Return Git output or an actionable capture error."""
    try:
        return subprocess.run(
            ["git", "-C", str(checkout), *arguments],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise BaselineError(
            f"cannot record Git state for {checkout}: {error}"
        ) from error


def reference_git_state(reference: Path, destination: Path) -> dict[str, object]:
    """Record the reference revision, branch, status, and dirty patch."""
    diff_path = destination / "dirty.diff"
    diff_path.write_text(
        git_output(reference, "diff", "--binary", "HEAD"), encoding="utf-8"
    )
    state = {
        "revision": git_output(reference, "rev-parse", "HEAD").strip(),
        "branch": git_output(reference, "branch", "--show-current").strip(),
        "status_porcelain": git_output(
            reference, "status", "--porcelain=v1", "--untracked-files=all"
        ),
        "dirty_diff_path": diff_path.name,
        "dirty_diff_sha256": sha256(diff_path),
    }
    write_json(destination / "git.json", state)
    return state


def source_candidates(reference: Path) -> list[Path]:
    """Select imported Python source and package metadata for preservation."""
    candidates = list((reference / "lyaforecast").rglob("*.py"))
    for name in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
        path = reference / name
        if path.is_file():
            candidates.append(path)
    egg_info = reference / "lyaforecast.egg-info"
    if egg_info.is_dir():
        candidates.extend(path for path in egg_info.rglob("*") if path.is_file())
    return sorted(set(path.resolve() for path in candidates))


def source_inventory(reference: Path) -> list[str]:
    """Return the bounded reference-source file inventory."""
    root = reference.resolve()
    return [path.relative_to(root).as_posix() for path in source_candidates(root)]


def directory_inventory(root: Path) -> list[dict[str, str]]:
    """Return a bounded recursive inventory of files and directories."""
    inventory = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            kind = "symlink"
        elif path.is_dir():
            kind = "directory"
        elif path.is_file():
            kind = "file"
        else:
            kind = "other"
        inventory.append({"path": path.relative_to(root).as_posix(), "kind": kind})
    return inventory


def snapshot_configurations(
    bundle: Path,
    examples: Path,
    case_names: list[str],
) -> tuple[dict[str, Path], dict[str, object]]:
    """Snapshot authoritative INIs before resolution or scientific execution."""
    destination = bundle / "reference" / "configurations"
    snapshot_root = destination / "snapshot"
    snapshot_root.mkdir(parents=True)
    inventory = sorted(path.name for path in examples.glob("*.ini"))
    records = []
    originals = {}
    for case_name in case_names:
        original = examples / case_name
        copied = snapshot_root / case_name
        shutil.copyfile(original, copied)
        digest = sha256(original)
        records.append(
            {
                "case": case_name,
                "reference_path": str(original.resolve()),
                "snapshot_path": relative(bundle, copied),
                "before_sha256": digest,
                "snapshot_sha256": sha256(copied),
                "after_sha256": None,
                "unchanged": None,
            }
        )
        originals[case_name] = copied
    manifest = {
        "scope_cases": case_names,
        "authoritative_inventory_before": inventory,
        "authoritative_inventory_after": None,
        "inventory_unchanged": None,
        "files": records,
        "all_unchanged": None,
    }
    write_json(destination / "configuration_manifest.json", manifest)
    return originals, manifest


def finalize_configuration_hashes(
    bundle: Path,
    examples: Path,
    config_manifest: dict[str, object],
) -> None:
    """Verify authoritative INI contents and inventory after capture."""
    inventory_after = sorted(path.name for path in examples.glob("*.ini"))
    config_manifest["authoritative_inventory_after"] = inventory_after
    config_manifest["inventory_unchanged"] = (
        inventory_after == config_manifest["authoritative_inventory_before"]
    )
    unchanged = config_manifest["inventory_unchanged"]
    for record in config_manifest["files"]:
        original = Path(record["reference_path"])
        after = sha256(original) if original.is_file() else None
        record["after_sha256"] = after
        record["unchanged"] = after == record["before_sha256"]
        unchanged = unchanged and record["unchanged"]
    config_manifest["all_unchanged"] = unchanged
    write_json(
        bundle / "reference" / "configurations" / "configuration_manifest.json",
        config_manifest,
    )


def snapshot_source(reference: Path, destination: Path) -> dict[str, object]:
    """Snapshot reference source/metadata and its pre-capture hashes."""
    records = []
    snapshot_root = destination / "snapshot"
    for original in source_candidates(reference):
        relpath = original.relative_to(reference.resolve())
        copied = snapshot_root / relpath
        copied.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, copied)
        digest = sha256(original)
        records.append(
            {
                "reference_path": relpath.as_posix(),
                "snapshot_path": (
                    Path("reference/source/snapshot") / relpath
                ).as_posix(),
                "before_sha256": digest,
                "snapshot_sha256": sha256(copied),
                "after_sha256": None,
                "unchanged": None,
            }
        )
    manifest = {
        "inventory_before": source_inventory(reference),
        "inventory_after": None,
        "added": None,
        "deleted": None,
        "inventory_unchanged": None,
        "files": records,
        "all_unchanged": None,
    }
    write_json(destination / "source_manifest.json", manifest)
    return manifest


def finalize_source_hashes(
    reference: Path, destination: Path, source_manifest: dict[str, object]
) -> None:
    """Record post-capture source hashes without changing the reference."""
    inventory_after = source_inventory(reference)
    before = source_manifest["inventory_before"]
    source_manifest["inventory_after"] = inventory_after
    source_manifest["added"] = sorted(set(inventory_after) - set(before))
    source_manifest["deleted"] = sorted(set(before) - set(inventory_after))
    source_manifest["inventory_unchanged"] = inventory_after == before
    unchanged = source_manifest["inventory_unchanged"]
    for record in source_manifest["files"]:
        original = reference / record["reference_path"]
        after = sha256(original) if original.is_file() else None
        record["after_sha256"] = after
        record["unchanged"] = after == record["before_sha256"]
        unchanged = unchanged and record["unchanged"]
    source_manifest["all_unchanged"] = unchanged
    write_json(destination / "source_manifest.json", source_manifest)


def copy_input(
    bundle: Path,
    reference: Path,
    original: Path,
    kind: str,
    input_id: str,
) -> dict[str, object]:
    """Copy one resolved input file/directory and record every entry."""
    original = original.resolve()
    try:
        relpath = original.relative_to(reference.resolve())
        snapshot = bundle / "inputs" / "reference" / relpath
    except ValueError:
        snapshot = bundle / "inputs" / "external" / f"{input_id}-{original.name}"

    if kind == "file":
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, snapshot)
        digest = sha256(original)
        return {
            "id": input_id,
            "kind": kind,
            "resolved_original_path": str(original),
            "snapshot_path": relative(bundle, snapshot),
            "entries": [
                {
                    "path": ".",
                    "before_sha256": digest,
                    "snapshot_sha256": sha256(snapshot),
                    "after_sha256": None,
                    "unchanged": None,
                }
            ],
            "directories": [],
            "immediate_directory_entries": [],
            "inventory_before": [{"path": ".", "kind": "file"}],
            "inventory_after": None,
            "inventory_unchanged": None,
            "references": [],
        }

    if kind != "directory" or not original.is_dir():
        raise BaselineError(f"unsupported or missing resolved input: {original}")
    snapshot.mkdir(parents=True, exist_ok=False)
    entries = []
    directories = []
    immediate = []
    for child in sorted(original.iterdir(), key=lambda item: item.name):
        immediate.append(
            {"name": child.name, "kind": "directory" if child.is_dir() else "file"}
        )
    for path in sorted(original.rglob("*")):
        relpath = path.relative_to(original)
        copied = snapshot / relpath
        if path.is_symlink():
            raise BaselineError(f"symlink input requires review: {path}")
        if path.is_dir():
            copied.mkdir(parents=True, exist_ok=True)
            directories.append(relpath.as_posix())
        elif path.is_file():
            copied.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, copied)
            digest = sha256(path)
            entries.append(
                {
                    "path": relpath.as_posix(),
                    "before_sha256": digest,
                    "snapshot_sha256": sha256(copied),
                    "after_sha256": None,
                    "unchanged": None,
                }
            )
        else:
            raise BaselineError(f"unsupported input directory entry: {path}")
    return {
        "id": input_id,
        "kind": kind,
        "resolved_original_path": str(original),
        "snapshot_path": relative(bundle, snapshot),
        "entries": entries,
        "directories": directories,
        "immediate_directory_entries": immediate,
        "inventory_before": directory_inventory(original),
        "inventory_after": None,
        "inventory_unchanged": None,
        "references": [],
    }


def finalize_input_hashes(bundle: Path, input_manifest: dict[str, object]) -> None:
    """Record post-capture input hashes and preservation results."""
    all_unchanged = True
    for record in input_manifest["inputs"]:
        original_root = Path(record["resolved_original_path"])
        inventory_after = (
            directory_inventory(original_root)
            if record["kind"] == "directory" and original_root.is_dir()
            else ([{"path": ".", "kind": "file"}] if original_root.is_file() else [])
        )
        record["inventory_after"] = inventory_after
        record["inventory_unchanged"] = inventory_after == record["inventory_before"]
        all_unchanged = all_unchanged and record["inventory_unchanged"]
        for entry in record["entries"]:
            original = (
                original_root if entry["path"] == "." else original_root / entry["path"]
            )
            after = sha256(original) if original.is_file() else None
            entry["after_sha256"] = after
            entry["unchanged"] = after == entry["before_sha256"]
            all_unchanged = all_unchanged and entry["unchanged"]
            snapshot_root = bundle / record["snapshot_path"]
            snapshot = (
                snapshot_root if entry["path"] == "." else snapshot_root / entry["path"]
            )
            entry["snapshot_sha256"] = sha256(snapshot)
    input_manifest["all_unchanged"] = all_unchanged
    write_json(bundle / "inputs" / "manifest.json", input_manifest)


def capture_tool_state(bundle: Path, worker: Path, invocation: list[str]) -> None:
    """Snapshot the capture tools and FishHighz Git state."""
    destination = bundle / "tool"
    destination.mkdir(parents=True)
    tools = []
    for source in (Path(__file__).resolve(), worker.resolve()):
        target = destination / source.name
        shutil.copyfile(source, target)
        tools.append(
            {
                "source_path": str(source),
                "snapshot_path": relative(bundle, target),
                "sha256": sha256(source),
            }
        )
    diff = destination / "git.diff"
    diff.write_text(git_output(ROOT, "diff", "--binary", "HEAD"), encoding="utf-8")
    state = {
        "invocation": invocation,
        "fishhighz_revision": git_output(ROOT, "rev-parse", "HEAD").strip(),
        "fishhighz_branch": git_output(ROOT, "branch", "--show-current").strip(),
        "fishhighz_status_porcelain": git_output(
            ROOT, "status", "--porcelain=v1", "--untracked-files=all"
        ),
        "git_diff_path": relative(bundle, diff),
        "git_diff_sha256": sha256(diff),
        "tools": tools,
        "thread_settings": {variable: "1" for variable in THREAD_VARIABLES},
    }
    write_json(destination / "tool.json", state)


def resolve_inputs(
    bundle: Path,
    reference: Path,
    python: Path,
    worker: Path,
    originals: dict[str, Path],
) -> dict[str, object]:
    """Resolve all path inputs via the reference utilities, then snapshot them."""
    resolutions = {}
    resolve_dir = bundle / "reference" / "resolution"
    resolve_dir.mkdir(parents=True)
    for case_name, original in originals.items():
        stem = Path(case_name).stem
        response = resolve_dir / f"{stem}.json"
        command = worker_command(
            python,
            worker,
            "resolve",
            reference,
            response,
            "--config",
            str(original),
        )
        completed = run_worker(
            command,
            resolve_dir,
            resolve_dir / f"{stem}.stdout.log",
            resolve_dir / f"{stem}.stderr.log",
        )
        if completed.returncode != 0:
            raise BaselineError(
                f"input resolution failed for {case_name}; see "
                f"{relative(bundle, resolve_dir / f'{stem}.stderr.log')}"
            )
        resolutions[case_name] = read_json(response)

    by_original = {}
    records = []
    for case_name in originals:
        for resolved in resolutions[case_name]["paths"]:
            original_path = str(Path(resolved["resolved_path"]).resolve())
            record = by_original.get(original_path)
            if record is None:
                input_id = f"input-{len(records) + 1:03d}"
                record = copy_input(
                    bundle,
                    reference,
                    Path(original_path),
                    resolved["kind"],
                    input_id,
                )
                records.append(record)
                by_original[original_path] = record
            reference_record = {
                "case": case_name,
                "section": resolved["section"],
                "option": resolved["option"],
                "original_setting": resolved["original_setting"],
            }
            if reference_record not in record["references"]:
                record["references"].append(reference_record)
            resolved["input_id"] = record["id"]
            resolved["snapshot_path"] = record["snapshot_path"]
            resolved["effective_path"] = str(
                (bundle / record["snapshot_path"]).resolve()
            )

    manifest = {
        "resolver": "lyaforecast.utils.get_file/get_dir",
        "inputs": records,
        "case_resolutions": resolutions,
        "all_unchanged": None,
    }
    write_json(bundle / "inputs" / "manifest.json", manifest)
    return manifest


def create_effective_config(
    bundle: Path,
    case_dir: Path,
    original: Path,
    resolutions: dict[str, object],
) -> tuple[Path, dict[str, object]]:
    """Create a path-only effective INI that runs entirely from snapshots."""
    effective = case_dir / "effective.ini"
    parser = parse_config(original)
    for resolved in resolutions["paths"]:
        snapshot = (bundle / resolved["snapshot_path"]).resolve()
        parser[resolved["section"]][resolved["option"]] = str(snapshot)
    parser["output"]["filename"] = str((case_dir / "generated" / "forecast").resolve())
    with effective.open("w", encoding="utf-8") as stream:
        parser.write(stream)
    comparison = compare_configs(original, effective)
    if comparison["unauthorized_changes"]:
        raise BaselineError(f"unauthorized effective configuration change: {original}")
    return effective, comparison


def case_artifacts(bundle: Path, case_dir: Path) -> dict[str, dict[str, str]]:
    """Hash all required artifacts for one completed worker."""
    names = {
        "original_config": case_dir / "original.ini",
        "effective_config": case_dir / "effective.ini",
        "configuration_comparison": case_dir / "configuration-comparison.json",
        "worker_request": case_dir / "worker-request.json",
        "stdout": case_dir / "stdout.log",
        "stderr": case_dir / "stderr.log",
        "status": case_dir / "status.json",
        "worker_response": case_dir / "worker-response.json",
        "result_pickle": case_dir / "result.pkl",
        "result_json": case_dir / "result.json",
        "metadata": case_dir / "metadata.json",
        "summary": case_dir / "summary.json",
        "forecast_log": case_dir / "generated" / "forecast.log",
    }
    return {name: artifact(bundle, path) for name, path in names.items()}


def run_case(
    bundle: Path,
    reference: Path,
    python: Path,
    worker: Path,
    case_name: str,
    original_source: Path,
    resolutions: dict[str, object],
    role: str,
    directory: Path,
) -> dict[str, object]:
    """Run one full configuration in a fresh subprocess and record evidence."""
    directory.mkdir(parents=True, exist_ok=False)
    original = directory / "original.ini"
    shutil.copyfile(original_source, original)
    effective, comparison = create_effective_config(
        bundle, directory, original, resolutions
    )
    comparison_path = directory / "configuration-comparison.json"
    write_json(comparison_path, comparison)
    response = directory / "worker-response.json"
    request = {
        "role": role,
        "case": case_name,
        "effective_config": str(effective.resolve()),
        "result_pickle": str((directory / "result.pkl").resolve()),
        "result_json": str((directory / "result.json").resolve()),
        "metadata": str((directory / "metadata.json").resolve()),
        "summary": str((directory / "summary.json").resolve()),
    }
    request_path = directory / "worker-request.json"
    write_json(request_path, request)
    command = worker_command(
        python,
        worker,
        "run",
        reference,
        response,
        "--request",
        str(request_path),
    )
    started_utc = utc_now()
    started = time.perf_counter()
    completed = run_worker(
        command, directory, directory / "stdout.log", directory / "stderr.log"
    )
    ended_utc = utc_now()
    status = {
        "role": role,
        "case": case_name,
        "state": "complete" if completed.returncode == 0 else "failed",
        "exit_code": completed.returncode,
        "started_utc": started_utc,
        "ended_utc": ended_utc,
        "wall_seconds": time.perf_counter() - started,
        "command": command,
        "worker_response": relative(bundle, response) if response.is_file() else None,
    }
    if response.is_file():
        status["worker"] = read_json(response)
    write_json(directory / "status.json", status)
    record = {
        "role": role,
        "case": case_name,
        "expected_selected_pair_count": PRIMARY_CASES[case_name],
        "directory": relative(bundle, directory),
        "state": status["state"],
        "exit_code": completed.returncode,
        "started_utc": started_utc,
        "ended_utc": ended_utc,
        "wall_seconds": status["wall_seconds"],
        "artifacts": {},
    }
    if completed.returncode == 0:
        record["artifacts"] = case_artifacts(bundle, directory)
    else:
        for name, path in {
            "original_config": original,
            "effective_config": effective,
            "configuration_comparison": comparison_path,
            "worker_request": request_path,
            "stdout": directory / "stdout.log",
            "stderr": directory / "stderr.log",
            "status": directory / "status.json",
            "worker_response": response,
        }.items():
            if path.is_file():
                record["artifacts"][name] = artifact(bundle, path)
    return record


def numeric_values(node: dict[str, object], location: str) -> list[float]:
    """Read a scalar or flat ndarray node as floats."""
    kind = node.get("kind")
    if kind == "ndarray":
        data = node.get("data")
        if not isinstance(data, list):
            raise BaselineError(f"{location}: ndarray data is not a list")
        return [decode_number(value, location) for value in data]
    if kind in {"python-scalar", "numpy-scalar"}:
        return [decode_number(node.get("value"), location)]
    raise BaselineError(f"{location}: expected numeric result, found {kind!r}")


def decode_number(value: object, location: str) -> float:
    """Decode a strict-JSON numeric value or explicit non-finite marker."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        if isinstance(value, str) and value in {"nan", "+inf", "-inf"}:
            raise BaselineError(f"{location}: non-finite value {value}")
        raise BaselineError(f"{location}: malformed numeric value {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise BaselineError(f"{location}: non-finite value")
    return number


def dict_node(node: object, location: str) -> tuple[list[str], dict[str, object]]:
    """Validate and unpack the ordered-dictionary JSON representation."""
    if not isinstance(node, dict) or node.get("kind") != "dict":
        raise BaselineError(f"{location}: expected dictionary node")
    items = node.get("items")
    if not isinstance(items, list):
        raise BaselineError(f"{location}: dictionary items are malformed")
    keys = []
    values = {}
    for index, item in enumerate(items):
        if not isinstance(item, dict) or set(item) != {"key", "value"}:
            raise BaselineError(f"{location}: malformed dictionary item {index}")
        key = item["key"]
        if not isinstance(key, str) or key in values:
            raise BaselineError(f"{location}: duplicate or non-string key {key!r}")
        keys.append(key)
        values[key] = item["value"]
    return keys, values


def array_node(
    node: object, location: str, expected_shape: list[int] | None = None
) -> list[float]:
    """Validate a numerical ndarray node and return its flattened values."""
    if not isinstance(node, dict) or node.get("kind") != "ndarray":
        raise BaselineError(f"{location}: expected ndarray")
    shape = node.get("shape")
    dtype = node.get("dtype")
    if (
        not isinstance(shape, list)
        or not all(isinstance(value, int) and value >= 0 for value in shape)
        or not isinstance(dtype, str)
    ):
        raise BaselineError(f"{location}: malformed ndarray metadata")
    if expected_shape is not None and shape != expected_shape:
        raise BaselineError(
            f"{location}: shape {shape} does not match expected {expected_shape}"
        )
    values = numeric_values(node, location)
    expected_size = math.prod(shape)
    if len(values) != expected_size:
        raise BaselineError(
            f"{location}: {len(values)} values do not match shape {shape}"
        )
    return values


def config_value(path: Path, section: str, option: str) -> str:
    """Return one explicit configuration value."""
    parser = parse_config(path)
    if section not in parser or option not in parser._sections[section]:  # noqa: SLF001
        raise BaselineError(f"{path}: missing [{section}] {option}")
    return parser._sections[section][option]  # noqa: SLF001


def expected_input_options(config: configparser.ConfigParser) -> list[tuple[str, str]]:
    """Return every reference path field consumed by the authoritative workflow."""
    options = [("cosmo", "filename")]
    for section in config.sections():
        if section.startswith("tracer "):
            options.append((section, "dn dz"))
            if "snr-file-dir" in config[section]:
                options.append((section, "snr-file-dir"))
    return options


def validate_config_evidence(
    bundle: Path,
    case: dict[str, object],
    input_manifest: dict[str, object],
) -> None:
    """Recompute configuration equivalence and verify path provenance."""
    artifacts = case["artifacts"]
    original = bundle / artifacts["original_config"]["path"]
    effective = bundle / artifacts["effective_config"]["path"]
    stored_path = bundle / artifacts["configuration_comparison"]["path"]
    recomputed = compare_configs(original, effective)
    if recomputed != read_json(stored_path):
        raise BaselineError(f"{case['case']}: configuration comparison is stale")
    if recomputed["unauthorized_changes"]:
        first = recomputed["unauthorized_changes"][0]
        raise BaselineError(
            f"{case['case']}: unauthorized configuration change: {first}"
        )
    if (
        not recomputed["section_order_matches"]
        or not recomputed["option_order_matches"]
    ):
        raise BaselineError(f"{case['case']}: configuration order changed")

    expected_changes = {("output", "filename")}
    original_parser = parse_config(original)
    resolutions = input_manifest["case_resolutions"].get(case["case"])
    if not isinstance(resolutions, dict):
        raise BaselineError(f"{case['case']}: missing input-resolution provenance")
    configuration_snapshot = (
        bundle / "reference" / "configurations" / "snapshot" / case["case"]
    )
    if configuration_snapshot.is_file():
        if sha256(original) != sha256(configuration_snapshot):
            raise BaselineError(
                f"{case['case']}: case original differs from pre-run snapshot"
            )
        recorded_config = resolutions.get("config")
        snapshot_suffix = relative(bundle, configuration_snapshot)
        if not isinstance(recorded_config, str) or not Path(
            recorded_config
        ).as_posix().endswith(snapshot_suffix):
            raise BaselineError(
                f"{case['case']}: resolver did not use the pre-run snapshot"
            )
    resolution_keys = [
        (item.get("section"), item.get("option"))
        for item in resolutions.get("paths", [])
    ]
    expected_resolution_keys = expected_input_options(original_parser)
    if resolution_keys != expected_resolution_keys:
        raise BaselineError(
            f"{case['case']}: input-resolution fields {resolution_keys} do not match "
            f"configuration fields {expected_resolution_keys}"
        )
    inputs_by_id = {item["id"]: item for item in input_manifest["inputs"]}
    for resolved in resolutions.get("paths", []):
        key = (resolved.get("section"), resolved.get("option"))
        expected_changes.add(key)
        input_record = inputs_by_id.get(resolved.get("input_id"))
        if input_record is None:
            raise BaselineError(
                f"{case['case']}: unresolved input provenance for {key}"
            )
        expected_effective = str((bundle / input_record["snapshot_path"]).resolve())
        recorded_effective = config_value(effective, *key)
        # Effective paths remain absolute records of the capture location after relocation.
        if recorded_effective != resolved.get("effective_path"):
            raise BaselineError(
                f"{case['case']}: effective input path record differs for {key}"
            )
        if (
            not Path(recorded_effective)
            .as_posix()
            .endswith(input_record["snapshot_path"])
            or Path(recorded_effective).name != Path(expected_effective).name
        ):
            raise BaselineError(
                f"{case['case']}: effective input target differs for {key}"
            )
        if config_value(original, *key) != resolved.get("original_setting"):
            raise BaselineError(
                f"{case['case']}: original input setting differs for {key}"
            )
    actual_changes = {
        (item["section"], item["option"]) for item in recomputed["allowed_changes"]
    }
    if actual_changes != expected_changes:
        raise BaselineError(
            f"{case['case']}: allowed path changes {sorted(actual_changes)} do not "
            f"match expected {sorted(expected_changes)}"
        )
    output = config_value(effective, "output", "filename")
    expected_suffix = f"{case['directory']}/generated/forecast"
    if not Path(output).as_posix().endswith(expected_suffix):
        raise BaselineError(f"{case['case']}: output path is outside its captured run")


def linear_grid(start: float, stop: float, count: int, location: str) -> list[float]:
    """Reproduce NumPy linspace for the supported reference configurations."""
    if count < 2 or not start < stop:
        raise BaselineError(
            f"{location}: unsupported grid requires count >= 2 and start < stop"
        )
    step = (stop - start) / (count - 1)
    values = [start + index * step for index in range(count)]
    values[-1] = stop
    return values


def assert_coordinates(
    actual: list[float], expected: list[float], location: str
) -> None:
    """Compare every stored coordinate with tight roundoff tolerance."""
    if len(actual) != len(expected):
        raise BaselineError(
            f"{location}: coordinate count {len(actual)} != {len(expected)}"
        )
    if any(left >= right for left, right in zip(actual, actual[1:])):
        raise BaselineError(f"{location}: coordinates are not strictly ordered")
    for index, (value, reference) in enumerate(zip(actual, expected)):
        if not math.isclose(
            value,
            reference,
            rel_tol=COORDINATE_RTOL,
            abs_tol=COORDINATE_ATOL,
        ):
            raise BaselineError(
                f"{location}: coordinate {index} is {value}, expected {reference}"
            )


def derived_configuration_identity(
    parser: configparser.ConfigParser, case_name: str
) -> tuple[list[str], list[str], list[str]]:
    """Derive tracer and pair identities using NewForecast's naming/order rules."""
    tracers = []
    for section in parser.sections():
        if not section.startswith("tracer "):
            continue
        simple_name = parser[section].get("tracer")
        if not simple_name:
            raise BaselineError(f"{case_name}: [{section}] has no tracer name")
        background = parser[section].get("background_tracer")
        tracers.append(f"{simple_name}({background})" if background else simple_name)
    if len(set(tracers)) != len(tracers):
        raise BaselineError(f"{case_name}: duplicate derived tracer identities")
    all_pairs = [
        f"{left}_{right}"
        for index, left in enumerate(tracers)
        for right in tracers[index:]
    ]
    configured = parser["control"].get("correlations", "all").split()
    if configured == ["all"]:
        selected = all_pairs
    else:
        accepted = set(configured)
        for pair in configured:
            parts = pair.split("_")
            if len(parts) != 2:
                raise BaselineError(f"{case_name}: invalid configured pair {pair!r}")
            if parts[0] != parts[1]:
                accepted.add(f"{parts[1]}_{parts[0]}")
        selected = [pair for pair in all_pairs if pair in accepted]
    return tracers, all_pairs, selected


def expected_redshifts(
    parser: configparser.ConfigParser, case_name: str
) -> tuple[list[float], list[float], list[float]]:
    """Reproduce Survey's redshift convention used by all seven DESI-2 INIs."""
    survey = parser["survey"]
    if survey.get("z bin edges") is not None or survey.get("z bin centres") is not None:
        raise BaselineError(
            f"{case_name}: unsupported explicit redshift-bin convention in checker"
        )
    count = survey.getint("num z bins", 1)
    edges = linear_grid(
        survey.getfloat("z bin min", 2.0),
        survey.getfloat("z bin max", 4.0),
        count + 1,
        f"{case_name}.redshift edges",
    )
    centers = [(left + right) / 2 for left, right in zip(edges, edges[1:])]
    edges_2d = [*edges[:-1], *edges[1:]]
    return centers, edges, edges_2d


def validate_result(case: dict[str, object], bundle: Path) -> None:
    """Inspect one stored result's identities, coordinates, and numerical values."""
    case_name = case["case"]
    metadata = read_json(bundle / case["artifacts"]["metadata"]["path"])
    result = read_json(bundle / case["artifacts"]["result_json"]["path"])
    keys, values = dict_node(result, f"{case_name} result")
    original = bundle / case["artifacts"]["original_config"]["path"]
    parser = parse_config(original)
    expected_tracers, expected_pairs, expected_selected = (
        derived_configuration_identity(parser, case_name)
    )
    tracer_order = metadata.get("tracer_order")
    all_pairs = metadata.get("all_pairs")
    selected_pairs = metadata.get("selected_pairs")
    if tracer_order != expected_tracers:
        raise BaselineError(f"{case_name}: tracer identities differ from configuration")
    if all_pairs != expected_pairs:
        raise BaselineError(f"{case_name}: pair identities differ from configuration")
    if selected_pairs != expected_selected:
        raise BaselineError(f"{case_name}: selected pairs differ from configuration")
    if len(selected_pairs) != case["expected_selected_pair_count"]:
        raise BaselineError(
            f"{case_name}: selected-pair count {len(selected_pairs)} does not match "
            f"expected {case['expected_selected_pair_count']}"
        )
    expected_keys = [
        "redshifts",
        "zedges",
        "fiducial redshift",
        *all_pairs,
        "sigma_at",
        "sigma_ap",
        "corr_coef",
    ]
    if keys != expected_keys:
        raise BaselineError(f"{case_name}: result key ordering or membership differs")

    expected_centers, expected_edges, expected_edges_2d = expected_redshifts(
        parser, case_name
    )
    nz = len(expected_centers)
    redshifts = array_node(values["redshifts"], f"{case_name}.redshifts", [nz])
    zedges = array_node(values["zedges"], f"{case_name}.zedges", [nz + 1])
    metadata_redshifts = array_node(
        metadata.get("redshift_centers"), f"{case_name}.metadata.redshifts", [nz]
    )
    metadata_edges = array_node(
        metadata.get("redshift_edges_2d"),
        f"{case_name}.metadata.redshift_edges_2d",
        [2, nz],
    )
    assert_coordinates(redshifts, expected_centers, f"{case_name}.redshifts")
    assert_coordinates(zedges, expected_edges, f"{case_name}.zedges")
    assert_coordinates(
        metadata_redshifts, expected_centers, f"{case_name}.metadata.redshifts"
    )
    for index, (value, expected) in enumerate(zip(metadata_edges, expected_edges_2d)):
        if not math.isclose(
            value,
            expected,
            rel_tol=COORDINATE_RTOL,
            abs_tol=COORDINATE_ATOL,
        ):
            raise BaselineError(
                f"{case_name}.metadata.redshift_edges_2d: coordinate {index} "
                f"is {value}, expected {expected}"
            )
    fiducial = numeric_values(
        values["fiducial redshift"], f"{case_name}.fiducial redshift"
    )
    expected_fiducial = parser["cosmo"].getfloat("z_ref")
    if len(fiducial) != 1 or not math.isclose(
        fiducial[0],
        expected_fiducial,
        rel_tol=COORDINATE_RTOL,
        abs_tol=COORDINATE_ATOL,
    ):
        raise BaselineError(
            f"{case_name}: fiducial redshift differs from configuration"
        )

    def validate_forecast(
        owner: str, forecast: dict[str, object], positive: bool
    ) -> None:
        metric_keys, metrics = dict_node(forecast, f"{case_name}.{owner}")
        if metric_keys != ["sigma_at", "sigma_ap", "corr_coef"]:
            raise BaselineError(f"{case_name}.{owner}: malformed forecast dictionary")
        sigma_at = array_node(
            metrics["sigma_at"], f"{case_name}.{owner}.sigma_at", [nz]
        )
        sigma_ap = array_node(
            metrics["sigma_ap"], f"{case_name}.{owner}.sigma_ap", [nz]
        )
        corr = array_node(metrics["corr_coef"], f"{case_name}.{owner}.corr_coef", [nz])
        if positive:
            if any(value <= 0 for value in sigma_at + sigma_ap):
                raise BaselineError(
                    f"{case_name}.{owner}: forecast sigma is not positive"
                )
            if any(abs(value) > 1 + 1e-12 for value in corr):
                raise BaselineError(f"{case_name}.{owner}: correlation exceeds unity")
        elif any(value != 0 for value in sigma_at + sigma_ap + corr):
            raise BaselineError(
                f"{case_name}.{owner}: unselected forecast is not zero-filled"
            )

    for pair in all_pairs:
        validate_forecast(pair, values[pair], pair in selected_pairs)
    combined = {
        "kind": "dict",
        "items": [
            {"key": metric, "value": values[metric]}
            for metric in ("sigma_at", "sigma_ap", "corr_coef")
        ],
    }
    validate_forecast("total", combined, True)

    grid = metadata.get("grid")
    if not isinstance(grid, dict):
        raise BaselineError(f"{case_name}: missing grid metadata")
    power = parser["power spectrum"]
    k_count = power.getint("num_k_bins")
    mu_count = power.getint("num_mu_bins")
    k_values = array_node(grid.get("k"), f"{case_name}.grid.k", [k_count])
    mu_values = array_node(grid.get("mu"), f"{case_name}.grid.mu", [mu_count])
    expected_k = linear_grid(
        power.getfloat("k_min_hmpc"),
        power.getfloat("k_max_hmpc"),
        k_count,
        f"{case_name}.grid.k",
    )
    mu_edges = linear_grid(
        power.getfloat("mu_min"),
        power.getfloat("mu_max"),
        mu_count + 1,
        f"{case_name}.grid.mu edges",
    )
    expected_mu = [
        left + (right - left) / 2 for left, right in zip(mu_edges, mu_edges[1:])
    ]
    assert_coordinates(k_values, expected_k, f"{case_name}.grid.k")
    assert_coordinates(mu_values, expected_mu, f"{case_name}.grid.mu")


def compare_result_nodes(
    primary: object,
    repeated: object,
    location: str = "result",
) -> dict[str, object]:
    """Compare two typed JSON results with fixed tolerances."""
    maximum_absolute = 0.0
    maximum_relative = 0.0
    relative_infinite = False
    count = 0

    def compare(left: object, right: object, where: str) -> None:
        nonlocal maximum_absolute, maximum_relative, relative_infinite, count
        if not isinstance(left, dict) or not isinstance(right, dict):
            raise BaselineError(f"repeat mismatch at {where}: malformed node")
        if left.get("kind") != right.get("kind"):
            raise BaselineError(f"repeat mismatch at {where}: node kinds differ")
        kind = left.get("kind")
        if kind == "dict":
            left_keys, left_values = dict_node(left, where)
            right_keys, right_values = dict_node(right, where)
            if left_keys != right_keys:
                raise BaselineError(f"repeat mismatch at {where}: key ordering differs")
            for key in left_keys:
                compare(left_values[key], right_values[key], f"{where}.{key}")
            return
        if kind == "ndarray":
            if left.get("dtype") != right.get("dtype"):
                raise BaselineError(f"repeat mismatch at {where}: dtypes differ")
            if left.get("shape") != right.get("shape"):
                raise BaselineError(f"repeat mismatch at {where}: shapes differ")
            left_numbers = numeric_values(left, where)
            right_numbers = numeric_values(right, where)
        elif kind in {"python-scalar", "numpy-scalar"}:
            if left.get("python_type") != right.get("python_type") or left.get(
                "dtype"
            ) != right.get("dtype"):
                raise BaselineError(f"repeat mismatch at {where}: scalar types differ")
            left_numbers = numeric_values(left, where)
            right_numbers = numeric_values(right, where)
        elif kind in {"string", "none", "bool"}:
            if left != right:
                raise BaselineError(f"repeat mismatch at {where}: values differ")
            return
        else:
            raise BaselineError(f"repeat mismatch at {where}: unsupported node kind")
        if len(left_numbers) != len(right_numbers):
            raise BaselineError(f"repeat mismatch at {where}: value counts differ")
        for left_value, right_value in zip(left_numbers, right_numbers):
            count += 1
            difference = abs(right_value - left_value)
            maximum_absolute = max(maximum_absolute, difference)
            if left_value == 0:
                if difference != 0:
                    relative_infinite = True
            else:
                maximum_relative = max(maximum_relative, difference / abs(left_value))
            if not math.isclose(right_value, left_value, rel_tol=1e-10, abs_tol=1e-12):
                raise BaselineError(
                    f"repeat mismatch at {where}: {right_value!r} != {left_value!r}"
                )

    compare(primary, repeated, location)
    return {
        "passed": True,
        "rtol": 1e-10,
        "atol": 1e-12,
        "number_count": count,
        "max_absolute_difference": maximum_absolute,
        "max_relative_difference": "infinity"
        if relative_infinite
        else maximum_relative,
        "zero_denominator_convention": (
            "relative difference is zero for 0 versus 0 and infinity for a nonzero "
            "repeat value versus a zero primary value"
        ),
        "structure": "identical keys, ordering, array shapes, and dtypes",
    }


def validate_artifact_record(bundle: Path, owner: str, record: object) -> None:
    """Validate one manifest artifact path and digest."""
    if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
        raise BaselineError(f"{owner}: malformed artifact record")
    path_value = record["path"]
    path = bundle_path(bundle, path_value, owner)
    if not path.is_file():
        raise BaselineError(f"{owner}: missing artifact {path_value}")
    actual = sha256(path)
    if actual != record["sha256"]:
        raise BaselineError(f"{owner}: corrupted artifact {path_value}")


def validate_preservation(
    bundle: Path, manifest: dict[str, object], schema_version: int
) -> tuple[dict[str, object], list[str]]:
    """Verify all shared snapshots and pre/post preservation records."""
    shared = manifest.get("shared_artifacts")
    if not isinstance(shared, dict):
        raise BaselineError("manifest has no shared artifacts")
    expected_shared = (
        LEGACY_SHARED_ARTIFACTS
        if schema_version == LEGACY_SCHEMA_VERSION
        else CURRENT_SHARED_ARTIFACTS
    )
    if set(shared) != expected_shared:
        raise BaselineError("shared artifact inventory differs from schema")
    for name, record in shared.items():
        validate_artifact_record(bundle, f"shared artifact {name}", record)

    tool_state = read_json(
        bundle_path(bundle, shared["tool_state"]["path"], "tool state")
    )
    tools = tool_state.get("tools") if isinstance(tool_state, dict) else None
    if not isinstance(tools, list) or len(tools) != 2:
        raise BaselineError("tool state must describe exactly two tool snapshots")
    tool_paths = []
    for index, record in enumerate(tools):
        if not isinstance(record, dict) or not {
            "source_path",
            "snapshot_path",
            "sha256",
        }.issubset(record):
            raise BaselineError(f"tool snapshot {index}: malformed record")
        snapshot = bundle_path(
            bundle, record["snapshot_path"], f"tool snapshot {index}"
        )
        if not snapshot.is_file() or sha256(snapshot) != record["sha256"]:
            raise BaselineError(f"tool snapshot {index}: corrupted snapshot")
        tool_paths.append(record["snapshot_path"])
    if len(set(tool_paths)) != 2:
        raise BaselineError("tool snapshots do not have distinct paths")
    expected_tool_files = sorted(["tool/git.diff", "tool/tool.json", *tool_paths])
    actual_tool_files = sorted(
        relative(bundle, path)
        for path in (bundle / "tool").rglob("*")
        if path.is_file()
    )
    if actual_tool_files != expected_tool_files:
        raise BaselineError("tool snapshot inventory differs")
    git_diff = bundle_path(bundle, tool_state.get("git_diff_path"), "tool Git diff")
    if tool_state.get("git_diff_sha256") != shared["tool_git_diff"][
        "sha256"
    ] or git_diff != bundle_path(
        bundle, shared["tool_git_diff"]["path"], "shared tool Git diff"
    ):
        raise BaselineError("tool Git diff metadata differs")

    source_manifest = read_json(
        bundle_path(bundle, shared["source_manifest"]["path"], "source manifest")
    )
    if not source_manifest.get("all_unchanged"):
        raise BaselineError("reference source was not recorded unchanged")
    source_records = source_manifest.get("files")
    if not isinstance(source_records, list) or not source_records:
        raise BaselineError("reference source manifest has no files")
    expected_source_paths = []
    for record in source_records:
        if (
            not record.get("unchanged")
            or record.get("before_sha256") != record.get("after_sha256")
            or record.get("before_sha256") != record.get("snapshot_sha256")
        ):
            raise BaselineError(
                f"reference source preservation failed for {record.get('reference_path')}"
            )
        snapshot = bundle_path(
            bundle, record.get("snapshot_path"), "reference source snapshot"
        )
        if not snapshot.is_file() or sha256(snapshot) != record["snapshot_sha256"]:
            raise BaselineError(f"corrupted source snapshot {record['snapshot_path']}")
        expected_source_paths.append(record.get("reference_path"))
    source_snapshot_root = bundle / "reference" / "source" / "snapshot"
    actual_source_paths = sorted(
        path.relative_to(source_snapshot_root).as_posix()
        for path in source_snapshot_root.rglob("*")
        if path.is_file()
    )
    if actual_source_paths != sorted(expected_source_paths):
        raise BaselineError("reference source snapshot inventory differs")

    input_manifest = read_json(
        bundle_path(bundle, shared["input_manifest"]["path"], "input manifest")
    )
    if not input_manifest.get("all_unchanged"):
        raise BaselineError("reference inputs were not recorded unchanged")
    input_records = input_manifest.get("inputs")
    if not isinstance(input_records, list) or not input_records:
        raise BaselineError("input manifest has no inputs")
    for record in input_records:
        snapshot_root = bundle_path(
            bundle, record.get("snapshot_path"), "input snapshot"
        )
        if record["kind"] == "directory":
            if not snapshot_root.is_dir():
                raise BaselineError(f"missing input snapshot {record['snapshot_path']}")
            expected_inventory = [
                *(
                    {"path": path, "kind": "directory"}
                    for path in record["directories"]
                ),
                *(
                    {"path": entry["path"], "kind": "file"}
                    for entry in record["entries"]
                ),
            ]
            if directory_inventory(snapshot_root) != sorted(
                expected_inventory, key=lambda item: item["path"]
            ):
                raise BaselineError(
                    f"input directory contents differ for {record['snapshot_path']}"
                )
        for entry in record.get("entries", []):
            if (
                not entry.get("unchanged")
                or entry.get("before_sha256") != entry.get("after_sha256")
                or entry.get("before_sha256") != entry.get("snapshot_sha256")
            ):
                raise BaselineError(
                    f"input preservation failed for {record['resolved_original_path']}"
                )
            snapshot = (
                snapshot_root if entry["path"] == "." else snapshot_root / entry["path"]
            )
            if not snapshot.is_file() or sha256(snapshot) != entry["snapshot_sha256"]:
                raise BaselineError(f"corrupted input snapshot {snapshot}")

    limitations = []
    if schema_version == LEGACY_SCHEMA_VERSION:
        limitations.append(
            "schema v1 records content preservation but not source, input, or "
            "configuration inventory preservation"
        )
    else:
        if (
            not source_manifest.get("inventory_unchanged")
            or source_manifest.get("inventory_before")
            != source_manifest.get("inventory_after")
            or source_manifest.get("added") != []
            or source_manifest.get("deleted") != []
            or len(source_manifest.get("inventory_before", []))
            != len(expected_source_paths)
            or set(source_manifest.get("inventory_before", []))
            != set(expected_source_paths)
        ):
            raise BaselineError("reference source inventory was not preserved")
        for record in input_records:
            if not record.get("inventory_unchanged") or record.get(
                "inventory_before"
            ) != record.get("inventory_after"):
                raise BaselineError(
                    f"input inventory was not preserved for {record['snapshot_path']}"
                )
        config_manifest = read_json(
            bundle_path(
                bundle,
                shared["configuration_manifest"]["path"],
                "configuration manifest",
            )
        )
        expected_inventory = sorted(PRIMARY_CASES)
        if (
            not config_manifest.get("all_unchanged")
            or not config_manifest.get("inventory_unchanged")
            or config_manifest.get("authoritative_inventory_before")
            != expected_inventory
            or config_manifest.get("authoritative_inventory_after")
            != expected_inventory
            or config_manifest.get("scope_cases") != SUITE_CASES[manifest["suite"]]
        ):
            raise BaselineError(
                "authoritative configuration inventory was not preserved"
            )
        config_records = config_manifest.get("files")
        if (
            not isinstance(config_records, list)
            or [record.get("case") for record in config_records]
            != SUITE_CASES[manifest["suite"]]
        ):
            raise BaselineError("configuration snapshot scope differs from suite")
        for record in config_records:
            if (
                not record.get("unchanged")
                or record.get("before_sha256") != record.get("after_sha256")
                or record.get("before_sha256") != record.get("snapshot_sha256")
            ):
                raise BaselineError(
                    f"configuration preservation failed for {record.get('case')}"
                )
            snapshot = bundle_path(
                bundle, record.get("snapshot_path"), "configuration snapshot"
            )
            if not snapshot.is_file() or sha256(snapshot) != record["snapshot_sha256"]:
                raise BaselineError("corrupted configuration snapshot")
    return input_manifest, limitations


def validate_case_evidence(
    bundle: Path,
    case: dict[str, object],
    input_manifest: dict[str, object],
    role: str,
    extra_artifacts: set[str] | None = None,
) -> None:
    """Validate common status, request, artifact, config, and result evidence."""
    case_name = case.get("case")
    if case_name not in PRIMARY_CASES:
        raise BaselineError(f"unknown case record {case_name!r}")
    if (
        case.get("role") != role
        or case.get("state") != "complete"
        or case.get("exit_code") != 0
    ):
        raise BaselineError(f"{case_name}: {role} worker did not complete successfully")
    if case.get("expected_selected_pair_count") != PRIMARY_CASES[case_name]:
        raise BaselineError(f"{case_name}: wrong expected selected-pair count")
    expected_artifacts = BASE_CASE_ARTIFACTS | (extra_artifacts or set())
    artifacts = case.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != expected_artifacts:
        raise BaselineError(f"{case_name}: required artifact inventory differs")
    artifact_paths = []
    for name, record in artifacts.items():
        validate_artifact_record(bundle, f"{case_name} {name}", record)
        artifact_paths.append(record["path"])
    if len(set(artifact_paths)) != len(artifact_paths):
        raise BaselineError(f"{case_name}: artifact paths are not distinct")
    directory = case.get("directory")
    if not isinstance(directory, str):
        raise BaselineError(f"{case_name}: malformed case directory")
    bundle_path(bundle, directory, f"{case_name} directory")
    if any(
        not Path(path).as_posix().startswith(f"{Path(directory).as_posix()}/")
        for path in artifact_paths
    ):
        raise BaselineError(f"{case_name}: artifact is outside its case directory")

    status = read_json(bundle_path(bundle, artifacts["status"]["path"], "status"))
    if (
        status.get("role") != role
        or status.get("case") != case_name
        or status.get("exit_code") != 0
        or status.get("state") != "complete"
    ):
        raise BaselineError(f"{case_name}: stored worker status is not successful")
    worker_state = status.get("worker")
    if (
        not isinstance(worker_state, dict)
        or worker_state.get("roundtrip_equal") is not True
    ):
        raise BaselineError(f"{case_name}: serialization round-trip was not verified")
    response = read_json(
        bundle_path(bundle, artifacts["worker_response"]["path"], "worker response")
    )
    if response != worker_state:
        raise BaselineError(f"{case_name}: worker response and status differ")
    if status.get("worker_response") != artifacts["worker_response"]["path"]:
        raise BaselineError(f"{case_name}: worker response path differs")
    if worker_state.get("thread_settings") != {
        variable: "1" for variable in THREAD_VARIABLES
    }:
        raise BaselineError(f"{case_name}: worker thread settings differ")

    request = read_json(
        bundle_path(bundle, artifacts["worker_request"]["path"], "worker request")
    )
    if request.get("role") != role or request.get("case") != case_name:
        raise BaselineError(f"{case_name}: worker request identity differs")
    request_artifacts = {
        "effective_config": "effective_config",
        "result_pickle": "result_pickle",
        "result_json": "result_json",
        "metadata": "metadata",
        "summary": "summary",
    }
    for request_key, artifact_key in request_artifacts.items():
        value = request.get(request_key)
        suffix = artifacts[artifact_key]["path"]
        if not isinstance(value, str) or not Path(value).as_posix().endswith(suffix):
            raise BaselineError(f"{case_name}: request path differs for {request_key}")
    validate_config_evidence(bundle, case, input_manifest)
    validate_result(case, bundle)


def suite_for_manifest(manifest: dict[str, object], schema_version: int) -> str:
    """Return the explicit v2 suite or legacy v1 full-suite interpretation."""
    if schema_version == LEGACY_SCHEMA_VERSION:
        return "full"
    suite = manifest.get("suite")
    if suite not in SUITE_CASES:
        raise BaselineError(f"unsupported capture suite {suite!r}")
    return suite


def full_case(manifest: dict[str, object], case_name: str) -> dict[str, object]:
    """Find one primary record in a validated full bundle."""
    return next(case for case in manifest["primary_cases"] if case["case"] == case_name)


def check_quick_baseline(
    bundle: Path, manifest: dict[str, object], quick_case: dict[str, object]
) -> None:
    """Validate embedded full-baseline compatibility and numerical comparison."""
    baseline = manifest.get("baseline")
    if (
        not isinstance(baseline, dict)
        or baseline.get("validated_at_capture") is not True
        or baseline.get("compatible") is not True
    ):
        raise BaselineError("quick suite has no validated full-baseline evidence")
    records = baseline.get("artifacts")
    if not isinstance(records, dict) or set(records) != {
        "manifest",
        "result",
        "compatibility",
        "comparison",
    }:
        raise BaselineError("quick baseline artifact inventory differs")
    for name, record in records.items():
        validate_artifact_record(bundle, f"quick baseline {name}", record)
    manifest_copy = bundle_path(
        bundle, records["manifest"]["path"], "baseline manifest"
    )
    if sha256(manifest_copy) != baseline.get("original_manifest_sha256"):
        raise BaselineError("embedded baseline manifest digest differs from source")
    copied_manifest = read_json(manifest_copy)
    copied_version = copied_manifest.get("schema_version")
    if (
        copied_manifest.get("schema") != SCHEMA
        or copied_version not in {LEGACY_SCHEMA_VERSION, SCHEMA_VERSION}
        or suite_for_manifest(copied_manifest, copied_version) != "full"
        or baseline.get("source_schema_version") != copied_version
    ):
        raise BaselineError("embedded baseline manifest is not a supported full suite")
    copied_case = full_case(copied_manifest, REPEAT_CASE)
    result_path = bundle_path(bundle, records["result"]["path"], "baseline result")
    if (
        sha256(result_path) != baseline.get("original_result_sha256")
        or baseline.get("original_result_sha256")
        != copied_case["artifacts"]["result_json"]["sha256"]
    ):
        raise BaselineError("embedded baseline result digest differs from source")
    compatibility = read_json(
        bundle_path(bundle, records["compatibility"]["path"], "compatibility")
    )
    if (
        not isinstance(compatibility, dict)
        or compatibility.get("compatible") is not True
        or not isinstance(compatibility.get("checks"), dict)
        or set(compatibility["checks"])
        != {"configuration", "source", "inputs", "environment"}
        or any(
            check.get("passed") is not True
            for check in compatibility["checks"].values()
            if isinstance(check, dict)
        )
        or not all(
            isinstance(check, dict) for check in compatibility["checks"].values()
        )
    ):
        raise BaselineError("quick baseline compatibility did not pass")
    comparisons = {
        "configuration": (
            compatibility["checks"]["configuration"].get("baseline_sha256"),
            compatibility["checks"]["configuration"].get("current_sha256"),
        ),
        **{
            name: (
                compatibility["checks"][name].get("baseline"),
                compatibility["checks"][name].get("current"),
            )
            for name in ("source", "inputs", "environment")
        },
    }
    if any(
        compatibility["checks"][name]["passed"] != (left == right)
        for name, (left, right) in comparisons.items()
    ):
        raise BaselineError("quick baseline compatibility flags are inconsistent")
    shared = manifest["shared_artifacts"]
    current_source = read_json(
        bundle_path(bundle, shared["source_manifest"]["path"], "source manifest")
    )
    current_inputs = read_json(
        bundle_path(bundle, shared["input_manifest"]["path"], "input manifest")
    )
    current_environment = read_json(
        bundle_path(bundle, shared["environment"]["path"], "environment")
    )
    expected_current = {
        "configuration": sha256(
            bundle_path(
                bundle,
                quick_case["artifacts"]["original_config"]["path"],
                "quick original configuration",
            )
        ),
        "source": [
            {"path": record["reference_path"], "sha256": record["before_sha256"]}
            for record in current_source["files"]
        ],
        "inputs": input_signature(current_inputs, REPEAT_CASE),
        "environment": environment_signature(current_environment),
    }
    recorded_current = {
        "configuration": compatibility["checks"]["configuration"]["current_sha256"],
        **{
            name: compatibility["checks"][name]["current"]
            for name in ("source", "inputs", "environment")
        },
    }
    if recorded_current != expected_current:
        raise BaselineError(
            "quick compatibility evidence differs from captured evidence"
        )
    if (
        compatibility["checks"]["configuration"]["baseline_sha256"]
        != copied_case["artifacts"]["original_config"]["sha256"]
    ):
        raise BaselineError(
            "quick baseline configuration identity differs from manifest"
        )
    comparison = read_json(
        bundle_path(bundle, records["comparison"]["path"], "quick comparison")
    )
    if comparison.get("passed") is not True:
        raise BaselineError("quick baseline numerical comparison did not pass")
    fresh = read_json(
        bundle_path(
            bundle, quick_case["artifacts"]["result_json"]["path"], "quick result"
        )
    )
    copied = read_json(result_path)
    if comparison != compare_result_nodes(copied, fresh, "quick versus full baseline"):
        raise BaselineError(
            "quick baseline comparison artifact is stale or inconsistent"
        )


def check_bundle(
    bundle: Path,
    require_complete: bool = True,
    require_suite: str | None = None,
) -> dict[str, object]:
    """Validate a captured bundle without the reference checkout or NumPy."""
    bundle = bundle.resolve()
    manifest = read_json(bundle / "manifest.json")
    if not isinstance(manifest, dict):
        raise BaselineError("manifest is not a JSON object")
    schema_version = manifest.get("schema_version")
    if manifest.get("schema") != SCHEMA or schema_version not in {
        LEGACY_SCHEMA_VERSION,
        SCHEMA_VERSION,
    }:
        raise BaselineError("unsupported baseline manifest schema")
    if require_complete and manifest.get("state") != "complete":
        raise BaselineError(f"bundle state is {manifest.get('state')!r}, not complete")
    suite = suite_for_manifest(manifest, schema_version)
    if require_suite is not None and suite != require_suite:
        raise BaselineError(
            f"bundle suite is {suite!r}, not required {require_suite!r}"
        )
    expected = SUITE_CASES[suite]
    if manifest.get("expected_primary_cases") != expected:
        raise BaselineError("manifest primary inventory differs from suite")
    cases = manifest.get("primary_cases")
    if not isinstance(cases, list):
        raise BaselineError("manifest primary cases are malformed")
    names = [case.get("case") for case in cases if isinstance(case, dict)]
    if (
        len(cases) != len(expected)
        or len(set(names)) != len(names)
        or names != expected
    ):
        raise BaselineError(
            f"expected exactly {len(expected)} unique ordered {suite} cases; found {names}"
        )
    if schema_version == SCHEMA_VERSION:
        expected_coverage = {
            "suite": suite,
            "primary_cases": expected,
            "primary_roles": ["quick" if suite == "quick" else "primary"]
            * len(expected),
            "repeat_case": REPEAT_CASE if suite == "full" else None,
            "scientific_run_count": 8 if suite == "full" else 1,
        }
        if manifest.get("coverage") != expected_coverage:
            raise BaselineError("manifest coverage differs from suite contract")
    input_manifest, limitations = validate_preservation(
        bundle, manifest, schema_version
    )
    role = "quick" if suite == "quick" else "primary"
    for case in cases:
        validate_case_evidence(bundle, case, input_manifest, role)

    repeated = manifest.get("repeat")
    if suite == "full":
        if schema_version == SCHEMA_VERSION and manifest.get("baseline") is not None:
            raise BaselineError(
                "full suite must not contain baseline comparison evidence"
            )
        if not isinstance(repeated, dict) or repeated.get("case") != REPEAT_CASE:
            raise BaselineError("missing separate 15x2pt repeat record")
        validate_case_evidence(
            bundle, repeated, input_manifest, "repeat", {"comparison"}
        )
        primary = full_case(manifest, REPEAT_CASE)
        if primary["directory"] == repeated["directory"]:
            raise BaselineError(
                "15x2pt primary and repeat directories are not distinct"
            )
        if set(record["path"] for record in primary["artifacts"].values()) & set(
            record["path"] for record in repeated["artifacts"].values()
        ):
            raise BaselineError("15x2pt primary and repeat artifacts overlap")
        comparison = read_json(
            bundle_path(
                bundle, repeated["artifacts"]["comparison"]["path"], "comparison"
            )
        )
        if comparison.get("passed") is not True:
            raise BaselineError("stored repeat comparison did not pass")
        primary_result = read_json(
            bundle_path(
                bundle, primary["artifacts"]["result_json"]["path"], "primary result"
            )
        )
        repeat_result = read_json(
            bundle_path(
                bundle, repeated["artifacts"]["result_json"]["path"], "repeat result"
            )
        )
        if comparison != compare_result_nodes(primary_result, repeat_result):
            raise BaselineError("repeat comparison artifact is stale or inconsistent")
    else:
        if repeated is not None:
            raise BaselineError("quick suite must not contain a repeat record")
        check_quick_baseline(bundle, manifest, cases[0])

    return {
        "bundle": str(bundle),
        "schema_version": schema_version,
        "suite": suite,
        "primary_cases": len(cases),
        "repeat_case": repeated["case"] if repeated is not None else None,
        "limitations": limitations,
        "state": "valid",
    }


def shared_artifacts(bundle: Path) -> dict[str, dict[str, str]]:
    """Hash the maintained shared evidence files."""
    paths = {
        "environment": bundle / "environment.json",
        "probe_stdout": bundle / "reference" / "probe.stdout.log",
        "probe_stderr": bundle / "reference" / "probe.stderr.log",
        "reference_git": bundle / "reference" / "git.json",
        "reference_dirty_diff": bundle / "reference" / "dirty.diff",
        "source_manifest": bundle / "reference" / "source" / "source_manifest.json",
        "input_manifest": bundle / "inputs" / "manifest.json",
        "tool_state": bundle / "tool" / "tool.json",
        "tool_git_diff": bundle / "tool" / "git.diff",
        "configuration_manifest": (
            bundle / "reference" / "configurations" / "configuration_manifest.json"
        ),
    }
    return {name: artifact(bundle, path) for name, path in paths.items()}


def input_signature(
    input_manifest: dict[str, object], case_name: str
) -> list[dict[str, object]]:
    """Return location-independent input bindings and content identities."""
    by_id = {record["id"]: record for record in input_manifest["inputs"]}
    signature = []
    for resolution in input_manifest["case_resolutions"][case_name]["paths"]:
        record = by_id[resolution["input_id"]]
        signature.append(
            {
                "section": resolution["section"],
                "option": resolution["option"],
                "kind": resolution["kind"],
                "original_setting": resolution["original_setting"],
                "entries": [
                    {"path": entry["path"], "sha256": entry["before_sha256"]}
                    for entry in record["entries"]
                ],
                "directories": record["directories"],
            }
        )
    return signature


def environment_signature(environment: dict[str, object]) -> dict[str, object]:
    """Return scientific environment fields, excluding location-only paths."""
    return {
        "python_version": environment.get("python_version"),
        "platform": environment.get("platform"),
        "distribution_versions": environment.get("distribution_versions"),
        "thread_settings": environment.get("thread_settings"),
        "dont_write_bytecode": environment.get("dont_write_bytecode"),
    }


def prepare_baseline_evidence(
    bundle: Path,
    baseline_bundle: Path,
    baseline_manifest: dict[str, object],
    config_manifest: dict[str, object],
    source_manifest: dict[str, object],
    input_manifest: dict[str, object],
    environment: dict[str, object],
) -> dict[str, object]:
    """Copy a full baseline and record all pre-run compatibility decisions."""
    destination = bundle / "baseline"
    destination.mkdir()
    baseline_case = full_case(baseline_manifest, REPEAT_CASE)
    manifest_source = baseline_bundle / "manifest.json"
    result_source = bundle_path(
        baseline_bundle,
        baseline_case["artifacts"]["result_json"]["path"],
        "full baseline result",
    )
    manifest_copy = destination / "manifest.json"
    result_copy = destination / "reference-result.json"
    shutil.copyfile(manifest_source, manifest_copy)
    shutil.copyfile(result_source, result_copy)

    baseline_shared = baseline_manifest["shared_artifacts"]
    baseline_source = read_json(
        bundle_path(
            baseline_bundle,
            baseline_shared["source_manifest"]["path"],
            "full baseline source manifest",
        )
    )
    baseline_inputs = read_json(
        bundle_path(
            baseline_bundle,
            baseline_shared["input_manifest"]["path"],
            "full baseline input manifest",
        )
    )
    baseline_environment = read_json(
        bundle_path(
            baseline_bundle,
            baseline_shared["environment"]["path"],
            "full baseline environment",
        )
    )
    baseline_source_hashes = [
        {"path": record["reference_path"], "sha256": record["before_sha256"]}
        for record in baseline_source["files"]
    ]
    current_source_hashes = [
        {"path": record["reference_path"], "sha256": record["before_sha256"]}
        for record in source_manifest["files"]
    ]
    current_config_hash = next(
        record["before_sha256"]
        for record in config_manifest["files"]
        if record["case"] == REPEAT_CASE
    )
    baseline_config_hash = baseline_case["artifacts"]["original_config"]["sha256"]
    baseline_input_signature = input_signature(baseline_inputs, REPEAT_CASE)
    current_input_signature = input_signature(input_manifest, REPEAT_CASE)
    baseline_environment_signature = environment_signature(baseline_environment)
    current_environment_signature = environment_signature(environment)
    checks = {
        "configuration": {
            "baseline_sha256": baseline_config_hash,
            "current_sha256": current_config_hash,
            "passed": baseline_config_hash == current_config_hash,
        },
        "source": {
            "baseline": baseline_source_hashes,
            "current": current_source_hashes,
            "passed": baseline_source_hashes == current_source_hashes,
        },
        "inputs": {
            "baseline": baseline_input_signature,
            "current": current_input_signature,
            "passed": baseline_input_signature == current_input_signature,
        },
        "environment": {
            "baseline": baseline_environment_signature,
            "current": current_environment_signature,
            "passed": baseline_environment_signature == current_environment_signature,
        },
    }
    compatibility = {
        "compatible": all(check["passed"] for check in checks.values()),
        "case": REPEAT_CASE,
        "path_policy": "absolute locations are excluded; logical bindings and hashes match",
        "checks": checks,
    }
    compatibility_path = destination / "compatibility.json"
    write_json(compatibility_path, compatibility)
    baseline_record = {
        "source_path": str(baseline_bundle),
        "validated_at_capture": True,
        "compatible": compatibility["compatible"],
        "source_schema_version": baseline_manifest["schema_version"],
        "original_manifest_sha256": sha256(manifest_source),
        "original_result_sha256": sha256(result_source),
        "artifacts": {
            "manifest": artifact(bundle, manifest_copy),
            "result": artifact(bundle, result_copy),
            "compatibility": artifact(bundle, compatibility_path),
        },
    }
    return baseline_record


def create_destination(requested: Path | None) -> Path:
    """Create a unique bundle destination and refuse every existing target."""
    if requested is None:
        parent = ROOT / ".validation" / "baseline"
        parent.mkdir(parents=True, exist_ok=True)
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        requested = parent / f"{stamp}-{uuid.uuid4().hex[:8]}"
    destination = requested.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        destination.mkdir(exist_ok=False)
    except FileExistsError as error:
        raise BaselineError(
            f"refusing existing capture destination: {destination}"
        ) from error
    return destination


def capture(args: argparse.Namespace) -> Path:
    """Perform a fresh quick or explicitly requested full reference capture."""
    reference = args.reference_checkout.expanduser().resolve()
    # Do not resolve a virtual-environment Python symlink: its lexical location
    # selects the pyvenv.cfg and therefore the requested scientific environment.
    python = Path(os.path.abspath(args.python.expanduser()))
    worker = args.worker_script.expanduser().resolve()
    if not (reference / "lyaforecast" / "forecast_new.py").is_file():
        raise BaselineError(f"not a lyaforecast source checkout: {reference}")
    if not python.is_file() or not os.access(python, os.X_OK):
        raise BaselineError(f"reference Python is not executable: {python}")
    if not worker.is_file():
        raise BaselineError(f"worker script does not exist: {worker}")
    suite = args.suite
    baseline_bundle = None
    baseline_manifest = None
    if suite == "quick":
        if args.baseline is None:
            raise BaselineError(
                "quick suite requires --baseline with a valid full bundle"
            )
        baseline_bundle = args.baseline.expanduser().resolve()
        check_bundle(baseline_bundle, require_suite="full")
        baseline_manifest = read_json(baseline_bundle / "manifest.json")
    elif args.baseline is not None:
        raise BaselineError("--baseline is only valid with the quick suite")

    bundle = create_destination(args.destination)
    case_names = SUITE_CASES[suite]
    manifest = {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "plan": "Step 02, revision 3",
        "created_utc": utc_now(),
        "state": "capturing",
        "suite": suite,
        "coverage": {
            "suite": suite,
            "primary_cases": case_names,
            "primary_roles": ["quick" if suite == "quick" else "primary"]
            * len(case_names),
            "repeat_case": REPEAT_CASE if suite == "full" else None,
            "scientific_run_count": 8 if suite == "full" else 1,
        },
        "expected_primary_cases": case_names,
        "primary_cases": [],
        "repeat": None,
        "baseline": None,
        "shared_artifacts": {},
        "failure": None,
    }
    write_json(bundle / "manifest.json", manifest)

    try:
        examples = reference / "examples" / "desi2"
        actual = sorted(path.name for path in examples.glob("*.ini"))
        expected_sorted = sorted(PRIMARY_CASES)
        if actual != expected_sorted:
            missing = sorted(set(expected_sorted) - set(actual))
            unexpected = sorted(set(actual) - set(expected_sorted))
            raise BaselineError(
                f"authoritative INI inventory differs; missing={missing}, "
                f"unexpected={unexpected}"
            )
        originals, config_manifest = snapshot_configurations(
            bundle, examples, case_names
        )
        reference_dir = bundle / "reference"
        source_dir = reference_dir / "source"
        source_dir.mkdir(parents=True)
        reference_git_state(reference, reference_dir)
        source_manifest = snapshot_source(reference, source_dir)
        capture_tool_state(bundle, worker, sys.argv)

        probe_response = bundle / "environment.json"
        probe_command = worker_command(
            python, worker, "probe", reference, probe_response
        )
        probe = run_worker(
            probe_command,
            reference_dir,
            reference_dir / "probe.stdout.log",
            reference_dir / "probe.stderr.log",
        )
        if probe.returncode != 0:
            raise BaselineError(
                "reference-environment probe failed; see reference/probe.stderr.log"
            )
        environment = read_json(probe_response)
        environment["requested_python"] = str(args.python)
        environment["resolved_python"] = str(python)
        environment["probe_command"] = probe_command
        write_json(probe_response, environment)

        input_manifest = resolve_inputs(bundle, reference, python, worker, originals)
        if suite == "quick":
            manifest["baseline"] = prepare_baseline_evidence(
                bundle,
                baseline_bundle,
                baseline_manifest,
                config_manifest,
                source_manifest,
                input_manifest,
                environment,
            )
            write_json(bundle / "manifest.json", manifest)
            if not manifest["baseline"]["compatible"]:
                compatibility = read_json(
                    bundle / manifest["baseline"]["artifacts"]["compatibility"]["path"]
                )
                failed = [
                    name
                    for name, check in compatibility["checks"].items()
                    if not check["passed"]
                ]
                finalize_configuration_hashes(bundle, examples, config_manifest)
                finalize_source_hashes(reference, source_dir, source_manifest)
                finalize_input_hashes(bundle, input_manifest)
                manifest["shared_artifacts"] = shared_artifacts(bundle)
                write_json(bundle / "manifest.json", manifest)
                raise BaselineError(
                    "quick suite is incompatible with full baseline: "
                    + ", ".join(failed)
                )

        for case_name, original in originals.items():
            case_dir = bundle / "cases" / Path(case_name).stem
            record = run_case(
                bundle,
                reference,
                python,
                worker,
                case_name,
                original,
                input_manifest["case_resolutions"][case_name],
                "quick" if suite == "quick" else "primary",
                case_dir,
            )
            manifest["primary_cases"].append(record)
            write_json(bundle / "manifest.json", manifest)

        primary_15 = manifest["primary_cases"][-1]
        if suite == "full":
            repeat_dir = bundle / "repeat" / Path(REPEAT_CASE).stem
            repeated = run_case(
                bundle,
                reference,
                python,
                worker,
                REPEAT_CASE,
                originals[REPEAT_CASE],
                input_manifest["case_resolutions"][REPEAT_CASE],
                "repeat",
                repeat_dir,
            )
            if primary_15["state"] == "complete" and repeated["state"] == "complete":
                primary_result = read_json(
                    bundle / primary_15["artifacts"]["result_json"]["path"]
                )
                repeat_result = read_json(
                    bundle / repeated["artifacts"]["result_json"]["path"]
                )
                comparison_path = repeat_dir / "comparison.json"
                try:
                    comparison = compare_result_nodes(primary_result, repeat_result)
                except BaselineError as error:
                    comparison = {
                        "passed": False,
                        "rtol": 1e-10,
                        "atol": 1e-12,
                        "error": str(error),
                        "zero_denominator_convention": (
                            "relative difference is zero for 0 versus 0 and infinity "
                            "for a nonzero repeat value versus a zero primary value"
                        ),
                    }
                write_json(comparison_path, comparison)
                repeated["artifacts"]["comparison"] = artifact(bundle, comparison_path)
            manifest["repeat"] = repeated
        elif primary_15["state"] == "complete":
            comparison_path = bundle / "baseline" / "comparison.json"
            baseline_result = read_json(
                bundle / manifest["baseline"]["artifacts"]["result"]["path"]
            )
            quick_result = read_json(
                bundle / primary_15["artifacts"]["result_json"]["path"]
            )
            try:
                comparison = compare_result_nodes(
                    baseline_result, quick_result, "quick versus full baseline"
                )
            except BaselineError as error:
                comparison = {
                    "passed": False,
                    "rtol": 1e-10,
                    "atol": 1e-12,
                    "error": str(error),
                    "zero_denominator_convention": (
                        "relative difference is zero for 0 versus 0 and infinity for "
                        "a nonzero quick value versus a zero baseline value"
                    ),
                }
            write_json(comparison_path, comparison)
            manifest["baseline"]["artifacts"]["comparison"] = artifact(
                bundle, comparison_path
            )

        finalize_configuration_hashes(bundle, examples, config_manifest)
        finalize_source_hashes(reference, source_dir, source_manifest)
        finalize_input_hashes(bundle, input_manifest)
        manifest["shared_artifacts"] = shared_artifacts(bundle)
        manifest["state"] = "complete"
        write_json(bundle / "manifest.json", manifest)
        try:
            check_bundle(bundle, require_suite=suite)
        except BaselineError as error:
            manifest["state"] = "incomplete"
            manifest["failure"] = str(error)
            write_json(bundle / "manifest.json", manifest)
            raise
    except Exception as error:
        current = read_json(bundle / "manifest.json")
        if isinstance(current, dict):
            current["state"] = "incomplete"
            current["failure"] = str(error)
            write_json(bundle / "manifest.json", current)
        (bundle / "capture-error.txt").write_text(f"{type(error).__name__}: {error}\n")
        raise
    return bundle


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    capture_parser = subparsers.add_parser(
        "capture", help="capture one quick case (default) or an explicit full suite"
    )
    capture_parser.add_argument(
        "--suite",
        choices=tuple(SUITE_CASES),
        default="quick",
        help="capture coverage; full is explicit opt-in (default: quick)",
    )
    capture_parser.add_argument(
        "--baseline",
        type=Path,
        help="validated full bundle required by the quick suite",
    )
    capture_parser.add_argument(
        "--reference-checkout", required=True, type=Path, help="lyaforecast checkout"
    )
    capture_parser.add_argument(
        "--python", required=True, type=Path, help="reference Python interpreter"
    )
    capture_parser.add_argument(
        "--destination", type=Path, help="new bundle path (must not exist)"
    )
    capture_parser.add_argument(
        "--worker-script",
        type=Path,
        default=DEFAULT_WORKER,
        help=argparse.SUPPRESS,
    )
    check_parser = subparsers.add_parser(
        "check", help="validate an existing bundle without rerunning forecasts"
    )
    check_parser.add_argument("bundle", type=Path)
    check_parser.add_argument(
        "--require-suite",
        choices=tuple(SUITE_CASES),
        help="reject a valid bundle with different coverage",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the capture/check CLI."""
    args = build_parser().parse_args(argv)
    try:
        if args.command == "capture":
            bundle = capture(args)
            print(f"Captured and validated {args.suite} baseline: {bundle}")
        else:
            result = check_bundle(args.bundle, require_suite=args.require_suite)
            detail = (
                f"{result['primary_cases']} primary cases and separate "
                f"{result['repeat_case']} repeat"
                if result["suite"] == "full"
                else "one 15x2pt quick case against embedded full-baseline evidence"
            )
            print(
                f"PASS: {result['suite']} suite, {detail} validated in {result['bundle']}"
            )
            for limitation in result["limitations"]:
                print(f"LIMITATION: {limitation}")
    except BaselineError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
