"""Portable tests for the lyaforecast baseline capture/check tool."""

from __future__ import annotations

import configparser
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL_PATH = ROOT / "scripts" / "lyaforecast_baseline.py"
STUB_WORKER = Path(__file__).parent / "fixtures" / "stub_lyaforecast_worker.py"
SPEC = importlib.util.spec_from_file_location("lyaforecast_baseline", TOOL_PATH)
assert SPEC is not None and SPEC.loader is not None
TOOL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TOOL)


def digest(path: Path) -> str:
    """Hash one test artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    """Write one strict test JSON file."""
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_tool(arguments, cwd, environment=None):
    """Run the public CLI from an unrelated directory."""
    env = os.environ.copy()
    if environment:
        env.update(environment)
    return subprocess.run(
        [sys.executable, str(TOOL_PATH), *map(str, arguments)],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def config_text(case_name: str, pair_names: list[str], selected_count: int) -> str:
    """Build a full-size-inventory but numerically synthetic reference INI."""
    correlations = (
        "all"
        if selected_count == len(pair_names)
        else " ".join(pair_names[:selected_count])
    )
    if case_name.startswith("lbg_lae"):
        tracers = [("lbg", None), ("lae", None)]
    elif case_name.startswith("lya_lbg"):
        tracers = [("lya", "lbg"), ("lbg", None), ("lae", None)]
    elif case_name == "lya_qso_2x2pt.ini":
        tracers = [("lya", "qso"), ("qso", None)]
    else:
        tracers = [
            ("lya", "qso"),
            ("qso", None),
            ("lbg", None),
            ("lae", None),
            ("lya", "lbg"),
        ]
    sections = [
        "[control]",
        "measurement type = bao",
        f"correlations = {correlations}",
        "",
        "[cosmo]",
        "filename = Planck18.ini",
        "z_ref = 2.3",
        "",
        "[output]",
        "overwrite = False",
        f"filename = old/{Path(case_name).stem}/forecast",
        "",
        "[survey]",
        "z bin min = 2.0",
        "z bin max = 2.4",
        "num z bins = 2",
        "",
        "[power spectrum]",
        "k_min_hmpc = 0.01",
        "k_max_hmpc = 0.5",
        "num_k_bins = 4",
        "mu_min = 0",
        "mu_max = 1",
        "num_mu_bins = 4",
        "linear power = True",
    ]
    for index, (tracer, background) in enumerate(tracers, start=1):
        sections.extend(
            [
                "",
                f"[tracer {index}]",
                f"dn dz = density-{index}.dat",
                f"tracer = {tracer}",
                "tracer_type = discrete",
                "target density = 1",
            ]
        )
        if background:
            sections.extend(["snr-file-dir = SNR", f"background_tracer = {background}"])
    sections.extend(["", "[biasing]", ""])
    return "\n".join(sections)


def make_reference(path: Path) -> None:
    """Create a standalone synthetic checkout with the exact seven-case inventory."""
    examples = path / "examples" / "desi2"
    package = path / "lyaforecast"
    resources = path / "resources"
    examples.mkdir(parents=True)
    package.mkdir()
    resources.mkdir()
    (package / "forecast_new.py").write_text("class NewForecast:\n    pass\n")
    (package / "utils.py").write_text("# fixture resolver marker\n")
    egg_info = path / "lyaforecast.egg-info"
    egg_info.mkdir()
    (egg_info / "PKG-INFO").write_text("Name: lyaforecast\nVersion: 0\n")
    (path / "pyproject.toml").write_text("[project]\nname='lyaforecast'\nversion='0'\n")
    (resources / "Planck18.ini").write_text("hubble = 67.36\n")
    (resources / "SNR").mkdir()
    (resources / "SNR" / "header-and-data.dat").write_text("# fixture\n1 2\n")
    for index in range(1, 6):
        (resources / f"density-{index}.dat").write_text(f"2.0 {index}.0 1.0\n")

    for case_name, selected_count in TOOL.PRIMARY_CASES.items():
        if case_name.startswith("lbg_lae"):
            pairs = ["lbg_lbg", "lbg_lae", "lae_lae"]
        elif case_name.startswith("lya_lbg"):
            pairs = [
                "lya(lbg)_lya(lbg)",
                "lya(lbg)_lbg",
                "lya(lbg)_lae",
                "lbg_lbg",
                "lbg_lae",
                "lae_lae",
            ]
        elif case_name == "lya_qso_2x2pt.ini":
            pairs = ["lya(qso)_lya(qso)", "lya(qso)_qso", "qso_qso"]
        else:
            tracers = ["lya(qso)", "qso", "lbg", "lae", "lya(lbg)"]
            pairs = [
                f"{left}_{right}"
                for index, left in enumerate(tracers)
                for right in tracers[index:]
            ]
        (examples / case_name).write_text(config_text(case_name, pairs, selected_count))

    subprocess.run(["git", "init", "-q", str(path)], check=True)
    subprocess.run(["git", "-C", str(path), "add", "."], check=True)
    environment = os.environ.copy()
    environment.update(
        {
            "GIT_AUTHOR_NAME": "FishHighz Test",
            "GIT_AUTHOR_EMAIL": "test@example.invalid",
            "GIT_COMMITTER_NAME": "FishHighz Test",
            "GIT_COMMITTER_EMAIL": "test@example.invalid",
        }
    )
    subprocess.run(
        [
            "git",
            "-C",
            str(path),
            "-c",
            "commit.gpgsign=false",
            "commit",
            "-q",
            "-m",
            "fixture",
        ],
        env=environment,
        check=True,
    )


def capture_stub(
    base: Path,
    destination_name: str,
    environment: dict[str, str] | None = None,
    suite: str | None = "full",
    baseline: Path | None = None,
    reference_name: str = "reference",
) -> tuple[Path, subprocess.CompletedProcess[str]]:
    """Capture one complete synthetic bundle with real subprocess isolation."""
    reference = base / reference_name
    if not reference.exists():
        make_reference(reference)
    outside = base / "outside"
    outside.mkdir(exist_ok=True)
    destination = base / destination_name
    arguments = [
        "capture",
        "--reference-checkout",
        reference,
        "--python",
        sys.executable,
        "--worker-script",
        STUB_WORKER,
        "--destination",
        destination,
    ]
    if suite is not None:
        arguments[1:1] = ["--suite", suite]
    if baseline is not None:
        arguments.extend(["--baseline", baseline])
    result = run_tool(
        arguments,
        outside,
        environment,
    )
    return destination, result


@pytest.fixture(scope="module")
def captured(tmp_path_factory):
    """Return one valid synthetic bundle created from outside the checkout."""
    base = tmp_path_factory.mktemp("baseline-tool")
    action_log = base / "full-actions.jsonl"
    bundle, completed = capture_stub(
        base, "valid", {"FISHHIGHZ_STUB_ACTION_LOG": str(action_log)}
    )
    assert completed.returncode == 0, completed.stderr
    return base, bundle


@pytest.fixture(scope="module")
def quick_captured(captured):
    """Return one valid default-suite quick bundle and its worker action log."""
    base, full = captured
    action_log = base / "quick-actions.jsonl"
    bundle, completed = capture_stub(
        base,
        "quick-valid",
        {"FISHHIGHZ_STUB_ACTION_LOG": str(action_log)},
        suite=None,
        baseline=full,
    )
    assert completed.returncode == 0, completed.stderr
    return bundle, action_log


def mutable_bundle(captured, tmp_path, name="bundle"):
    """Copy the valid ignored evidence for one destructive checker test."""
    _, source = captured
    target = tmp_path / name
    shutil.copytree(source, target)
    return target


def manifest_for(bundle: Path) -> dict[str, object]:
    """Load a bundle manifest."""
    return json.loads((bundle / "manifest.json").read_text())


def update_artifact(bundle: Path, owner: dict[str, object], artifact_name: str) -> None:
    """Refresh one intentionally modified artifact digest in the manifest."""
    record = owner["artifacts"][artifact_name]
    record["sha256"] = digest(bundle / record["path"])


def update_shared(bundle: Path, manifest: dict[str, object], name: str) -> None:
    """Refresh one intentionally modified shared artifact digest."""
    record = manifest["shared_artifacts"][name]
    record["sha256"] = digest(bundle / record["path"])


def actions(path: Path) -> list[dict[str, object]]:
    """Read subprocess actions recorded by the portable worker."""
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_valid_bundle_and_relocated_cli_check(captured, tmp_path):
    """A valid bundle, including zero-filled pairs, survives relocation."""
    _, bundle = captured
    checked = TOOL.check_bundle(bundle)
    assert checked["primary_cases"] == 7
    manifest = manifest_for(bundle)
    first_status = json.loads(
        (
            bundle / manifest["primary_cases"][0]["artifacts"]["status"]["path"]
        ).read_text()
    )
    assert Path(first_status["command"][0]).resolve() == Path(sys.executable).resolve()
    relocated = tmp_path / "relocated"
    shutil.copytree(bundle, relocated)
    outside = tmp_path / "unrelated-working-directory"
    outside.mkdir()
    completed = run_tool(["check", relocated], outside)
    assert completed.returncode == 0, completed.stderr
    assert "7 primary cases" in completed.stdout


@pytest.mark.parametrize("mode", ["missing", "duplicate"])
def test_missing_and_duplicate_primary_cases(captured, tmp_path, mode):
    """Primary inventory must contain exactly seven unique cases."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    if mode == "missing":
        manifest["primary_cases"].pop()
    else:
        manifest["primary_cases"][-1] = manifest["primary_cases"][0]
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="7 unique ordered full"):
        TOOL.check_bundle(bundle)


def test_failed_worker_leaves_incomplete_inspectable_bundle(tmp_path):
    """A subprocess failure is retained and makes capture exit nonzero."""
    bundle, completed = capture_stub(
        tmp_path,
        "failed",
        {"FISHHIGHZ_STUB_FAIL_CASE": "lya_qso_2x2pt.ini"},
    )
    assert completed.returncode != 0
    manifest = manifest_for(bundle)
    assert manifest["state"] == "incomplete"
    failed = next(
        case
        for case in manifest["primary_cases"]
        if case["case"] == "lya_qso_2x2pt.ini"
    )
    assert failed["state"] == "failed"
    assert (bundle / failed["artifacts"]["stderr"]["path"]).is_file()


@pytest.mark.parametrize("mode", ["missing", "corrupted"])
def test_missing_and_corrupted_artifacts(captured, tmp_path, mode):
    """All recorded artifacts must exist with their captured digest."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    record = manifest["primary_cases"][0]["artifacts"]["summary"]
    path = bundle / record["path"]
    if mode == "missing":
        path.unlink()
        match = "missing artifact"
    else:
        path.write_text("corrupted\n")
        match = "corrupted artifact"
    with pytest.raises(TOOL.BaselineError, match=match):
        TOOL.check_bundle(bundle)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [("malformed", "shape"), ("nonfinite", "non-finite")],
)
def test_malformed_and_nonfinite_results(captured, tmp_path, mutation, message):
    """The checker inspects numerical JSON rather than trusting worker status."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    case = manifest["primary_cases"][0]
    result_path = bundle / case["artifacts"]["result_json"]["path"]
    result = json.loads(result_path.read_text())
    redshifts = result["items"][0]["value"]
    if mutation == "malformed":
        redshifts["shape"] = [3]
    else:
        redshifts["data"][0] = "nan"
    write_json(result_path, result)
    update_artifact(bundle, case, "result_json")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match=message):
        TOOL.check_bundle(bundle)


def test_unauthorized_scientific_configuration_change(captured, tmp_path):
    """A changed scientific setting fails even with refreshed artifact hashes."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    case = manifest["primary_cases"][0]
    original = bundle / case["artifacts"]["original_config"]["path"]
    effective = bundle / case["artifacts"]["effective_config"]["path"]
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser.read(effective)
    parser["survey"]["num z bins"] = "99"
    with effective.open("w") as stream:
        parser.write(stream)
    comparison_path = bundle / case["artifacts"]["configuration_comparison"]["path"]
    write_json(comparison_path, TOOL.compare_configs(original, effective))
    update_artifact(bundle, case, "effective_config")
    update_artifact(bundle, case, "configuration_comparison")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="unauthorized configuration"):
        TOOL.check_bundle(bundle)


def test_destination_refusal_and_consecutive_captures_are_fresh(tmp_path):
    """Captures refuse reuse and separate consecutive evidence trees."""
    first, first_run = capture_stub(tmp_path, "first")
    assert first_run.returncode == 0, first_run.stderr
    _, refused = capture_stub(tmp_path, "first")
    assert refused.returncode != 0
    assert "refusing existing capture destination" in refused.stderr
    second, second_run = capture_stub(tmp_path, "second")
    assert second_run.returncode == 0, second_run.stderr
    first_manifest = manifest_for(first)
    second_manifest = manifest_for(second)
    first_commands = [case["directory"] for case in first_manifest["primary_cases"]]
    second_commands = [case["directory"] for case in second_manifest["primary_cases"]]
    assert first.resolve() != second.resolve()
    assert first_commands == second_commands
    first_effective = (
        first
        / first_manifest["primary_cases"][0]["artifacts"]["effective_config"]["path"]
    ).read_text()
    second_effective = (
        second
        / second_manifest["primary_cases"][0]["artifacts"]["effective_config"]["path"]
    ).read_text()
    assert str(first) in first_effective
    assert str(second) in second_effective
    assert str(second) not in first_effective


@pytest.mark.parametrize("mode", ["number", "structure"])
def test_repeat_comparison_failures(tmp_path, mode):
    """Changed repeat numbers and structures both make capture incomplete."""
    bundle, completed = capture_stub(
        tmp_path,
        f"repeat-{mode}",
        {"FISHHIGHZ_STUB_REPEAT_MODE": mode},
    )
    assert completed.returncode != 0
    manifest = manifest_for(bundle)
    assert manifest["state"] == "incomplete"
    comparison_path = bundle / manifest["repeat"]["artifacts"]["comparison"]["path"]
    comparison = json.loads(comparison_path.read_text())
    assert comparison["passed"] is False
    assert "repeat mismatch" in comparison["error"]


def test_explicit_full_and_default_quick_action_counts(captured, quick_captured):
    """Suite selection controls scientific runs and resolution-worker counts."""
    base, full = captured
    full_actions = actions(base / "full-actions.jsonl")
    assert [item["action"] for item in full_actions].count("run") == 8
    assert [item["action"] for item in full_actions].count("resolve") == 7
    assert manifest_for(full)["coverage"]["scientific_run_count"] == 8

    quick, quick_log = quick_captured
    quick_actions = actions(quick_log)
    assert [item["action"] for item in quick_actions].count("run") == 1
    assert [item["action"] for item in quick_actions].count("resolve") == 1
    assert [item["case"] for item in quick_actions if item["action"] == "run"] == [
        TOOL.REPEAT_CASE
    ]
    quick_manifest = manifest_for(quick)
    assert quick_manifest["suite"] == "quick"
    assert quick_manifest["repeat"] is None


def test_quick_bundle_suite_requirement_and_relocation(quick_captured, tmp_path):
    """Quick evidence is self-contained after relocation but cannot satisfy full."""
    quick, _ = quick_captured
    checked = TOOL.check_bundle(quick)
    assert checked["suite"] == "quick"
    with pytest.raises(TOOL.BaselineError, match="not required 'full'"):
        TOOL.check_bundle(quick, require_suite="full")

    relocated = tmp_path / "relocated-quick"
    shutil.copytree(quick, relocated)
    result = run_tool(["check", relocated], tmp_path)
    assert result.returncode == 0, result.stderr
    assert "embedded full-baseline evidence" in result.stdout


def test_quick_embedded_compatibility_is_recomputed(quick_captured, tmp_path):
    """Relocated checks do not trust stale compatibility pass flags."""
    quick, _ = quick_captured
    bundle = tmp_path / "stale-compatibility"
    shutil.copytree(quick, bundle)
    manifest = manifest_for(bundle)
    record = manifest["baseline"]["artifacts"]["compatibility"]
    path = bundle / record["path"]
    compatibility = json.loads(path.read_text())
    compatibility["checks"]["source"]["current"][0]["sha256"] = "0" * 64
    write_json(path, compatibility)
    record["sha256"] = digest(path)
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="flags are inconsistent"):
        TOOL.check_bundle(bundle)


@pytest.mark.parametrize("mode", ["number", "structure"])
def test_quick_numerical_and_structural_mismatch(captured, tmp_path, mode):
    """Quick results must match copied full results numerically and structurally."""
    _, full = captured
    bundle, completed = capture_stub(
        tmp_path,
        f"quick-{mode}",
        {"FISHHIGHZ_STUB_QUICK_MODE": mode},
        suite="quick",
        baseline=full,
    )
    assert completed.returncode != 0
    manifest = manifest_for(bundle)
    assert manifest["state"] == "incomplete"
    comparison = json.loads(
        (bundle / manifest["baseline"]["artifacts"]["comparison"]["path"]).read_text()
    )
    assert comparison["passed"] is False
    assert "quick versus full baseline" in comparison["error"]


def test_quick_missing_baseline_fails_without_destination(tmp_path):
    """The default quick suite requires an explicit full baseline."""
    bundle, completed = capture_stub(tmp_path, "missing-baseline", suite=None)
    assert completed.returncode != 0
    assert "requires --baseline" in completed.stderr
    assert not bundle.exists()


@pytest.mark.parametrize("kind", ["quick", "incomplete", "corrupt"])
def test_invalid_baseline_rejected_before_capture(
    captured, quick_captured, tmp_path, kind
):
    """Only a complete, valid full bundle can authorize a quick forecast."""
    _, full = captured
    quick, _ = quick_captured
    baseline = quick if kind == "quick" else tmp_path / "invalid-full"
    if kind != "quick":
        shutil.copytree(full, baseline)
        manifest = manifest_for(baseline)
        if kind == "incomplete":
            manifest["state"] = "incomplete"
            write_json(baseline / "manifest.json", manifest)
        else:
            record = manifest["primary_cases"][0]["artifacts"]["summary"]
            (baseline / record["path"]).write_text("corrupt\n")
    action_log = tmp_path / f"{kind}-actions.jsonl"
    bundle, completed = capture_stub(
        tmp_path,
        f"rejected-{kind}",
        {"FISHHIGHZ_STUB_ACTION_LOG": str(action_log)},
        suite="quick",
        baseline=baseline,
    )
    assert completed.returncode != 0
    assert not bundle.exists()
    assert not action_log.exists()


@pytest.mark.parametrize(
    "mismatch", ["configuration", "source", "input", "environment"]
)
def test_quick_science_or_environment_mismatch_prevents_run(
    captured, tmp_path, mismatch
):
    """Compatibility mismatches retain evidence and launch no forecast worker."""
    base, full = captured
    reference = tmp_path / "reference"
    shutil.copytree(base / "reference", reference)
    if mismatch == "configuration":
        path = reference / "examples" / "desi2" / TOOL.REPEAT_CASE
        path.write_text(path.read_text().replace("z_ref = 2.3", "z_ref = 2.31"))
    elif mismatch == "source":
        path = reference / "lyaforecast" / "forecast_new.py"
        path.write_text(path.read_text() + "# source mismatch\n")
    elif mismatch == "input":
        path = reference / "resources" / "SNR" / "header-and-data.dat"
        path.write_text(path.read_text() + "3 4\n")
    environment = {"FISHHIGHZ_STUB_ACTION_LOG": str(tmp_path / "actions.jsonl")}
    if mismatch == "environment":
        environment["FISHHIGHZ_STUB_NUMPY_VERSION"] = "different"
    bundle, completed = capture_stub(
        tmp_path,
        f"incompatible-{mismatch}",
        environment,
        suite="quick",
        baseline=full,
    )
    assert completed.returncode != 0
    assert [item["action"] for item in actions(tmp_path / "actions.jsonl")].count(
        "run"
    ) == 0
    manifest = manifest_for(bundle)
    assert manifest["state"] == "incomplete"
    compatibility = json.loads(
        (
            bundle / manifest["baseline"]["artifacts"]["compatibility"]["path"]
        ).read_text()
    )
    compatibility_key = "inputs" if mismatch == "input" else mismatch
    assert compatibility["checks"][compatibility_key]["passed"] is False


def test_location_only_compatibility(captured, tmp_path):
    """A relocated but byte-identical reference remains compatible."""
    base, full = captured
    shutil.copytree(base / "reference", tmp_path / "moved-reference")
    bundle, completed = capture_stub(
        tmp_path,
        "location-only",
        suite="quick",
        baseline=full,
        reference_name="moved-reference",
    )
    assert completed.returncode == 0, completed.stderr
    assert TOOL.check_bundle(bundle)["suite"] == "quick"


def test_manifest_coverage_and_partial_capture_rejected(captured, tmp_path):
    """Suite labels cannot relabel full coverage or conceal partial captures."""
    relabeled = mutable_bundle(captured, tmp_path, "relabeled")
    manifest = manifest_for(relabeled)
    manifest["suite"] = "quick"
    write_json(relabeled / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="primary inventory"):
        TOOL.check_bundle(relabeled)

    partial = mutable_bundle(captured, tmp_path, "partial")
    manifest = manifest_for(partial)
    manifest["primary_cases"].pop()
    write_json(partial / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="unique ordered full"):
        TOOL.check_bundle(partial)


@pytest.mark.parametrize("mode", ["missing", "corrupt", "escape"])
def test_nested_tool_snapshot_validation(captured, tmp_path, mode):
    """Nested tool paths and bytes are independently checked."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    tool_record = manifest["shared_artifacts"]["tool_state"]
    tool_path = bundle / tool_record["path"]
    tool_state = json.loads(tool_path.read_text())
    snapshot = bundle / tool_state["tools"][0]["snapshot_path"]
    if mode == "missing":
        snapshot.unlink()
        match = "corrupted snapshot"
    elif mode == "corrupt":
        snapshot.write_text(snapshot.read_text() + "# corrupt\n")
        match = "corrupted snapshot"
    else:
        tool_state["tools"][0]["snapshot_path"] = "../../outside.py"
        write_json(tool_path, tool_state)
        update_shared(bundle, manifest, "tool_state")
        write_json(bundle / "manifest.json", manifest)
        match = "escapes bundle"
    with pytest.raises(TOOL.BaselineError, match=match):
        TOOL.check_bundle(bundle)


def test_required_shared_artifact_inventory(captured, tmp_path):
    """Removing shared evidence cannot be hidden by changing the manifest."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    manifest["shared_artifacts"].pop("tool_state")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="shared artifact inventory"):
        TOOL.check_bundle(bundle)


@pytest.mark.parametrize("location", ["manifest", "status"])
def test_failed_repeat_status_is_semantic_failure(captured, tmp_path, location):
    """Neither manifest nor stored repeat failure can be masked by agreement."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    repeated = manifest["repeat"]
    if location == "manifest":
        repeated["state"] = "failed"
        repeated["exit_code"] = 17
    else:
        status_path = bundle / repeated["artifacts"]["status"]["path"]
        status = json.loads(status_path.read_text())
        status["state"] = "failed"
        status["exit_code"] = 17
        write_json(status_path, status)
        update_artifact(bundle, repeated, "status")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="did not complete|not successful"):
        TOOL.check_bundle(bundle)


def test_repeat_roundtrip_and_artifact_aliases_rejected(captured, tmp_path):
    """Repeat evidence must verify serialization and remain independent."""
    bundle = mutable_bundle(captured, tmp_path, "roundtrip")
    manifest = manifest_for(bundle)
    repeated = manifest["repeat"]
    status_path = bundle / repeated["artifacts"]["status"]["path"]
    response_path = bundle / repeated["artifacts"]["worker_response"]["path"]
    status = json.loads(status_path.read_text())
    response = json.loads(response_path.read_text())
    status["worker"]["roundtrip_equal"] = False
    response["roundtrip_equal"] = False
    write_json(status_path, status)
    write_json(response_path, response)
    update_artifact(bundle, repeated, "status")
    update_artifact(bundle, repeated, "worker_response")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="round-trip"):
        TOOL.check_bundle(bundle)

    aliased = mutable_bundle(captured, tmp_path, "aliased")
    manifest = manifest_for(aliased)
    repeated = manifest["repeat"]
    primary = manifest["primary_cases"][-1]
    for name in TOOL.BASE_CASE_ARTIFACTS:
        repeated["artifacts"][name] = primary["artifacts"][name]
    write_json(aliased / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="outside its case directory|overlap"):
        TOOL.check_bundle(aliased)


@pytest.mark.parametrize("coordinate", ["k", "mu"])
def test_interior_grid_coordinate_mutation(captured, tmp_path, coordinate):
    """Every interior k and mu coordinate is checked against the INI formula."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    case = manifest["primary_cases"][0]
    metadata_path = bundle / case["artifacts"]["metadata"]["path"]
    metadata = json.loads(metadata_path.read_text())
    metadata["grid"][coordinate]["data"][1] = -1000
    write_json(metadata_path, metadata)
    update_artifact(bundle, case, "metadata")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(
        TOOL.BaselineError, match="coordinates are not strictly ordered"
    ):
        TOOL.check_bundle(bundle)


@pytest.mark.parametrize("location", ["result", "metadata", "bin-count"])
def test_redshift_configuration_linkage(captured, tmp_path, location):
    """Result and metadata redshifts must both match the configured binning."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    case = manifest["primary_cases"][0]
    result_path = bundle / case["artifacts"]["result_json"]["path"]
    metadata_path = bundle / case["artifacts"]["metadata"]["path"]
    if location == "result":
        result = json.loads(result_path.read_text())
        result["items"][0]["value"]["data"][0] += 0.001
        write_json(result_path, result)
        update_artifact(bundle, case, "result_json")
        match = "expected"
    else:
        metadata = json.loads(metadata_path.read_text())
        if location == "metadata":
            metadata["redshift_centers"]["data"][0] += 0.001
            match = "expected"
        else:
            metadata["redshift_centers"]["shape"] = [3]
            metadata["redshift_centers"]["data"].append(2.5)
            match = "shape"
        write_json(metadata_path, metadata)
        update_artifact(bundle, case, "metadata")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match=match):
        TOOL.check_bundle(bundle)


@pytest.mark.parametrize("identity", ["tracer_order", "all_pairs", "selected_pairs"])
def test_tracer_and_pair_identities_come_from_configuration(
    captured, tmp_path, identity
):
    """Metadata identities cannot be self-consistent inventions."""
    bundle = mutable_bundle(captured, tmp_path)
    manifest = manifest_for(bundle)
    case = manifest["primary_cases"][-1]
    metadata_path = bundle / case["artifacts"]["metadata"]["path"]
    metadata = json.loads(metadata_path.read_text())
    metadata[identity][0] = "invented"
    write_json(metadata_path, metadata)
    update_artifact(bundle, case, "metadata")
    write_json(bundle / "manifest.json", manifest)
    with pytest.raises(TOOL.BaselineError, match="identities|selected pairs"):
        TOOL.check_bundle(bundle)


@pytest.mark.parametrize(
    "mutation",
    [
        "source-change",
        "source-add",
        "source-delete",
        "input-change",
        "input-add",
        "input-delete",
        "config-change",
        "config-add",
        "config-delete",
    ],
)
def test_capture_detects_reference_mutation(captured, tmp_path, mutation):
    """Content and bounded inventory changes invalidate preservation evidence."""
    base, full = captured
    shutil.copytree(base / "reference", tmp_path / "reference")
    bundle, completed = capture_stub(
        tmp_path,
        f"mutated-{mutation}",
        {"FISHHIGHZ_STUB_MUTATE": mutation},
        suite="quick",
        baseline=full,
    )
    assert completed.returncode != 0
    manifest = manifest_for(bundle)
    assert manifest["state"] == "incomplete"
    assert manifest["primary_cases"][0]["state"] == "complete"
    assert "preserv" in manifest["failure"] or "unchanged" in manifest["failure"]
