"""Focused controls for the W01 signed legacy weighting diagnostic."""

import copy
import importlib.util
import json
import sys
from dataclasses import replace
from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts/diagnose_legacy_forest_iterations.py"
)
SPEC = importlib.util.spec_from_file_location("legacy_forest_iterations", SCRIPT)
TOOL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = TOOL
SPEC.loader.exec_module(TOOL)


def legacy_inputs(density, variance, *, length=2.0, pixel=1.0, signal=1.0, alias=1.0):
    density = np.asarray(density, dtype=np.float64)
    return TOOL.LegacyInputs(
        magnitudes=np.arange(len(density), dtype=np.float64),
        density=density,
        variance=np.asarray(variance, dtype=np.float64),
        magnitude_step=1.0,
        forest_length=float(length),
        pixel_width=float(pixel),
        auxiliary_signal=float(signal),
        auxiliary_p1d=float(alias),
    )


@pytest.mark.parametrize("d", [0.5, 1.0, 2.0])
def test_one_cell_coefficients_do_not_follow_weight_amplitude(d):
    length, pixel, signal, variance = 2.0, 1.0, 1.0, 1.0
    mass = d * pixel * variance / (length * signal)
    inputs = legacy_inputs(
        [mass], [variance], length=length, pixel=pixel, signal=signal
    )
    expected_a = 1 / (length * mass)
    expected_p = pixel * variance / (length * mass)
    states = [TOOL.float_snapshot(inputs, count) for count in (0, 3, 6)]
    for state in states:
        assert state["a"] == pytest.approx(expected_a, rel=2e-15)
        assert state["p_pixel"] == pytest.approx(expected_p, rel=2e-15)
    assert not TOOL.coefficient_comparison(states[-1], states[1])["triggered"]


@pytest.mark.parametrize(
    "density,variance",
    [
        ([0.2, 0.7], [0.3, 2.0]),
        ([0.2, -0.01, 0.7], [0.3, 1.1, 2.0]),
    ],
)
@pytest.mark.parametrize("count", [0, 3, 6])
def test_vectorized_literal_order_matches_scalar_decimal(density, variance, count):
    inputs = legacy_inputs(
        density, variance, length=5.0, pixel=0.7, signal=1.3, alias=2.1
    )
    floating = TOOL.float_snapshot(inputs, count)
    decimal = TOOL.decimal_snapshot(inputs, count, 100)
    np.testing.assert_allclose(
        floating["weights"],
        [float(value) for value in decimal["weights"]],
        rtol=5e-13,
        atol=0,
    )
    for name in ("i1", "i2", "i3", "a", "p_pixel"):
        assert floating[name] == pytest.approx(float(decimal[name]), rel=5e-13)


def test_signed_cancellation_is_classified_as_a_singularity():
    inputs = legacy_inputs([1.0, -1.0], [1.0, 1.0])
    previous = np.geterr()
    with pytest.raises(
        TOOL.DiagnosticArithmeticError, match="zero cumulative effective density"
    ) as caught:
        TOOL.float_snapshot(inputs, 1)
    assert caught.value.kind == "singular_signed_recurrence"
    assert caught.value.index == 1
    assert np.geterr() == previous


def test_coefficient_trigger_is_not_hidden_by_small_absolute_weights():
    reference = {"a": 1.0, "p_pixel": 2.0, "weights": np.array([1e-200])}
    changed = {"a": 1.002, "p_pixel": 2.0, "weights": np.array([2e-200])}
    comparison = TOOL.coefficient_comparison(changed, reference)
    assert comparison["triggered"]
    stable = {"a": 1.0, "p_pixel": 2.0, "weights": np.array([1e-250])}
    assert not TOOL.coefficient_comparison(stable, reference)["triggered"]


def test_scalar_fisher_normalization_and_analytic_inverse():
    jacobian = np.array([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
    covariance = np.array([2.0, 4.0, 5.0])
    fisher, errors = TOOL.scalar_fisher_oracle(jacobian, covariance)
    expected = sum(
        np.outer(row, row) / variance for row, variance in zip(jacobian, covariance)
    )
    np.testing.assert_allclose(fisher, expected, rtol=0, atol=0)
    np.testing.assert_allclose(errors, np.sqrt(np.diag(np.linalg.inv(expected))))
    summary = TOOL._rank_aware_summary(fisher)
    np.testing.assert_allclose(summary["errors"], errors, rtol=5e-15, atol=0)


def test_null_fisher_direction_has_no_invented_error_bar():
    summary = TOOL._rank_aware_summary(np.diag([4.0, 0.0]))
    assert summary["rank"] == 1
    assert summary["errors"][0] == pytest.approx(0.5)
    assert np.isnan(summary["errors"][1])
    np.testing.assert_array_equal(summary["constrained"], [1, 0])


def test_field_and_pair_permutations_resolve_columns_by_identity():
    fields = ("qso", "lya(qso)", "lbg")
    selected = np.array([[0, 0], [1, 2], [1, 1], [0, 1]])
    required = np.array([[1, 2], [0, 0], [0, 1], [1, 1]])
    assert TOOL._resolve_auto_columns(fields, selected, required) == (2, 3)
    bad = required.copy()
    bad[3] = [2, 2]
    with pytest.raises(ValueError, match="expected one"):
        TOOL._resolve_auto_columns(fields, selected, bad)


def synthetic_source():
    inputs = legacy_inputs(
        [0.2, -0.01, 0.7],
        [0.3, 1.1, 2.0],
        length=5.0,
        pixel=0.7,
        signal=1.3,
        alias=2.1,
    )
    captured = TOOL.float_snapshot(inputs, 3)
    pairs = np.array([[0, 0], [0, 1], [1, 1]], dtype=np.int64)
    projection = TOOL.ProjectionInputs(
        k=np.array([0.05, 0.08, 0.12, 0.17]),
        mu=np.array([0.1, 0.4, 0.7, 0.9]),
        modes=np.array([20.0, 30.0, 40.0, 50.0]),
        baseline_total=np.array([12.0, 11.0, 10.0, 9.0]),
        observed_j=np.array([[1.0, 0.2], [0.3, 1.1], [0.8, 0.4], [0.2, 0.9]]),
        saved_reference_fisher=np.zeros((2, 2)),
        saved_pair_fisher=np.zeros((2, 2)),
        distance_to_velocity=100.0,
        angle_to_distance=60.0,
        covariance_redshift=2.1,
        resolution_width=110.0,
        selected_pairs=pairs,
        required_pairs=pairs,
        field_ids=("qso", "lya(qso)"),
        selected_column=2,
        required_column=2,
    )
    source = TOOL.SourceInputs(
        legacy=inputs,
        projection=projection,
        captured_weights=captured["weights"],
        captured_a=captured["a"],
        captured_p_pixel=captured["p_pixel"],
        provenance={
            "fixture": "self-contained",
            "report": {"path": "/synthetic/report.json", "sha256": "a" * 64},
        },
    )
    reference = TOOL.project_one_spectrum(source, captured, captured)["fisher"]
    return replace(
        source,
        projection=replace(
            projection,
            saved_reference_fisher=reference,
            saved_pair_fisher=reference,
        ),
    )


@pytest.fixture(scope="module")
def payload():
    source = synthetic_source()
    result = TOOL.run_diagnostic(source)
    summary = TOOL._summary(source, result, "synthetic focused test")
    arrays = TOOL._payload_arrays(source, result)
    assert TOOL.validate_payload(
        summary,
        arrays,
        authoritative_source=source,
    )
    return source, summary, arrays


@pytest.mark.parametrize("mutation", ["scalar", "density", "population", "pair"])
def test_source_and_pair_binding_rejects_mutations(payload, mutation):
    source, original_summary, original_arrays = payload
    summary = copy.deepcopy(original_summary)
    arrays = {name: value.copy() for name, value in original_arrays.items()}
    if mutation == "scalar":
        summary["inputs"]["scalars"]["auxiliary_signal"] *= 1.01
    elif mutation == "density":
        arrays["source_density"][1] *= 1.01
    elif mutation == "population":
        summary["identity"]["population"] = "lya(lbg)"
    else:
        summary["identity"]["selected_column"] = 0
    with pytest.raises((ValueError, AssertionError)):
        TOOL.validate_payload(
            summary,
            arrays,
            authoritative_source=source,
        )


def test_changed_count_with_unchanged_trajectory_is_rejected(payload):
    source, summary, original_arrays = payload
    arrays = {name: value.copy() for name, value in original_arrays.items()}
    arrays["state_counts"][-1] = 12
    with pytest.raises((ValueError, AssertionError)):
        TOOL.validate_payload(
            summary,
            arrays,
            authoritative_source=source,
        )


def test_attempt_inventory_preserves_failures_caps_and_early_stop():
    failure = [
        {"count": 0, "status": "completed"},
        {"count": 3, "status": "completed"},
        {"count": 6, "status": "arithmetic_failure"},
        {"count": 12, "status": "not_attempted_after_failure"},
        {"count": 24, "status": "not_attempted_after_failure"},
    ]
    assert TOOL.validate_attempt_inventory(failure, "arithmetic_failure")
    with pytest.raises(ValueError, match="every prescribed count"):
        TOOL.validate_attempt_inventory(failure[:-1], "arithmetic_failure")
    capped = [{"count": count, "status": "completed"} for count in TOOL.CHECKPOINTS]
    assert TOOL.validate_attempt_inventory(capped, "no_sensitivity_detected_within_cap")
    with pytest.raises(ValueError, match="every prescribed count"):
        TOOL.validate_attempt_inventory(
            capped[:-1], "no_sensitivity_detected_within_cap"
        )


def test_first_discrepancy_stops_later_counts_and_broad_calls():
    evaluated = []
    confirmed = []
    broad_calls = []

    def state(count):
        evaluated.append(count)
        return {
            "a": 1.002 if count == 6 else 1.0,
            "p_pixel": 2.0,
            "weights": np.array([1e-200 / (count + 1)]),
        }

    def confirm(count, error):
        confirmed.append((count, error))
        return {"confirmed": True}

    states, _, attempts, terminal = TOOL.select_bounded_attempts(state, confirm)
    if terminal == "coefficient_sensitive":
        for count in (3, 6):
            broad_calls.append(("one_spectrum", count))
    assert sorted(states) == [0, 3, 6]
    assert evaluated == [0, 3, 6]
    assert [row["status"] for row in attempts[-2:]] == [
        "not_attempted_after_identification",
        "not_attempted_after_identification",
    ]
    assert broad_calls == [("one_spectrum", 3), ("one_spectrum", 6)]
    assert all(name == "one_spectrum" for name, _ in broad_calls)


def test_failed_baseline_gate_stops_before_later_count():
    evaluated = []

    def state(count):
        evaluated.append(count)
        return {"a": 1.0, "p_pixel": 2.0, "weights": np.ones(1)}

    def reject(_state):
        raise ValueError("baseline identity")

    with pytest.raises(ValueError, match="baseline identity"):
        TOOL.select_bounded_attempts(state, lambda count, error: {}, reject)
    assert evaluated == [0, 3]


def test_output_directory_is_exclusive_and_cli_paths_are_explicit(tmp_path, payload):
    source, _, _ = payload
    result = TOOL.run_diagnostic(source)
    instructions = tmp_path / "instructions.md"
    instructions.write_text("approved synthetic instructions\n")
    existing = tmp_path / "already-exists"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        TOOL.write_evidence(existing, instructions, source, result, "test")
    args = TOOL.parse_args(
        [
            "--input",
            "saved-input",
            "--output",
            "new-output",
            "--instructions",
            "instructions.md",
        ]
    )
    assert (args.input, args.output, args.instructions) == (
        "saved-input",
        "new-output",
        "instructions.md",
    )


@pytest.mark.parametrize("quantity", ["a", "p_pixel"])
@pytest.mark.parametrize(
    "value,reference,passed",
    [(0.0, 0.0, True), (1.0, 0.0, False), (0.0, 1.0, False)],
)
def test_scalar_baseline_zero_matching_is_exact(quantity, value, reference, passed):
    source = synthetic_source()
    state = TOOL.float_snapshot(source.legacy, TOOL.BASELINE_COUNT)
    state[quantity] = value
    if quantity == "a":
        source = replace(source, captured_a=reference)
    else:
        source = replace(source, captured_p_pixel=reference)
    check = TOOL.baseline_gate(source, state)[quantity]
    assert check["passed"] is passed
    assert check["exact_zero_match"] is (value == reference == 0)


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def _finalized_bundle(tmp_path, source):
    output = tmp_path / "evidence"
    instructions = tmp_path / "instructions.md"
    handoff = tmp_path / "handoff.md"
    checks = tmp_path / "checks.json"
    instructions.write_text("approved W01 revision 2 instructions\n")
    handoff.write_text("W01 revision 2 synthetic handoff\n")
    _write_json(
        checks,
        {
            "schema": "fishhighz-weighting-diagnostic-checks-v1",
            "checks": [{"command": "synthetic focused check", "exit_code": 0}],
        },
    )
    result = TOOL.run_diagnostic(source)
    TOOL.write_evidence(output, instructions, source, result, "synthetic CLI")
    TOOL.finalize_evidence(output, "synthetic-input", handoff, checks)
    return output, instructions, handoff


def _refresh_manifest_hash(output, name):
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["files"][name] = TOOL.sha256_file(output / name)
    _write_json(manifest_path, manifest)


@pytest.mark.parametrize(
    "mutation",
    [
        "terminal_verdict",
        "trigger_operands",
        "table_coefficients_errors",
        "forecast_comparison",
        "projected_errors",
        "projected_total",
        "projected_covariance",
        "p1d",
        "response",
        "decimal_coefficients",
        "script_identity",
        "source_report_path",
        "source_report_hash",
        "manifest_artifact_hash",
    ],
)
def test_full_check_only_rejects_bound_claim_corruption(
    tmp_path, monkeypatch, mutation
):
    source = synthetic_source()
    monkeypatch.setattr(TOOL, "load_source", lambda _path: source)
    output, instructions, handoff = _finalized_bundle(tmp_path, source)
    summary_path = output / "summary.json"
    summary = json.loads(summary_path.read_text())

    if mutation == "terminal_verdict":
        summary["terminal"] = "arithmetic_failure"
        summary["verdicts"] = ["arithmetic_failure"]
    elif mutation == "trigger_operands":
        summary["attempts"][2]["relative_to_t3"]["a_relative"] = 0.0
        summary["attempts"][2]["relative_to_t3"]["p_pixel_relative"] = 0.0
    elif mutation == "table_coefficients_errors":
        summary["table"][1]["a"] *= 1.1
        summary["table"][1]["errors"] = [1.0, 1.0]
    elif mutation == "forecast_comparison":
        summary["forecast_comparison"]["fisher_relative"] = 0.0
    elif mutation == "decimal_coefficients":
        summary["decimal_confirmations"]["3"]["results"][0]["a"] = "0"
    elif mutation == "script_identity":
        summary["artifacts"]["diagnose_legacy_forest_iterations.py"] = "0" * 64
    elif mutation == "source_report_path":
        summary["provenance"]["report"]["path"] = "/unrelated/report.json"
    elif mutation == "source_report_hash":
        summary["provenance"]["report"]["sha256"] = "0" * 64
    elif mutation == "manifest_artifact_hash":
        manifest_path = output / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"]["arrays.npz"] = "0" * 64
        _write_json(manifest_path, manifest)
    else:
        with np.load(output / "arrays.npz", allow_pickle=False) as saved:
            arrays = {name: saved[name].copy() for name in saved.files}
        key = mutation
        arrays[key].flat[0] *= 1.1
        np.savez_compressed(output / "arrays.npz", **arrays)
        summary["artifacts"]["arrays.npz"] = TOOL.sha256_file(output / "arrays.npz")

    if mutation != "manifest_artifact_hash":
        _write_json(summary_path, summary)
        _refresh_manifest_hash(output, "summary.json")
    if mutation in {
        "projected_errors",
        "projected_total",
        "projected_covariance",
        "p1d",
        "response",
    }:
        _refresh_manifest_hash(output, "arrays.npz")
    with pytest.raises((ValueError, AssertionError)):
        TOOL.main(
            [
                "--input",
                "synthetic-input",
                "--output",
                str(output),
                "--instructions",
                str(instructions),
                "--handoff",
                str(handoff),
                "--check-only",
            ]
        )


@pytest.mark.parametrize(
    "route,expected_terminal",
    [
        ("initialization", "singular_signed_recurrence"),
        ("three_update", "arithmetic_failure"),
        ("decimal_disagreement", "arithmetic_failure"),
        ("baseline_mismatch", "source_discrepancy"),
    ],
)
def test_failed_baseline_cli_writes_and_rechecks(
    tmp_path, monkeypatch, route, expected_terminal
):
    source = synthetic_source()
    original_float = TOOL.float_snapshot
    original_decimal = TOOL.decimal_snapshot
    if route == "initialization":
        legacy = replace(
            source.legacy,
            pixel_width=1.0,
            auxiliary_p1d=3.0,
            variance=np.array([-3.0, 1.1, 2.0]),
        )
        source = replace(source, legacy=legacy)
    elif route == "three_update":

        def fail_three_update(inputs, count):
            if count == 3:
                raise TOOL.DiagnosticArithmeticError(
                    "arithmetic_failure", "update 3", "synthetic float64 range"
                )
            return original_float(inputs, count)

        monkeypatch.setattr(TOOL, "float_snapshot", fail_three_update)
    elif route == "decimal_disagreement":

        def disagree(inputs, count, precision):
            state = original_decimal(inputs, count, precision)
            if count == 3 and precision == 160:
                state = dict(state)
                state["a"] *= Decimal("1.1")
            return state

        monkeypatch.setattr(TOOL, "decimal_snapshot", disagree)
    else:
        source = replace(source, captured_a=source.captured_a * 1.1)

    monkeypatch.setattr(TOOL, "load_source", lambda _path: source)
    output = tmp_path / "evidence"
    instructions = tmp_path / "instructions.md"
    handoff = tmp_path / "handoff.md"
    checks = tmp_path / "checks.json"
    instructions.write_text("approved W01 revision 2 instructions\n")
    handoff.write_text("bounded failed-baseline handoff\n")
    _write_json(
        checks,
        {
            "schema": "fishhighz-weighting-diagnostic-checks-v1",
            "checks": [{"command": "failure-route CLI", "exit_code": 0}],
        },
    )
    common = [
        "--input",
        "synthetic-input",
        "--output",
        str(output),
        "--instructions",
        str(instructions),
    ]
    assert TOOL.main(common) == 0
    assert (
        TOOL.main(
            common + ["--finalize", "--handoff", str(handoff), "--checks", str(checks)]
        )
        == 0
    )
    assert TOOL.main(common + ["--check-only", "--handoff", str(handoff)]) == 0

    summary = json.loads((output / "summary.json").read_text())
    assert summary["terminal"] == expected_terminal
    assert "reproduced_legacy" not in summary["verdicts"]
    assert summary["projected_counts"] == []
    assert all(
        row["status"] == "not_attempted_after_failure"
        for row in summary["attempts"]
        if row["count"] > summary["terminal_event"]["count"]
    )
    with np.load(output / "arrays.npz", allow_pickle=False) as arrays:
        assert arrays["projected_counts"].size == 0
        assert arrays["p1d"].size == 0


def test_float_failure_with_successful_decimal_is_not_algebraic_singularity(
    monkeypatch,
):
    source = synthetic_source()
    original_float = TOOL.float_snapshot

    def fail_float(inputs, count):
        if count == 0:
            raise TOOL.DiagnosticArithmeticError(
                "singular_signed_recurrence", "initialization", "float cancellation"
            )
        return original_float(inputs, count)

    monkeypatch.setattr(TOOL, "float_snapshot", fail_float)
    result = TOOL.run_diagnostic(source)
    assert result["terminal"] == "arithmetic_failure"
    assert result["terminal_event"]["decimal_failures"] == [None, None]
