"""Per-spectrum trial binding across sixteen decades of information."""

import copy
import json
import runpy
from pathlib import Path

import numpy as np
import pytest

from fishhighz.validation import schema, trials
from fishhighz.validation.cases import selection
from fishhighz.validation.evidence import check, execute
from fishhighz.validation.numerics import relative
from fishhighz.validation.study import DEFAULT, study

CASE = "lya_qso_lbg_lae_15x2pt"


def provenance():
    return dict(
        fishhighz=dict(origin="fixture", module_hashes={"a": "a" * 64}),
        wheel=dict(origin="fixture", modules={"a": "a" * 64}, sha256="b" * 64),
        reference=dict(
            reference_origin="/fixture",
            sources={"/fixture/a": "c" * 64},
            versions={"fixture": "1"},
        ),
        resources={"fixture": "d" * 64},
    )


@pytest.fixture(autouse=True)
def tiny_nodes(monkeypatch):
    # Only quadrature construction is replaced, with four deterministic nodes.
    # The trial controller retains every real bounded control and its identity.
    monkeypatch.setattr(
        schema,
        "grid_nodes",
        lambda _: (
            np.array([0.02, 0.1, 0.2, 0.4]),
            np.array([0.1, 0.3, 0.6, 0.9]),
            np.ones(4),
        ),
    )


def payload(ratio=1e-16, pair=0, *, converged=False, task=None):
    task = task or schema.request(CASE, 0, "accuracy")
    settings = dict(
        profile=task["profile"],
        parameters=task["parameters"],
        bounds=task["bounds"],
        fields=[f["id"] for f in task["fields"]],
        controls=dict(DEFAULT),
        grid=dict(
            kind="gauss_legendre",
            volume=1000.0,
            h_fid=0.7,
            k_intervals=128,
            k_order=4,
            mu_order=32,
        ),
    )
    k, mu, _ = schema.grid_nodes(settings)
    pairs = np.array(task["required_pairs"])
    total = np.tile(np.eye(5)[pairs[:, 0], pairs[:, 1]], (len(k), 1))
    j = np.ones((len(k), 15, 2))
    j[:, :, 1] = mu[:, None] ** 2
    j[:, pair] *= np.sqrt(ratio)
    base, report = schema.assemble(task, total, j, settings)

    class WeakStudy:
        selection = selection(CASE)
        _prepared = {}

        def evaluate(self, task, controls):
            a, r = copy.deepcopy((base, report))
            factor = 1 if converged else (1 + 1 / controls["iterations"]) / (1 + 1 / 24)
            a["pair_fisher"][pair] *= factor
            a["fisher"] += a["pair_fisher"][pair] - base["pair_fisher"][pair]
            return a, r

    a, r = study(WeakStudy(), task)
    r["settings"]["controls"] = r["final_controls"].copy()
    r["effective_hash"] = schema.canonical(r["settings"])
    a["effective_token"] = schema.token(r["settings"])
    r["provenance"] = provenance()
    return task, a, r


def rejected_everywhere(tmp_path, task, arrays, report, reason, record_property):
    with pytest.raises(ValueError, match=reason) as error:
        schema.validate_payload(task, arrays, report, require_pass=False)
    record_property("validator_rejection", str(error.value))
    out = tmp_path / "writer"
    m = execute(
        out,
        suite="full",
        cases=[CASE],
        profiles=("accuracy",),
        bin_indices=[0],
        worker=lambda _: (arrays, report),
    )
    row = m["records"][0]
    assert row["status"] == "failed" and reason in row["error"]
    record_property("writer_rejection", row["error"])
    # Writer saved the mutated arrays/report and their actual hashes. Promote
    # only outer flags to reach the offline reader's semantic validation.
    row.update(status="completed", scientific_passed=True)
    m["complete"] = True
    (out / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match=reason) as error:
        check(out, verify_sources=False)
    record_property("offline_rejection", str(error.value))


@pytest.mark.parametrize("ratio", [1.0, 1e-8, 1e-16])
@pytest.mark.parametrize("pair", [0, 1, 8, 14])
@pytest.mark.parametrize("family", ["weights", "combined", "all"])
def test_weak_lower_operand(tmp_path, ratio, pair, family, record_property):
    task, a, r = payload(ratio, pair)
    assert schema.validate_payload(task, a, r, require_pass=False)
    i = r["metric_names"].index("weights" if family == "all" else family)
    # Analytic inverse of s*F gives sigma/sqrt(s), independent of amplitude.
    expected = 1 - np.sqrt((1 + 1 / 24) / (1 + 1 / 12))
    assert r["metrics"][i]["pair_error_relative"] == pytest.approx(expected, abs=2e-15)
    original = a["study_pair_fisher"].copy()
    indices = range(len(r["metric_names"])) if family == "all" else [i]
    for index in indices:
        a["metric_pair_fisher"][index, 0, pair] = a["metric_pair_fisher"][
            index, 1, pair
        ]
    r["metrics"] = schema.metric_values(a, r["metric_names"])
    r["passed"] = True
    assert np.array_equal(original, a["study_pair_fisher"])
    rejected_everywhere(tmp_path, task, a, r, "metric/trial operand", record_property)


@pytest.mark.parametrize(
    "mutation",
    [
        "primary",
        "final",
        "replay",
        "metric",
        "verdict",
        "false_pass",
        "pair_fisher",
        "pair_errors",
        "pair_covariance",
        "pair_correlation",
    ],
)
def test_weak_other_bindings(tmp_path, mutation, record_property):
    t, a, r = payload(converged=True)
    if mutation == "primary":
        a["pair_fisher"][0] *= 1.04
        with pytest.raises(ValueError, match="primary/final pair trial"):
            trials.validate(a, r)
        reason = "pair_fisher block"
    elif mutation == "final":
        i = r["trial_contract"]["successful_ids"].index(r["trial_contract"]["final_id"])
        a["study_pair_fisher"][i, 0] *= 1.04
        reason = "primary/final pair trial"
    elif mutation == "replay":
        # Failed combined trial has no successful operand binding. Its saved
        # placeholder must nevertheless agree with the controller replay.
        refs = r["trial_contract"]["metric_trials"][-1]
        identity = refs[0]
        i = r["trial_contract"]["successful_ids"].index(identity)
        row = next(o for o in r["trial_contract"]["outcomes"] if o["id"] == identity)
        row.update(outcome="failed", error="unrepresentable arithmetic")
        del row["array_index"]
        r["trial_contract"]["successful_ids"].pop(i)
        r["study_controls"].pop(i)
        for key in ("study_fisher", "study_pair_fisher", "study_volume"):
            a[key] = np.delete(a[key], i, axis=0)
        for row in r["trial_contract"]["outcomes"]:
            if row.get("array_index", -1) > i:
                row["array_index"] -= 1
        r["unresolved_controls"] = [
            dict(control="combined", error="unrepresentable arithmetic")
        ]
        r["passed"] = False
        assert schema.validate_payload(t, a, r, require_pass=False)
        a["metric_pair_fisher"][-1, 0, 0] *= 1.04
        r["metrics"] = schema.metric_values(a, r["metric_names"])
        reason = "replayed trial operands metric_pair_fisher"
    elif mutation == "metric":
        r["metrics"][0]["pair_error_relative"] = 0.0001
        reason = "replayed trial metrics differ"
    elif mutation == "false_pass":
        t, a, r = payload()
        r["passed"] = True
        reason = "replayed trial convergence verdict differs"
    elif mutation == "verdict":
        r["passed"] = False
        reason = "replayed trial convergence verdict differs"
    else:
        a[mutation][0] *= 1.04
        reason = mutation + " block"
    rejected_everywhere(tmp_path, t, a, r, reason, record_property)


@pytest.mark.parametrize("ratio", [1.0, 1e-8, 1e-16])
def test_valid_weak_converged_roundoff_and_unconverged(tmp_path, ratio):
    t, a, r = payload(ratio, converged=True)
    a["study_pair_fisher"][:, 0] *= 1 + np.finfo(float).eps
    assert schema.validate_payload(t, a, r)
    m = execute(
        tmp_path / "valid",
        suite="full",
        cases=[CASE],
        profiles=("accuracy",),
        bin_indices=[0],
        worker=lambda _: (a, r),
    )
    assert m["complete"] and check(tmp_path / "valid", verify_sources=False)["complete"]
    t, a, r = payload(ratio)
    assert not r["passed"] and schema.validate_payload(t, a, r, require_pass=False)
    with pytest.raises(ValueError, match="scientific validation failed"):
        schema.validate_payload(t, a, r)


@pytest.mark.parametrize("scale", [1e-300, 1e-150, 1.0, 1e150, 1e300])
def test_range_safe_block_comparison(scale):
    a = np.diag([scale, scale / 2])
    schema._close_blocks(a * (1 + 1e-13), a, "range")
    with pytest.raises(ValueError, match="range block"):
        schema._close_blocks(a * 1.04, a, "range")
    assert relative(a, np.zeros_like(a)) == float("inf")
    assert relative(np.zeros_like(a), np.zeros_like(a)) == 0


@pytest.mark.parametrize("value", [np.inf, np.nan])
def test_nonfinite_block_rejected(value):
    with pytest.raises(ValueError, match="finite block"):
        schema._close_blocks(np.array([[value]]), np.ones((1, 1)), "finite")


def test_null_and_partially_constrained_trials():
    for null in (np.zeros((2, 2)), np.diag([1e-250, 0.0])):
        t, a, r = payload(converged=True)
        a["pair_fisher"][0] = null
        a["study_pair_fisher"][:, 0] = null
        a["metric_pair_fisher"][:, :, 0] = null
        r["metrics"] = schema.metric_values(a, r["metric_names"])
        assert trials.validate(a, r)


def test_full_scoped_inventory(tmp_path, monkeypatch):
    module = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts/check_desi2_15x2pt.py")
    )
    diagnostics = [
        schema.request(CASE, i, "accuracy", kind="diagnostic", diagnostic_id=p)
        for i in range(6)
        for p in module["POLICIES"]
    ]

    def worker(t):
        _, a, r = payload(converged=True, task=t)
        if t["profile"] == "compatibility":
            from fishhighz.validation.numerics import change

            a["reference_fisher"] = a["fisher"].copy()
            a["reference_pair_fisher"] = a["pair_fisher"].copy()
            r["comparison"] = change(
                a["fisher"], a["fisher"], a["pair_fisher"], a["pair_fisher"], 1, 1
            )
        return a, r

    out = tmp_path / "scope"
    m = execute(
        out,
        suite="full",
        cases=[CASE],
        profiles=("compatibility", "accuracy"),
        diagnostic_requests=diagnostics,
        worker=worker,
    )
    assert m["complete"] and len(m["records"]) == 12 and len(m["diagnostics"]) == 72
    assert module["scoped_gate"](out)["complete"]
    for kind in ("omission", "extra", "duplicate", "kind", "diagnostic"):
        bad = copy.deepcopy(m)
        if kind == "omission":
            bad["records"].pop()
        elif kind == "extra":
            bad["records"].append(bad["records"][0])
        elif kind == "duplicate":
            bad["records"][1] = bad["records"][0]
        elif kind == "kind":
            bad["requested"][0]["kind"] = "synthetic_bao"
        else:
            bad["diagnostics"].pop()
        (out / "manifest.json").write_text(json.dumps(bad))
        with pytest.raises(ValueError, match="12 primary|missing/extra"):
            module["scoped_gate"](out)


@pytest.mark.parametrize("offset", [1e-16, 1e-4])
def test_compatibility_metric_roundoff(offset):
    task = schema.request(CASE, 0, "compatibility")
    _, a, r = payload(converged=True, task=task)
    a["reference_fisher"] = a["fisher"].copy()
    a["reference_pair_fisher"] = a["pair_fisher"].copy()
    r["comparison"] = dict(
        fisher_relative=offset,
        error_relative=0.0,
        pair_error_relative=0.0,
        volume_relative=0.0,
    )
    if offset < 1e-15:
        assert schema.validate_payload(task, a, r)
    else:
        with pytest.raises(ValueError, match="incorrect reference comparison metrics"):
            schema.validate_payload(task, a, r, require_pass=False)
