"""Actual-control and failed-trial mutations reach schema-3 semantics."""

import copy
import json

import numpy as np
import pytest

from fishhighz.validation.evidence import (
    check,
    execute,
)
from fishhighz.validation.schema import metric_values, request, validate_payload
from fishhighz.validation.synthetic import convergence_payload
from fishhighz.validation.trials import validate


@pytest.fixture(scope="module")
def original():
    task = request("lbg_lae_3x2pt", 0, "accuracy")
    arrays, report = convergence_payload(task)
    report["provenance"] = dict(
        fishhighz=dict(origin="fixture", module_hashes={"a": "a" * 64}),
        wheel=dict(origin="fixture", modules={"a": "a" * 64}, sha256="b" * 64),
        reference=dict(
            reference_origin="/fixture",
            sources={"/fixture/a": "c" * 64},
            versions={"fixture": "1"},
        ),
        resources={"fixture": "d" * 64},
    )
    return task, arrays, report


def test_distinct_controls_identical_information(original):
    t, a, r = original
    assert validate_payload(t, a, r)
    assert all(m["fisher_relative"] == 0 for m in r["metrics"])


@pytest.mark.parametrize(
    "mutation",
    [
        "duplicate_operands",
        "wrong_ref",
        "repeated_ref",
        "missing_trial",
        "reordered_trials",
        "controls",
        "primary",
        "pair",
        "volume",
        "lower",
        "missing_contract",
        "missing_outcome",
    ],
)
def test_trial_mutations_writer_and_checker(tmp_path, original, mutation):
    t, a, r = copy.deepcopy(original)
    # Begin with a failed metric and retain the contradictory actual trial.
    if mutation == "duplicate_operands":
        index = r["trial_contract"]["successful_ids"].index(
            r["trial_contract"]["metric_trials"][0][0]
        )
        a["study_fisher"][index] *= 2
        a["study_pair_fisher"][index] *= 2
    elif mutation == "wrong_ref":
        r["trial_contract"]["metric_trials"][0][0] = r["trial_contract"][
            "metric_trials"
        ][1][0]
    elif mutation == "repeated_ref":
        r["trial_contract"]["metric_trials"][0][0] = r["trial_contract"][
            "metric_trials"
        ][0][1]
    elif mutation == "missing_trial":
        a["study_fisher"] = a["study_fisher"][:-1]
    elif mutation == "reordered_trials":
        r["study_controls"].reverse()
    elif mutation == "controls":
        r["study_controls"][0]["iterations"] = 24
    elif mutation == "primary":
        a["study_fisher"][
            r["trial_contract"]["successful_ids"].index(r["trial_contract"]["final_id"])
        ] *= 2
    elif mutation == "pair":
        a["study_pair_fisher"][
            r["trial_contract"]["successful_ids"].index(
                r["trial_contract"]["metric_trials"][0][0]
            )
        ] *= 2
    elif mutation == "volume":
        a["study_volume"][
            r["trial_contract"]["successful_ids"].index(
                r["trial_contract"]["metric_trials"][0][0]
            )
        ] *= 2
    elif mutation == "lower":
        r["combined_lower_controls"]["step"] *= 2
    elif mutation == "missing_contract":
        del r["trial_contract"]
    elif mutation == "missing_outcome":
        r["trial_contract"]["outcomes"].pop()
    r["metrics"] = metric_values(a, r["metric_names"])
    with pytest.raises((ValueError, KeyError)):
        validate_payload(t, a, r)
    out = tmp_path / "writer"
    m = execute(
        out, suite="full", cases=[t["case"]], bin_indices=[0], worker=lambda _: (a, r)
    )
    assert not m["complete"] and m["records"][0]["status"] == "failed"
    # Consistently rehash and force all manifest flags to reach offline semantics.
    row = m["records"][0]
    row.update(status="completed", scientific_passed=True)
    m["complete"] = True
    (out / "manifest.json").write_text(json.dumps(m))
    with pytest.raises((ValueError, KeyError)):
        check(out, verify_sources=False)


def test_failed_attempt_cannot_be_suppressed():
    from types import SimpleNamespace

    from fishhighz.validation.study import study

    class FailedStudy:
        selection = SimpleNamespace(
            fields=[SimpleNamespace(kind="forest")], selected_pairs=np.array([[0, 0]])
        )
        _prepared = {}

        def evaluate(self, task, controls):
            if controls["iterations"] == 24:
                raise ValueError("first unrepresentable operation")
            f = np.eye(2) * controls["iterations"]
            return dict(fisher=f, pair_fisher=f[None]), dict(
                settings=dict(grid=dict(volume=1.0))
            )

    a, r = study(FailedStudy(), {})
    # Replay itself must require the attempted 24-update failure even if both
    # the outcomes and unbound legacy failure list were maliciously cleared.
    assert any(o["outcome"] == "failed" for o in r["trial_contract"]["outcomes"])
    r["trial_contract"]["outcomes"] = [
        o for o in r["trial_contract"]["outcomes"] if o["outcome"] == "success"
    ]
    r["unresolved_controls"] = []
    r["settings"]["controls"] = r["final_controls"]
    r["settings"]["grid"].update(k_intervals=128, mu_order=32, k_order=4)
    with pytest.raises(ValueError, match="replay"):
        validate(a, r)


def test_exact_scoped_dispatch(tmp_path):
    from fishhighz.validation.synthetic import evidence_payload

    seen = []

    def worker(t):
        seen.append(t)
        return evidence_payload(t)

    m = execute(
        tmp_path / "b",
        suite="full",
        cases=["lya_qso_lbg_lae_15x2pt"],
        profiles=("compatibility", "accuracy"),
        kind="synthetic_bao",
        worker=worker,
    )
    assert len(seen) == 12 and sum(len(t["selected_pairs"]) for t in seen) == 180
    assert {t["case"] for t in seen} == {"lya_qso_lbg_lae_15x2pt"}
    assert m["complete"]
    assert len(m["not_run_cases"]) == 6


def test_cumulative_first_cell_fixed_point_and_refinement():
    from fishhighz.kernels.weights import _iterate

    # At the first cell the map is w -> c*w/(c*w+variance/S).
    # Subdivision reduces c; a positive fixed point disappears below threshold.
    for mass in (0.25, 0.5, 2.0):
        w, _ = _iterate(np.array([mass]), np.ones(1), 1, 1, 1, 1, 100)
        np.testing.assert_allclose(w, [max(0, 1 - 1 / mass)], atol=1e-25)


def test_guard_diagnosis_identifies_mixed_product():
    from fishhighz.validation.weight_diagnosis import first_underflow

    result = first_underflow(np.array([1e-100]), np.array([1e-110]), np.array([1e100]))
    assert result["operation"] == "r*w*w"
    assert result["index"] == 0 and "underflow" in result["error"]
    assert (
        first_underflow(np.array([0.0, 1.0]), np.array([1.0, 1.0]), np.ones(2)) is None
    )


@pytest.mark.parametrize("null", [False, True])
def test_rank_and_roundoff_trial_controls(null):
    task = request("lbg_lae_3x2pt", 0, "accuracy", kind="synthetic_bao")
    a, r = convergence_payload(task, null=null)
    a["study_fisher"][:, 0, 0] *= 1 + np.finfo(float).eps
    assert validate(a, r)


def test_scoped_diagnostic_inventory_and_gate(tmp_path):
    import runpy
    from pathlib import Path

    from fishhighz.validation.schema import request
    from fishhighz.validation.synthetic import evidence_payload

    script = Path(__file__).resolve().parents[1] / "scripts/check_desi2_15x2pt.py"
    module = runpy.run_path(str(script))
    diagnostics = [
        request(
            "lya_qso_lbg_lae_15x2pt", i, "accuracy", kind="diagnostic", diagnostic_id=p
        )
        for i in range(6)
        for p in module["POLICIES"]
    ]
    m = execute(
        tmp_path / "b",
        suite="full",
        cases=["lya_qso_lbg_lae_15x2pt"],
        profiles=("compatibility", "accuracy"),
        kind="synthetic_bao",
        diagnostic_requests=diagnostics,
        worker=evidence_payload,
    )
    assert len(m["records"]) == 12 and len(m["diagnostics"]) == 72 and m["complete"]
    # Synthetic coverage is never relabeled real scientific acceptance.
    with pytest.raises(ValueError, match="12 primary"):
        module["scoped_gate"](tmp_path / "b")
    m["schema"] = 2
    (tmp_path / "b/manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match="schema 3"):
        module["scoped_gate"](tmp_path / "b")
