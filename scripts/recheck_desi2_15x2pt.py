"""Derive six-bin trial bindings from immutable revision-2 numerical evidence."""

import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from fishhighz.validation.evidence import digest, execute, record_report
from fishhighz.validation.schema import (
    canonical,
    metric_values,
    token,
    validate_payload,
)
from fishhighz.validation.study import study
from fishhighz.validation.trials import trial_id

CASE = "lya_qso_lbg_lae_15x2pt"


def derive(arrays, report):
    """Replay historical attempts; ambiguous/missing failure reconstruction rejects."""
    old = copy.deepcopy(report)
    controls = old["study_controls"]
    by_id = {trial_id(c): i for i, c in enumerate(controls)}
    failures = iter(old["unresolved_controls"])

    class HistoricalTrials:
        selection = SimpleNamespace(
            fields=[SimpleNamespace(kind="forest")], selected_pairs=np.array([[0, 0]])
        )
        _prepared = {}

        def evaluate(self, task, c):
            identity = trial_id(c)
            if identity not in by_id:
                try:
                    failure = next(failures)
                except StopIteration as error:
                    raise KeyError("unrecorded historical trial failure") from error
                raise ValueError(failure["error"])
            i = by_id[identity]
            return dict(
                fisher=arrays["study_fisher"][i],
                pair_fisher=arrays["study_pair_fisher"][i],
            ), dict(settings=dict(grid=dict(volume=arrays["study_volume"][i])))

    saved, replay = study(HistoricalTrials(), {})
    for name in (
        "study_controls",
        "actual_levels",
        "final_controls",
        "combined_lower_controls",
        "unresolved_controls",
        "metric_names",
    ):
        if replay[name] != old[name]:
            raise ValueError(f"historical replay mismatch: {name}")
    for name in (
        "study_fisher",
        "study_pair_fisher",
        "study_volume",
        "metric_fisher",
        "metric_pair_fisher",
        "metric_volume",
    ):
        np.testing.assert_allclose(arrays[name], saved[name], rtol=5e-12, atol=0)
    report["trial_contract"] = replay["trial_contract"]
    # Cached preparation could retain another derivative-step label. Bind the
    # effective report to the saved, independently verified final actual controls.
    report["settings"]["controls"] = dict(report["final_controls"])
    report["effective_hash"] = canonical(report["settings"])
    arrays["effective_token"] = token(report["settings"])
    report["derived_evidence"] = dict(
        historical_schema=2,
        failure_controls="inferred by exact bounded schedule replay; original errors and successful arrays matched",
        checker_sha256=digest(Path(__file__)),
        original_settings_hash=old["effective_hash"],
    )
    return arrays, report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    root = Path(args.source).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    records = {
        canonical(r["task"]): r
        for g in ("records", "diagnostics")
        for r in manifest[g]
        if r["task"]["case"] == CASE
    }
    diagnostics = [t for t in manifest["diagnostics_requested"] if t["case"] == CASE]
    checks = []

    def worker(task):
        rec = records[canonical(task)]
        if "arrays" not in rec:
            raise ValueError("Historical unavailable diagnostic: " + rec["error"])
        path = root / rec["arrays"]
        assert path.parent == root and digest(path) == rec["sha256"]
        report = record_report(root, rec)
        with np.load(path, allow_pickle=False) as data:
            arrays = {n: data[n] for n in data.files}
        validate_payload(task, arrays, report, require_pass=False, schema=2)
        if task["kind"] == "real_bao" and task["profile"] == "accuracy":
            arrays, report = derive(arrays, report)
            validate_payload(task, arrays, report, require_pass=False)
            bad = {k: v.copy() for k, v in arrays.items()}
            changed = copy.deepcopy(report)
            for name in ("metric_fisher", "metric_pair_fisher", "metric_volume"):
                bad[name][:, 0] = bad[name][:, 1]
            changed.update(passed=True, unresolved_controls=[])
            changed["metrics"] = metric_values(bad, changed["metric_names"])
            try:
                validate_payload(task, bad, changed)
            except ValueError as error:
                checks.append(
                    dict(
                        bin=task["bin"],
                        original_passed=report["passed"],
                        duplicate_operands_rejected=str(error),
                    )
                )
            else:
                raise AssertionError("historical R2 mutation accepted")
        report["historical_source"] = dict(
            path=str(path),
            sha256=rec["sha256"],
            manifest_sha256=digest(root / "manifest.json"),
        )
        return arrays, report

    result = execute(
        args.output,
        suite="full",
        cases=[CASE],
        profiles=("compatibility", "accuracy"),
        worker=worker,
        inputs=[root / "manifest.json"],
        diagnostic_requests=diagnostics,
    )
    (Path(args.output) / "historical-binding-check.json").write_text(
        json.dumps(checks, indent=2)
    )
    print(
        json.dumps(
            dict(
                primary=len(result["records"]),
                diagnostics=len(result["diagnostics"]),
                checked_accuracy=len(checks),
                statuses=[r["status"] for r in result["records"]],
                complete=result["complete"],
            )
        )
    )


if __name__ == "__main__":
    main()
