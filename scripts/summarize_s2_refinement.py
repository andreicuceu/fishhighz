"""Summarize actual S2 per-spectrum/combined finite-refinement operands."""

import argparse
import json
from pathlib import Path

import numpy as np

from fishhighz.validation.evidence import record_report
from fishhighz.validation.numerics import information
from fishhighz.validation.schema import LIMITS


def summarize(bundle):
    bundle = Path(bundle)
    manifest = json.loads((bundle / "manifest.json").read_text())
    result = []
    for record in manifest["records"]:
        if record["task"]["profile"] != "accuracy":
            continue
        row = dict(
            bin=record["task"]["bin"],
            status=record.get("numerical_status"),
            execution_status=record["status"],
        )
        if record["status"] != "completed":
            row["error"] = record.get("error")
            result.append(row)
            continue
        report = record_report(bundle, record)
        task = record["task"]
        labels = [
            "_".join(task["fields"][i]["id"] for i in p) for p in task["selected_pairs"]
        ]
        row.update(
            final_controls=report["final_controls"],
            actual_levels=report["actual_levels"],
            unresolved=report["unresolved_controls"],
            families=[],
            forests={},
        )
        with np.load(bundle / record["arrays"]) as a:
            row["errors"] = a["errors"].tolist()
            row["pair_errors"] = {
                name: errors.tolist() for name, errors in zip(labels, a["pair_errors"])
            }
            for i, (name, metrics) in enumerate(
                zip(report["metric_names"], report["metrics"])
            ):
                joint = np.array(
                    [information(f)["errors"] for f in a["metric_fisher"][i]]
                )
                pairs = np.array(
                    [
                        [information(f)["errors"] for f in side]
                        for side in a["metric_pair_fisher"][i]
                    ]
                )
                joint_change = joint[1] / joint[0] - 1
                pair_change = pairs[1] / pairs[0] - 1
                maximum = np.unravel_index(
                    np.argmax(np.abs(pair_change)), pair_change.shape
                )
                row["families"].append(
                    dict(
                        name=name,
                        passed=all(metrics[k] <= LIMITS[k] for k in LIMITS),
                        metrics=metrics,
                        joint_fractional_change=joint_change.tolist(),
                        pair_fractional_change={
                            label: values.tolist()
                            for label, values in zip(labels, pair_change)
                        },
                        largest_pair=dict(
                            spectrum=labels[maximum[0]],
                            component=["parallel", "transverse"][maximum[1]],
                            fractional_change=float(pair_change[maximum]),
                        ),
                    )
                )
        for outcome in report["trial_contract"]["outcomes"]:
            if outcome["outcome"] != "success":
                continue
            for name, field in (
                outcome.get("forest_weighting", {}).get("forests", {}).items()
            ):
                row["forests"].setdefault(name, []).append(
                    dict(
                        controls=outcome["controls"],
                        status=field["result"]["status"],
                        updates=field["result"]["updates"],
                        candidate=field["result"]["candidate"],
                        last_step=field["result"]["last_step"],
                        confirmation=field["result"]["confirmation"],
                        A=field["A"],
                        P_pixel=field["P_pixel"],
                        auxiliary=field["auxiliary"],
                    )
                )
        result.append(row)
    return dict(
        execution_finished=manifest["execution_finished"],
        complete=manifest["complete"],
        bins=result,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    summary = summarize(args.bundle)
    Path(args.output).write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            [
                {
                    k: v
                    for k, v in row.items()
                    if k not in ("forests", "pair_errors", "families")
                }
                for row in summary["bins"]
            ],
            indent=2,
        )
    )
