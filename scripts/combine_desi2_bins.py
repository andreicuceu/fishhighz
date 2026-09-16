"""Check independent-bin combinations from saved per-bin Fisher information."""

import argparse
import json
from itertools import islice
from pathlib import Path

import numpy as np

from fishhighz.validation.cases import CASE_IDS, bins
from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import information, relative
from fishhighz.validation.plots import records
from fishhighz.validation.schema import LIMITS


def block(matrices):
    result = np.zeros((2 * len(matrices), 2 * len(matrices)))
    for i, matrix in enumerate(matrices):
        result[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = matrix
    return result


def combine(bundle, output, cases=CASE_IDS):
    root = Path(bundle)
    source = json.loads((root / "manifest.json").read_text())
    groups = {}
    for record, arrays in islice(records(root), len(source["records"])):
        task = record["task"]
        compact = (
            None
            if arrays is None
            else {
                key: arrays[key]
                for key in (
                    "fisher",
                    "errors",
                    "reference_fisher",
                    "metric_fisher",
                    "metric_volume",
                )
                if key in arrays
            }
        )
        metadata = dict(
            task=task,
            scientific_passed=record.get("scientific_passed", False),
            report={
                key: record.get("report", {}).get(key)
                for key in ("metric_names", "metrics")
            },
        )
        groups.setdefault((task["case"], task["profile"]), []).append(
            (metadata, compact)
        )
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    reports, payload = [], {}
    for case in cases:
        for profile in ("compatibility", "accuracy"):
            rows = groups[(case, profile)]
            if [r["task"]["bin"] for r, _ in rows] != list(range(len(bins(case)))):
                raise ValueError("combined-bin check requires every original bin")
            if any(a is None for _, a in rows):
                raise ValueError("combined-bin information is unavailable")
            fisher = block([a["fisher"] for _, a in rows])
            info = information(fisher)
            np.testing.assert_allclose(
                info["errors"],
                np.concatenate([a["errors"] for _, a in rows]),
                rtol=5e-12,
                atol=0,
            )
            row = dict(
                case=case,
                profile=profile,
                parameters=[p for r, _ in rows for p in r["task"]["parameters"]],
                primary_passed=all(r["scientific_passed"] for r, _ in rows),
                metrics={},
            )
            prefix = case + "__" + profile + "__"
            row["array_prefix"] = prefix
            payload[prefix + "fisher"] = fisher
            payload.update({prefix + k: np.asarray(v) for k, v in info.items()})
            if profile == "compatibility":
                reference = block([a["reference_fisher"] for _, a in rows])
                original = information(reference)
                row["comparison"] = dict(
                    fisher_relative=relative(fisher, reference),
                    error_relative=float(
                        np.max(abs(info["errors"] / original["errors"] - 1))
                    ),
                    correlation_absolute=float(
                        np.max(abs(info["correlation"] - original["correlation"]))
                    ),
                )
                payload[prefix + "reference_fisher"] = reference
                assert row["comparison"]["fisher_relative"] <= 5e-12
                assert row["comparison"]["error_relative"] <= 1e-6
            else:
                for control in (
                    "k",
                    "mu",
                    "magnitude",
                    "volume",
                    "step",
                    "weights",
                    "combined",
                ):
                    indices = [
                        r["report"]["metric_names"].index(control) for r, _ in rows
                    ]
                    paired = [
                        block(
                            [
                                a["metric_fisher"][i, j]
                                for (_, a), i in zip(rows, indices)
                            ]
                        )
                        for j in (0, 1)
                    ]
                    volume = [
                        sum(
                            a["metric_volume"][i, j] for (_, a), i in zip(rows, indices)
                        )
                        for j in (0, 1)
                    ]
                    lower, upper = map(information, paired)
                    metric = dict(
                        fisher_relative=relative(*paired),
                        error_relative=float(
                            np.max(abs(lower["errors"] / upper["errors"] - 1))
                        ),
                        pair_error_relative=max(
                            r["report"]["metrics"][i]["pair_error_relative"]
                            for (r, _), i in zip(rows, indices)
                        ),
                        volume_relative=abs(float(volume[0] / volume[1]) - 1),
                    )
                    row["metrics"][control] = dict(
                        values=metric,
                        passed=all(metric[k] <= LIMITS[k] for k in LIMITS),
                    )
                    payload[prefix + control + "_paired_fisher"] = np.asarray(paired)
                    payload[prefix + control + "_volume"] = np.asarray(volume)
            reports.append(row)
    arrays_path = out / "matrices.npz"
    np.savez_compressed(arrays_path, **payload)
    report = dict(
        kind="derived-independent-bin-comparison",
        source_schema=source["schema"],
        source_bundle=str(root.resolve()),
        interpretation="Each case combines its original independent redshift bins with distinct ap/at parameters and zero cross-bin information. Per-bin marginalized errors must be preserved. Accuracy control comparisons combine the recorded final-two levels of each bin; unavailable/stability findings remain primary failures.",
        source_sha256=digest(root / "manifest.json"),
        script_sha256=digest(__file__),
        arrays=arrays_path.name,
        arrays_sha256=digest(arrays_path),
        records=reports,
    )
    (out / "manifest.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(f"Checked {len(reports)} case/profile independent-bin combinations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cases", nargs="+", choices=CASE_IDS, default=CASE_IDS)
    args = parser.parse_args()
    combine(args.bundle, args.output, args.cases)
