"""Bounded offline stage-1 replay of six saved compatibility bins; no new forecast."""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from fishhighz.validation.compatibility_weights import reassemble
from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import contract, relative, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    results = []
    for index in range(6):
        stem = args.input / f"records-{2 * index:03d}"
        report_path = stem.with_suffix(".report.json")
        array_path = stem.with_suffix(".npz")
        report = json.loads(report_path.read_text())
        with np.load(array_path, allow_pickle=False) as data:
            arrays = {name: data[name] for name in data.files}
        task = report["context"]
        if task["profile"] != "compatibility" or task["bin"] != index:
            raise ValueError("unexpected saved baseline identity")
        rows = report["settings"]["pair_inputs"]
        total, covariance, states, noise = reassemble(arrays, task, rows)
        checks = {}
        for name, state in states.items():
            row = rows[name]
            checks[name] = dict(
                weights_relative=relative(state["weights"], row["_w_lya"]),
                coefficients_relative=relative(
                    state["coefficients"][-2:],
                    [row["_aliasing_weights"][-1], row["_effective_noise_power"][-1]],
                ),
                negative_density_nodes=int(
                    np.count_nonzero(np.asarray(row["density"]) < 0)
                ),
                auxiliary_signal=row["auxiliary_signal"],
            )
        decomposition = {}
        for column, (i, j) in enumerate(task["required_pairs"]):
            if task["fields"][i]["kind"] != "forest":
                continue
            name = task["fields"][i]["id"] + "_" + task["fields"][j]["id"]
            selected = task["selected_pairs"].index([i, j])
            intrinsic = (
                noise[name]["intrinsic"] if i == j else arrays["total"][:, column]
            )
            decomposition[name] = relative(intrinsic, arrays["mean"][:, selected])
        f, pair_f = contract(covariance, arrays["observed_j"])
        recovered = summaries(f, pair_f)
        metrics = {
            name: relative(recovered[name], arrays[name])
            for name in ("fisher", "pair_fisher", "errors", "pair_errors")
        }
        metrics.update(
            total=relative(total, arrays["total"]),
            covariance=relative(covariance, arrays["selected_covariance"]),
        )
        if (
            max(
                v
                for row in checks.values()
                for k, v in row.items()
                if k.endswith("relative")
            )
            > 5e-12
        ):
            raise ValueError("saved weights/coefficient mismatch")
        if max(metrics.values()) > 1e-10:
            raise ValueError(
                f"baseline decomposition/forecast mismatch: {decomposition}, {metrics}"
            )
        results.append(
            dict(
                bin=index,
                context_count=len(states),
                contexts=checks,
                intrinsic_mean_relative=decomposition,
                mean_z=report["settings"]["mean_z"],
                covariance_z=report["settings"]["covariance_z"],
                metrics=metrics,
                source_report=str(report_path.resolve()),
                report_sha256=digest(report_path),
                arrays_sha256=digest(array_path),
            )
        )
    result = dict(
        passed=True,
        seconds=time.monotonic() - start,
        bins=results,
        scope="saved three-update compatibility control only; no variant BAO forecasts",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            dict(
                passed=True,
                contexts=sum(r["context_count"] for r in results),
                seconds=result["seconds"],
            )
        )
    )


if __name__ == "__main__":
    main()
