"""S2 selected bin-1 rectangular refinement from immutable stopped states."""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from fishhighz.validation.compatibility_bao import forecast
from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import contract, relative
from fishhighz.validation.profile_definitions import REVISION
from fishhighz.validation.schema import request


def check(root):
    start = time.monotonic()
    root = Path(root)
    stem = root / "step12-r5-20260914T191855Z/profiles-checked/records-000"
    source = json.loads(stem.with_suffix(".report.json").read_text())
    with np.load(stem.with_suffix(".npz")) as data:
        original = {k: data[k] for k in data.files}
    task = request(
        source["context"]["case"], 0, "fixed-compatibility", recipe_revision=REVISION
    )
    old_task = source["context"]
    selected = [old_task["selected_pairs"].index(p) for p in task["selected_pairs"]]
    required = [old_task["required_pairs"].index(p) for p in task["required_pairs"]]
    arrays = {
        **original,
        "total": original["total"][:, required],
        "observed_j": original["observed_j"][:, selected],
    }
    stage2_path = root / "compatibility-weighting-stage2/bin-0.json"
    stage3_path = root / "compatibility-weighting-stage3/bin-0.json"
    stage2 = json.loads(stage2_path.read_text())
    stage3 = json.loads(stage3_path.read_text())
    assert (
        stage2["report_sha256"]
        == stage3["source_sha256"]
        == digest(stem.with_suffix(".report.json"))
    )
    assert stage3["stage2_json_sha256"] == digest(stage2_path)
    assert stage3["stage2_npz_sha256"] == digest(stage2_path.with_suffix(".npz"))
    name = "lya(qso)_lya(qso)"
    records = []
    preparations = []
    checks = []
    with np.load(stage2_path.with_suffix(".npz")) as trajectories:
        for nodes in (107, 213, 425):
            stopped = next(
                r
                for r in stage3["records"]
                if r["nodes"] == nodes
                and r["context"] == name
                and r["variant"] == "sum_historical"
                and r["rtol"] == 1e-4
            )
            assert stopped["status"] == "converged"
            coefficients = trajectories[f"{name}/{nodes}/sum_historical/coefficients"][
                stopped["updates"]
            ]
            np.testing.assert_array_equal(coefficients, stopped["coefficients"])
            singles, joint, total, cov, valid = forecast(
                arrays, task, source["settings"]["pair_inputs"], {name: coefficients}
            )
            assert np.all(valid) and joint is not None and cov.shape[1:] == (3, 3)
            direct = contract(cov, arrays["observed_j"], independent=True)[0]
            difference = relative(direct, joint["fisher"])
            assert difference < 5e-12
            # Literal scalar Wick contractions for all cells; three retained pairs.
            lookup = {tuple(p): i for i, p in enumerate(task["required_pairs"])}

            def power(a, b):
                return total[:, lookup[tuple(sorted((a, b)))]]

            scalar = np.empty_like(cov)
            for s, (a, b) in enumerate(task["selected_pairs"]):
                for t, (c, d) in enumerate(task["selected_pairs"]):
                    scalar[:, s, t] = (
                        power(a, c) * power(b, d) + power(a, d) * power(b, c)
                    ) / arrays["modes"]
            np.testing.assert_array_equal(scalar, cov)
            for p, item in enumerate(singles + [joint]):
                label = (
                    "joint"
                    if p == 3
                    else "_".join(
                        task["fields"][i]["id"] for i in task["selected_pairs"][p]
                    )
                )
                records.append(
                    dict(
                        nodes=nodes,
                        spectrum=label,
                        errors=item["errors"].tolist(),
                        fisher=item["fisher"].tolist(),
                    )
                )
                if p < 3 and all(
                    task["fields"][i]["kind"] == "galaxy"
                    for i in task["selected_pairs"][p]
                ):
                    np.testing.assert_array_equal(
                        item["fisher"], original["pair_fisher"][selected[p]]
                    )
            preparations.append(
                dict(
                    nodes=nodes,
                    field="lya(qso)",
                    status=stopped["status"],
                    updates=stopped["updates"],
                    coefficients=coefficients.tolist(),
                )
            )
            checks.append(
                dict(
                    nodes=nodes,
                    direct_joint_relative=difference,
                    scalar_covariance_exact=True,
                )
            )
    changes = []
    for spectrum in [r["spectrum"] for r in records if r["nodes"] == 107]:
        for low, high in ((107, 213), (213, 425), (107, 425)):
            values = [
                next(
                    r for r in records if r["nodes"] == n and r["spectrum"] == spectrum
                )["errors"]
                for n in (low, high)
            ]
            changes.append(
                dict(
                    spectrum=spectrum,
                    low=low,
                    high=high,
                    fractional_change=(np.array(values[1]) / values[0] - 1).tolist(),
                )
            )
    reuse = []
    for index in range(1, 6):
        path = root / f"early-lyaforecast-grid-bao/bin-{index}.json"
        prior = json.loads(path.read_text())
        checked = {}
        for source_path, sha in prior["sources"].items():
            source_path = Path(source_path)
            # Numerical operands and prior reports; old code hashes stay historical.
            if source_path.parts[0] == ".validation":
                actual = root.joinpath(*source_path.parts[1:])
                assert digest(actual) == sha
                checked[str(actual)] = sha
        assert checked and prior["checks"]["original_grid_stage4_bitwise"]
        reuse.append(
            dict(
                bin=index,
                path=str(path),
                sha256=digest(path),
                unchanged_operands=checked,
                reason="Identical all-15 selection and immutable numerical operands; S1 fixed original-grid replay bitwise",
            )
        )
    return dict(
        passed=True,
        bin=0,
        task=task,
        records=records,
        preparations=preparations,
        changes=changes,
        checks=checks,
        reused_bins=reuse,
        sources={
            str(p): digest(p)
            for p in (
                stem.with_suffix(".report.json"),
                stem.with_suffix(".npz"),
                stage2_path,
                stage3_path,
                stage2_path.with_suffix(".npz"),
            )
        },
        seconds=time.monotonic() - start,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-root", default=".validation")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = check(args.validation_root)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            dict(
                passed=result["passed"],
                seconds=result["seconds"],
                changes=result["changes"],
            ),
            indent=2,
        )
    )
