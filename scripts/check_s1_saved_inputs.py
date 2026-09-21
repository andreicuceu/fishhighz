"""S1 saved-input controls, without raw interpolation or new model forecasts."""

import argparse
import json
from pathlib import Path

import numpy as np

from fishhighz.kernels.full_sum_weights import solve
from fishhighz.validation.compatibility_weights import WeightInputs
from fishhighz.validation.numerics import contract, relative
from fishhighz.validation.profile_definitions import REVISION
from fishhighz.validation.revised_compatibility import run
from fishhighz.validation.schema import request


def check(root):
    root = Path(root)
    bundle = root / "step12-r5-20260914T191855Z/profiles-checked"
    rows = []
    for index in range(6):
        stem = bundle / f"records-{2 * index:03d}"
        report = json.loads(stem.with_suffix(".report.json").read_text())
        stage3 = json.loads(
            (root / f"compatibility-weighting-stage3/bin-{index}.json").read_text()
        )
        stage4 = json.loads(
            (root / f"compatibility-weighting-stage4/bin-{index}.json").read_text()
        )
        case = report["context"]["case"]
        task = request(case, index, "fixed-compatibility", recipe_revision=REVISION)
        controls = []
        for i, j in task["required_pairs"]:
            field = task["fields"][i]
            if i != j or field["kind"] != "forest":
                continue
            name = field["id"] + "_" + field["id"]
            inputs = WeightInputs.from_pair(report["settings"]["pair_inputs"][name])
            for variant in ("sum_historical", "sum_aliasing"):
                state = solve(inputs, variant)
                old = next(
                    r
                    for r in stage3["records"]
                    if r["nodes"] == 107
                    and r["context"] == name
                    and r["variant"] == variant
                    and r["rtol"] == 1e-4
                )
                assert (
                    state["status"] == old["status"]
                    and state["updates"] == old["updates"]
                )
                np.testing.assert_array_equal(
                    state["coefficients"], old["coefficients"]
                )
                with np.load(
                    root / f"compatibility-weighting-stage2/bin-{index}.npz"
                ) as trajectory:
                    np.testing.assert_array_equal(
                        state["weights"],
                        trajectory[f"{name}/107/{variant}/weights"][
                            state["state_updates"]
                        ],
                    )
                controls.append(
                    dict(
                        field=field["id"],
                        variant=variant,
                        status=state["status"],
                        updates=state["updates"],
                    )
                )
        a, _ = run(task, bundle, {})
        direct = contract(a["selected_covariance"], a["observed_j"], independent=True)[
            0
        ]
        assert relative(direct, a["fisher"]) < 5e-12
        full, _ = run(
            request(case, index, "full-compatibility", recipe_revision=REVISION),
            bundle,
            {},
        )
        with np.load(stem.with_suffix(".npz")) as old:
            for j, pair in enumerate(task["selected_pairs"]):
                old_j = report["context"]["selected_pairs"].index(pair)
                np.testing.assert_array_equal(
                    full["pair_fisher"][j], old["pair_fisher"][old_j]
                )
                if all(task["fields"][i]["kind"] == "galaxy" for i in pair):
                    np.testing.assert_array_equal(
                        a["pair_fisher"][j], old["pair_fisher"][old_j]
                    )
                name = "_".join(task["fields"][i]["id"] for i in pair)
                saved = next(
                    r
                    for r in stage4["records"]
                    if r["variant"] == "sum_historical"
                    and r["mode"] == "converged_1e-4"
                    and r["spectrum"] == name
                )
                np.testing.assert_array_equal(a["pair_fisher"][j], saved["fisher"])
        if index:
            saved = next(
                r
                for r in stage4["records"]
                if r["variant"] == "sum_historical"
                and r["mode"] == "converged_1e-4"
                and r["spectrum"] == "joint"
            )
            np.testing.assert_array_equal(a["fisher"], saved["fisher"])
        rows.append(
            dict(
                bin=index,
                selected=len(task["selected_pairs"]),
                weights=controls,
                joint_solve_relative=relative(direct, a["fisher"]),
                errors=a["errors"].tolist(),
            )
        )
    return dict(passed=True, bins=rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--validation-root", default=".validation")
    p.add_argument("--output", required=True)
    args = p.parse_args()
    result = check(args.validation_root)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))
