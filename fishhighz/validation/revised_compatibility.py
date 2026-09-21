"""Selected compatibility forecasts from immutable captured powers/Jacobians."""

import json
from pathlib import Path

import numpy as np

from ..adapters.legacy_compat import plain
from ..kernels.full_sum_weights import METHODS, solve
from .compatibility_weights import WeightInputs, forest_noise
from .evidence import digest, record_report
from .profile_definitions import STOPPING, identity
from .schema import assemble


def run(task, bundle, provenance):
    """Keep the literal 107-node measure; replace only required forest auto noise."""
    root = Path(bundle)
    manifest = json.loads((root / "manifest.json").read_text())
    record = next(
        row
        for row in manifest["records"]
        if row["task"]["case"] == task["case"]
        and row["task"]["bin"] == task["bin"]
        and row["task"]["profile"] == "compatibility"
    )
    source = record_report(root, record)
    path = root / record["arrays"]
    if digest(path) != record["sha256"]:
        raise ValueError("captured compatibility array hash mismatch")
    with np.load(path, allow_pickle=False) as data:
        old = {key: data[key] for key in data.files}
    old_task = source["context"]
    required = [old_task["required_pairs"].index(p) for p in task["required_pairs"]]
    selected = [old_task["selected_pairs"].index(p) for p in task["selected_pairs"]]
    total = old["total"][:, required].copy()
    jacobian = old["observed_j"][:, selected]
    settings = dict(source["settings"])
    settings.update(profile=task["profile"], recipe_identity=identity(task["profile"]))
    states = {}
    fixed = task["profile"] == "fixed-compatibility"
    if not fixed:
        for i, j in task["required_pairs"]:
            field = task["fields"][i]
            if i == j and field["kind"] == "forest":
                row = settings["pair_inputs"][field["id"] + "_" + field["id"]]
                states[field["id"]] = dict(
                    status="fixed_count",
                    updates=3,
                    weights=row["_w_lya"],
                    A=row["_aliasing_weights"][-1],
                    P_pixel=row["_effective_noise_power"][-1],
                )
    if fixed:
        rows = settings["pair_inputs"]
        for column, (i, j) in enumerate(task["required_pairs"]):
            field = task["fields"][i]
            if i != j or field["kind"] != "forest":
                continue
            name = field["id"] + "_" + field["id"]
            row = rows[name]
            inputs = WeightInputs.from_pair(row)
            if len(inputs.magnitudes) != 107:
                raise ValueError("primary fixed-compatibility requires 107 nodes")
            state = solve(inputs, METHODS["early_lyaforecast"], **STOPPING)
            states[field["id"]] = plain(state)
            parallel = 0.00035 * row["_distance_to_velocity"]
            k = float(np.hypot(parallel, 2.4 / row["_angle_to_distance"]))
            states[field["id"]]["auxiliary"] = dict(
                k=k,
                mu=parallel / k,
                P=inputs.signal,
                B=inputs.p1d,
                k_t_deg=2.4,
                k_p_velocity=0.00035,
            )
            if state["status"] != "converged":
                raise ValueError(f"{name}: {state['status']}; no qualified forecast")
            baseline = forest_noise(
                row,
                old["k"],
                old["mu"],
                row["_aliasing_weights"][-1],
                row["_effective_noise_power"][-1],
            )
            changed = forest_noise(
                row, old["k"], old["mu"], *state["coefficients"][-2:]
            )
            total[:, column] = total[:, column] - baseline + changed
    settings["forest_weighting"] = dict(
        method="early_lyaforecast" if fixed else "legacy",
        stopping=STOPPING if fixed else None,
        iterations=None if fixed else 3,
        forests=states,
        quadrature="original rectangular 107; full endpoints",
        input_policy="literal signed compatibility",
    )
    arrays, report = assemble(task, total, jacobian, settings)
    report.update(
        provenance=provenance,
        captured_input=dict(
            path=str(path.resolve()),
            sha256=digest(path),
            report_hash=record["effective_hash"],
        ),
        selection_exclusions=[
            p for p in old_task["selected_pairs"] if p not in task["selected_pairs"]
        ],
        interpretation="Selected Wick covariance recomputed; historical joint Fisher is not reused",
    )
    return arrays, report
