"""Schema-3 convergence operands bound to ordered, actual refinement trials."""

import numpy as np

from .schema import LIMITS, _close, _close_blocks, canonical

CONTROL_KEYS = dict(
    k="k_intervals",
    mu="mu_order",
    magnitude="magnitude_order",
    volume="z_order",
    step="step",
    weights="iterations",
)
LEVELS = dict(
    k=[32, 64, 128, 256, 512],
    mu=[8, 16, 32, 64, 128],
    magnitude=[4, 8, 16, 32, 64],
    volume=[8, 16, 32, 64, 128],
    step=[1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5],
    weights=[3, 6, 12, 24],
)
FIXED_CONTROL_KEYS = {
    name: key for name, key in CONTROL_KEYS.items() if name != "weights"
}
FIXED_LEVELS = {name: values for name, values in LEVELS.items() if name != "weights"}
ADAPTIVE_CONTROL_KEYS = dict(FIXED_CONTROL_KEYS, weights="weight_rtol")
ADAPTIVE_LEVELS = dict(FIXED_LEVELS, weights=[1e-4, 1e-5])
FIXED_REFERENCE = dict(
    convention="intrinsic_p1d_times_field_response_squared",
    q_star=0.00035,
    q_star_units="s/km",
)


def trial_id(controls):
    return canonical(controls)


def bind(arrays, report, outcomes, *, method="legacy"):
    """Persist references without substituting or recomputing numerical operands."""
    if method == "legacy":
        version = 1
        control_keys = CONTROL_KEYS
        levels = LEVELS
    elif method == "inverse_variance":
        version = 2
        control_keys = FIXED_CONTROL_KEYS
        levels = FIXED_LEVELS
    elif method in ("early_lyaforecast", "mcdonald"):
        version = 3
        control_keys = ADAPTIVE_CONTROL_KEYS
        levels = ADAPTIVE_LEVELS
    else:
        raise ValueError("unknown accuracy forest-weight method")
    final = report["final_controls"]
    references = []
    for name in report["metric_names"]:
        if name == "combined":
            controls = [report["combined_lower_controls"], final]
        else:
            key = control_keys[name]
            values = report["actual_levels"][name]
            controls = [{**final, key: value} for value in values[-2:]]
        references.append([trial_id(c) for c in controls])
    report["trial_contract"] = dict(
        version=version,
        allowed_levels=levels,
        outcomes=outcomes,
        successful_ids=[trial_id(c) for c in report["study_controls"]],
        final_id=trial_id(final),
        metric_trials=references,
    )
    if version == 2:
        report["trial_contract"].update(
            method=method,
            fixed_reference=FIXED_REFERENCE,
        )
    if version == 3:
        from .profile_definitions import REFERENCE, REVISION, STOPPING

        report["trial_contract"].update(
            method=method,
            reference=REFERENCE,
            stopping=STOPPING,
            recipe_revision=REVISION,
        )


def validate(arrays, report):
    """Check identities and controls before any convergence metric is evaluated."""
    contract = report.get("trial_contract", {})
    version = contract.get("version")
    if version == 1:
        control_keys = CONTROL_KEYS
        levels = LEVELS
        method = "legacy"
        weighting = report.get("settings", {}).get("forest_weighting")
        if weighting is not None and weighting.get("method") != method:
            raise ValueError("version-1 trials require explicit legacy weighting")
    elif version == 2:
        control_keys = FIXED_CONTROL_KEYS
        levels = FIXED_LEVELS
        method = "inverse_variance"
        weighting = report.get("settings", {}).get("forest_weighting", {})
        if (
            contract.get("method") != method
            or contract.get("fixed_reference") != FIXED_REFERENCE
            or weighting.get("method") != method
            or weighting.get("reference") != FIXED_REFERENCE
            or weighting.get("iterations")
            != {"applicable": False, "status": "inapplicable"}
        ):
            raise ValueError("version-2 trials require the fixed-reference convention")
    elif version == 3:
        from .profile_definitions import REFERENCE, REVISION, STOPPING

        control_keys, levels = ADAPTIVE_CONTROL_KEYS, ADAPTIVE_LEVELS
        method = contract.get("method")
        weighting = report.get("settings", {}).get("forest_weighting", {})
        if (
            method not in ("early_lyaforecast", "mcdonald")
            or weighting.get("method") != method
            or contract.get("reference") != REFERENCE
            or contract.get("stopping") != STOPPING
            or contract.get("recipe_revision") != REVISION
        ):
            raise ValueError("version-3 trials require the adaptive recipe")
        for row in weighting.get("forests", {}).values():
            if row.get("result", {}).get("status") != "converged":
                raise ValueError(
                    "unconverged adaptive weights cannot supply a forecast"
                )
    else:
        raise ValueError("missing or changed predeclared trial contract")
    if contract.get("allowed_levels") != levels:
        raise ValueError("missing or changed predeclared trial contract")
    names = list(control_keys) + ["combined"]
    if report["metric_names"] != names:
        raise ValueError("wrong ordered refinement families")
    outcomes = contract.get("outcomes")
    if not isinstance(outcomes, list) or not outcomes:
        raise ValueError("missing attempted trial outcomes")
    by_id = {}
    successful = []
    failed = []
    for row in outcomes:
        controls = row.get("controls", {})
        if set(controls) != set(control_keys.values()):
            raise ValueError("wrong actual trial controls")
        for family, key in control_keys.items():
            if isinstance(controls[key], bool) or controls[key] not in levels[family]:
                raise ValueError("trial outside bounded refinement levels")
        identity = trial_id(controls)
        if row.get("id") != identity or identity in by_id:
            raise ValueError("wrong/repeated trial identity")
        by_id[identity] = row
        if row.get("outcome") == "success":
            if row.get("array_index") != len(successful):
                raise ValueError("missing/reordered successful trial")
            successful.append(identity)
        elif (
            row.get("outcome") == "failed"
            and isinstance(row.get("error"), str)
            and row["error"]
        ):
            if "array_index" in row:
                raise ValueError("failed trial cannot supply numerical operands")
            failed.append(identity)
        else:
            raise ValueError("invalid trial outcome")
    if successful != contract.get("successful_ids"):
        raise ValueError("suppressed/reordered trial inventory")
    controls = report.get("study_controls", [])
    if [by_id[i]["controls"] for i in successful] != controls:
        raise ValueError("trial controls/array ordering mismatch")
    n = len(successful)
    for name, shape in [
        ("study_fisher", (n, *arrays["fisher"].shape)),
        ("study_pair_fisher", (n, *arrays["pair_fisher"].shape)),
        ("study_volume", (n,)),
    ]:
        if name not in arrays or arrays[name].shape != shape:
            raise ValueError("wrong trial numerical dimensions")
    if np.any(arrays["study_volume"] <= 0):
        raise ValueError("nonpositive trial volume")
    from .numerics import information

    for f in arrays["study_fisher"]:
        information(f)
    for f in arrays["study_pair_fisher"].reshape(-1, *arrays["fisher"].shape):
        information(f)
    final = report.get("final_controls")
    final_id = trial_id(final)
    if contract.get("final_id") != final_id or final_id not in successful:
        raise ValueError("primary final trial missing or failed")
    if report["settings"].get("controls") != final:
        raise ValueError("primary settings differ from final controls")
    grid = report["settings"]["grid"]
    if (
        grid["k_intervals"] != final["k_intervals"]
        or grid["mu_order"] != final["mu_order"]
        or grid["k_order"] != 4
    ):
        raise ValueError("primary grid differs from final controls")
    index = by_id[final_id]["array_index"]
    for a, b in [
        (arrays["fisher"], arrays["study_fisher"][index]),
        (grid["volume"], arrays["study_volume"][index]),
    ]:
        _close(a, b, "primary/final trial")
    _close_blocks(
        arrays["pair_fisher"],
        arrays["study_pair_fisher"][index],
        "primary/final pair trial",
    )
    refs = contract.get("metric_trials")
    if not isinstance(refs, list) or len(refs) != len(names):
        raise ValueError("missing metric-to-trial references")
    lower = dict(final)
    for name, key in control_keys.items():
        values = report.get("actual_levels", {}).get(name, [])
        if len(values) < 2 or values != [v for v in levels[name] if v in values]:
            raise ValueError("missing/repeated/reordered refinement levels")
        if values[-1] != final[key]:
            raise ValueError("metric upper control is not final")
        lower[key] = values[-2]
        expected_values = [
            v for v in levels[name] if trial_id({**final, key: v}) in successful
        ]
        if values != expected_values:
            raise ValueError("actual levels do not match successful isolated trials")
    if report.get("combined_lower_controls") != lower:
        raise ValueError("combined lower controls not bound to isolated lower trials")
    for i, name in enumerate(names):
        expected_controls = (
            [lower, final]
            if name == "combined"
            else [
                {**final, control_keys[name]: v}
                for v in report["actual_levels"][name][-2:]
            ]
        )
        expected_ids = [trial_id(c) for c in expected_controls]
        if refs[i] != expected_ids or len(set(refs[i])) != 2:
            raise ValueError("wrong/repeated/unrelated metric trial selection")
        for side, identity in enumerate(refs[i]):
            if identity not in by_id:
                raise ValueError("missing metric trial")
            row = by_id[identity]
            if row["outcome"] != "success":
                if report["passed"]:
                    raise ValueError("failed required trial disguised as convergence")
                continue
            j = row["array_index"]
            for metric, study in [
                ("metric_fisher", "study_fisher"),
                ("metric_pair_fisher", "study_pair_fisher"),
                ("metric_volume", "study_volume"),
            ]:
                _close_blocks(
                    arrays[metric][i, side],
                    arrays[study][j],
                    f"metric/trial operand {name} {side} {metric}",
                    ndim=0 if metric == "metric_volume" else 2,
                )
    if failed and (report["passed"] or not report.get("unresolved_controls")):
        raise ValueError("suppressed failed trial outcomes")
    # Replay the controller on saved summaries only. Every lookup must name an
    # actual recorded trial; replay cannot run a model or invent a missing result.
    from types import SimpleNamespace

    from .study import study

    class SavedTrials:
        selection = SimpleNamespace(
            fields=[SimpleNamespace(kind="forest")], selected_pairs=np.array([[0, 0]])
        )
        _prepared = {}
        weight_method = method

        def evaluate(self, task, controls):
            identity = trial_id(controls)
            if identity not in by_id:
                raise KeyError("missing attempted refinement trial")
            row = by_id[identity]
            if row["outcome"] == "failed":
                raise ValueError(row["error"])
            j = row["array_index"]
            return dict(
                fisher=arrays["study_fisher"][j],
                pair_fisher=arrays["study_pair_fisher"][j],
            ), dict(
                settings=dict(
                    grid=dict(volume=arrays["study_volume"][j]),
                    forest_weighting=row.get("forest_weighting", {}),
                )
            )

    saved = SavedTrials()
    context = report.get("context", {})
    if context.get("fields") and not any(
        f["kind"] == "forest" for f in context["fields"]
    ):
        saved.selection = SimpleNamespace(
            fields=[SimpleNamespace(kind="galaxy")], selected_pairs=np.array([[0, 0]])
        )
    try:
        replay_arrays, replay_report = study(saved, {})
    except (KeyError, ValueError, FloatingPointError) as error:
        raise ValueError(f"trial schedule replay failed: {error}") from error
    for name in (
        "trial_contract",
        "study_controls",
        "actual_levels",
        "final_controls",
        "combined_lower_controls",
        "unresolved_controls",
        "metric_names",
    ):
        if report.get(name) != replay_report[name]:
            raise ValueError(f"trial schedule replay differs: {name}")
    for name in ("metric_fisher", "metric_pair_fisher", "metric_volume"):
        _close_blocks(
            arrays[name],
            replay_arrays[name],
            f"replayed trial operands {name}",
            ndim=0 if name == "metric_volume" else 2,
        )
    # Replay constructs these metrics from the referenced study summaries,
    # independently of the detached metric operands supplied by the report.
    for supplied, reconstructed in zip(report["metrics"], replay_report["metrics"]):
        if set(supplied) != set(LIMITS) or any(
            not np.isfinite(supplied[key])
            or not np.isclose(supplied[key], reconstructed[key], rtol=5e-12, atol=1e-15)
            for key in LIMITS
        ):
            raise ValueError("replayed trial metrics differ")
    if report["passed"] is not replay_report["passed"]:
        raise ValueError("replayed trial convergence verdict differs")
    return True
