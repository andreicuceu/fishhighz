"""Stage 4: bounded per-bin BAO comparison from immutable saved inputs."""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from fishhighz.validation.compatibility_bao import dependencies, forecast
from fishhighz.validation.compatibility_weights import (
    VARIANTS,
    coefficients,
    fixed_weights,
    prepare_contexts,
)
from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import contract, relative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", type=int, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    root = Path(".validation")
    stem = (
        root
        / "step12-r5-20260914T191855Z/profiles-checked"
        / f"records-{2 * args.bin:03d}"
    )
    report_path = stem.with_suffix(".report.json")
    array_path = stem.with_suffix(".npz")
    report = json.loads(report_path.read_text())
    with np.load(array_path) as data:
        arrays = {key: data[key] for key in data.files}
    task, rows = report["context"], report["settings"]["pair_inputs"]
    assert task["profile"] == "compatibility" and task["bin"] == args.bin
    stage2_path = root / "compatibility-weighting-stage2" / f"bin-{args.bin}.json"
    stage3_path = root / "compatibility-weighting-stage3" / f"bin-{args.bin}.json"
    stage2, stage3 = (
        json.loads(stage2_path.read_text()),
        json.loads(stage3_path.read_text()),
    )
    assert stage2["report_sha256"] == digest(report_path) == stage3["source_sha256"]
    assert stage3["stage2_json_sha256"] == digest(stage2_path)
    assert stage3["stage2_npz_sha256"] == digest(stage2_path.with_suffix(".npz"))
    trajectories = np.load(stage2_path.with_suffix(".npz"))
    contexts = prepare_contexts(rows)
    records, preparations, checks = [], [], {}
    for variant in VARIANTS:
        for mode, tol in [
            ("three", None),
            ("converged_1e-3", 1e-3),
            ("converged_1e-4", 1e-4),
        ]:
            autos, states = {}, {}
            for name, inputs in contexts.items():
                auto = name in dependencies(task, range(len(task["fields"])))
                state = dict(context=name, auto=auto, variant=variant, mode=mode)
                if tol is None or not auto:
                    w = fixed_weights(inputs, variant, 3)
                    c = coefficients(inputs, w)
                    state.update(
                        status="exactly_three"
                        if tol is None
                        else "nonauto_legacy_three_not_noise",
                        updates=3,
                        coefficients=c.tolist(),
                    )
                else:
                    saved = next(
                        r
                        for r in stage2["records"]
                        if r["nodes"] == 107
                        and r["variant"] == variant
                        and r["context"] == name
                    )
                    evidence = saved["results"][str(tol)]
                    state.update(
                        status=evidence["status"], updates=None, coefficients=None
                    )
                    if evidence["status"] == "finite_nonzero_convergence":
                        t = evidence["confirmed_at"]
                        if variant.startswith("sum_"):
                            adaptive = next(
                                r
                                for r in stage3["records"]
                                if r["nodes"] == 107
                                and r["variant"] == variant
                                and r["context"] == name
                                and r["rtol"] == tol
                            )
                            assert (
                                adaptive["status"] == "converged"
                                and adaptive["updates"] == t
                            )
                        c = trajectories[f"{name}/107/{variant}/coefficients"][t]
                        state.update(
                            status="stage2_tail_tested"
                            if variant.startswith("prefix_")
                            else "adaptive_tail_tested",
                            updates=t,
                            coefficients=c.tolist(),
                            candidate=evidence["candidate"],
                        )
                    else:
                        state["reason"] = (
                            str(saved["failure"])
                            if saved["failure"]
                            else evidence["status"]
                        )
                states[name] = state
                preparations.append(state)
                if auto and state["coefficients"] is not None:
                    autos[name] = state["coefficients"]
            singles, joint, total, covariance, valid = forecast(
                arrays, task, rows, autos
            )
            for index, result in enumerate(singles + [joint]):
                pair = (
                    task["selected_pairs"][index]
                    if index < len(singles)
                    else range(len(task["fields"]))
                )
                needed = dependencies(task, pair)
                missing = [n for n in needed if n not in autos]
                errors = (
                    None
                    if result is None
                    else [
                        float(x) if ok else None
                        for x, ok in zip(
                            result["errors"], result["constrained"], strict=True
                        )
                    ]
                )
                reason = (
                    ("missing converged auto: " + ", ".join(missing))
                    if missing
                    else (
                        "singular Fisher"
                        if errors is not None and None in errors
                        else None
                    )
                )
                omitted = reason or ("sigma > 0.2" if max(errors) > 0.2 else None)
                records.append(
                    dict(
                        bin=args.bin,
                        z=report["settings"]["mean_z"],
                        variant=variant,
                        mode=mode,
                        spectrum="joint"
                        if index == len(singles)
                        else "_".join(task["fields"][i]["id"] for i in pair),
                        dependencies={n: states[n] for n in needed},
                        errors=errors,
                        reason=reason,
                        plot_omission=omitted,
                        fisher=None if result is None else result["fisher"].tolist(),
                    )
                )
                if result is not None and not needed:
                    np.testing.assert_array_equal(
                        result["fisher"], arrays["pair_fisher"][index]
                    )
                    np.testing.assert_array_equal(
                        result["errors"], arrays["pair_errors"][index]
                    )
            if mode == "three" and variant == "prefix_intrinsic":
                checks["baseline"] = dict(
                    total=relative(total, arrays["total"]),
                    covariance=relative(covariance, arrays["selected_covariance"]),
                    fisher=relative(joint["fisher"], arrays["fisher"]),
                    errors=relative(joint["errors"], arrays["errors"]),
                    pair_errors=relative(
                        [s["errors"] for s in singles], arrays["pair_errors"]
                    ),
                )
                assert max(checks["baseline"].values()) < 1e-10
            if mode == "three" and variant == "sum_aliasing":
                # Direct scalar Wick construction and per-cell solves, independent
                # of production factor kernels and vectorized Wick indexing.
                cells = [17, 101, 777]
                direct_f = np.zeros((2, 2))
                lookup = {tuple(p): i for i, p in enumerate(task["required_pairs"])}

                def power(cell, a, b):
                    return total[cell, lookup[tuple(sorted((a, b)))]]

                direct_c = []
                for cell in cells:
                    c = np.array(
                        [
                            [
                                (
                                    power(cell, a, c) * power(cell, b, d)
                                    + power(cell, a, d) * power(cell, b, c)
                                )
                                / arrays["modes"][cell]
                                for c, d in task["selected_pairs"]
                            ]
                            for a, b in task["selected_pairs"]
                        ]
                    )
                    j = arrays["observed_j"][cell]
                    direct_f += j.T @ np.linalg.solve(c, j)
                    direct_c.append(c)
                checks["independent_covariance"] = relative(direct_c, covariance[cells])
                checks["independent_fisher"] = relative(
                    direct_f,
                    contract(covariance[cells], arrays["observed_j"][cells])[0],
                )
                assert (
                    checks["independent_covariance"] < 1e-12
                    and checks["independent_fisher"] < 1e-10
                )
    tolerance_changes = []
    for tight in [
        r for r in records if r["mode"] == "converged_1e-4" and r["reason"] is None
    ]:
        loose = next(
            r
            for r in records
            if r["mode"] == "converged_1e-3"
            and r["variant"] == tight["variant"]
            and r["spectrum"] == tight["spectrum"]
        )
        assert loose["reason"] is None
        delta = (np.array(tight["errors"]) / loose["errors"] - 1).tolist()
        tolerance_changes.append(
            dict(
                variant=tight["variant"],
                spectrum=tight["spectrum"],
                relative=delta,
                hidden=tight["plot_omission"] is not None
                or loose["plot_omission"] is not None,
            )
        )
    checks["tolerance_max"] = max(
        abs(x) for r in tolerance_changes for x in r["relative"]
    )
    assert checks["tolerance_max"] < 0.001
    checks["galaxy_only_bitwise_unchanged"] = True
    output = dict(
        bin=args.bin,
        seconds=time.monotonic() - start,
        checks=checks,
        records=records,
        preparations=preparations,
        tolerance_changes=tolerance_changes,
        sources={
            str(p): digest(p)
            for p in [
                report_path,
                array_path,
                stage2_path,
                stage2_path.with_suffix(".npz"),
                stage3_path,
            ]
        },
    )
    out = root / "compatibility-weighting-stage4"
    out.mkdir(exist_ok=True)
    (out / f"bin-{args.bin}.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(dict(bin=args.bin, seconds=output["seconds"], checks=checks)))


if __name__ == "__main__":
    main()
