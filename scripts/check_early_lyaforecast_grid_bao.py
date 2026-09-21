"""Bounded historical-sum magnitude-grid BAO diagnostic, user bins 2--6."""

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from fishhighz.validation.compatibility_bao import dependencies, forecast, plot_ratio
from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import contract, relative

ROOT = Path(".validation")
OUT = ROOT / "early-lyaforecast-grid-bao"
GRIDS = (107, 213, 425)


def run(bin_index):
    start = time.monotonic()
    stem = (
        ROOT
        / "step12-r5-20260914T191855Z/profiles-checked"
        / f"records-{2 * bin_index:03d}"
    )
    report_path = stem.with_suffix(".report.json")
    report = json.loads(report_path.read_text())
    with np.load(stem.with_suffix(".npz")) as source:
        arrays = {k: source[k] for k in source.files}
    paths = [
        ROOT / f"compatibility-weighting-stage{s}/bin-{bin_index}.json"
        for s in (2, 3, 4)
    ]
    stage2, stage3, stage4 = [json.loads(p.read_text()) for p in paths]
    assert all(s["bin"] == bin_index for s in (stage2, stage3, stage4))
    assert stage2["report_sha256"] == stage3["source_sha256"] == digest(report_path)
    assert stage3["stage2_json_sha256"] == digest(paths[0])
    assert stage3["stage2_npz_sha256"] == digest(paths[0].with_suffix(".npz"))
    task, rows = report["context"], report["settings"]["pair_inputs"]
    assert task["bin"] == bin_index and task["profile"] == "compatibility"
    records, preparations, checks = [], [], {}
    with np.load(paths[0].with_suffix(".npz")) as trajectories:
        for nodes in GRIDS:
            autos = {}
            for name in dependencies(task, range(len(task["fields"]))):
                saved = next(
                    r
                    for r in stage2["records"]
                    if r["nodes"] == nodes
                    and r["variant"] == "sum_historical"
                    and r["context"] == name
                )
                state = saved["results"]["0.0001"]
                stopped = next(
                    r
                    for r in stage3["records"]
                    if r["nodes"] == nodes
                    and r["variant"] == "sum_historical"
                    and r["context"] == name
                    and r["rtol"] == 1e-4
                )
                assert state["status"] == "finite_nonzero_convergence"
                assert stopped["status"] == "converged"
                count = state["confirmed_at"]
                assert (
                    stopped["updates"] == count
                    and stopped["candidate"] == state["candidate"]
                )
                assert max(stopped["tail_max"].values()) < 1e-4
                coeff = trajectories[f"{name}/{nodes}/sum_historical/coefficients"][
                    count
                ]
                np.testing.assert_array_equal(coeff, stopped["coefficients"])
                autos[name] = coeff
                preparations.append(
                    dict(
                        nodes=nodes,
                        context=name,
                        updates=count,
                        candidate=state["candidate"],
                        coefficients=coeff.tolist(),
                        tail_max=stopped["tail_max"],
                    )
                )
            singles, joint, total, covariance, valid = forecast(
                arrays, task, rows, autos
            )
            assert np.all(valid) and joint is not None
            for index, result in enumerate(singles + [joint]):
                pair = (
                    task["selected_pairs"][index]
                    if index < len(singles)
                    else range(len(task["fields"]))
                )
                name = (
                    "_".join(task["fields"][i]["id"] for i in pair)
                    if index < len(singles)
                    else "joint"
                )
                assert np.all(result["constrained"]) and np.all(
                    np.isfinite(result["errors"])
                )
                records.append(
                    dict(
                        bin=bin_index,
                        user_bin=bin_index + 1,
                        z=report["settings"]["mean_z"],
                        nodes=nodes,
                        spectrum=name,
                        errors=result["errors"].tolist(),
                        fisher=result["fisher"].tolist(),
                        plot_omission="sigma > 0.2"
                        if max(result["errors"]) > 0.2
                        else None,
                    )
                )
                if not dependencies(task, pair):
                    np.testing.assert_array_equal(
                        result["fisher"], arrays["pair_fisher"][index]
                    )
                    np.testing.assert_array_equal(
                        result["errors"], arrays["pair_errors"][index]
                    )
                if nodes == 107:
                    original = next(
                        r
                        for r in stage4["records"]
                        if r["variant"] == "sum_historical"
                        and r["mode"] == "converged_1e-4"
                        and r["spectrum"] == name
                    )
                    np.testing.assert_array_equal(result["errors"], original["errors"])
                    np.testing.assert_array_equal(result["fisher"], original["fisher"])
            if bin_index == 1 and nodes == 425:
                cells = [17, 101, 777]
                lookup = {tuple(p): i for i, p in enumerate(task["required_pairs"])}
                direct_c, direct_f = [], np.zeros((2, 2))
                for cell in cells:

                    def power(a, b):
                        return total[cell, lookup[tuple(sorted((a, b)))]]

                    c = np.array(
                        [
                            [
                                (power(a, c) * power(b, d) + power(a, d) * power(b, c))
                                / arrays["modes"][cell]
                                for c, d in task["selected_pairs"]
                            ]
                            for a, b in task["selected_pairs"]
                        ]
                    )
                    j = arrays["observed_j"][cell]
                    direct_c.append(c)
                    direct_f += j.T @ np.linalg.solve(c, j)
                checks["independent_covariance"] = relative(direct_c, covariance[cells])
                checks["independent_fisher"] = relative(
                    direct_f,
                    contract(covariance[cells], arrays["observed_j"][cells])[0],
                )
                assert (
                    checks["independent_covariance"] < 1e-12
                    and checks["independent_fisher"] < 1e-10
                )
    checks.update(original_grid_stage4_bitwise=True, galaxy_only_bitwise=True)
    sources = [
        report_path,
        stem.with_suffix(".npz"),
        *paths,
        paths[0].with_suffix(".npz"),
        Path(__file__),
        Path("fishhighz/validation/compatibility_bao.py"),
        Path("fishhighz/validation/compatibility_weights.py"),
        Path("fishhighz/validation/numerics.py"),
    ]
    output = dict(
        bin=bin_index,
        seconds=time.monotonic() - start,
        checks=checks,
        records=records,
        preparations=preparations,
        sources={str(p): digest(p) for p in sources},
    )
    OUT.mkdir(exist_ok=True)
    (OUT / f"bin-{bin_index}.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(dict(bin=bin_index, seconds=output["seconds"], checks=checks)))


def summarize():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = [json.loads((OUT / f"bin-{b}.json").read_text()) for b in range(1, 6)]
    records = [r for d in data for r in d["records"]]
    names = [r["spectrum"] for r in records if r["bin"] == 1 and r["nodes"] == 107]
    comparisons = [(107, 213), (213, 425), (107, 425)]
    changes = []
    for b in range(1, 6):
        for name in names:
            for low, high in comparisons:
                a, c = [
                    next(
                        r
                        for r in records
                        if r["bin"] == b and r["nodes"] == n and r["spectrum"] == name
                    )
                    for n in (low, high)
                ]
                delta = np.array(c["errors"]) / a["errors"] - 1
                for component in range(2):
                    changes.append(
                        dict(
                            bin=b,
                            user_bin=b + 1,
                            z=a["z"],
                            spectrum=name,
                            component=["parallel", "transverse"][component],
                            coarse=low,
                            fine=high,
                            sigma_coarse=a["errors"][component],
                            sigma_fine=c["errors"][component],
                            fraction=float(delta[component]),
                            hidden=bool(a["plot_omission"] or c["plot_omission"]),
                        )
                    )
    maxima = []
    for low, high in comparisons:
        for b in [None, 1, 2, 3, 4, 5]:
            for group in ["individual", "joint"]:
                chosen = [
                    r
                    for r in changes
                    if r["coarse"] == low
                    and r["fine"] == high
                    and (b is None or r["bin"] == b)
                    and ((r["spectrum"] == "joint") == (group == "joint"))
                ]
                maxima.append(
                    dict(
                        bin=b,
                        group=group,
                        coarse=low,
                        fine=high,
                        maximum=max(chosen, key=lambda r: abs(r["fraction"])),
                        below_0p1_percent=all(
                            abs(r["fraction"]) < 0.001 for r in chosen
                        ),
                    )
                )
    (OUT / "summary.json").write_text(
        json.dumps(dict(maxima=maxima, changes=changes), indent=2, allow_nan=False)
        + "\n"
    )
    with (OUT / "fractional_changes.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(changes[0]))
        writer.writeheader()
        writer.writerows(changes)
    for component, label in enumerate(["parallel", "transverse"]):
        fig, axes = plt.subplots(4, 4, figsize=(14, 11), sharex=True)
        for ax, name in zip(axes.flat, names, strict=True):
            for low, high in comparisons:
                operands = [
                    [
                        next(
                            r
                            for r in records
                            if r["bin"] == b
                            and r["nodes"] == n
                            and r["spectrum"] == name
                        )
                        for b in range(1, 6)
                    ]
                    for n in (low, high)
                ]
                delta = plot_ratio(
                    [r["errors"] for r in operands[1]],
                    [r["errors"] for r in operands[0]],
                )[:, component]
                ax.plot(
                    [r["z"] for r in operands[0]],
                    100 * delta,
                    "o-",
                    markersize=3,
                    label=f"{low} → {high}",
                )
            ax.axhline(0, color="black", lw=0.5)
            ax.axhline(0.1, color="gray", lw=0.5, ls=":")
            ax.axhline(-0.1, color="gray", lw=0.5, ls=":")
            ax.set_title(name, fontsize=9)
            ax.set_ylabel("Δσ / σ [%]")
            ax.set_xlabel("Mean redshift")
        axes.flat[0].legend(fontsize=8)
        fig.suptitle(f"Historical sum: {label} BAO uncertainty; user bins 2–6")
        fig.tight_layout()
        fig.savefig(OUT / f"{label}-refinement.png", dpi=160)
        plt.close(fig)
    print(json.dumps([m for m in maxima if m["bin"] is None], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", type=int, choices=range(1, 6))
    parser.add_argument("--summarize", action="store_true")
    args = parser.parse_args()
    if args.summarize:
        summarize()
    elif args.bin is not None:
        run(args.bin)
    else:
        parser.error("choose --bin or --summarize")
