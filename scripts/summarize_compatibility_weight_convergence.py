"""Tables and plots of bounded compatibility trajectories (no new updates)."""

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fishhighz.validation.compatibility_weights import VARIANTS
from fishhighz.validation.weight_convergence import changes, classify, relative_change


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--bin", type=int)
    args = parser.parse_args()
    if args.bin is not None:
        path = args.directory / f"bin-{args.bin}.json"
        report = json.loads(path.read_text())
        with np.load(path.with_suffix(".npz")) as arrays:
            names = list(dict.fromkeys(r["context"] for r in report["records"]))
            fw, aw = plt.subplots(9, 5, figsize=(16, 22))
            fc, ac = plt.subplots(9, 5, figsize=(16, 22))
            for r in report["records"]:
                prefix = f"{r['context']}/{r['nodes']}/{r['variant']}"
                data = {
                    k: arrays[prefix + "/" + k]
                    for k in (
                        "weights",
                        "coefficients",
                        "changes",
                        "residual",
                        "amplitude",
                    )
                }
                data["failure"] = r["failure"]
                r["results"] = {str(tol): classify(data, tol) for tol in (1e-3, 1e-4)}
                r["checkpoints"] = {
                    str(t): dict(
                        amplitude=float(data["amplitude"][t]),
                        coefficients=data["coefficients"][t].tolist(),
                    )
                    for t in (0, 3, 6, 12, 24, 48, 96)
                    if t < len(data["weights"])
                }
                if r["nodes"] != 107:
                    continue
                i, j = names.index(r["context"]), VARIANTS.index(r["variant"])
                for t in (0, 3, 6, 12, 24, 48, 96):
                    if t < len(data["weights"]):
                        aw[i, j].plot(
                            arrays[f"{r['context']}/107/magnitudes"],
                            data["weights"][t],
                            label=str(t),
                            lw=0.8,
                        )
                aw[i, j].set_yscale("symlog", linthresh=1e-6)
                d = data["changes"]
                for k, label in enumerate(("amplitude", "shape", "A", "Ppixel")):
                    ac[i, j].semilogy(
                        np.arange(1, len(d) + 1), d[:, k], label=label, lw=0.8
                    )
                ac[i, j].semilogy(
                    np.arange(len(data["residual"])),
                    data["residual"],
                    label="residual",
                    ls=":",
                    lw=0.8,
                )
                ac[i, j].axhline(1e-3, color="grey", lw=0.5)
                ac[i, j].axhline(1e-4, color="grey", lw=0.5, ls="--")
                for ax in (aw[i, j], ac[i, j]):
                    ax.set_title(r["context"] + "\n" + r["variant"], fontsize=8)
                    ax.tick_params(labelsize=6)
            aw[0, 1].legend(fontsize=6, ncol=2)
            fw.subplots_adjust(
                hspace=0.55, wspace=0.25, bottom=0.04, top=0.97, left=0.06, right=0.98
            )
            fc.subplots_adjust(
                hspace=0.55, wspace=0.25, bottom=0.04, top=0.97, left=0.06, right=0.98
            )
            ac[0, 0].legend(fontsize=6)
            fw.supxlabel("Magnitude (signed weights; symlog scale)")
            fc.supxlabel(
                "Update count (relative changes and forward fixed-point residual)"
            )
            fw.savefig(args.directory / f"weights-bin-{args.bin}.png", dpi=100)
            fc.savefig(args.directory / f"convergence-bin-{args.bin}.png", dpi=100)
            plt.close("all")
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print("plotted/reclassified bin", args.bin)
        return
    reports = [
        json.loads((args.directory / f"bin-{i}.json").read_text()) for i in range(6)
    ]
    rows = sum([r["records"] for r in reports], [])
    summary = []
    comparisons = []
    archives = {i: np.load(args.directory / f"bin-{i}.npz") for i in range(6)}
    for nodes in (107, 213, 425):
        for variant in VARIANTS:
            for tol in ("0.001", "0.0001"):
                for auto in (False, True):
                    selected = [
                        r
                        for r in rows
                        if r["nodes"] == nodes
                        and r["variant"] == variant
                        and (r["auto"] or not auto)
                    ]
                    counts = [
                        r["results"][tol]["confirmed_at"]
                        for r in selected
                        if r["results"][tol]["confirmed_at"] is not None
                    ]
                    candidates = [
                        r["results"][tol]["candidate"]
                        for r in selected
                        if r["results"][tol]["candidate"] is not None
                    ]
                    summary.append(
                        dict(
                            nodes=nodes,
                            variant=variant,
                            tolerance=tol,
                            scope="autos12" if auto else "all54",
                            statuses=dict(
                                Counter(r["results"][tol]["status"] for r in selected)
                            ),
                            candidate_range=[min(candidates), max(candidates)]
                            if candidates
                            else None,
                            confirmed_range=[min(counts), max(counts)]
                            if counts
                            else None,
                        )
                    )
    for r in rows:
        if r["nodes"] != 107:
            continue
        a = archives[r["bin"]]
        key = f"{r['context']}/107/{r['variant']}"
        w, c = a[key + "/weights"], a[key + "/coefficients"]
        endpoints = [r["results"][tol]["confirmed_at"] for tol in ("0.001", "0.0001")]
        tolerance_changes = None
        if all(t is not None for t in endpoints):
            t, u = endpoints
            tolerance_changes = changes(w[u], w[t], c[u], c[t]).tolist()
        item = dict(
            bin=r["bin"],
            context=r["context"],
            variant=r["variant"],
            auto=r["auto"],
            tolerance_changes=tolerance_changes,
            grid=[],
        )
        for low, high in ((107, 213), (213, 425)):
            lr = next(
                x
                for x in rows
                if x["bin"] == r["bin"]
                and x["context"] == r["context"]
                and x["variant"] == r["variant"]
                and x["nodes"] == low
            )
            hr = next(
                x
                for x in rows
                if x["bin"] == r["bin"]
                and x["context"] == r["context"]
                and x["variant"] == r["variant"]
                and x["nodes"] == high
            )
            lc = a[f"{r['context']}/{low}/{r['variant']}/coefficients"]
            hc = a[f"{r['context']}/{high}/{r['variant']}/coefficients"]
            t, u = (
                lr["results"]["0.0001"]["confirmed_at"],
                hr["results"]["0.0001"]["confirmed_at"],
            )
            item["grid"].append(
                dict(
                    low=low,
                    high=high,
                    three_update=[
                        relative_change(hc[3, k : k + 1], lc[3, k : k + 1])
                        for k in (3, 4)
                    ],
                    confirmed=[
                        relative_change(hc[u, k : k + 1], lc[t, k : k + 1])
                        for k in (3, 4)
                    ]
                    if t is not None and u is not None
                    else None,
                )
            )
        comparisons.append(item)
    result = dict(
        summary=summary,
        comparisons=comparisons,
        reader_max={
            k: max(
                r["reader_checks"][name][k]
                for r in reports
                for name in r["reader_checks"]
            )
            for k in ("density", "variance")
        },
    )
    (args.directory / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(result["reader_max"]))


if __name__ == "__main__":
    main()
