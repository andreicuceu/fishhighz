"""Offline, failure-preserving numerical tables and scientific comparison plots."""

import csv
import json
from pathlib import Path

import numpy as np

from .attribution import chain
from .cases import CASE_IDS, bins, selection
from .evidence import canonical, digest, modern_requests, record_report
from .numerics import information, relative
from .schema import validate_payload


def records(root):
    """Read exact inventory and validate mathematics, retaining scientific failures."""
    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    expected = modern_requests(
        manifest["suite"],
        manifest["cases"],
        manifest["bin_indices"],
        profiles=manifest["profiles"],
        kind=manifest["payload_kind"],
    )
    if (
        manifest["schema"] not in (2, 3)
        or manifest["requested"] != expected
        or [r["task"] for r in manifest["records"]] != expected
    ):
        raise ValueError("wrong primary plotting inventory")
    if manifest.get("execution_finished") is not True:
        raise ValueError("execution unfinished")
    for group in ("records", "diagnostics"):
        if (
            group == "diagnostics"
            and [r["task"] for r in manifest[group]]
            != manifest["diagnostics_requested"]
        ):
            raise ValueError("diagnostic inventory")
        for row in manifest[group]:
            if "arrays" not in row:
                yield row, None
                continue
            row = {**row, "report": record_report(root, row)}
            path = (root / row["arrays"]).resolve()
            if path.parent != root or digest(path) != row["sha256"]:
                raise ValueError("plot source hash/path")
            if canonical(row["report"]) != row["effective_hash"]:
                raise ValueError("plot report hash")
            with np.load(path, allow_pickle=False) as f:
                arrays = {k: f[k] for k in f.files}
            if {k: list(v.shape) for k, v in arrays.items()} != row["inventory"]:
                raise ValueError("plot array inventory")
            validate_payload(
                row["task"],
                arrays,
                row["report"],
                require_pass=False,
                schema=manifest["schema"],
            )
            yield row, arrays


def values(fisher):
    info = information(fisher)
    return [
        float(info["errors"][i]) if info["constrained"][i] else None for i in range(2)
    ] + [float(info["correlation"][0, 1]) if np.all(info["constrained"]) else None]


def difference(a, b):
    return (
        [
            (x / y - 1) * 100
            if i < 2 and y
            else x - y
            if i == 2 and x is not None and y is not None
            else None
            for i, (x, y) in enumerate(zip(a, b))
        ]
        if all(x is not None for x in a + b)
        else [None] * 3
    )


def tables(bundle, output):
    """Every primary and selected pair survives even if execution/science failed."""
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    rows = {}
    diagnostics = []
    attributions = []
    pending = {}
    for record, a in records(bundle):
        task = record["task"]
        key = (task["case"], task["bin"])
        if task["kind"] == "diagnostic":
            diagnostics.append(
                dict(
                    case=key[0],
                    bin=key[1],
                    policy=task["diagnostic_id"],
                    status=record["status"],
                    passed=record.get("scientific_passed", False),
                    sensitivity=record.get("report", {}).get("sensitivity"),
                    error=record.get("error"),
                )
            )
            continue
        if key not in rows:
            sel = selection(key[0])
            lo, hi = bins(key[0])[key[1]]
            rows[key] = dict(
                case=key[0],
                bin=key[1],
                bounds=[lo, hi],
                centre=(lo + hi) / 2,
                z_eval=float(np.sqrt((1 + lo) * (1 + hi)) - 1),
                profiles={},
                pairs=[
                    dict(
                        pair=[int(i), int(j)],
                        label=f"{sel.fields[i].id} x {sel.fields[j].id}",
                        profiles={},
                    )
                    for i, j in sel.selected_pairs
                ],
            )
        row = rows[key]
        p = task["profile"]
        row["profiles"][p] = dict(
            values=values(a["fisher"]) if a is not None else [None] * 3,
            passed=record.get("scientific_passed", False),
            metrics=record.get("report", {}).get("metrics"),
            controls=record.get("report", {}).get("final_controls"),
            error=record.get("scientific_error", record.get("error")),
        )
        for i, pair in enumerate(row["pairs"]):
            pair["profiles"][p] = (
                values(a["pair_fisher"][i]) if a is not None else [None] * 3
            )
        if a is not None and p == "compatibility":
            row["profiles"]["reference"] = dict(
                values=values(a["reference_fisher"]), passed=True
            )
            for i, pair in enumerate(row["pairs"]):
                pair["profiles"]["reference"] = values(a["reference_pair_fisher"][i])
            pending[key] = (a, record["report"]["settings"]["grid"]["volume"])
        elif a is not None and key in pending:
            legacy, volume = pending.pop(key)
            try:
                arrays, report = chain(
                    legacy,
                    a,
                    task,
                    volume,
                    record["report"]["settings"]["grid"]["volume"],
                )
                name = f"attribution-{key[0]}-{key[1]}.npz"
                np.savez_compressed(out / name, **arrays)
                report.update(
                    array=name,
                    sha256=digest(out / name),
                    values=[values(arrays[f"fisher_{i}"]) for i in range(5)],
                )
                attributions.append(report)
            except (ValueError, np.linalg.LinAlgError) as error:
                attributions.append(
                    dict(context=task, error=str(error), complete=False)
                )
            row["fisher_accuracy_compatibility"] = relative(
                a["fisher"], legacy["fisher"]
            )
        if p == "accuracy":
            for entry in [row, *row["pairs"]]:
                profiles = entry["profiles"]
                v = {
                    name: (item["values"] if entry is row else item)
                    for name, item in profiles.items()
                }
                entry["differences"] = {
                    f"{x}/{y}": difference(v.get(x, [None] * 3), v.get(y, [None] * 3))
                    for x, y in [
                        ("compatibility", "reference"),
                        ("accuracy", "reference"),
                        ("accuracy", "compatibility"),
                    ]
                }
    table = dict(
        schema=2,
        source_bundle=str(Path(bundle).resolve()),
        source_sha256=digest(Path(bundle) / "manifest.json"),
        columns=["sigma_ap", "sigma_at", "correlation"],
        difference_units=["percent", "percent", "absolute"],
        rows=list(rows.values()),
        diagnostics=diagnostics,
        attributions=attributions,
        caption="Maximum accuracy denotes the prescribed existing-model profile; failed convergence is explicitly marked. Raw-policy perturbations are diagnostics, not calibrated systematic errors.",
    )
    (out / "tables.json").write_text(
        json.dumps(table, indent=2, allow_nan=False) + "\n"
    )
    with (out / "tables.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case",
                "bin",
                "pair",
                "profile",
                "sigma_ap",
                "sigma_at",
                "correlation",
                "passed",
                "ap_percent_reference",
                "at_percent_reference",
                "correlation_difference_reference",
            ]
        )
        for row in table["rows"]:
            for entry in [row, *row["pairs"]]:
                for p, v in entry["profiles"].items():
                    value = v["values"] if entry is row else v
                    diff = entry.get("differences", {}).get(
                        p + "/reference", [0, 0, 0] if p == "reference" else [None] * 3
                    )
                    w.writerow(
                        [
                            row["case"],
                            row["bin"],
                            entry.get("label", "combined"),
                            p,
                            *value,
                            row["profiles"].get(p, {}).get("passed", False),
                            *diff,
                        ]
                    )
    return table


def draw_case(plt, case, rows, pair=None):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    styles = dict(
        reference=("black", "o", "actual lyaforecast"),
        compatibility=("tab:blue", "D", "FishHighz maximum compatibility"),
        accuracy=("tab:orange", "x", "FishHighz maximum accuracy"),
    )
    z = [r["centre"] for r in rows]
    width = [(r["bounds"][1] - r["bounds"][0]) / 2 for r in rows]
    entries = [r if pair is None else r["pairs"][pair] for r in rows]
    for p, (color, marker, label) in styles.items():
        val = [e["profiles"].get(p, {"values": [None] * 3}) for e in entries]
        val = [v["values"] if isinstance(v, dict) else v for v in val]
        for j in range(3):
            y = [np.nan if v[j] is None else v[j] for v in val]
            axes[0, j].errorbar(
                z,
                y,
                xerr=width,
                color=color,
                marker=marker,
                mfc="none",
                label=label,
                lw=1,
                ms=5,
                capsize=2,
            )
            if p != "reference":
                d = [
                    e.get("differences", {}).get(p + "/reference", [None] * 3)[j]
                    for e in entries
                ]
                axes[1, j].plot(
                    z,
                    [np.nan if x is None else x for x in d],
                    color=color,
                    marker=marker,
                    mfc="none",
                )
    for j, title in enumerate(
        (r"$\sigma_{a_\parallel}$", r"$\sigma_{a_\perp}$", "correlation")
    ):
        axes[0, j].set_title(title)
        axes[1, j].set_ylabel(
            "% change to legacy" if j < 2 else "absolute change to legacy"
        )
        axes[1, j].axhline(0, color="gray", lw=0.7)
        axes[1, j].set_xlabel("arithmetic bin centre z")
        for ax in axes[:, j]:
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7)
    failed = [
        r["bin"] for r in rows if not r["profiles"].get("accuracy", {}).get("passed")
    ]
    title = case + (" | " + rows[0]["pairs"][pair]["label"] if pair is not None else "")
    fig.suptitle(
        title + "\nAccuracy convergence unresolved in bins " + str(failed)
        if failed
        else title,
        fontsize=11,
    )
    fig.tight_layout()
    return fig


def plot(table, output):
    """Render only checked table values; optional Matplotlib is imported lazily."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    out = Path(output)
    files = []

    def save(fig, name):
        for suffix in ("png", "pdf"):
            path = out / f"{name}.{suffix}"
            fig.savefig(path, dpi=140)
            files.append(path)
        plt.close(fig)

    for case in CASE_IDS:
        rows = [r for r in table["rows"] if r["case"] == case]
        if not rows:
            continue
        save(draw_case(plt, case, rows), case)
        path = out / f"{case}-pairs.pdf"
        with PdfPages(path) as pdf:
            for pair in range(len(rows[0]["pairs"])):
                fig = draw_case(plt, case, rows, pair)
                pdf.savefig(fig)
                plt.close(fig)
        files.append(path)
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    for p, color in [("compatibility", "tab:blue"), ("accuracy", "tab:orange")]:
        for j, ax in enumerate(axes):
            vals = [
                r.get("differences", {}).get(p + "/reference", [None] * 3)[j]
                for r in table["rows"]
            ]
            ax.plot(
                range(len(vals)),
                [np.nan if x is None else x for x in vals],
                marker="o",
                ms=4,
                label=p,
                color=color,
            )
            bad = [
                i
                for i, r in enumerate(table["rows"])
                if not r["profiles"].get(p, {}).get("passed")
            ]
            ax.scatter(
                bad,
                [vals[i] for i in bad],
                marker="x",
                s=60,
                color="red",
                label="unresolved" if p == "accuracy" else None,
            )
            ax.set_ylabel(("ap" if j == 0 else "at") + " change to legacy (%)")
            ax.grid(alpha=0.2)
    axes[0].legend()
    axes[1].set_xticks(
        range(len(table["rows"])),
        [
            r["case"].replace("lya_qso_lbg_lae_", "") + ":" + str(r["bin"])
            for r in table["rows"]
        ],
        rotation=90,
        fontsize=7,
    )
    fig.tight_layout()
    save(fig, "all-cases")
    for case in CASE_IDS:
        diag = [
            d
            for d in table["diagnostics"]
            if d["case"] == case and d.get("sensitivity")
        ]
        attr = [
            a
            for a in table["attributions"]
            if a["context"]["case"] == case and a.get("complete")
        ]
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for policy in sorted({d["policy"] for d in diag}):
            ds = [d for d in diag if d["policy"] == policy]
            axes[0].plot(
                [d["bin"] for d in ds],
                [d["sensitivity"]["change"]["error_relative"] * 100 for d in ds],
                marker=".",
                label=policy,
            )
        axes[0].set_ylabel("max uncertainty change (%)")
        axes[0].set_xlabel("bin")
        axes[0].legend(fontsize=6)
        for a in attr:
            v = np.array(a["values"])
            axes[1].plot(
                range(5),
                100 * (v[:, 0] / v[0, 0] - 1),
                marker=".",
                label="bin " + str(a["context"]["bin"]),
            )
        axes[1].set_xticks(
            range(5), ["legacy", "grid", "volume", "total P", "full J"], rotation=30
        )
        axes[1].set_ylabel("ap change from legacy (%)")
        axes[1].legend(fontsize=7)
        fig.suptitle(
            case
            + " | diagnostic swaps; grouped physical attribution remains unresolved"
        )
        fig.tight_layout()
        save(fig, case + "-diagnostics")
    manifest = dict(
        table_sha256=digest(out / "tables.json"),
        source_sha256=table["source_sha256"],
        files={p.name: digest(p) for p in files},
    )
    (out / "plots-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
