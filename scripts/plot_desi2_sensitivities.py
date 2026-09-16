"""Plot verified sensitivity effects with gaps for unavailable calculations."""

import argparse
import json
import textwrap
from pathlib import Path

import numpy as np

from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import change


def fisher_arrays(root, record):
    path = (root / record["arrays"]).resolve()
    if path.parent != root.resolve():
        raise ValueError("invalid source array path")
    with np.load(path, allow_pickle=False) as arrays:
        return arrays["fisher"], arrays["pair_fisher"]


def plot(plots, output):
    source = Path(plots)
    table = json.loads((source / "tables.json").read_text())
    proof = json.loads((source / "plots-manifest.json").read_text())
    if digest(source / "tables.json") != proof["table_sha256"]:
        raise ValueError("stale checked table")
    root = Path(table["source_bundle"])
    if digest(root / "manifest.json") != table["source_sha256"]:
        raise ValueError("stale source bundle")
    bundle = json.loads((root / "manifest.json").read_text())
    primary = {
        (r["task"]["case"], r["task"]["bin"]): fisher_arrays(root, r)
        for r in bundle["records"]
        if r["task"]["profile"] == "accuracy"
    }
    diagnostic = {
        (r["task"]["case"], r["task"]["bin"], r["task"]["diagnostic_id"]): r
        for r in bundle["diagnostics"]
    }
    checked = 0
    for row in table["diagnostics"]:
        key = (row["case"], row["bin"], row["policy"])
        record = diagnostic[key]
        if row["status"] != record["status"]:
            raise ValueError("wrong diagnostic status")
        if record["status"] != "completed":
            if row.get("sensitivity") is not None:
                raise ValueError("failed diagnostic has a finite plotted result")
            continue
        fisher, pairs = fisher_arrays(root, record)
        base, base_pairs = primary[key[:2]]
        actual = change(fisher, base, pairs, base_pairs, 1, 1)
        if set(actual) != set(row["sensitivity"]["change"]):
            raise ValueError("incorrect sensitivity metrics")
        for name, value in actual.items():
            np.testing.assert_allclose(
                row["sensitivity"]["change"][name], value, rtol=5e-12, atol=1e-15
            )
        checked += 1
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = {}
    for case in dict.fromkeys(row["case"] for row in table["rows"]):
        primary_rows = [r for r in table["rows"] if r["case"] == case]
        indices = [r["bin"] for r in primary_rows]
        rows = [r for r in table["diagnostics"] if r["case"] == case]
        policies = sorted({r["policy"] for r in rows})
        lookup = {(r["bin"], r["policy"]): r for r in rows}
        fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
        unavailable = []
        for policy in policies:
            values = []
            missing = []
            for index in indices:
                row = lookup[index, policy]
                sensitivity = row.get("sensitivity")
                if sensitivity is None:
                    missing.append(index)
                    values.append([np.nan, np.nan])
                else:
                    values.append(
                        [
                            100 * sensitivity["change"]["error_relative"],
                            100 * sensitivity["change"]["fisher_relative"],
                        ]
                    )
            if missing:
                unavailable.append(policy + ": " + ",".join(map(str, missing)))
            label = policy + (" (unavailable)" if len(missing) == len(indices) else "")
            for j, ax in enumerate(axes):
                (line,) = ax.plot(
                    indices, np.asarray(values)[:, j], marker="o", ms=4, label=label
                )
                assert np.array_equal(
                    np.isnan(line.get_ydata()), np.isin(indices, missing)
                )
        for ax, label in zip(
            axes, ["max combined uncertainty change (%)", "Fisher Frobenius change (%)"]
        ):
            ax.set_xlabel("original bin index")
            ax.set_ylabel(label)
            ax.set_xticks(indices)
            ax.grid(alpha=0.2)
            ax.legend(fontsize=7)
        failed = [
            r["bin"] for r in primary_rows if not r["profiles"]["accuracy"]["passed"]
        ]
        fig.suptitle(
            case
            + " | sensitivity at recorded controls\nPrimary accuracy convergence unresolved in bins "
            + str(failed),
            fontsize=11,
        )
        footer = "Unavailable diagnostic bins: " + ("; ".join(unavailable) or "none")
        fig.text(0.5, 0.02, textwrap.fill(footer, 125), ha="center", fontsize=8)
        fig.tight_layout(rect=(0, 0.08, 1, 0.93))
        for suffix in ("png", "pdf"):
            path = out / f"{case}.{suffix}"
            fig.savefig(path, dpi=140)
            files[path.name] = digest(path)
        plt.close(fig)
    result = dict(
        kind="failure-aware-sensitivity-plots",
        source_table=str((source / "tables.json").resolve()),
        table_sha256=digest(source / "tables.json"),
        source_sha256=table["source_sha256"],
        script_sha256=digest(__file__),
        numerical_diagnostics_checked=checked,
        unavailable_diagnostics=len(table["diagnostics"]) - checked,
        files=files,
    )
    (out / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Verified {checked} finite diagnostic effects and explicit gaps for {len(table['diagnostics']) - checked} unavailable variants"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plots", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    plot(args.plots, args.output)
