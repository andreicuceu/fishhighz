"""Offline plots of fixed-weight quadrature and coupled finite-update diagnostics."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fishhighz.validation.evidence import digest


def plot(source, output):
    root = Path(source)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    d = json.loads((root / "diagnosis.json").read_text())
    assert [b["bin"] for b in d["bins"]] == list(range(6))
    files = {}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    for b, ax in zip(d["bins"], axes.flat):
        for count in (0, 3, 6, 12, 24):
            rows = [r for r in b["comparisons"] if r["iterations"] == count]
            if rows:
                ax.semilogy(
                    [r["orders"][1] for r in rows],
                    [max(r["pair_error_relative"], 1e-16) for r in rows],
                    marker="o",
                    ms=3,
                    label="fixed initial w(m)" if count == 0 else f"{count} updates",
                )
        ax.axhline(1e-3, c="black", ls="--", lw=1, label="required tolerance")
        unavailable = []
        for count in (3, 6, 12, 24):
            orders = [
                r["order"]
                for r in b["rows"]
                if r["iterations"] == count and not r["available"]
            ]
            if orders:
                label = (
                    "all orders"
                    if len(orders) == 5
                    else "order " + ", ".join(map(str, orders))
                )
                unavailable.append(f"{count} updates: {label}")
        ax.text(
            0.04,
            0.25,
            "Unavailable:\n" + "\n".join(unavailable),
            transform=ax.transAxes,
            fontsize=8,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"),
        )
        ax.set_title(f"15x2pt bin {b['bin']}")
        ax.set_xlabel("magnitude quadrature order")
        ax.set_ylabel("max relative single-spectrum\nerror change")
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Fixed initial weights versus cumulative updates | consecutive available orders\nGalaxy noise, Fourier nodes, volume and Jacobians held fixed; zero changes displayed at 1e-16"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for suffix in ("png", "pdf"):
        p = out / f"magnitude-convergence.{suffix}"
        fig.savefig(p, dpi=150)
        files[p.name] = digest(p)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
    for b, ax in zip(d["bins"], axes.flat):
        assert digest(root / b["arrays"]) == b["sha256"]
        with np.load(root / b["arrays"], allow_pickle=False) as arrays:
            for fi, ls in ((0, "-"), (4, "--")):
                prefix = f"order16_field{fi}_"
                if prefix + "magnitudes" not in arrays:
                    continue
                for count in (3, 6, 12, 24):
                    key = prefix + f"weights_{count}"
                    if key in arrays:
                        w = arrays[key]
                        ax.semilogy(
                            arrays[prefix + "magnitudes"],
                            np.where(w > 0, w, np.nan),
                            ls=ls,
                            color=f"C{(3, 6, 12, 24).index(count)}",
                            label=f"{'F_Q' if fi == 0 else 'F_L'}, {count} updates",
                        )
        ax.set_title(f"15x2pt bin {b['bin']}, order 16")
        ax.set_xlabel("source magnitude")
        ax.set_ylabel("unnormalized cumulative weight")
        ax.set_ylim(1e-285, 1)
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle(
        "Finite cumulative updates | F_Q solid, F_L dashed\nRecorded weights include trials whose coefficient integrals fail the unchanged arithmetic guards"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for suffix in ("png", "pdf"):
        p = out / f"weight-decay.{suffix}"
        fig.savefig(p, dpi=150)
        files[p.name] = digest(p)
    plt.close(fig)
    (out / "manifest.json").write_text(
        json.dumps(
            dict(source_sha256=digest(root / "diagnosis.json"), files=files), indent=2
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    plot(a.source, a.output)
