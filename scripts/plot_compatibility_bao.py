"""Stage-4 matched panels and full-covariance joint curves from saved tables."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fishhighz.validation.compatibility_bao import plot_errors, plot_ratio
from fishhighz.validation.compatibility_weights import VARIANTS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-only", action="store_true")
    args = parser.parse_args()
    root = Path(".validation/compatibility-weighting-stage4")
    data = [json.loads((root / f"bin-{i}.json").read_text()) for i in range(6)]
    rows = [r for d in data for r in d["records"]]
    spectra = list(dict.fromkeys(r["spectrum"] for r in rows))
    z = [d["records"][0]["z"] for d in data]

    def values(variant, mode, spectrum):
        return np.array(
            [
                next(
                    r["errors"] or [np.nan, np.nan]
                    for r in d["records"]
                    if r["variant"] == variant
                    and r["mode"] == mode
                    and r["spectrum"] == spectrum
                )
                for d in data
            ],
            dtype=float,
        )

    caption = "Both components omitted if either sigma > 0.2, unavailable or nonfinite; ratios require both operands. Gaps are retained."
    if not args.tables_only:
        for mode in ["three", "converged_1e-4"]:
            for component, label in enumerate(["parallel", "perp"]):
                for kind in ["errors", "baseline_fraction"]:
                    fig, axes = plt.subplots(3, 5, figsize=(19, 10), sharex=True)
                    for ax, spectrum in zip(axes.flat, spectra[:-1], strict=True):
                        for variant in VARIANTS:
                            a = values(variant, mode, spectrum)
                            b = values("prefix_intrinsic", "three", spectrum)
                            v = plot_errors(a) if kind == "errors" else plot_ratio(a, b)
                            ax.plot(z, v[:, component], ".-", label=variant)
                        ax.set_title(spectrum, fontsize=10)
                        ax.set_xlabel("z")
                        ax.grid(alpha=0.2)
                        ax.set_ylabel(
                            f"sigma(alpha_{label})"
                            if kind == "errors"
                            else "sigma / baseline - 1"
                        )
                    handles, labels = axes.flat[0].get_legend_handles_labels()
                    fig.legend(handles, labels, loc="upper center", ncol=5)
                    fig.suptitle(f"{mode}: {label}, {kind}", y=0.96)
                    fig.text(0.5, 0.01, caption, ha="center", fontsize=9)
                    fig.tight_layout(rect=(0, 0.035, 1, 0.94))
                    fig.savefig(root / f"{mode}-{label}-{kind}.png", dpi=130)
                    plt.close(fig)
        for component, label in enumerate(["parallel", "perp"]):
            fig, axes = plt.subplots(3, 5, figsize=(19, 10), sharex=True)
            for ax, spectrum in zip(axes.flat, spectra[:-1], strict=True):
                for variant in VARIANTS:
                    v = plot_ratio(
                        values(variant, "converged_1e-4", spectrum),
                        values(variant, "three", spectrum),
                    )
                    ax.plot(z, v[:, component], ".-", label=variant)
                ax.set_title(spectrum, fontsize=10)
                ax.set_xlabel("z")
                ax.grid(alpha=0.2)
                ax.set_ylabel("sigma(converged) / sigma(3) - 1")
            fig.legend(
                *axes.flat[0].get_legend_handles_labels(), loc="upper center", ncol=5
            )
            fig.suptitle(f"Within-variant iteration change: {label}", y=0.96)
            fig.text(0.5, 0.01, caption, ha="center", fontsize=9)
            fig.tight_layout(rect=(0, 0.035, 1, 0.94))
            fig.savefig(root / f"convergence-{label}-fraction.png", dpi=130)
            plt.close(fig)
        fig, axes = plt.subplots(2, 5, figsize=(22, 8), sharex=True)
        for component, label in enumerate(["parallel", "perp"]):
            for variant in VARIANTS:
                three = values(variant, "three", "joint")
                conv = values(variant, "converged_1e-4", "joint")
                baseline = values("prefix_intrinsic", "three", "joint")
                series = [
                    plot_errors(three),
                    plot_errors(conv),
                    plot_ratio(three, baseline),
                    plot_ratio(conv, baseline),
                    plot_ratio(conv, three),
                ]
                for col, v in enumerate(series):
                    axes[component, col].plot(z, v[:, component], ".-", label=variant)
            for col, title in enumerate(
                [
                    "three",
                    "converged",
                    "three / baseline - 1",
                    "converged / baseline - 1",
                    "converged / three - 1",
                ]
            ):
                axes[component, col].set_title(title)
                axes[component, col].set_xlabel("z")
                axes[component, col].set_ylabel(label)
                axes[component, col].grid(alpha=0.2)
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper center", ncol=5)
        fig.suptitle("Joint 15×2pt: full inter-spectrum covariance", y=0.94)
        fig.text(0.5, 0.01, caption, ha="center", fontsize=9)
        fig.tight_layout(rect=(0, 0.04, 1, 0.92))
        fig.savefig(root / "joint.png", dpi=130)
        plt.close(fig)
    ratios = []
    for row in rows:
        if row["mode"] == "converged_1e-3":
            continue
        for comparison in ["baseline", "within_variant"]:
            if comparison == "within_variant" and row["mode"] == "three":
                continue
            denominator = next(
                r
                for r in rows
                if r["bin"] == row["bin"]
                and r["spectrum"] == row["spectrum"]
                and r["mode"] == "three"
                and r["variant"]
                == ("prefix_intrinsic" if comparison == "baseline" else row["variant"])
            )
            change = None
            if row["reason"] is None and denominator["reason"] is None:
                change = (np.array(row["errors"]) / denominator["errors"] - 1).tolist()
            ratios.append(
                dict(
                    bin=row["bin"],
                    spectrum=row["spectrum"],
                    variant=row["variant"],
                    mode=row["mode"],
                    comparison=comparison,
                    fractional_change=change,
                    numerator_omission=row["plot_omission"],
                    denominator_omission=denominator["plot_omission"],
                )
            )
    (root / "fractional_changes.json").write_text(
        json.dumps(ratios, indent=2, allow_nan=False) + "\n"
    )
    summary = {}
    for variant in VARIANTS:
        summary[variant] = {}
        for mode in ["three", "converged_1e-4"]:
            rr = [r for r in rows if r["variant"] == variant and r["mode"] == mode]
            changes = []
            within = []
            for r in rr:
                if r["errors"] is None:
                    continue
                base = values("prefix_intrinsic", "three", r["spectrum"])[r["bin"]]
                changes.append(np.array(r["errors"]) / base - 1)
                if mode != "three":
                    within.append(
                        np.array(r["errors"])
                        / values(variant, "three", r["spectrum"])[r["bin"]]
                        - 1
                    )
            summary[variant][mode] = dict(
                individual_available=sum(
                    r["reason"] is None and r["spectrum"] != "joint" for r in rr
                ),
                joint_available=sum(
                    r["reason"] is None and r["spectrum"] == "joint" for r in rr
                ),
                omitted=sum(r["plot_omission"] is not None for r in rr),
                fractional_baseline_range=[
                    np.min(changes, axis=0).tolist(),
                    np.max(changes, axis=0).tolist(),
                ],
                within_range=None
                if not within
                else [np.min(within, axis=0).tolist(), np.max(within, axis=0).tolist()],
            )
    output = dict(
        variants=summary,
        tolerance_max=max(d["checks"]["tolerance_max"] for d in data),
        forecast_seconds=sum(d["seconds"] for d in data),
    )
    (root / "summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output))


if __name__ == "__main__":
    main()
