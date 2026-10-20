"""Create scientific plots, tables and validation summaries from saved runs."""

import csv
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path.cwd()
REPORT = ROOT / "report"
REPORT.mkdir(exist_ok=True)
PLOTS = REPORT / "plots"
PLOTS.mkdir(exist_ok=True)
spec = importlib.util.spec_from_file_location(
    "study", ROOT / "fishhighz/examples/eft_3x2_study.py"
)
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def read(name):
    return json.loads((ROOT / "results" / name / "forecasts.json").read_text())


runs = {n: read(n) for n in ["main", "refined", "step_half", "quadrature", "fft128"]}
rows = runs["refined"]
metadata = json.loads((ROOT / "results/refined/metadata.json").read_text())
arrays = np.load(ROOT / "results/refined/arrays.npz")
labels = [
    r"$\alpha_{\rm wiggle}$",
    r"$\varphi_{\rm wiggle}$",
    r"$\alpha_{\rm smooth}$",
    r"$\varphi_{\rm smooth}$",
    r"$f\sigma_8$",
]
colors = {"FF": "#287f8e", "GG": "#c77228", "3x2": "#6255a5"}
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
    }
)


def save(fig, name):
    fig.savefig(PLOTS / (name + ".pdf"))
    fig.savefig(PLOTS / (name + ".png"), dpi=180)
    plt.close(fig)


def selected(data, order, prior, source=rows):
    return [
        r
        for r in source
        if r["data"] == data and r["order"] == order and r["prior"] == prior
    ]


def relative_errors(r):
    a = np.array(r["errors"][:5], float)
    a[4] *= metadata["sigma8"] / metadata["fsigma8"]
    return a * 100


# Target forecast curves, including explicit unverified high-k region.
fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
for i, ax in enumerate(axes.flat[:5]):
    for data in colors:
        for order, prior, style in [
            ("tree", None, ":"),
            ("loop", None, "--"),
            ("loop", 1.0, "-"),
        ]:
            rr = selected(data, order, prior)
            yy = [relative_errors(r)[i] for r in rr]
            ax.plot([r["kmax"] for r in rr], yy, style, color=colors[data], lw=1.7)
    ax.axvspan(0.2, 1.0, color="#eeeeee", zorder=-10)
    ax.set(
        xscale="log",
        yscale="log",
        title=labels[i],
        xlabel=r"$k_{\max}\ [h/{\rm Mpc}]$",
        ylabel=r"Marginalized error [%]",
    )
    ax.grid(alpha=0.15)
axes.flat[4].text(
    0.04,
    0.08,
    "FF tree growth is unidentifiable",
    transform=axes.flat[4].transAxes,
    fontsize=8,
)
axes.flat[5].axis("off")
for j, (data, color) in enumerate(colors.items()):
    axes.flat[5].plot([], [], color=color, label=data)
for style, text in [
    (":", "Tree: linear biases free"),
    ("--", "Loop: all EFT free"),
    ("-", "Loop: unit bias/ct priors"),
]:
    axes.flat[5].plot([], [], style, color="black", label=text)
axes.flat[5].legend(loc="center", frameon=False)
axes.flat[5].text(
    0.04,
    0.08,
    "Shaded: high-k extrapolation;\nphysical validity is unverified.",
    transform=axes.flat[5].transAxes,
    fontsize=10,
)
fig.suptitle("Idealized z=2.3, V=20 (Gpc/h)^3; common fixed loop covariance")
fig.tight_layout()
save(fig, "targets")

# Prior response for joint analysis.
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))
for ax, i in zip(axes, [2, 3, 4]):
    for cut, style in [(0.2, "-"), (0.5, "--"), (1.0, ":")]:
        values = [
            next(
                r
                for r in rows
                if r["data"] == "3x2"
                and r["order"] == "loop"
                and r["kmax"] == cut
                and r["prior"] == p
            )
            for p in [0.3, 1.0, 3.0]
        ]
        no = next(
            r
            for r in rows
            if r["data"] == "3x2"
            and r["order"] == "loop"
            and r["kmax"] == cut
            and r["prior"] is None
        )
        ax.plot(
            [0.3, 1.0, 3.0],
            [r["errors"][i] / no["errors"][i] for r in values],
            style + "o",
            label=f"kmax={cut}",
        )
    ax.set(
        xscale="log",
        ylim=(0, 1.05),
        title=labels[i],
        xlabel=r"Prior width $\sigma_{\rm EFT}$",
        ylabel="Error / no-prior error",
    )
    ax.grid(alpha=0.15)
axes[0].legend(frameon=False)
fig.tight_layout()
save(fig, "priors")

# Higher EFT principal modes after marginalizing cosmology and linear biases.
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5))
pca_counts = []
for ax, data in zip(axes, colors):
    for cut in [0.08, 0.2, 0.5, 1.0]:
        r = next(r for r in selected(data, "loop", None) if r["kmax"] == cut)
        sig = np.array(
            [m["sigma"] if m["sigma"] is not None else np.inf for m in r["pca"]]
        )
        ax.plot(np.arange(1, len(sig) + 1), sig, label=f"{cut}")
        pca_counts.append(
            {
                "data": data,
                "kmax": cut,
                "modes_sigma_lt_1": int(np.sum(sig < 1)),
                "n_modes": len(sig),
            }
        )
    ax.axhline(1, color="black", lw=0.8)
    ax.set(
        yscale="log",
        ylim=(1e-4, 1e6),
        title=data,
        xlabel="Higher EFT mode (best first)",
        ylabel="Dimensionless mode error",
    )
    ax.grid(alpha=0.15)
axes[0].legend(title="kmax", frameon=False)
fig.tight_layout()
save(fig, "pca")

# Model correction ratios.
k = arrays["k"]
mu = arrays["mu"]
uniq_k = np.unique(k)
uniq_mu = np.unique(mu)
tree = arrays["tree_power"].reshape(len(uniq_k), len(uniq_mu), 3)
loop = arrays["loop_power"].reshape(len(uniq_k), len(uniq_mu), 3)
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.3))
loop_ratios = {}
for p, ax in enumerate(axes):
    for muval in [0.1, 0.6, 0.95]:
        im = np.argmin(abs(uniq_mu - muval))
        ratio = np.abs((loop[:, im, p] - tree[:, im, p]) / tree[:, im, p])
        ax.plot(uniq_k, ratio, label=f"mu={uniq_mu[im]:.2f}")
    ax.axhline(0.3, color="black", lw=0.8, ls="--")
    ax.set(
        xscale="log",
        yscale="log",
        title=["FF", "FG", "GG"][p],
        xlabel=r"$k\ [h/{\rm Mpc}]$",
        ylabel=r"$|P_{\rm loop\ model}-P_{\rm tree}|/|P_{\rm tree}|$",
    )
    ax.grid(alpha=0.15)
    loop_ratios[["FF", "FG", "GG"][p]] = {
        str(c): float(np.max(abs((loop - tree) / tree)[uniq_k <= c, :, p]))
        for c in [0.2, 0.5, 1.0]
    }
axes[0].legend(frameon=False)
fig.tight_layout()
save(fig, "loop_size")

# Isolated numerical comparisons.
diag = {"comparisons": {}, "pca_counts": pca_counts, "loop_ratios": loop_ratios}


def signature(r):
    return r["order"], r["data"], r["kmax"], r["prior"]


base = {signature(r): r for r in runs["main"]}
fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
for variation in ["step_half", "quadrature", "refined", "fft128"]:
    comp = []
    for rr in runs[variation]:
        bb = base[signature(rr)]
        errors = np.array(rr["errors"][:5], float)
        b = np.array(bb["errors"][:5], float)
        rel = np.abs(errors / b - 1)
        comp.append(
            {"signature": signature(rr), "max_target_change": float(np.nanmax(rel))}
        )
    diag["comparisons"][variation] = comp
    for prior, ax in zip([None, 1.0], axes):
        rr = [
            r
            for r in runs[variation]
            if r["data"] == "3x2" and r["order"] == "loop" and r["prior"] == prior
        ]
        y = [
            max(
                c["max_target_change"]
                for c in comp
                if tuple(c["signature"]) == signature(r)
            )
            for r in rr
        ]
        ax.plot([r["kmax"] for r in rr], 100 * np.array(y), "o-", label=variation)
for ax, pr in zip(axes, ["No priors", "Unit bias/ct priors"]):
    ax.set(
        xscale="log",
        yscale="log",
        title="3x2: " + pr,
        xlabel=r"$k_{\max}\ [h/{\rm Mpc}]$",
        ylabel="Max change among five target errors [%]",
    )
    ax.grid(alpha=0.15)
axes[0].legend(frameon=False)
fig.tight_layout()
save(fig, "convergence")

ja = np.load(ROOT / "results/main/arrays.npz")
jb = np.load(ROOT / "results/step_half/arrays.npz")
diag["derivatives"] = {}
for order in ["tree", "loop"]:
    a = ja[order + "_jacobian"]
    b = jb[order + "_jacobian"]
    ratio = np.linalg.norm(a - b, axis=(0, 1)) / np.maximum(
        np.linalg.norm(b, axis=(0, 1)), 1e-100
    )
    diag["derivatives"][order] = {
        "max_relative_column_norm": float(np.max(ratio)),
        "columns": ratio.tolist(),
    }

# Independent loop sums and quadrature precision.
sums = json.loads((ROOT / "results/halo_loop_sums.json").read_text())
diag["loop_sum_max_linear_fraction"] = max(abs(r["linear_fraction"]) for r in sums)
diag["loop_sum_max_relative"] = max(abs(r["relative_loop_sum"]) for r in sums)
fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
for bias, label in [
    ([1.0, 0.0, 0.0], "Unbiased matter limit"),
    ([2.4, 0.5, -0.3], "Biased halo"),
]:
    for muv in [0.0, 0.6, 1.0]:
        rr = [r for r in sums if r["f"] > 0 and r["bias"] == bias and r["mu"] == muv]
        axes[0].plot(
            [r["k"] for r in rr],
            [100 * r["linear_fraction"] for r in rr],
            "o-",
            label=f"{label}, mu={muv}",
        )
        axes[1].plot(
            [r["k"] for r in rr], [100 * r["relative_loop_sum"] for r in rr], "o-"
        )
axes[0].set(ylabel="Loop-sum residual / linear power [%]", xlabel=r"$k\ [h/{\rm Mpc}]$")
axes[1].set(
    ylabel="Loop-sum residual / reference loop sum [%]", xlabel=r"$k\ [h/{\rm Mpc}]$"
)
axes[1].set_yscale("symlog", linthresh=0.01)
for ax in axes:
    ax.axhline(0, color="black", lw=0.8)
    ax.grid(alpha=0.15)
axes[0].legend(frameon=False, fontsize=7)
fig.tight_layout()
save(fig, "physical_validation")


# Tables and flat CSV.
def esc(s):
    return s.replace("_", r"\_")


def num(x):
    return r"$\infty$" if x is None else f"{x:.4g}"


table = []
for cut in [0.2, 0.5, 1.0]:
    for data in colors:
        for order, prior in [("tree", None), ("loop", None), ("loop", 1.0)]:
            r = next(
                r
                for r in rows
                if r["kmax"] == cut
                and r["data"] == data
                and r["order"] == order
                and r["prior"] == prior
            )
            e = r["errors"][:4] + [r["fsigma8_error"]]
            table.append(
                f"{cut:g} & {data} & "
                + (
                    "tree"
                    if order == "tree"
                    else "loop, " + ("free" if prior is None else "prior 1")
                )
                + " & "
                + " & ".join(num(x) for x in e)
                + r" \\"
            )
(REPORT / "forecast_table.tex").write_text("\n".join(table))
with (ROOT / "results/forecast_summary.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(
        [
            "order",
            "data",
            "kmax",
            "prior",
            "sigma_alpha_w",
            "sigma_phi_w",
            "sigma_alpha_s",
            "sigma_phi_s",
            "sigma_f_fixed_sigma8",
            "sigma_fsigma8",
            "rank",
            "nparameters",
        ]
    )
    for r in rows:
        writer.writerow(
            [r["order"], r["data"], r["kmax"], r["prior"]]
            + r["errors"][:5]
            + [r["fsigma8_error"], r["rank"], r["nparameters"]]
        )

mode_table = []
for data in colors:
    for cut in [0.2, 0.5, 1.0]:
        r = next(r for r in selected(data, "loop", None) if r["kmax"] == cut)
        for i, m in enumerate(r["pca"][:3]):
            terms = " + ".join(
                f"{v:+.2f} " + r"\texttt{" + esc(n) + "}" for n, v in m["terms"].items()
            )
            mode_table.append(
                f"{data} & {cut:g} & {i + 1} & {num(m['sigma'])} & " + terms + r" \\"
            )
(REPORT / "mode_table.tex").write_text("\n".join(mode_table))

param_table = []
for name, value in metadata["defaults"].items():
    scale = metadata["scales"].get(name, 1.0)
    prior = (
        "none"
        if name in study.TARGETS or name.startswith("bias") or name.startswith("stoch")
        else r"$\sigma_{\rm EFT}$"
    )
    param_table.append(
        r"\texttt{"
        + esc(name)
        + "} & "
        + num(value)
        + " & "
        + num(scale)
        + " & "
        + num(value * scale)
        + " & "
        + prior
        + r" \\"
    )
(REPORT / "parameters.tex").write_text("\n".join(param_table))

fixed = []
for r in rows:
    if r["order"] != "loop" or r["prior"] is not None:
        continue
    f = arrays[f"loop_{r['data']}_{r['kmax']}_fisher"]
    n = 5 + sum(name in ("bias_F", "bias_eta_F", "bias_G") for name in r["names"])
    cov = np.linalg.inv(f[:n, :n])
    fixed.append(
        {
            "data": r["data"],
            "kmax": r["kmax"],
            "errors": np.sqrt(np.diag(cov)).tolist(),
            "fsigma8_error": float(np.sqrt(cov[4, 4]) * metadata["sigma8"]),
        }
    )
diag["loop_higher_fixed"] = fixed
(ROOT / "results/diagnostics.json").write_text(json.dumps(diag, indent=2))

physical_table = []
for r in sums:
    if r["f"] > 0 and r["mu"] in [0.0, 1.0]:
        physical_table.append(
            ("matter" if r["bias"][0] == 1 else "halo")
            + f" & {r['k']:g} & {r['mu']:g} & {r['fft']:.5g} & {r['reference']:.5g} & {100 * r['linear_fraction']:.3f} & {100 * r['relative_loop_sum']:.2f}"
            + r" \\"
        )
(REPORT / "physical_table.tex").write_text("\n".join(physical_table))

macros = {
    "LoopLinearResidual": 100 * diag["loop_sum_max_linear_fraction"],
    "LoopRelativeResidual": 100 * diag["loop_sum_max_relative"],
    "DerivativeResidual": 100 * diag["derivatives"]["loop"]["max_relative_column_norm"],
}
for name in ["step_half", "quadrature", "refined", "fft128"]:
    cases = [
        c["max_target_change"]
        for c in diag["comparisons"][name]
        if c["signature"][0] == "loop" and c["signature"][1] == "3x2"
    ]
    macro_name = (
        "Changefftresolution" if name == "fft128" else "Change" + name.replace("_", "")
    )
    macros[macro_name] = 100 * max(cases)
(REPORT / "numbers.tex").write_text(
    "\n".join(
        "\\newcommand{\\" + n + "}{" + f"{v:.3g}" + "}" for n, v in macros.items()
    )
)
print(
    json.dumps(
        {"macros": macros, "pca_counts": pca_counts, "loop_ratios": loop_ratios},
        indent=2,
    )
)

# Keep the final rule inside each input file: a rule following LaTeX's
# input macro can otherwise start a new alignment row after its trailing token.
for table_name in [
    "forecast_table.tex",
    "mode_table.tex",
    "parameters.tex",
    "physical_table.tex",
]:
    table_path = REPORT / table_name
    table_path.write_text(table_path.read_text() + "\n\\bottomrule\n")
