"""Saved three-profile BAO results; plot cuts never alter numerical tables."""

import json
from pathlib import Path

import numpy as np

PROFILES = ("full-compatibility", "fixed-compatibility", "accuracy")
COMPARISONS = (
    (PROFILES[1], PROFILES[0]),
    (PROFILES[2], PROFILES[0]),
    (PROFILES[2], PROFILES[1]),
)
REVISION = "early-lyaforecast-2026-09-18"


def build_tables(records):
    """Extract saved individual and joint contractions, preserving missing results.

    Parameters
    ----------
    records : iterable of (record, arrays)
        Evidence records with task/report metadata and saved arrays, or None on
        failure. Numerical qualification is explicit metadata, never inferred.
    """
    rows, weights = [], []
    seen = set()
    for record, arrays in records:
        task = record["task"]
        profile, index = task["profile"], task["bin"]
        if profile not in PROFILES or task.get("recipe_revision") != REVISION:
            raise ValueError(
                "comparison requires explicit current three-profile recipe"
            )
        key = (profile, index)
        if key in seen:
            raise ValueError("duplicate profile/bin")
        seen.add(key)
        fields = [f["id"] for f in task["fields"]]
        pairs = [tuple(p) for p in task["selected_pairs"]]
        if len(fields) != 5 or len(pairs) != (3 if index == 0 else 15):
            raise ValueError("comparison requires six-bin DESI-2 selection")
        report = record.get("report", {})
        settings = report.get("settings", {})
        status = record.get("numerical_status", "unqualified")
        if report.get("unresolved_controls"):
            status = "unresolved: " + str(report["unresolved_controls"])
        forest_rows = []
        for field, state in (
            (settings.get("forest_weighting") or {}).get("forests", {}).items()
        ):
            coefficients = state.get("coefficients", [])
            forest_rows.append(
                dict(
                    field=field,
                    magnitudes=settings.get(
                        "magnitude_nodes",
                        settings.get("pair_inputs", {})
                        .get(field + "_" + field, {})
                        .get("magnitudes"),
                    ),
                    weights=state.get("weights"),
                    A=state.get(
                        "A", coefficients[-2] if len(coefficients) >= 2 else None
                    ),
                    P_pixel=state.get(
                        "P_pixel", coefficients[-1] if coefficients else None
                    ),
                    iteration_status=state.get("result", state).get(
                        "status", "unreported"
                    ),
                )
            )
        weights.append(
            dict(
                profile=profile,
                bin=index,
                numerical_status=status,
                recipe_identity=settings.get("recipe_identity"),
                fields=forest_rows,
                qualification={
                    key: report.get(key)
                    for key in (
                        "passed",
                        "unresolved_controls",
                        "final_controls",
                        "metric_names",
                        "metrics",
                        "trial_contract",
                    )
                },
                forest_weighting=settings.get("forest_weighting"),
            )
        )
        for pair in [(i, j) for i in range(5) for j in range(i, 5)] + [None]:
            excluded = pair is not None and pair not in pairs
            prefix = "" if pair is None else "pair_"
            position = (
                None if pair is None else pairs.index(pair) if not excluded else None
            )

            def get(name):
                value = np.asarray(arrays[prefix + name])
                return value if pair is None else value[position]

            available = not excluded and arrays is not None
            rank = int(get("rank").reshape(-1)[0]) if available else None
            errors = get("errors") if available else None
            available = (
                available
                and rank == 2
                and np.all(np.isfinite(errors))
                and np.all(errors > 0)
            )
            reason = (
                "excluded_by_selection"
                if excluded
                else None
                if available
                else "failed_or_unconstrained"
            )
            fisher = get("fisher") if not excluded and arrays is not None else None
            rho = float(get("correlation")[0, 1]) if available else None
            area = (
                float(np.sqrt(np.linalg.det(get("covariance")))) if available else None
            )
            rows.append(
                dict(
                    profile=profile,
                    bin=index,
                    z=float(np.mean(task["bounds"])),
                    spectrum="joint"
                    if pair is None
                    else f"{fields[pair[0]]} × {fields[pair[1]]}",
                    selected_spectra=len(pairs),
                    availability=reason or "available",
                    numerical_status=status,
                    scientific_status="diagnostic; scientific acceptance not inferred",
                    sigma_parallel=float(errors[0]) if available else None,
                    sigma_transverse=float(errors[1]) if available else None,
                    correlation=rho,
                    rank=rank,
                    ellipse_area_over_pi=area,
                    fisher=fisher.tolist()
                    if fisher is not None and np.all(np.isfinite(fisher))
                    else None,
                    plot_omission=reason
                    or ("sigma_above_0.2" if np.any(errors > 0.2) else None),
                )
            )
    if seen != {(p, i) for p in PROFILES for i in range(6)}:
        raise ValueError("comparison requires exactly three profiles and six bins")
    lookup = {(r["profile"], r["bin"], r["spectrum"]): r for r in rows}
    ratios = []
    for numerator, denominator in COMPARISONS:
        for a in [r for r in rows if r["profile"] == numerator]:
            b = lookup[denominator, a["bin"], a["spectrum"]]
            valid = a["availability"] == b["availability"] == "available"
            ratios.append(
                dict(
                    comparison=f"{numerator}/{denominator}",
                    bin=a["bin"],
                    z=a["z"],
                    spectrum=a["spectrum"],
                    sigma_parallel_fraction=a["sigma_parallel"] / b["sigma_parallel"]
                    - 1
                    if valid
                    else None,
                    sigma_transverse_fraction=a["sigma_transverse"]
                    / b["sigma_transverse"]
                    - 1
                    if valid
                    else None,
                    ellipse_area_ratio=a["ellipse_area_over_pi"]
                    / b["ellipse_area_over_pi"]
                    if valid
                    else None,
                    numerator_status=a["numerical_status"],
                    denominator_status=b["numerical_status"],
                    plot_omission=a["plot_omission"] or b["plot_omission"],
                )
            )
    return dict(rows=rows, ratios=ratios, weighting=weights)


def render(table, output):
    """Save full JSON tables and matched panels with paired-component masking."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "three-profile-tables.json").write_text(
        json.dumps(table, indent=2, allow_nan=False) + "\n"
    )
    spectra = list(
        dict.fromkeys(r["spectrum"] for r in table["rows"] if r["spectrum"] != "joint")
    )
    statuses = "; ".join(
        f"{p}: "
        + ", ".join(
            sorted(
                {str(r["numerical_status"]) for r in table["rows"] if r["profile"] == p}
            )
        )
        for p in PROFILES
    )
    caption = (
        "Paired sigma > 0.2 cut; ratios require both operands. Bin 1 exclusions are gaps. Numerical qualification: "
        + statuses
    )
    for joint in (False, True):
        for ratio in (False, True):
            for component in ("parallel", "transverse"):
                fig, axes = plt.subplots(
                    1 if joint else 3,
                    1 if joint else 5,
                    figsize=(10, 5) if joint else (20, 11),
                    squeeze=False,
                )
                source = table["ratios"] if ratio else table["rows"]
                groups = [f"{a}/{b}" for a, b in COMPARISONS] if ratio else PROFILES
                for ax, spectrum in zip(
                    axes.flat, ["joint"] if joint else spectra, strict=True
                ):
                    for group in groups:
                        series = sorted(
                            (
                                r
                                for r in source
                                if r["spectrum"] == spectrum
                                and r["comparison" if ratio else "profile"] == group
                            ),
                            key=lambda r: r["bin"],
                        )
                        quantity = f"sigma_{component}" + ("_fraction" if ratio else "")
                        ax.plot(
                            [r["z"] for r in series],
                            [
                                r[quantity] if r["plot_omission"] is None else np.nan
                                for r in series
                            ],
                            ".-",
                            label=group,
                        )
                    ax.set_title(
                        "Joint: 3 spectra in bin 1; 15 in bins 2–6"
                        if joint
                        else spectrum
                    )
                    ax.set_xlabel("z")
                    ax.set_ylabel(
                        f"Δ sigma({component}) / sigma"
                        if ratio
                        else f"sigma(alpha_{component})"
                    )
                    ax.grid(alpha=0.2)
                fig.legend(
                    *axes.flat[0].get_legend_handles_labels(),
                    loc="upper center",
                    ncol=3,
                    fontsize=8,
                )
                fig.text(0.5, 0.01, caption, ha="center", fontsize=6, wrap=True)
                fig.tight_layout(rect=(0, 0.06, 1, 0.94))
                fig.savefig(
                    output
                    / f"{'joint' if joint else 'individual'}-{component}-{'fractions' if ratio else 'errors'}.png",
                    dpi=120,
                )
                plt.close(fig)

    field_names = sorted({f["field"] for r in table["weighting"] for f in r["fields"]})
    for field in field_names:
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        for index, ax in enumerate(axes.flat):
            for row in table["weighting"]:
                if row["bin"] != index:
                    continue
                for state in row["fields"]:
                    if (
                        state["field"] != field
                        or state["magnitudes"] is None
                        or state["weights"] is None
                    ):
                        continue
                    ax.plot(
                        state["magnitudes"],
                        state["weights"],
                        label=f"{row['profile']}: {state['iteration_status']}",
                    )
            ax.set_title(f"{field}, bin {index + 1}")
            ax.set_xlabel("magnitude")
            ax.set_ylabel("forest weight")
            if ax.lines:
                ax.legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(output / f"weights-{field}.png", dpi=120)
        plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        for ax, quantity in zip(axes, ("A", "P_pixel"), strict=True):
            for profile in PROFILES:
                points = [
                    (row["bin"] + 1, state[quantity])
                    for row in table["weighting"]
                    if row["profile"] == profile
                    for state in row["fields"]
                    if state["field"] == field and state[quantity] is not None
                ]
                if points:
                    points.sort()
                    ax.plot(*zip(*points, strict=True), ".-", label=profile)
            ax.set_xlabel("redshift bin")
            ax.set_ylabel(quantity)
            if ax.lines:
                ax.legend(fontsize=8)
        fig.suptitle(
            f"{field}: saved coefficients; iteration stability is not numerical qualification"
        )
        fig.tight_layout()
        fig.savefig(output / f"coefficients-{field}.png", dpi=120)
        plt.close(fig)
