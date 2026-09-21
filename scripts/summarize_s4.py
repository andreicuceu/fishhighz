"""Summarize saved S4 controls without evaluating a model or loading cell arrays."""

import argparse
import json
from pathlib import Path

import numpy as np

FAMILIES = {
    "magnitude": "107-node full-endpoint rectangular → composite Gauss–Legendre",
    "negative": "Signed interpolated density → negative-density floor 1e−20",
    "fourier": "Legacy rectangular k/midpoint mu → Gauss–Legendre",
    "volume": "Centre-volume approximation → integrated volume",
    "redshift": "Arithmetic mean/damping redshift → geometric centre",
    "response": "Resolution sigma=c/R → FWHM conversion",
    "growth": "EdS power scaling → CAMB sigma8 power scaling",
    "power": "CAMB interpolated linear spectrum → supplied Vega template",
    "covariance_damping": "Undamped covariance power → damped template wiggle",
    "cross_damping": "Mixed forest–galaxy forest width → mean auto-width squares",
    "peak": "Observed-power polynomial residual → template PK−PKSB wiggle",
    "operator": "Backward k derivative → inverse AP derivative with fixed observed response",
    "field_response": "Pair response/forest roundtrip → per-field response/exact angle",
    "constants": "c=299800 → 299792.458 in velocity/noise/response, volume fixed",
    "stopping": "Weight convergence rtol 1e−4 → 1e−5",
}


def vector(row):
    return np.r_[row["errors"], np.ravel(row["pair_errors"])]


def difference(new, old):
    return (100 * (np.asarray(new) / np.asarray(old) - 1)).tolist()


def small_arrays(root, row):
    keys = (
        "fisher",
        "rank",
        "pair_rank",
        "correlation",
        "pair_correlation",
        "covariance",
        "pair_covariance",
        "selected_pairs",
    )
    with np.load(root / row["arrays"], allow_pickle=False) as data:
        return {key: data[key] for key in keys}


def geometry_change(new, old):
    """Ellipse area is proportional to sqrt(det(parameter covariance))."""
    return {
        "rank": np.asarray(new["rank"]).tolist(),
        "rank_changed": bool(np.any(new["rank"] != old["rank"])),
        "pair_rank": new["pair_rank"].tolist(),
        "pair_rank_changed": bool(np.any(new["pair_rank"] != old["pair_rank"])),
        "correlation": float(new["correlation"][0, 1]),
        "correlation_change": float(
            new["correlation"][0, 1] - old["correlation"][0, 1]
        ),
        "pair_correlation_change": (
            new["pair_correlation"][:, 0, 1] - old["pair_correlation"][:, 0, 1]
        ).tolist(),
        "ellipse_area_change_percent": float(
            100
            * (
                np.sqrt(
                    np.linalg.det(new["covariance"]) / np.linalg.det(old["covariance"])
                )
                - 1
            )
        ),
        "pair_ellipse_area_change_percent": (
            100
            * (
                np.sqrt(
                    np.linalg.det(new["pair_covariance"])
                    / np.linalg.det(old["pair_covariance"])
                )
                - 1
            )
        ).tolist(),
    }


def summarize(root):
    records = json.loads((root / "summary.json").read_text())
    by_bin = {i: {r["label"]: r for r in records if r["bin"] == i} for i in range(6)}
    output = {
        "complete": True,
        "directions": {
            "forward": "100*(sigma(fixed with one accuracy setting)/sigma(fixed)-1)",
            "reverse": "100*(sigma(accuracy with one legacy setting)/sigma(accuracy)-1)",
            "interaction": "ln[sigma_AB*sigma_0/(sigma_A*sigma_B)]; no sign reversal",
        },
        "families": FAMILIES,
        "bins": [],
    }
    for index, rows in by_bin.items():
        required = ["fixed", "accuracy", "accuracy-refined"]
        required += [f"{d}-{k}" for d in ("forward", "reverse") for k in FAMILIES]
        required += ["forward-rect425", "reverse-rect425"]
        missing = sorted(set(required) - set(rows))
        # Each GL single-switch control and the accuracy endpoint must have its
        # matching numerical refinement; two-switch interactions were not refined.
        required_refined = (
            ["forward-fourier-refined"]
            + [f"reverse-{k}-refined" for k in FAMILIES if k != "fourier"]
            + ["reverse-rect425-refined"]
        )
        missing += sorted(set(required_refined) - set(rows))
        if missing:
            raise ValueError(f"Bin {index + 1} incomplete: {missing}")
        values = {key: small_arrays(root, row) for key, row in rows.items()}
        item = {
            "bin": index + 1,
            "selected_pairs": values["fixed"]["selected_pairs"].tolist(),
            "endpoint_change_percent": difference(
                rows["accuracy"]["errors"], rows["fixed"]["errors"]
            ),
            "endpoint_pair_change_percent": difference(
                rows["accuracy"]["pair_errors"], rows["fixed"]["pair_errors"]
            ),
            "endpoint_geometry": geometry_change(values["accuracy"], values["fixed"]),
            "closure": {
                key: {
                    name: rows[key][name] for name in ("closure", "closure_max_error")
                }
                for key in ("fixed", "accuracy")
            },
            "switches": {},
            "refinements": {},
            "interactions": {},
        }
        for key in list(FAMILIES) + ["rect425"]:
            item["switches"][key] = {}
            for direction, base in [("forward", "fixed"), ("reverse", "accuracy")]:
                label = direction + "-" + key
                item["switches"][key][direction] = {
                    "joint_percent": difference(
                        rows[label]["errors"], rows[base]["errors"]
                    ),
                    "individual_percent": difference(
                        rows[label]["pair_errors"], rows[base]["pair_errors"]
                    ),
                    "geometry": geometry_change(values[label], values[base]),
                }
        for label, row in rows.items():
            if label.endswith("-refined"):
                base = label.removesuffix("-refined")
                item["refinements"][base] = {
                    "joint_fisher_relative": float(
                        np.linalg.norm(values[label]["fisher"] - values[base]["fisher"])
                        / np.linalg.norm(values[base]["fisher"])
                    ),
                    "joint_percent": difference(row["errors"], rows[base]["errors"]),
                    "individual_percent": difference(
                        row["pair_errors"], rows[base]["pair_errors"]
                    ),
                    "max_all_percent": float(
                        np.max(abs(100 * (vector(row) / vector(rows[base]) - 1)))
                    ),
                }
            elif "+" in label:
                direction, pair = label.split("-", 1)
                a, b = pair.split("+")
                base = "fixed" if direction == "forward" else "accuracy"
                residual = np.log(
                    vector(row)
                    * vector(rows[base])
                    / (
                        vector(rows[f"{direction}-{a}"])
                        * vector(rows[f"{direction}-{b}"])
                    )
                )
                item["interactions"][label] = {
                    "joint_log": residual[:2].tolist(),
                    "individual_log": residual[2:].reshape(-1, 2).tolist(),
                    "max_abs_log": float(np.max(abs(residual))),
                }
        if not item["interactions"]:
            raise ValueError(f"Bin {index + 1} lacks interaction controls")
        output["bins"].append(item)
    output["closure_max_error_fraction"] = max(
        c["closure_max_error"] for b in output["bins"] for c in b["closure"].values()
    )
    output["refinement_max_error_percent"] = max(
        r["max_all_percent"] for b in output["bins"] for r in b["refinements"].values()
    )
    output["interaction_max_abs_log"] = max(
        r["max_abs_log"] for b in output["bins"] for r in b["interactions"].values()
    )
    output["refinement_max_joint_fisher_relative"] = max(
        r["joint_fisher_relative"]
        for b in output["bins"]
        for r in b["refinements"].values()
    )
    output["refinement_component_max_percent"] = {
        kind: np.max(
            np.abs(
                np.concatenate(
                    [
                        np.asarray(r[kind]).reshape(-1, 2)
                        for b in output["bins"]
                        for r in b["refinements"].values()
                    ]
                )
            ),
            axis=0,
        ).tolist()
        for kind in ("joint_percent", "individual_percent")
    }
    geometries = [
        d["geometry"]
        for b in output["bins"]
        for switch in b["switches"].values()
        for d in switch.values()
    ]
    output["any_rank_change"] = any(
        g["rank_changed"] or g["pair_rank_changed"] for g in geometries
    )
    return output


def interval(values):
    return f"{min(values):+.6g} to {max(values):+.6g}"


def report(result):
    lines = [
        "# S4 joint BAO effects from saved direct controls",
        "",
        "Each range spans bins 1–6; components are radial and transverse. "
        "Forward changes introduce one accuracy setting into fixed-compatibility. "
        "Reverse changes restore one legacy setting in accuracy; their signs are "
        "**not reversed**. Percentages are relative to that direction’s own endpoint. "
        "They are conditional effects, not additive shares of the endpoint difference.",
        "",
        "| Setting (legacy → accuracy) | Forward radial (%) | Forward transverse (%) | Reverse radial (%) | Reverse transverse (%) |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for key, description in FAMILIES.items():
        ranges = [
            interval(
                [b["switches"][key][d]["joint_percent"][c] for b in result["bins"]]
            )
            for d in ("forward", "reverse")
            for c in range(2)
        ]
        lines.append("| " + " | ".join([description] + ranges) + " |")
    lines += [
        "",
        "Per-bin values, individual-spectrum effects, exact tested interactions, "
        "rank, correlation and ellipse-area changes are retained in "
        "[attribution-summary.json](../.validation/s4/attribution-summary.json). "
        "Source definitions and unchanged assumptions are in "
        "[s4-profile-inventory.md](s4-profile-inventory.md).",
        "",
        f"Maximum endpoint error closure residual: {result['closure_max_error_fraction']:.6g} fraction. "
        f"Maximum tested refinement change across joint and individual errors: {result['refinement_max_error_percent']:.6g}%. "
        f"Maximum joint Fisher Frobenius relative refinement: {result['refinement_max_joint_fisher_relative']:.6g}. "
        "These refinements cover the Gauss–Legendre endpoint and single-switch controls "
        "(including the accuracy-side rectangular-425 check), not all two-switch interactions "
        "or the literal rectangular Fourier estimator.",
        "",
        "The interaction statistic is ln[sigma_AB sigma_0/(sigma_A sigma_B)], "
        "computed separately in each direction, bin, spectrum and component. "
        f"The maximum absolute value across retained controls is {result['interaction_max_abs_log']:.6g}. "
        "A nonzero value measures failure of multiplicative single-switch attribution; "
        "it is not a statistical significance.",
        "",
        "| Bin | Accuracy/fixed radial (%) | Accuracy/fixed transverse (%) | Correlation change | Ellipse-area change (%) | Rect425 forward radial (%) | Rect425 forward transverse (%) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for b in result["bins"]:
        row = [str(b["bin"])] + [f"{x:+.7g}" for x in b["endpoint_change_percent"]]
        row += [
            f"{b['endpoint_geometry']['correlation_change']:+.7g}",
            f"{b['endpoint_geometry']['ellipse_area_change_percent']:+.7g}",
        ]
        row += [
            f"{x:+.7g}" for x in b["switches"]["rect425"]["forward"]["joint_percent"]
        ]
        lines.append("| " + " | ".join(row) + " |")
    lines += [
        "",
        f"Any joint or individual rank change among single-switch controls: {result['any_rank_change']}.",
        "Maximum absolute refined joint changes (radial, transverse), percent: "
        + str(result["refinement_component_max_percent"]["joint_percent"])
        + ". "
        "Corresponding individual maxima: "
        + str(result["refinement_component_max_percent"]["individual_percent"])
        + ".",
        "",
        "| Interaction | Tested bins | Joint radial log range | Joint transverse log range | Largest individual absolute log |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    labels = sorted({label for b in result["bins"] for label in b["interactions"]})
    for label in labels:
        tested = [
            (b["bin"], b["interactions"][label])
            for b in result["bins"]
            if label in b["interactions"]
        ]
        lines.append(
            "| "
            + " | ".join(
                [
                    label,
                    ", ".join(str(i) for i, _ in tested),
                    interval([r["joint_log"][0] for _, r in tested]),
                    interval([r["joint_log"][1] for _, r in tested]),
                    f"{max(np.max(abs(np.asarray(r['individual_log']))) for _, r in tested):.6g}",
                ]
            )
            + " |"
        )
    lines += [
        "",
        "The rectangular-425 reverse check restores a finer rectangular measure on accuracy: "
        + "; ".join(
            f"bin {b['bin']}: "
            + ", ".join(
                f"{x:+.6g}%"
                for x in b["switches"]["rect425"]["reverse"]["joint_percent"]
            )
            for b in result["bins"]
        )
        + ".",
        "",
    ]
    lines += [
        "",
        "Shared prescriptions have no independent parameter change here: "
        "early-lyaforecast recurrence; representative angular/velocity mode; "
        "density normalization and support; raw first-spacing density-width convention; "
        "SNR bright clamp, sentinel and floor; P1D prescription; independent-sampling "
        "noise; forecast selection; covariance-derivative exclusion. Pair versus field "
        "instrument response agrees for the equal forest pixel widths of these inputs. "
        "The field-response switch additionally removes the legacy forest-auto "
        "mu=k_parallel/(k+1e−10) roundtrip and therefore can retain a tiny nonzero effect.",
        "",
        "The redshift switch changes the derivative mean and damping evaluation; "
        "covariance, source queries and geometry already use geometric redshift. "
        "The direct reconstruction uses the newly prepared CAMB redshift samples for "
        "galaxy f in both hybrid endpoints. The captured legacy forecast interpolated "
        "f between its arithmetic-bin samples. That small interpolation-versus-exact "
        "difference is not separately switched by redshift: it is included in the "
        "reported endpoint reconstruction residual, with no independently measured "
        "BAO contribution assigned. "
        "The covariance-damping hybrid adds the template-wiggle damping correction "
        "to the selected linear power. The peak switch adopts the template wiggle "
        "even when the covariance power source is legacy. These hybrid definitions "
        "must be retained when interpreting forward/reverse asymmetry. The constants "
        "switch changes remaining velocity, noise and response inputs with volume held fixed. The separate volume switch replaces saved centre volume with the physical integral and therefore already includes the c convention within volume. "
        "Refined rectangular integration, floor/width changes, fixed-comoving "
        "representative modes, McDonald and W12 are not new adopted profiles.",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(".validation/s4/direct-r1"))
    parser.add_argument(
        "--output", type=Path, default=Path(".validation/s4/attribution-summary.json")
    )
    parser.add_argument(
        "--report", type=Path, default=Path("reviews/s4-impact-table.md")
    )
    args = parser.parse_args()
    result = summarize(args.input)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    args.report.write_text(report(result))
    print(
        json.dumps(
            {
                k: v
                for k, v in result.items()
                if k.endswith(("fraction", "percent", "log"))
            }
        )
    )


if __name__ == "__main__":
    main()
