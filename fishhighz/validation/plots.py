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
    """Read exact inventory and validate mathematics, retaining scientific failures.

    Parameters
    ----------
    root : str or pathlib.Path
        Root directory of the input checkout or saved evidence bundle.

    Returns
    -------
    records : iterator of tuple
        Manifest records paired with verified numerical arrays, or None for
        records without numerical output.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
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
            with np.load(path, allow_pickle=False) as archive:
                arrays = {k: archive[k] for k in archive.files}
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
    """Extract rank-aware two-parameter BAO uncertainties and correlation.

    Parameters
    ----------
    fisher : array_like, shape (n_parameter, n_parameter)
        Fisher information in inverse products of parameter units.

    Returns
    -------
    values : list of float or None
        Parallel error, transverse error and correlation; unconstrained
        quantities are None.
    """
    info = information(fisher)
    return [
        float(info["errors"][i]) if info["constrained"][i] else None for i in range(2)
    ] + [float(info["correlation"][0, 1]) if np.all(info["constrained"]) else None]


def difference(a, b):
    """Compare available BAO errors and correlation in their reported units.

    Parameters
    ----------
    a : list of float or None
        Comparison parallel error, transverse error and correlation.
    b : list of float or None
        Reference parallel error, transverse error and correlation.

    Returns
    -------
    differences : list of float or None
        Percentage changes in both errors and absolute correlation change; all
        None if either result is unavailable.
    """
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
    """Every primary and selected pair survives even if execution/science failed.

    Parameters
    ----------
    bundle : str or pathlib.Path
        Saved validation evidence directory.
    output : str or pathlib.Path
        Destination directory for the generated evidence or figures.

    Returns
    -------
    table : dict
        All primary, individual-spectrum, sensitivity and sequential-attribution
        rows with explicit failures.

    Notes
    -----
    Creates a new directory and writes JSON/CSV tables plus any successful attribution arrays.
    """
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    rows = {}
    diagnostics = []
    attributions = []
    pending = {}
    for record, record_arrays in records(bundle):
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
            pair_selection = selection(key[0])
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
                        label=f"{pair_selection.fields[i].id} x {pair_selection.fields[j].id}",
                        profiles={},
                    )
                    for i, j in pair_selection.selected_pairs
                ],
            )
        row = rows[key]
        profile = task["profile"]
        row["profiles"][profile] = dict(
            values=values(record_arrays["fisher"])
            if record_arrays is not None
            else [None] * 3,
            passed=record.get("scientific_passed", False),
            metrics=record.get("report", {}).get("metrics"),
            controls=record.get("report", {}).get("final_controls"),
            error=record.get("scientific_error", record.get("error")),
        )
        for i, pair in enumerate(row["pairs"]):
            pair["profiles"][profile] = (
                values(record_arrays["pair_fisher"][i])
                if record_arrays is not None
                else [None] * 3
            )
        if record_arrays is not None and profile == "compatibility":
            row["profiles"]["reference"] = dict(
                values=values(record_arrays["reference_fisher"]), passed=True
            )
            for i, pair in enumerate(row["pairs"]):
                pair["profiles"]["reference"] = values(
                    record_arrays["reference_pair_fisher"][i]
                )
            pending[key] = (
                record_arrays,
                record["report"]["settings"]["grid"]["volume"],
            )
        elif record_arrays is not None and key in pending:
            legacy, volume = pending.pop(key)
            try:
                arrays, report = chain(
                    legacy,
                    record_arrays,
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
                record_arrays["fisher"], legacy["fisher"]
            )
        if profile == "accuracy":
            for entry in [row, *row["pairs"]]:
                profiles = entry["profiles"]
                profile_values = {
                    name: (item["values"] if entry is row else item)
                    for name, item in profiles.items()
                }
                entry["differences"] = {
                    f"{x}/{y}": difference(
                        profile_values.get(x, [None] * 3),
                        profile_values.get(y, [None] * 3),
                    )
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
    with (out / "tables.csv").open("w") as table_stream:
        table_writer = csv.writer(table_stream)
        table_writer.writerow(
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
                for profile, profile_values in entry["profiles"].items():
                    value = profile_values["values"] if entry is row else profile_values
                    diff = entry.get("differences", {}).get(
                        profile + "/reference",
                        [0, 0, 0] if profile == "reference" else [None] * 3,
                    )
                    table_writer.writerow(
                        [
                            row["case"],
                            row["bin"],
                            entry.get("label", "combined"),
                            profile,
                            *value,
                            row["profiles"].get(profile, {}).get("passed", False),
                            *diff,
                        ]
                    )
    return table


def draw_case(plt, case, rows, pair=None):
    """Draw one case using saved bin centres and availability flags.

    Parameters
    ----------
    plt : module
        Imported matplotlib.pyplot module.
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.
    rows : list of dict
        Checked table rows for this case, containing bounds and profile
        comparisons.
    pair : int or None
        Individual-spectrum index; None draws the joint result. Default is
        ``None``.

    Returns
    -------
    figure : matplotlib.figure.Figure
        Joint or individual-spectrum BAO errors, correlation and reference
        differences.
    """
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    styles = dict(
        reference=("black", "o", "actual lyaforecast"),
        compatibility=("tab:blue", "D", "FishHighz maximum compatibility"),
        accuracy=("tab:orange", "x", "FishHighz maximum accuracy"),
    )
    redshift_centres = [r["centre"] for r in rows]
    width = [(r["bounds"][1] - r["bounds"][0]) / 2 for r in rows]
    entries = [r if pair is None else r["pairs"][pair] for r in rows]
    for profile, (color, marker, label) in styles.items():
        val = [e["profiles"].get(profile, {"values": [None] * 3}) for e in entries]
        val = [
            profile_values["values"]
            if isinstance(profile_values, dict)
            else profile_values
            for profile_values in val
        ]
        for j in range(3):
            plotted_values = [
                np.nan if profile_values[j] is None else profile_values[j]
                for profile_values in val
            ]
            axes[0, j].errorbar(
                redshift_centres,
                plotted_values,
                xerr=width,
                color=color,
                marker=marker,
                mfc="none",
                label=label,
                lw=1,
                ms=5,
                capsize=2,
            )
            if profile != "reference":
                profile_differences = [
                    e.get("differences", {}).get(profile + "/reference", [None] * 3)[j]
                    for e in entries
                ]
                axes[1, j].plot(
                    redshift_centres,
                    [np.nan if x is None else x for x in profile_differences],
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
    """Render only checked table values; optional Matplotlib is imported lazily.

    Parameters
    ----------
    table : dict
        Checked comparison table returned by tables.
    output : str or pathlib.Path
        Destination directory for the generated evidence or figures.

    Returns
    -------
    manifest : dict
        Table/source hashes and hashes of every generated figure.

    Notes
    -----
    Selects the noninteractive Agg backend, writes PNG/PDF figures and a figure manifest, and closes created figures.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    out = Path(output)
    files = []

    def save(fig, name):
        """Write both figure formats and retain their paths.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Completed scientific comparison figure.
        name : str
            Quantity or record label used in diagnostics.

        Notes
        -----
        Writes PNG and PDF files, appends their paths to the enclosing list, then closes the figure.
        """
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
    for profile, color in [("compatibility", "tab:blue"), ("accuracy", "tab:orange")]:
        for j, ax in enumerate(axes):
            vals = [
                r.get("differences", {}).get(profile + "/reference", [None] * 3)[j]
                for r in table["rows"]
            ]
            ax.plot(
                range(len(vals)),
                [np.nan if x is None else x for x in vals],
                marker="o",
                ms=4,
                label=profile,
                color=color,
            )
            bad = [
                i
                for i, r in enumerate(table["rows"])
                if not r["profiles"].get(profile, {}).get("passed")
            ]
            ax.scatter(
                bad,
                [vals[i] for i in bad],
                marker="x",
                s=60,
                color="red",
                label="unresolved" if profile == "accuracy" else None,
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
            stage_values = np.array(a["values"])
            axes[1].plot(
                range(5),
                100 * (stage_values[:, 0] / stage_values[0, 0] - 1),
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
        files={path.name: digest(path) for path in files},
    )
    (out / "plots-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
