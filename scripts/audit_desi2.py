"""Offline per-bin condition/convergence and policy audit, including failures."""

import argparse
import json
from itertools import islice
from pathlib import Path

import numpy as np

from fishhighz.validation.evidence import digest
from fishhighz.validation.plots import records


def audit(bundle):
    result = []
    primary_count = len(
        json.loads((Path(bundle) / "manifest.json").read_text())["records"]
    )
    # This report covers primary bins. The plot reader checks all diagnostics too.
    for record, arrays in islice(records(bundle), primary_count):
        task = record["task"]
        if task["kind"] == "diagnostic":
            continue
        row = dict(
            case=task["case"],
            bin=task["bin"],
            profile=task["profile"],
            status=record["status"],
            passed=record.get("scientific_passed", False),
            error=record.get("error", record.get("scientific_error")),
        )
        if arrays is not None:
            report = record["report"]
            settings = report["settings"]
            eigen = np.linalg.eigvalsh(arrays["selected_covariance"])
            from fishhighz.validation.numerics import field_matrix

            active = np.unique(arrays["selected_pairs"])
            total = field_matrix(
                arrays["total"], arrays["required_pairs"], len(task["fields"])
            )[:, active][:, :, active]
            scale = np.sqrt(np.diagonal(total, axis1=1, axis2=2))
            normalized = total / scale[:, :, None] / scale[:, None, :]
            field_eigen = np.linalg.eigvalsh(normalized)[:, 0]
            row["field_psd_material_violations"] = int(
                np.count_nonzero(field_eigen < -64 * np.finfo(float).eps * len(active))
            )
            row["field_min_normalized_eigenvalue"] = float(np.min(field_eigen))

            row.update(
                fisher_condition=float(np.linalg.cond(arrays["fisher"])),
                covariance_condition_max=float(np.max(eigen[:, -1] / eigen[:, 0])),
                field_negative_nodes=int(
                    np.count_nonzero(arrays["field_min_eigenvalue"] < 0)
                ),
                convergence=dict(zip(report["metric_names"], report["metrics"])),
                controls=report.get("final_controls"),
                levels=report.get("actual_levels"),
                unresolved=report.get("unresolved_controls"),
                comparison=report.get("comparison"),
                field_policies={},
            )
            for name, sample in settings.get("samples", {}).items():
                density = np.asarray(sample["density"])
                measure = np.asarray(sample["quadrature"])
                masks = sample["density_diagnostics"]["masks"]
                contributions = {
                    key: float(np.sum(density * np.asarray(mask) * measure))
                    for key, mask in masks.items()
                }
                snr = sample.get("snr_diagnostics", {})
                row["field_policies"][name] = dict(
                    density_counts=sample["density_diagnostics"]["counts"],
                    snr_counts=snr.get("counts", {}),
                    integrated_density=float(density @ measure),
                    density_contributions=contributions,
                    snr_density_contributions={
                        key: float(np.sum(density * np.asarray(mask) * measure))
                        for key, mask in snr.get("masks", {}).items()
                    },
                )
        result.append(row)
    return dict(source_sha256=digest(Path(bundle) / "manifest.json"), rows=result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    out = Path(a.output)
    if out.exists():
        raise FileExistsError(out)
    out.write_text(json.dumps(audit(a.bundle), indent=2, allow_nan=False) + "\n")
