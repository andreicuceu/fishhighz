"""Read-only scientific table extraction and independent S4 operand checks."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from fishhighz.validation.numerics import contract, relative, wick


def vector(row):
    return np.r_[row["errors"], np.ravel(row["pair_errors"])]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root
    rows = json.loads((root / "summary.json").read_text())
    by = {(r["bin"], r["label"]): r for r in rows}
    effects = []
    interactions = []
    checks = []
    tables = []
    for row in rows:
        index, label = row["bin"], row["label"]
        a = dict(np.load(root / row["arrays"]))
        pairs = a["required_pairs"]
        selected = a["selected_pairs"]
        for n, (total, modes, cov, jac) in enumerate(
            zip(
                a["sample_total"],
                a["sample_modes"],
                a["sample_selected_covariance"],
                a["sample_observed_j"],
            )
        ):
            entries = {tuple(p): v for p, v in zip(pairs, total)}

            def p(i, j):
                return entries[tuple(sorted((int(i), int(j))))]

            scalar = np.array(
                [
                    [
                        (p(i, k) * p(j, ell) + p(i, ell) * p(j, k)) / modes
                        for k, ell in selected
                    ]
                    for i, j in selected
                ]
            )
            if not np.allclose(scalar, cov, rtol=2e-14, atol=0):
                raise ValueError("scalar Wick")
            np.linalg.cholesky(scalar)
            if not np.all(np.isfinite(jac.T @ np.linalg.solve(scalar, jac))):
                raise ValueError("sample solve")
        check = dict(bin=index, label=label, scalar_wick_cells=len(a["sample_total"]))
        if label in ("fixed", "accuracy"):
            c = wick(a["total"], a["modes"], pairs, selected, 5)
            f, individual = contract(c, a["observed_j"], independent=True)
            check.update(
                fisher_relative=relative(f, a["fisher"]),
                pair_fisher_relative=relative(individual, a["pair_fisher"]),
            )
            if max(check["fisher_relative"], check["pair_fisher_relative"]) > 5e-12:
                raise ValueError("endpoint independent Fisher")
        checks.append(check)
        for s in range(-1, len(selected)):
            f = a["fisher"] if s == -1 else a["pair_fisher"][s]
            cov = np.linalg.inv(f)
            errors = np.sqrt(np.diag(cov))
            corr = cov[0, 1] / np.prod(errors)
            error = a["errors"] if s == -1 else a["pair_errors"][s]
            if not np.allclose(errors, error, rtol=5e-12, atol=0):
                raise ValueError("error inversion")
            tables.append(
                dict(
                    bin=index + 1,
                    label=label,
                    spectrum="joint" if s == -1 else str(tuple(selected[s])),
                    parallel=errors[0],
                    transverse=errors[1],
                    correlation=corr,
                    ellipse_area=np.sqrt(np.linalg.det(cov)),
                    rank=np.linalg.matrix_rank(f),
                )
            )
        if label.startswith(("forward-", "reverse-")) and not label.endswith(
            "-refined"
        ):
            direction, change = label.split("-", 1)
            base = by[index, "fixed" if direction == "forward" else "accuracy"]
            delta = 100 * (vector(row) / vector(base) - 1)
            effects.append(
                dict(
                    bin=index + 1,
                    direction=direction,
                    change=change,
                    joint=delta[:2].tolist(),
                    individual=delta[2:].reshape(-1, 2).tolist(),
                )
            )
            if "+" in change:
                x, y = change.split("+")
                one = by[index, direction + "-" + x]
                two = by[index, direction + "-" + y]
                interaction = np.log(
                    vector(row) * vector(base) / (vector(one) * vector(two))
                )
                interactions.append(
                    dict(
                        bin=index + 1,
                        direction=direction,
                        change=change,
                        joint=interaction[:2].tolist(),
                        individual=interaction[2:].reshape(-1, 2).tolist(),
                    )
                )
    volume_checks = []
    for index in sorted({r["bin"] for r in rows}):
        for direction, endpoint in [("forward", "fixed"), ("reverse", "accuracy")]:
            base = by[index, endpoint]
            change = by[index, direction + "-volume"]
            b = json.loads((root / base["report"]).read_text())
            c = json.loads((root / change["report"]).read_text())
            ratio = np.sqrt(b["volume"] / c["volume"])
            err = float(np.max(abs(vector(change) / vector(base) / ratio - 1)))
            if err > 5e-12:
                raise ValueError("volume scaling")
            volume_checks.append(
                dict(
                    bin=index + 1,
                    direction=direction,
                    expected_ratio=ratio,
                    max_relative=err,
                )
            )
    output = dict(
        effects=effects,
        interactions=interactions,
        interaction_statistic="natural log of sigma_AB sigma_0 / (sigma_A sigma_B)",
        checks=checks,
        volume_checks=volume_checks,
        endpoint_closures=[
            {
                k: v
                for k, v in r.items()
                if k in ("bin", "label", "closure", "closure_max_error")
            }
            for r in rows
            if r["label"] in ("fixed", "accuracy")
        ],
        refinements=[
            dict(
                bin=r["bin"] + 1,
                label=r["label"],
                max_error_refinement=r["max_error_refinement"],
            )
            for r in rows
            if "max_error_refinement" in r
        ],
    )
    (root / "analysis.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n"
    )
    with (root / "forecasts.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(tables[0]))
        writer.writeheader()
        writer.writerows(tables)
    print(
        json.dumps(
            dict(
                trials=len(rows),
                forecasts=len(tables),
                scalar_cells=sum(c["scalar_wick_cells"] for c in checks),
                max_endpoint_closure=max(r.get("closure_max_error", 0) for r in rows),
                max_refinement=max(r.get("max_error_refinement", 0) for r in rows),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
