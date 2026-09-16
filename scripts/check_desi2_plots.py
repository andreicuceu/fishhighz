"""Verify all full-comparison plot tables, references, differences and artist data."""

import argparse
import json
from pathlib import Path

import numpy as np

from fishhighz.validation.cases import CASE_IDS, bins, selection
from fishhighz.validation.evidence import digest
from fishhighz.validation.plots import difference, draw_case, values


def check(plots, reference, cases=CASE_IDS):
    plots = Path(plots)
    reference = Path(reference)
    table = json.loads((plots / "tables.json").read_text())
    manifest = json.loads((plots / "plots-manifest.json").read_text())
    assert manifest["table_sha256"] == digest(plots / "tables.json")
    for name, sha in manifest["files"].items():
        assert Path(name).name == name and digest(plots / name) == sha
    expected = [(case, i) for case in cases for i in range(len(bins(case)))]
    assert [(r["case"], r["bin"]) for r in table["rows"]] == expected
    bundle = Path(table["source_bundle"])
    source = json.loads((bundle / "manifest.json").read_text())
    assert table["source_sha256"] == digest(bundle / "manifest.json")
    original = {
        case: json.loads((reference / case / "result.json").read_text())
        for case in cases
    }
    meta = {
        case: json.loads((reference / case / "metadata.json").read_text())
        for case in cases
    }
    primary = {
        (r["task"]["case"], r["task"]["bin"], r["task"]["profile"]): r
        for r in source["records"]
    }
    for row in table["rows"]:
        case, index = row["case"], row["bin"]
        sel = selection(case)
        assert (
            row["bounds"] == list(bins(case)[index])
            and [r["pair"] for r in row["pairs"]] == sel.selected_pairs.tolist()
        )
        for profile in ("compatibility", "accuracy"):
            record = primary[(case, index, profile)]
            if "arrays" not in record:
                assert row["profiles"][profile]["values"] == [None] * 3
                continue
            with np.load(bundle / record["arrays"], allow_pickle=False) as a:
                assert row["profiles"][profile]["values"] == values(a["fisher"])
                assert [r["profiles"][profile] for r in row["pairs"]] == [
                    values(f) for f in a["pair_fisher"]
                ]
        ref = original[case]
        np.testing.assert_allclose(
            row["profiles"]["reference"]["values"],
            [ref[n][index] for n in ("sigma_ap", "sigma_at", "corr_coef")],
            rtol=5e-12,
            atol=0,
        )
        for pair, label in zip(
            row["pairs"], meta[case]["records"][index]["selected_labels"]
        ):
            np.testing.assert_allclose(
                pair["profiles"]["reference"],
                [ref[label][n][index] for n in ("sigma_ap", "sigma_at", "corr_coef")],
                rtol=5e-12,
                atol=0,
            )
        for entry in [row, *row["pairs"]]:
            v = {
                p: (x["values"] if entry is row else x)
                for p, x in entry["profiles"].items()
            }
            for x, y in [
                ("compatibility", "reference"),
                ("accuracy", "reference"),
                ("accuracy", "compatibility"),
            ]:
                assert entry["differences"][x + "/" + y] == difference(v[x], v[y])
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    artists = 0
    for case in cases:
        rows = [r for r in table["rows"] if r["case"] == case]
        for pair in [None, *range(len(rows[0]["pairs"]))]:
            fig = draw_case(plt, case, rows, pair)
            entries = [r if pair is None else r["pairs"][pair] for r in rows]
            for j, ax in enumerate(fig.axes[:3]):
                for container, profile in zip(
                    ax.containers, ("reference", "compatibility", "accuracy")
                ):
                    expected_y = [
                        e["profiles"][profile]["values"][j]
                        if pair is None
                        else e["profiles"][profile][j]
                        for e in entries
                    ]
                    np.testing.assert_array_equal(
                        container.lines[0].get_xdata(), [r["centre"] for r in rows]
                    )
                    np.testing.assert_allclose(
                        container.lines[0].get_ydata(orig=False),
                        [np.nan if x is None else x for x in expected_y],
                        rtol=0,
                        atol=0,
                        equal_nan=True,
                    )
                    artists += 1
            for j, ax in enumerate(fig.axes[3:]):
                for line, profile in zip(ax.lines[:2], ("compatibility", "accuracy")):
                    expected_y = [
                        e["differences"][profile + "/reference"][j] for e in entries
                    ]
                    np.testing.assert_allclose(
                        line.get_ydata(orig=False),
                        [np.nan if x is None else x for x in expected_y],
                        rtol=0,
                        atol=0,
                        equal_nan=True,
                    )
                    artists += 1
            plt.close(fig)
    return dict(
        case_bins=len(expected),
        profiles=2 * len(expected),
        checked_artists=artists,
        files=len(manifest["files"]),
        passed=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plots", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--cases", nargs="+", choices=CASE_IDS, default=CASE_IDS)
    a = p.parse_args()
    print(json.dumps(check(a.plots, a.reference, a.cases), indent=2))
