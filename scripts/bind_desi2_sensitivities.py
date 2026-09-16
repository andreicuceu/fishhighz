"""Rebind resumed diagnostic comparisons to their verified final primary arrays.

Derive an exclusive evidence directory; numerical NPZ files remain byte-identical
and are hard-linked as immutable inputs. Prior reports and bundles are untouched.
"""

import argparse
import copy
import json
import os
from pathlib import Path

import numpy as np

from fishhighz.validation.evidence import canonical, digest, record_report
from fishhighz.validation.numerics import change


def bind(source, output):
    source, out = Path(source).resolve(), Path(output).resolve()
    original = json.loads((source / "manifest.json").read_text())
    if original["schema"] != 3 or not original["execution_finished"]:
        raise ValueError("require finished schema-3 evidence")
    if {r["task"]["case"] for r in original["records"]} != {"lya_qso_lbg_lae_15x2pt"}:
        raise ValueError("15x2pt-only comparison repair")
    out.mkdir(parents=True, exist_ok=False)
    manifest = copy.deepcopy(original)
    primaries = {}
    for record in manifest["records"]:
        if record["task"]["profile"] == "accuracy":
            path = source / record["arrays"]
            if path.parent != source or digest(path) != record["sha256"]:
                raise ValueError("primary path/hash mismatch")
            with np.load(path, allow_pickle=False) as data:
                primaries[record["task"]["bin"]] = (
                    data["fisher"],
                    data["pair_fisher"],
                    record["sha256"],
                )
    if set(primaries) != set(range(6)):
        raise ValueError("all six final primary results required")
    for group in ("records", "diagnostics"):
        for record in manifest[group]:
            if "arrays" not in record:
                continue
            path = source / record["arrays"]
            if path.parent != source or digest(path) != record["sha256"]:
                raise ValueError("diagnostic path/hash mismatch")
            os.link(path, out / path.name)
            report = record_report(source, record)
            if group == "diagnostics":
                with np.load(path, allow_pickle=False) as data:
                    f, pf = data["fisher"], data["pair_fisher"]
                primary, pairs, sha = primaries[record["task"]["bin"]]
                report["comparison_rebound"] = dict(
                    source_report_sha256=record["report_sha256"],
                    original_change=report["sensitivity"]["change"],
                    final_primary_sha256=sha,
                    interpretation="Comparison recomputed after roundoff-level primary C/F reassembly; diagnostic numerical arrays unchanged",
                )
                report["sensitivity"]["change"] = change(f, primary, pf, pairs, 1, 1)
            target = out / record["report_file"]
            if target.parent != out:
                raise ValueError("invalid report path")
            target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            record["report_sha256"] = digest(target)
            record["effective_hash"] = canonical(report)
    manifest["inputs"].update(
        {
            str(source / "manifest.json"): digest(source / "manifest.json"),
            str(Path(__file__).resolve()): digest(__file__),
        }
    )
    manifest["comparison_derivation"] = dict(
        source=str(source), numerical_arrays="byte-identical immutable hard links"
    )
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    bind(args.source, args.output)
