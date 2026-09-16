"""Schema-3 scientific gate for six 15x2pt bins, both profiles and 72 diagnostics."""

import argparse
import json
from pathlib import Path

from fishhighz.validation.evidence import check, modern_requests
from fishhighz.validation.schema import request

POLICIES = (
    "floor_low",
    "floor_high",
    "remove_bright",
    "remove_sentinel",
    "width_minus",
    "width_plus",
    "growth_eds",
    "legacy_resolution",
    "rectangular_magnitude",
    "no_reconstruction",
    "arithmetic_mean",
    "three_weights",
)


def scoped_gate(bundle):
    root = Path(bundle)
    m = json.loads((root / "manifest.json").read_text())
    tasks = modern_requests(
        "full", ["lya_qso_lbg_lae_15x2pt"], profiles=("compatibility", "accuracy")
    )
    diagnostics = [
        request(t["case"], t["bin"], "accuracy", kind="diagnostic", diagnostic_id=p)
        for t in tasks
        if t["profile"] == "accuracy"
        for p in POLICIES
    ]
    if m.get("schema") != 3:
        raise ValueError("15x2pt acceptance requires schema 3 trial bindings")
    if m.get("requested") != tasks or m.get("diagnostics_requested") != diagnostics:
        raise ValueError(
            "15x2pt acceptance requires exactly 12 primary and 72 diagnostic requests"
        )
    if [r["task"] for r in m["records"]] != tasks or [
        r["task"] for r in m["diagnostics"]
    ] != diagnostics:
        raise ValueError("missing/extra/reordered 15x2pt records")
    return check(root)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle")
    args = parser.parse_args()
    try:
        result = scoped_gate(args.bundle)
    except (ValueError, KeyError) as error:
        print(json.dumps(dict(passed=False, error=str(error))))
        raise SystemExit(1)
    print(
        json.dumps(
            dict(passed=True, primary=12, selected_pair_results=180, diagnostics=72)
        )
    )
