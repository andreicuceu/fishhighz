"""Strict full-assignment gate: exact 7/39/78 and all declared policy swaps."""

import argparse
import json
from pathlib import Path

from fishhighz.validation.evidence import check, modern_requests
from fishhighz.validation.schema import request

# Bound full-assignment diagnostic inventory, without importing model setup.
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


def full_gate(bundle):
    root = Path(bundle)
    m = json.loads((root / "manifest.json").read_text())
    tasks = modern_requests("full", profiles=("compatibility", "accuracy"))
    diagnostics = [
        request(t["case"], t["bin"], "accuracy", kind="diagnostic", diagnostic_id=p)
        for t in tasks
        if t["profile"] == "accuracy"
        for p in POLICIES
    ]
    if m["requested"] != tasks or m["diagnostics_requested"] != diagnostics:
        raise ValueError(
            "full acceptance requires exactly 78 primary and 468 diagnostic requests"
        )
    if [r["task"] for r in m["records"]] != tasks or [
        r["task"] for r in m["diagnostics"]
    ] != diagnostics:
        raise ValueError("missing/extra/reordered full assignment records")
    return check(root)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bundle")
    a = p.parse_args()
    try:
        m = full_gate(a.bundle)
    except (ValueError, KeyError) as error:
        print(json.dumps(dict(passed=False, error=str(error)), indent=2))
        raise SystemExit(1)
    print(
        json.dumps(
            dict(
                passed=True,
                primary=len(m["records"]),
                diagnostics=len(m["diagnostics"]),
            ),
            indent=2,
        )
    )
