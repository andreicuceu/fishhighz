"""Validation-only CLI; full real execution is explicitly opt-in."""

import argparse
import json
from pathlib import Path

from fishhighz.validation.cases import verify_inventory
from fishhighz.validation.evidence import check, requests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["preflight", "run", "check"])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--suite", choices=["quick", "full"], default="quick")
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--bins", nargs="+", type=int)
    parser.add_argument("--negative-policy", choices=["reject", "floor_negative"])
    parser.add_argument("--allow-real-full", action="store_true")
    args = parser.parse_args()
    if args.mode == "check":
        print(json.dumps(dict(complete=check(args.output)["complete"])))
        return
    if args.reference is None:
        parser.error("--reference required")
    inventory = verify_inventory(args.reference / "examples/desi2")
    work = requests(args.suite, args.cases, args.bins)
    if args.mode == "preflight":
        print(json.dumps(dict(inventory=inventory, requested=work), indent=2))
        return
    parser.error(
        "Revision-2 real runs use scripts/desi2_profiles.py with an exact wheel and fresh reference bundle; historical r1 bundles require inspect_legacy"
    )


if __name__ == "__main__":
    main()
