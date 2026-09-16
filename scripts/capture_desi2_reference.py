"""Explicit fresh seven-case actual reference capture with upstream arrays."""

import argparse

from fishhighz.validation.reference_capture import capture

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--suite", choices=["full"], required=True)
    args = p.parse_args()
    if not capture(args.reference, args.output)["complete"]:
        raise SystemExit(1)
