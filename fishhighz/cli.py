"""Command-line entry point for the native FishHighz forecast facade."""

from __future__ import annotations

import argparse

from .public import Forecast


def main(argv=None):
    """Run one native INI forecast and save its result directory."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ini", nargs="?", default=None, help="native FishHighz INI")
    parser.add_argument("--output", required=True, help="new result directory")
    args = parser.parse_args(argv)
    Forecast(args.ini).run().save(args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through python -m
    raise SystemExit(main())
