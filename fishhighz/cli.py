"""Command-line entry point for the native FishHighz forecast facade."""

from __future__ import annotations

import argparse

from .public import Forecast


def main(argv=None):
    """Run one native INI forecast and save its result directory."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ini", nargs="?", default=None, help="native FishHighz INI")
    parser.add_argument("--output", required=True, help="new result directory")
    parser.add_argument(
        "--joint-only", action="store_true", help="skip individual spectra"
    )
    args = parser.parse_args(argv)
    forecast = Forecast(args.ini)
    result = forecast.run(individuals=False) if args.joint_only else forecast.run()
    result.save(args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through python -m
    raise SystemExit(main())
