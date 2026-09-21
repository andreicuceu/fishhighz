"""Render saved three-profile DESI-2 evidence without running forecasts."""

import argparse
import json
from pathlib import Path

import numpy as np

from fishhighz.validation.evidence import digest, record_report
from fishhighz.validation.schema import validate_payload
from fishhighz.validation.three_profile_plots import build_tables, render


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    records = []
    for bundle in args.bundle:
        root = Path(bundle)
        manifest = json.loads((root / "manifest.json").read_text())
        for record in manifest["records"]:
            arrays = None
            if "arrays" in record:
                path = root / record["arrays"]
                if digest(path) != record["sha256"]:
                    raise ValueError("saved array hash mismatch")
                record = {**record, "report": record_report(root, record)}
                with np.load(path, allow_pickle=False) as saved:
                    arrays = dict(saved)
                validate_payload(
                    record["task"],
                    arrays,
                    record["report"],
                    require_pass=False,
                    schema=manifest["schema"],
                )
            records.append((record, arrays))
    render(build_tables(records), args.output)


if __name__ == "__main__":
    main()
