"""Shared result output for the three DESI Run-2 examples."""

import json
from pathlib import Path

import numpy as np

from fishhighz.adapters.legacy_compat import plain


def create_output(path):
    """Create and return a new output directory without replacing prior results."""
    output = Path(path).resolve()
    output.mkdir(parents=True, exist_ok=False)
    return output


def print_results(records):
    """Print individual and joint AP constraints in bin order."""
    print("bin  result                         sigma(ap)   sigma(at)   corr(ap,at)")
    for row in records:
        label = "joint" if row["pair"] is None else " x ".join(row["pair"])
        if row["status"] != "available":
            print(f"{row['bin_index'] + 1:>3}  {label:<29} {row['status']}")
            continue
        print(
            f"{row['bin_index'] + 1:>3}  {label:<29} "
            f"{row['sigma_ap']:>10.6g} {row['sigma_at']:>11.6g} "
            f"{row['correlation']:>13.6g}"
        )


def write_results(output, *, settings, records):
    """Save JSON identities/status and compact numeric NPZ forecast arrays."""
    output = Path(output).resolve()
    manifest = []
    for index, row in enumerate(records):
        manifest.append(
            {
                "array_index": index,
                "bin_index": row["bin_index"],
                "bounds": list(row["bounds"]),
                "redshift": row["redshift"],
                "kind": row["kind"],
                "pair": None if row["pair"] is None else list(row["pair"]),
                "parameter_ids": row.get(
                    "parameter_ids", settings.get("parameter_order", [])
                ),
                "status": row["status"],
                "sigma_ap": row["sigma_ap"],
                "sigma_at": row["sigma_at"],
                "correlation": row["correlation"],
            }
        )
    document = plain(settings)
    document["results"] = manifest
    (output / "settings.json").write_text(
        json.dumps(document, indent=2, allow_nan=False) + "\n"
    )
    np.savez_compressed(
        output / "results.npz",
        bin_index=np.asarray([row["bin_index"] for row in records], dtype=np.int64),
        bin_bounds=np.asarray([row["bounds"] for row in records], dtype=float),
        redshift=np.asarray([row["redshift"] for row in records], dtype=float),
        kind=np.asarray([row["kind"] == "joint" for row in records], dtype=np.int8),
        pair_indices=np.asarray(
            [
                (-1, -1) if row["pair_indices"] is None else row["pair_indices"]
                for row in records
            ],
            dtype=np.int64,
        ),
        available=np.asarray(
            [row["status"] == "available" for row in records], dtype=np.int8
        ),
        fisher=np.asarray([row["fisher"] for row in records], dtype=float),
        sigma_ap=np.asarray(
            [np.nan if row["sigma_ap"] is None else row["sigma_ap"] for row in records]
        ),
        sigma_at=np.asarray(
            [np.nan if row["sigma_at"] is None else row["sigma_at"] for row in records]
        ),
        correlation=np.asarray(
            [
                np.nan if row["correlation"] is None else row["correlation"]
                for row in records
            ]
        ),
    )
    return output / "settings.json", output / "results.npz"
