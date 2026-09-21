"""Replay stopping on saved stage-2 inputs and compare every remaining state."""

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from fishhighz.validation.adaptive_weights import METRICS, adaptive_weights, residuals
from fishhighz.validation.compatibility_weights import WeightInputs, fixed_weights


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--bin", type=int, required=True)
    args = parser.parse_args()
    path = args.input / f"bin-{args.bin}.json"
    data = json.loads(path.read_text())
    source = Path(data["source"])
    assert digest(source) == data["report_sha256"]
    report = json.loads(source.read_text())
    arrays_path = path.with_suffix(".npz")
    records = []
    with np.load(arrays_path) as arrays:
        for row in data["records"]:
            context, nodes, variant = row["context"], row["nodes"], row["variant"]
            prefix = f"{context}/{nodes}"
            original = WeightInputs.from_pair(
                report["settings"]["pair_inputs"][context]
            )
            m = arrays[prefix + "/magnitudes"]
            inputs = replace(
                original,
                magnitudes=m,
                density=arrays[prefix + "/density"],
                variance=arrays[prefix + "/variance"],
                quadrature=np.full(m.shape, m[1] - m[0]),
            )
            key = prefix + "/" + variant
            w, c = arrays[key + "/weights"], arrays[key + "/coefficients"]
            np.testing.assert_array_equal(fixed_weights(inputs, variant, 3), w[3])
            for tol in (1e-3, 1e-4):
                result = adaptive_weights(inputs, variant, context=context, rtol=tol)
                expected_eligible = row["auto"] and variant.startswith("sum_")
                if not expected_eligible:
                    assert result["status"] == "ineligible"
                    records.append(
                        dict(
                            context=context,
                            nodes=nodes,
                            variant=variant,
                            rtol=tol,
                            status="ineligible",
                            updates=0,
                        )
                    )
                    continue
                t = result["updates"]
                assert result["state_updates"] == t
                np.testing.assert_array_equal(result["weights"], w[t])
                np.testing.assert_array_equal(result["coefficients"], c[t])
                expected = row["results"][str(tol)]["status"]
                assert (result["status"] == "converged") == (
                    expected == "finite_nonzero_convergence"
                )
                tail = None
                if result["status"] == "converged":
                    assert t == 2 * result["candidate"]
                    assert t == row["results"][str(tol)]["confirmed_at"]
                    tail = np.max(
                        [residuals(w[j], w[t], c[j], c[t]) for j in range(t, len(w))],
                        axis=0,
                    )
                    assert np.all(tail <= tol), (
                        args.bin,
                        context,
                        nodes,
                        variant,
                        tail,
                    )
                else:
                    assert result["status"] == "capped" and t == 96
                records.append(
                    dict(
                        context=context,
                        nodes=nodes,
                        variant=variant,
                        rtol=tol,
                        status=result["status"],
                        updates=t,
                        candidate=result["candidate"],
                        last_step=result["last_step"],
                        forward_residual=result["forward_residual"],
                        coefficients=result["coefficients"].tolist(),
                        tail_max=None
                        if tail is None
                        else dict(zip(METRICS, tail.tolist(), strict=True)),
                    )
                )
    args.output.mkdir(parents=True, exist_ok=True)
    output = dict(
        bin=args.bin,
        source_sha256=digest(source),
        stage2_json_sha256=digest(path),
        stage2_npz_sha256=digest(arrays_path),
        records=records,
    )
    (args.output / path.name).write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n"
    )
    print(
        f"bin {args.bin}: {len(records)} stopping checks and 135 fixed-three checks passed"
    )


if __name__ == "__main__":
    main()
