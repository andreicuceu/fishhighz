#!/usr/bin/env python3
"""W05 r1: four saved-input coefficient comparisons, with independent arithmetic."""

import hashlib
import json
import os
import platform
import signal
import time
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / ".validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json"
)
W04 = ROOT / ".validation/forest-weight-diagnostics/w04-r1-20260915T225726Z"
COMMAND = (
    "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
    ".venv/bin/python -B scripts/compare_inverse_variance_forest_weights.py"
)
NAMES = ("I1", "I2", "I3", "A", "P_pixel", "AB", "Q")


def agree(actual, expected):
    actual, expected = (
        np.asarray(actual, dtype=float),
        np.asarray(expected, dtype=float),
    )
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    np.testing.assert_allclose(actual, expected, rtol=5e-12, atol=0)
    mask = expected != 0
    return float(
        np.max(np.abs((actual[mask] - expected[mask]) / expected[mask]), initial=0)
    )


def moments(w, rho, v, dm, length, pixel, b):
    i1 = np.cumsum((rho * w) * dm)[-1]
    i2 = np.cumsum((rho * w**2) * dm)[-1]
    i3 = np.cumsum(((rho * w**2) * v) * dm)[-1]
    assert i1 != 0, "undefined normalization"
    denominator = length * i1**2
    a, p = i2 / denominator, pixel * i3 / denominator
    return np.array([i1, i2, i3, a, p, a * b, a * b + p])


def scalar_coefficients(w, r, v, length, pixel, b):
    # Independent direct scalar sums; no float moment helper or cumulative sum.
    i1 = sum(x * y for x, y in zip(r, w, strict=True))
    i2 = sum(x * y * y for x, y in zip(r, w, strict=True))
    i3 = sum(x * y * y * z for x, y, z in zip(r, w, v, strict=True))
    denominator = length * i1 * i1
    a, p = i2 / denominator, pixel * i3 / denominator
    direct = (
        sum(x * (b + pixel * z) * y * y for x, y, z in zip(r, w, v, strict=True))
        / denominator
    )
    return [i1, i2, i3, a, p, a * b, a * b + p], direct


def rational_control():
    f = Fraction
    v, r = [f(1), f(4)], [f(1), f(1)]
    nu = [1 / (1 + x) for x in v]
    assert nu == [f(1, 2), f(1, 5)]
    k = sum(x / (1 + z) for x, z in zip(r, v, strict=True))
    assert k == sum(nu) == f(7, 10)
    values, direct = scalar_coefficients(nu, r, v, f(1), f(1), f(1))
    assert values[3:5] == [f(29, 49), f(41, 49)]
    assert values[-1] == direct == 1 / k == f(10, 7)
    scaled, _ = scalar_coefficients([3 * x for x in nu], r, v, f(1), f(1), f(1))
    uniform, _ = scalar_coefficients(r, r, v, f(1), f(1), f(1))
    assert scaled[3:] == values[3:]
    assert uniform[-1] == f(7, 4) > direct
    for w, expected in ((nu, values), ([3 * x for x in nu], scaled), (r, uniform)):
        agree(
            moments(
                np.array(w, dtype=float),
                np.ones(2),
                np.array(v, dtype=float),
                1.0,
                1.0,
                1.0,
                1.0,
            ),
            expected,
        )
    return {
        "nu": list(map(str, nu)),
        "K": str(k),
        "n_eff": str(k),
        "values": list(map(str, values)),
        "uniform_Q": str(uniform[-1]),
        "normalization_by_3": "PASS",
    }


def timeout_handler(_signum, _frame):
    raise TimeoutError("30-second cap; no extension")


def main():
    assert all(
        os.environ.get(k) == "1"
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    )
    started = time.monotonic()
    out = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / ("w05-r1-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    )
    out.mkdir(exist_ok=False)
    arrays, result = {}, {"command": COMMAND, "states": {}}
    sources = [
        REPORT,
        W04 / "arrays.npz",
        W04 / "summary.json",
        Path(__file__).resolve(),
        ROOT / "WEIGHTING_DIAGNOSTIC_STEP.md",
        ROOT / "IMPLEMENTATION_STEP.md",
        ROOT / "scripts/compare_forest_weight_signals.py",
        ROOT / "fishhighz/kernels/weights.py",
        ROOT.parent / "lyaforecast/lyaforecast/weights.py",
        ROOT.parent / "lyaforecast/lyaforecast/covariance.py",
    ]
    hashes = {
        os.path.relpath(p, ROOT): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sources
    }
    result["source_sha256"] = hashes
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(30)
    try:
        result["rational_control"] = rational_control()
        assert hashes[os.path.relpath(REPORT, ROOT)] == (
            "3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6"
        )
        report = json.loads(REPORT.read_text())
        context = report["context"]
        assert context["profile"] == "compatibility" and context["bin"] == 0
        assert context["bounds"] == [2.0, 2.235]
        assert report["settings"]["profile"] == "compatibility"
        row = report["settings"]["pair_inputs"]["lya(qso)_lya(qso)"]
        old = np.load(W04 / "arrays.npz", allow_pickle=False)
        history = json.loads((W04 / "summary.json").read_text())
        for key in ("magnitudes", "density", "variance"):
            arrays[key] = np.array(row[key], dtype=np.float64)
            assert arrays[key].shape == (107,) and np.all(np.isfinite(arrays[key]))
            np.testing.assert_array_equal(arrays[key], old[key])
        mag, rho, v = (arrays[k] for k in ("magnitudes", "density", "variance"))
        assert (mag[0], mag[-1]) == (16.1, 26.75)
        assert np.all(np.diff(mag) > 0) and np.count_nonzero(rho < 0) == 9
        dm = float(mag[1] - mag[0])
        length, pixel, p0, b = (
            float(row[k])
            for k in ("forest_length", "_pix_kms", "auxiliary_signal", "auxiliary_p1d")
        )
        held = dict(
            zip(("dm", "L", "lp", "P0", "B"), (dm, length, pixel, p0, b), strict=True)
        )
        assert held == history["held_scalars"]
        assert (
            np.all(np.isfinite(list(held.values()))) and min(dm, length, pixel, b) > 0
        )
        assert np.all(v >= 0)
        result.update(
            held_scalars=held,
            identity={
                "profile": "compatibility",
                "bin": 0,
                "bounds": context["bounds"],
                "pair": "lya(qso)_lya(qso)",
                "negative_density_indices": np.flatnonzero(rho < 0).tolist(),
            },
        )
        with np.errstate(all="raise"), localcontext() as ctx:
            ctx.prec = 80
            ddm, dl, dp, db = map(Decimal.from_float, (dm, length, pixel, b))
            dr = [Decimal.from_float(float(x)) * ddm for x in rho]
            dv = [Decimal.from_float(float(x)) for x in v]
            dden = [db + dp * x for x in dv]
            denominator = b + pixel * v
            result["D_range"] = [float(denominator.min()), float(denominator.max())]
            assert np.all(denominator > 0) and all(x > 0 for x in dden), "nonpositive D"
            nu = b / denominator
            seed = (b / pixel) / (b / pixel + v)
            np.testing.assert_array_equal(seed, old["seed"])
            result["seed_relative_error"] = agree(nu, seed)
            dnu = [db / x for x in dden]
            result["decimal_nu_relative_error"] = agree(nu, dnu)
            k = np.cumsum((rho / denominator) * dm)[-1]
            dk = sum(x / y for x, y in zip(dr, dden, strict=True))
            result["K"] = float(k)
            result["decimal_K"] = str(dk)
            agree(k, dk)
            assert k > 0 and dk > 0, "nonpositive signed K; algebraic domain stop"
            result["formal_n_eff"] = float(length * b * k)
            result["K_cancellation_ratio"] = float(
                sum(abs(x / y) for x, y in zip(dr, dden, strict=True)) / abs(dk)
            )
            weights = {
                "nu": nu,
                "prefix_3": np.array(row["_w_lya"]),
                "fixed_24": old["fixed_24_weights"],
                "refreshed_24": old["refreshed_24_weights"],
            }
            decimal_values = {}
            for name, w in weights.items():
                assert w.shape == (107,) and np.all(np.isfinite(w))
                assert np.max(np.abs(w)) > 0, "undefined shape"
                arrays[name + "_weights"] = w
                arrays[name + "_shape"] = w / np.max(np.abs(w))
                values = moments(w, rho, v, dm, length, pixel, b)
                dw = [Decimal.from_float(float(x)) for x in w]
                dec, direct = scalar_coefficients(dw, dr, dv, dl, dp, db)
                error = agree(values, dec)
                agree(values[-1], direct)
                agree(dec[-1], direct)
                decimal_values[name] = dec
                state = dict(zip(NAMES, map(float, values), strict=True))
                state.update(
                    decimal_values=list(map(str, dec)),
                    decimal_direct_Q=str(direct),
                    max_relative_error=error,
                    negative_weights=int(np.count_nonzero(w < 0)),
                    shape_difference=float(
                        np.max(
                            np.abs(arrays[name + "_shape"] - nu / np.max(np.abs(nu)))
                        )
                    ),
                    moment_cancellation_ratios=[
                        float(sum(abs(x) for x in terms) / abs(total))
                        for terms, total in (
                            ([x * y for x, y in zip(dr, dw, strict=True)], dec[0]),
                            ([x * y * y for x, y in zip(dr, dw, strict=True)], dec[1]),
                            (
                                [
                                    x * y * y * z
                                    for x, y, z in zip(dr, dw, dv, strict=True)
                                ],
                                dec[2],
                            ),
                        )
                    ],
                )
                result["states"][name] = state
                assert np.all(np.isfinite(values)) and all(x.is_finite() for x in dec)
                assert min(values[3:]) > 0 and min(dec[3:]) > 0, (
                    "nonpositive signed coefficient"
                )
                if name == "prefix_3":
                    agree(
                        values[3:5],
                        [
                            row["_aliasing_weights"][-1],
                            row["_effective_noise_power"][-1],
                        ],
                    )
                elif name != "nu":
                    agree(values[:5], old[name + "_moments_coefficients"])
                    agree(values[:5], history["states"][name]["values"])
                else:
                    reference, reference_direct = scalar_coefficients(
                        dnu, dr, dv, dl, dp, db
                    )
                    agree(values, reference)
                    agree(values[-1], 1 / (length * k))
                    agree(values[-1], b / result["formal_n_eff"])
                    agree(reference_direct, 1 / (dl * dk))
                    agree(
                        reference[-1],
                        db / (dl * sum(x * y for x, y in zip(dr, dnu, strict=True))),
                    )
                    state["recomputed_decimal_nu_values"] = list(map(str, reference))
                    state["decimal_Q_from_K"] = str(1 / (dl * dk))
            for name, state in result["states"].items():
                contrasts = {}
                for key, index in (("A", 3), ("P_pixel", 4), ("Q", 6)):
                    contrast = state[key] / result["states"]["nu"][key] - 1
                    independent = (
                        decimal_values[name][index] / decimal_values["nu"][index] - 1
                    )
                    assert abs(contrast - float(independent)) <= 5e-12
                    contrasts[key] = contrast
                state["fractional_change_from_nu"] = contrasts
            contrast = (
                result["states"]["refreshed_24"]["Q"]
                / result["states"]["fixed_24"]["Q"]
                - 1
            )
            independent = (
                decimal_values["refreshed_24"][-1] / decimal_values["fixed_24"][-1] - 1
            )
            assert abs(contrast - float(independent)) <= 5e-12
            result["refreshed_vs_fixed_Q_fraction"] = contrast
        result["checks"] = (
            "PASS: exact control, normalization, seed, historical coefficients, 80-digit sums and contrasts"
        )
    except Exception as error:
        result["obstruction"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        signal.alarm(0)
        result.update(
            seconds=time.monotonic() - started,
            python=platform.python_version(),
            numpy=np.__version__,
        )
        np.savez(out / "arrays.npz", **arrays)
        (out / "summary.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        print(out.relative_to(ROOT))
        print(result.get("checks", result.get("obstruction")))
        assert hashes == {
            os.path.relpath(p, ROOT): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sources
        }, "input/source changed during check"


if __name__ == "__main__":
    main()
