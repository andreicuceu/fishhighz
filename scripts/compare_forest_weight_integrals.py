#!/usr/bin/env python3
"""W03 r1: bounded, saved-input prefix/full-sample coefficient comparison."""

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
REPORT_HASH = "3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6"
RTOL = 5e-12


def agree(actual, expected):
    """Require finite, elementwise relative agreement, including exact zeros."""
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=0)
    nonzero = expected != 0
    return float(
        np.max(
            np.abs((actual[nonzero] - expected[nonzero]) / expected[nonzero]), initial=0
        )
    )


def update(w, rho, v, dm, length, pixel, s, full):
    """Change only selection of the feedback integral; retain literal order."""
    integral = np.cumsum((rho * w) * dm)
    feedback = integral[-1] if full else integral
    density = feedback * (length / pixel)
    assert np.all(density != 0), "zero effective density"
    denominator = s + v / density
    assert np.all(np.isfinite(denominator)) and np.all(denominator != 0)
    result = s / denominator
    assert np.all(np.isfinite(result))
    return result


def moments(w, rho, v, dm, length, pixel):
    """Use full-sample final cumulative elements in W01 product order."""
    i1 = np.cumsum((rho * w) * dm)[-1]
    squared = np.square(w)
    i2 = np.cumsum((rho * squared) * dm)[-1]
    i3 = np.cumsum(((rho * squared) * v) * dm)[-1]
    denominator = length * np.square(i1)
    assert denominator != 0, "undefined coefficient normalization"
    result = np.array([i1, i2, i3, i2 / denominator, pixel * i3 / denominator])
    assert np.all(np.isfinite(result))
    return result


def decimal_state(rho, v, dm, length, pixel, s, b, count, full):
    """Independent scalar recurrence and direct moment sums at 80 digits."""
    with localcontext() as ctx:
        ctx.prec = 80
        d = Decimal.from_float
        r, variance = [d(float(x)) for x in rho], [d(float(x)) for x in v]
        q, ell, lp, sig, power = map(d, (dm, length, pixel, s, b))
        w = [(power / lp) / (power / lp + x) for x in variance]
        for _ in range(count):
            total = sum((r[j] * w[j] * q for j in range(len(w))), Decimal(0))
            prefix = Decimal(0)
            new = []
            for j in range(len(w)):
                prefix += r[j] * w[j] * q
                n = (ell / lp) * (total if full else prefix)
                new.append(sig / (sig + variance[j] / n))
            w = new
        i1 = sum((r[j] * w[j] * q for j in range(len(w))), Decimal(0))
        i2 = sum((r[j] * w[j] ** 2 * q for j in range(len(w))), Decimal(0))
        i3 = sum(
            (r[j] * w[j] ** 2 * variance[j] * q for j in range(len(w))), Decimal(0)
        )
        values = [i1, i2, i3, i2 / (ell * i1**2), lp * i3 / (ell * i1**2)]
        return np.array([float(x) for x in w]), values


def changes(current, reference):
    """Signed fractional changes, with explicit zero-reference behavior."""
    return [
        {"difference": float(x - y), "fraction": float((x - y) / y) if y != 0 else None}
        for x, y in zip(current, reference, strict=True)
    ]


def tiny_check():
    q = Fraction
    w = [q(1, 2), q(1, 5)]
    prefix = [sum(w[: i + 1]) / (sum(w[: i + 1]) + v) for i, v in enumerate((1, 4))]
    full = [sum(w) / (sum(w) + v) for v in (1, 4)]
    assert prefix == [q(1, 3), q(7, 47)]
    assert full == [q(7, 17), q(7, 47)] and prefix[-1] == full[-1]
    for shared, expected in ((False, prefix), (True, full)):
        agree(
            update(
                np.array([0.5, 0.2]),
                np.ones(2),
                np.array([1.0, 4.0]),
                1.0,
                1.0,
                1.0,
                1.0,
                shared,
            ),
            [float(x) for x in expected],
        )


def main():
    started = time.monotonic()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        assert os.environ.get(name) == "1", name
    out = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w03-r1-%Y%m%dT%H%M%SZ")
    )
    out.mkdir(exist_ok=False)
    arrays, states = {}, {}
    summary = {"output": str(out.relative_to(ROOT)), "states": states}
    try:
        signal.alarm(30)
        tiny_check()
        assert hashlib.sha256(REPORT.read_bytes()).hexdigest() == REPORT_HASH
        report = json.loads(REPORT.read_text())
        context = report["context"]
        assert context["profile"] == "compatibility" and context["bin"] == 0
        assert context["bounds"] == [2.0, 2.235]
        assert report["settings"]["profile"] == "compatibility"
        assert "lya(qso)" in report["settings"]["fields"]
        row = report["settings"]["pair_inputs"]["lya(qso)_lya(qso)"]
        for key in (
            "magnitudes",
            "density",
            "variance",
            "_w_lya",
            "_aliasing_weights",
            "_effective_noise_power",
        ):
            arrays[key] = np.array(row[key], dtype=np.float64)
            assert arrays[key].shape == (107,) and np.all(np.isfinite(arrays[key]))
            arrays[key].flags.writeable = False
        mag, rho, v = (arrays[k] for k in ("magnitudes", "density", "variance"))
        assert (mag[0], mag[-1]) == (16.1, 26.75)
        assert np.all(np.diff(mag) > 0) and np.count_nonzero(rho < 0) == 9
        dm = float(mag[1] - mag[0])
        length, pixel, s, b = (
            float(row[k])
            for k in ("forest_length", "_pix_kms", "auxiliary_signal", "auxiliary_p1d")
        )
        assert np.all(np.isfinite([dm, length, pixel, s, b]))
        arrays.update(
            dm=dm,
            forest_length=length,
            pixel_width=pixel,
            auxiliary_signal=s,
            auxiliary_p1d=b,
        )
        summary.update(
            report=str(REPORT.relative_to(ROOT)),
            report_sha256=REPORT_HASH,
            identity={
                "profile": "compatibility",
                "bin": 0,
                "bounds": context["bounds"],
                "population": "lya(qso)",
            },
            scalars={
                k: float(arrays[k])
                for k in (
                    "dm",
                    "forest_length",
                    "pixel_width",
                    "auxiliary_signal",
                    "auxiliary_p1d",
                )
            },
            negative_density_nodes=np.flatnonzero(rho < 0).tolist(),
        )
        with np.errstate(all="raise"):
            seed = (b / pixel) / (b / pixel + v)
            arrays["initial_weights"] = seed.copy()
            weights = {"prefix": seed.copy(), "full": seed.copy()}
            assert not np.shares_memory(weights["prefix"], weights["full"])
            args = (rho, v, dm, length, pixel, s)

            def advance(branch, start, stop):
                for count in range(start + 1, stop + 1):
                    weights[branch] = update(weights[branch], *args, branch == "full")
                    arrays[f"{branch}_last_finite_weights"] = weights[branch].copy()
                    summary[f"{branch}_last_finite_count"] = count
                w = weights[branch]
                values = moments(w, rho, v, dm, length, pixel)
                dw, dv = decimal_state(*args, b, stop, branch == "full")
                error_w = agree(w, dw)
                error_m = agree(values, [float(x) for x in dv])
                label = f"{branch}_{stop}"
                amplitude = float(np.max(np.abs(w)))
                assert amplitude != 0, "undefined all-zero shape"
                arrays[label + "_weights"] = w.copy()
                arrays[label + "_shape"] = w / amplitude
                arrays[label + "_moments_coefficients"] = values
                states[label] = dict(
                    zip(
                        ("I1", "I2", "I3", "A", "P_pixel"), values.tolist(), strict=True
                    )
                )
                states[label].update(
                    amplitude=amplitude,
                    negative_weights=int(np.count_nonzero(w < 0)),
                    decimal_weights_max_relative=error_w,
                    decimal_moments_coefficients_max_relative=error_m,
                    decimal_values=[str(x) for x in dv],
                )
                return values[-2:]

            p3 = advance("prefix", 0, 3)
            agree(p3, [0.022118395206768268, 0.5480269062497443])
            agree(weights["prefix"], arrays["_w_lya"])
            agree(
                p3,
                [arrays["_aliasing_weights"][-1], arrays["_effective_noise_power"][-1]],
            )
            p6 = advance("prefix", 3, 6)
            agree(p6, [0.022308973469159228, 0.5544833486755246])
            summary["prefix_replay"] = "passed"
            f3 = advance("full", 0, 3)
            f6 = advance("full", 3, 6)
            summary["changes"] = {
                "prefix_6_vs_3": changes(p6, p3),
                "full_6_vs_3": changes(f6, f3),
                "full_vs_prefix_3": changes(f3, p3),
                "full_vs_prefix_6": changes(f6, p6),
            }

            def sensitive(items):
                return any(
                    abs(x["fraction"]) > 1e-3
                    if x["fraction"] is not None
                    else x["difference"] != 0
                    for x in items
                )

            if sensitive(summary["changes"]["full_6_vs_3"]):
                summary["stop"] = (
                    "confirmed full-sample sensitivity at t=3/6; t=12 not attempted"
                )
            else:
                f12 = advance("full", 6, 12)
                summary["changes"]["full_12_vs_6"] = changes(f12, f6)
                summary["changes"]["full_12_vs_3"] = changes(f12, f3)
                renewed = any(
                    sensitive(summary["changes"][k])
                    for k in ("full_12_vs_6", "full_12_vs_3")
                )
                summary["stop"] = (
                    "renewed sensitivity at t=12"
                    if renewed
                    else "bounded coefficient stability at t=3/6/12"
                )
            summary["checks"] = (
                "PASS: exact two-bin, prefix replay, independent 80-digit Decimal at every checkpoint"
            )
    except Exception as error:
        summary["obstruction"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        signal.alarm(0)
        summary.update(
            seconds=time.monotonic() - started,
            python=platform.python_version(),
            numpy=np.__version__,
        )
        summary["source_sha256"] = {
            str(p.relative_to(ROOT))
            if p.is_relative_to(ROOT)
            else str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                Path(__file__).resolve(),
                ROOT / "WEIGHTING_DIAGNOSTIC_STEP.md",
                ROOT / "scripts/diagnose_legacy_forest_iterations.py",
                ROOT.parent / "lyaforecast/lyaforecast/weights.py",
            )
        }
        np.savez(out / "arrays.npz", **arrays)
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\n"
        )
        print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
