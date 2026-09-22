#!/usr/bin/env python3
"""W04 r1: two bounded full-sample trajectories on saved signed inputs."""

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
W03 = ROOT / ".validation/forest-weight-diagnostics/w03-r1-20260915T224027Z"
REPORT_HASH = "3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6"
COMMAND = (
    "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
    ".venv/bin/python -B scripts/compare_forest_weight_signals.py"
)


def agree(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    np.testing.assert_allclose(actual, expected, rtol=5e-12, atol=0)
    mask = expected != 0
    return float(
        np.max(np.abs((actual[mask] - expected[mask]) / expected[mask]), initial=0)
    )


def moments(w, rho, v, dm, length, pixel):
    i1 = np.cumsum((rho * w) * dm)[-1]
    squared = np.square(w)
    i2 = np.cumsum((rho * squared) * dm)[-1]
    i3 = np.cumsum(((rho * squared) * v) * dm)[-1]
    denominator = length * np.square(i1)
    assert denominator != 0, "undefined moment normalization"
    values = np.array([i1, i2, i3, i2 / denominator, pixel * i3 / denominator])
    assert np.all(np.isfinite(values)), "nonfinite moments/coefficients"
    assert i1 * (length / pixel) > 0, "nonpositive effective density"
    assert np.all(values[-2:] >= 0), "negative coefficient"
    return values


def update(w, rho, v, dm, length, pixel, p0, b, refresh):
    values = moments(w, rho, v, dm, length, pixel)
    s = p0 + values[3] * b if refresh else p0
    assert np.isfinite(s) and s > 0, "nonpositive/nonfinite signal"
    n = values[0] * (length / pixel)
    denominator = s + v / n
    assert np.all(np.isfinite(denominator)) and np.all(denominator != 0)
    result = s / denominator
    assert np.all(np.isfinite(result)), "nonfinite updated weights"
    return result, s


def tiny_checks():
    # One bin: A=1, N=w; fixed f(w)=w/(1+w), refreshed f(w)=2w/(1+2w).
    for refresh, expected in ((False, 1 / 3), (True, 1 / 2)):
        w, s = update(
            np.array([0.5]), np.ones(1), np.ones(1), 1.0, 1.0, 1.0, 1.0, 1.0, refresh
        )
        agree(w, [expected])
        assert s == (2 if refresh else 1)
    q = Fraction
    w = [q(1, 2), q(1, 5)]
    first = [q(39, 74), q(39, 179)]
    # At first iterate, A1=(179^2+74^2)/253^2; S1=101526/64009.
    a1 = q(179**2 + 74**2, 253**2)
    s1 = 1 + a1
    n1 = sum(first)
    second = [s1 * n1 / (s1 * n1 + x) for x in (1, 4)]
    assert sum(x * x for x in w) / sum(w) ** 2 == q(29, 49)
    assert [(q(78, 49) * sum(w)) / (q(78, 49) * sum(w) + x) for x in (1, 4)] == first
    frozen = [q(78, 49) * n1 / (q(78, 49) * n1 + x) for x in (1, 4)]
    assert second != frozen
    actual = np.array([float(x) for x in w])
    for expected, sig in ((first, q(78, 49)), (second, s1)):
        actual, s = update(
            actual, np.ones(2), np.array([1.0, 4.0]), 1.0, 1.0, 1.0, 1.0, 1.0, True
        )
        agree(actual, [float(x) for x in expected])
        agree([s], [float(sig)])
    return {
        "one_bin": "passed",
        "two_bin": "passed",
        "second_weights": list(map(str, second)),
        "second_signal": str(s1),
        "second_A": str(a1),
    }


def decimal_trajectory(rho, v, dm, length, pixel, p0, b, count, refresh):
    """Independent scalar sums and recurrence; no float moments/update calls."""
    with localcontext() as ctx:
        ctx.prec = 80
        d = Decimal.from_float
        r, var = ([d(float(x)) for x in a] for a in (rho, v))
        q, ell, lp, intrinsic, power = map(d, (dm, length, pixel, p0, b))
        w = [(power / lp) / (power / lp + x) for x in var]
        states, signals = {}, []
        for t in range(count + 1):
            i1 = sum((rj * x * q for rj, x in zip(r, w, strict=True)), Decimal(0))
            i2 = sum((rj * x * x * q for rj, x in zip(r, w, strict=True)), Decimal(0))
            i3 = sum(
                (rj * x * x * y * q for rj, x, y in zip(r, w, var, strict=True)),
                Decimal(0),
            )
            assert i1 != 0, "Decimal undefined normalization"
            a, p = i2 / (ell * i1 * i1), lp * i3 / (ell * i1 * i1)
            n = (ell / lp) * i1
            s = intrinsic + a * power if refresh else intrinsic
            assert n > 0 and s > 0 and a >= 0 and p >= 0, (
                "Decimal physical-domain obstruction"
            )
            if t in (0, 3, 6, 12, 24):
                states[t] = (list(w), [i1, i2, i3, a, p])
            if t == count:
                break
            signals.append(s)
            denominators = [s + y / n for y in var]
            assert all(x != 0 for x in denominators), "Decimal zero update denominator"
            w = [s / x for x in denominators]
        return states, signals


def changes(current, reference):
    return [
        {
            "difference": float(x - y),
            "fraction": float((x - y) / y) if y != 0 else None,
            "unchanged": bool(x == y),
        }
        for x, y in zip(current, reference, strict=True)
    ]


def bounded_stability(states, branch, triple):
    a, b, c = triple
    comparisons = {
        f"{j}_vs_{i}": changes(
            states[f"{branch}_{j}"]["values"][-2:],
            states[f"{branch}_{i}"]["values"][-2:],
        )
        for i, j in ((a, b), (b, c), (a, c))
    }
    stable = all(
        abs(x["fraction"]) <= 1e-3 if x["fraction"] is not None else x["unchanged"]
        for pair in comparisons.values()
        for x in pair
    )
    return {"stable": stable, "comparisons": comparisons}


def timeout_handler(signum, frame):
    raise TimeoutError("30-second invocation cap")


def main():
    started = time.monotonic()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        assert os.environ.get(name) == "1"
    out = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w04-r1-%Y%m%dT%H%M%SZ")
    )
    out.mkdir(exist_ok=False)
    arrays, states = {}, {}
    summary = {
        "command": COMMAND,
        "states": states,
        "output": str(out.relative_to(ROOT)),
        "signals": {"fixed": [], "refreshed": []},
        "attempts": [],
        "stability": {},
    }
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(30)
    try:
        summary["analytic_checks"] = tiny_checks()
        assert hashlib.sha256(REPORT.read_bytes()).hexdigest() == REPORT_HASH
        report = json.loads(REPORT.read_text())
        context = report["context"]
        assert context["profile"] == "compatibility" and context["bin"] == 0
        assert context["bounds"] == [2.0, 2.235]
        assert report["settings"]["profile"] == "compatibility"
        assert "lya(qso)" in report["settings"]["fields"]
        row = report["settings"]["pair_inputs"]["lya(qso)_lya(qso)"]
        old = np.load(W03 / "arrays.npz")
        for key in ("magnitudes", "density", "variance"):
            arrays[key] = np.array(row[key], dtype=np.float64)
            assert arrays[key].shape == (107,) and np.all(np.isfinite(arrays[key]))
            np.testing.assert_array_equal(arrays[key], old[key])
            arrays[key].flags.writeable = False
        mag, rho, v = (arrays[k] for k in ("magnitudes", "density", "variance"))
        assert (mag[0], mag[-1]) == (16.1, 26.75)
        assert np.all(np.diff(mag) > 0) and np.count_nonzero(rho < 0) == 9
        dm = float(mag[1] - mag[0])
        length, pixel, p0, b = (
            float(row[k])
            for k in ("forest_length", "_pix_kms", "auxiliary_signal", "auxiliary_p1d")
        )
        assert np.all(np.isfinite([dm, length, pixel, p0, b]))
        args = (rho, v, dm, length, pixel, p0, b)
        summary["held_scalars"] = dict(
            zip(("dm", "L", "lp", "P0", "B"), (dm, length, pixel, p0, b), strict=True)
        )
        summary["identity"] = {
            "profile": "compatibility",
            "bin": 0,
            "bounds": context["bounds"],
            "pair": "lya(qso)_lya(qso)",
            "negative_density_indices": np.flatnonzero(rho < 0).tolist(),
        }
        with np.errstate(all="raise"):
            seed = (b / pixel) / (b / pixel + v)
            np.testing.assert_array_equal(seed, old["initial_weights"])
            arrays["seed"] = seed.copy()
            weights = {branch: seed.copy() for branch in ("fixed", "refreshed")}
            assert not np.shares_memory(*weights.values())
            previous = 0
            for count in (3, 6, 12, 24):
                for branch in weights:
                    refresh = branch == "refreshed"
                    for t in range(previous, count):
                        summary["attempts"].append(
                            {"branch": branch, "from_t": t, "to_t": t + 1}
                        )
                        w, s = update(weights[branch], *args, refresh)
                        # Check new state before accepting it; preserve previous finite state on failure.
                        moments(w, rho, v, dm, length, pixel)
                        weights[branch] = w
                        arrays[f"{branch}_last_finite_weights"] = w.copy()
                        summary[f"{branch}_last_finite_count"] = t + 1
                        summary["signals"][branch].append({"t": t, "S": float(s)})
                    values = moments(w, rho, v, dm, length, pixel)
                    label = f"{branch}_{count}"
                    arrays[label + "_weights"] = w.copy()
                    arrays[label + "_shape"] = w / np.max(np.abs(w))
                    arrays[label + "_moments_coefficients"] = values
                    states[label] = {
                        "values": values.tolist(),
                        "max_abs_w": float(np.max(np.abs(w))),
                        "negative_weights": int(np.count_nonzero(w < 0)),
                        "A_B_over_P0": float(values[3] * b / p0),
                    }
                    if branch == "fixed" and count in (3, 6):
                        agree(w, old[f"full_{count}_weights"])
                        agree(values, old[f"full_{count}_moments_coefficients"])
                        expected = {
                            3: [0.022021379679306282, 0.5298546380249203],
                            6: [0.021818865253617570, 0.5405425825641871],
                        }
                        agree(values[-2:], expected[count])
                if count >= 12:
                    triple = (3, 6, 12) if count == 12 else (6, 12, 24)
                    assessment = {
                        branch: bounded_stability(states, branch, triple)
                        for branch in weights
                    }
                    summary["stability"][str(triple)] = assessment
                    if all(x["stable"] for x in assessment.values()):
                        summary["stop"] = f"both branches stable over {triple}"
                        break
                previous = count
            else:
                summary["stop"] = "24-update cap; not both branches stable"
            summary["decimal_checks"] = {}
            for branch in weights:
                independent, signals = decimal_trajectory(
                    *args, count, branch == "refreshed"
                )
                errors = []
                for t, (dw, dv) in independent.items():
                    if t == 0:
                        errors.append(agree(seed, [float(x) for x in dw]))
                        continue
                    label = f"{branch}_{t}"
                    errors.extend(
                        [
                            agree(arrays[label + "_weights"], [float(x) for x in dw]),
                            agree(states[label]["values"], [float(x) for x in dv]),
                        ]
                    )
                    states[label]["decimal_values"] = list(map(str, dv))
                    arrays[label + "_decimal_weights"] = np.array(
                        [float(x) for x in dw]
                    )
                errors.append(
                    agree(
                        [x["S"] for x in summary["signals"][branch]],
                        [float(x) for x in signals],
                    )
                )
                summary["decimal_checks"][branch] = {
                    "max_relative_error": max(errors),
                    "signals_80_digit": list(map(str, signals)),
                }
            summary["matched_count_changes"] = {
                str(t): changes(
                    states[f"refreshed_{t}"]["values"][-2:],
                    states[f"fixed_{t}"]["values"][-2:],
                )
                for t in (3, 6, 12, 24)
                if f"fixed_{t}" in states
            }
            summary["checks"] = (
                "PASS: analytic controls, W03 replay, independent 80-digit checkpoints and signals"
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
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                REPORT,
                W03 / "arrays.npz",
                W03 / "summary.json",
                Path(__file__).resolve(),
                ROOT / "docs/archive/notes/WEIGHTING_DIAGNOSTIC_STEP.md",
                ROOT / "docs/archive/notes/IMPLEMENTATION_STEP.md",
                ROOT / "scripts/compare_forest_weight_integrals.py",
            )
        }
        np.savez(out / "arrays.npz", **arrays)
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\n"
        )
        print(out.relative_to(ROOT))
        print(summary.get("stop", summary.get("obstruction")))


if __name__ == "__main__":
    main()
