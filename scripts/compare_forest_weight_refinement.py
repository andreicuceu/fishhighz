"""W06: nine saved-input coefficient comparisons; no trajectory or provider calls."""

import hashlib
import json
import os
import platform
import signal
import time
from datetime import datetime, timezone
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / ".validation/step12-r5-20260914T191855Z"
INPUTS = {
    "weight-diagnosis/diagnosis.json": "069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008",
    "weight-diagnosis/bin-0.npz": "d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73",
    "profiles-checked/records-001.report.json": "9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666",
}
NAMES = ["I1", "I2", "I3", "A", "P_pixel", "AB", "Q"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def agree(actual, expected):
    a, e = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(e))
    np.testing.assert_allclose(a, e, rtol=5e-12, atol=0)
    nz = e != 0
    return float(np.max(np.abs((a[nz] - e[nz]) / e[nz]), initial=0))


def moments(r, w, v, length, pixel, b):
    rw = r * w
    rww = rw * w
    rwwv = rww * v
    i1, i2, i3 = [np.cumsum(x)[-1] for x in (rw, rww, rwwv)]
    assert i1 != 0
    a = i2 / i1 / i1 / length
    p = i3 / i1 / i1 * pixel / length
    return np.array([i1, i2, i3, a, p, a * b, a * b + p])


def scalar_sums(r, w, v, length, pixel, b):
    # Independent scalar arithmetic, used with Fraction or 80-digit Decimal.
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


def control():
    f = Fraction
    q, rho, v = [f(1, 2), f(2)], [f(2), f(1, 2)], [f(1), f(4)]
    r = [x * y for x, y in zip(q, rho, strict=True)]
    assert r == [1, 1]
    nu = [1 / (1 + z) for z in v]
    assert nu == [f(1, 2), f(1, 5)]
    values, direct = scalar_sums(r, nu, v, f(1), f(1), f(1))
    neff = sum(x / (1 + z) for x, z in zip(r, v, strict=True))
    assert values[3:5] == [f(29, 49), f(41, 49)]
    assert values[-1] == direct == 1 / neff == f(10, 7)
    assert neff == f(7, 10)
    agree(
        moments(np.array(r, float), np.array(nu, float), np.array(v, float), 1, 1, 1),
        values,
    )
    return dict(
        q=list(map(str, q)),
        rho=list(map(str, rho)),
        values=list(map(str, values)),
        n_eff=str(neff),
    )


def timeout_handler(_signum, _frame):
    raise TimeoutError("30-second cap; no extension")


def compare(result):
    for name, expected in INPUTS.items():
        assert digest(BASE / name) == expected, name
    sources = [
        "scripts/compare_forest_weight_refinement.py",
        "scripts/compare_inverse_variance_forest_weights.py",
        "scripts/diagnose_desi2_weights.py",
        "scripts/diagnose_weight_limit.py",
        "fishhighz/kernels/weights.py",
        "IMPLEMENTATION_STEP.md",
        "WEIGHTING_DIAGNOSTIC_STEP.md",
        "../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md",
        "reviews/step-12-r5.md",
        "reviews/step-12-review-r5.md",
        "reviews/weighting-diagnostics-w05-r1.md",
        "reviews/weighting-diagnostics-w05-review-r1.md",
    ]
    result["source_sha256"] = {p: digest(ROOT / p) for p in sources}
    result["input_sha256"] = INPUTS
    diagnosis = json.loads((BASE / "weight-diagnosis/diagnosis.json").read_text())[
        "bins"
    ][0]
    report = json.loads((BASE / "profiles-checked/records-001.report.json").read_text())
    s, context = report["settings"], report["context"]
    assert diagnosis["bin"] == context["bin"] == 0
    assert context["profile"] == s["profile"] == "accuracy"
    assert context["bounds"] == s["bounds"] == [2.0, 2.235]
    assert context["fields"][0] == dict(
        id="lya(qso)", kind="forest", physical="lya", background="qso"
    )
    assert s["fields"][0] == "lya(qso)"
    assert s["controls"]["magnitude_order"] == 32
    assert diagnosis["arrays"] == "bin-0.npz"
    assert diagnosis["sha256"] == INPUTS["weight-diagnosis/bin-0.npz"]
    partition = np.array(diagnosis["partition"])
    assert np.array_equal(partition, s["partition"])
    assert np.all(np.diff(partition) > 0) and tuple(partition[[0, -1]]) == (16.1, 26.75)
    sample = s["samples"]["lya(qso)"]
    policies = {
        key: {
            k: v
            for k, v in sample[key + "_diagnostics"]["policies"].items()
            if not isinstance(v, (list, dict))
        }
        for key in ("density", "snr")
    }
    assert policies["density"]["negative"] == "floor_negative"
    assert policies["snr"]["snr"] == "legacy_floor_clamp"
    length, zsource, zeval = (
        sample["length_velocity"],
        sample["z_source"],
        s["model"]["z_eval"],
    )
    assert (length, zsource, zeval) == (
        44147.09373976799,
        2.3830099582078716,
        2.1152848986890427,
    )
    assert policies["snr"]["pixel_width_angstrom"] == 0.8
    pixel = (
        299792.458 * policies["snr"]["pixel_width_angstrom"] / (1215.67 * (1 + zeval))
    )
    result.update(
        policies=policies,
        partition=partition.tolist(),
        scalars=dict(L=length, l_p=pixel, z_source=zsource, z_eval=zeval),
        source_width_policy=s["source_width_policy"],
        control=control(),
    )
    result["coefficients"], result["domains"] = {}, {}
    float_values, decimal_values = {}, {}
    discrepancies = []
    getcontext().prec = 80
    dec = Decimal.from_float
    ld, pd = dec(length), dec(pixel)
    with np.load(BASE / "weight-diagnosis/bin-0.npz", allow_pickle=False) as saved:
        for order in (16, 32, 64):
            candidates = [row for row in diagnosis["rows"] if row["order"] == order]
            metadata = next(
                f for f in candidates[0]["fields"] if f["field"] == "lya(qso)"
            )
            assert candidates[0]["fields"][0]["field"] == "lya(qso)"
            assert all(
                next(f for f in row["fields"] if f["field"] == "lya(qso)") == metadata
                for row in candidates
            )
            b, p0 = metadata["alias"], metadata["signal"]
            assert (b, p0) == (16.356748967852234, 1.7182791088399005)
            result["scalars"].update(B=b, P0=p0)
            prefix = f"order{order}_field0_"
            m, q, r, v = [
                saved[prefix + key]
                for key in ("magnitudes", "measure", "masses", "variance")
            ]
            assert all(
                x.shape == (162 * order,) and np.all(np.isfinite(x))
                for x in (m, q, r, v)
            )
            assert np.all(np.diff(m) > 0) and m[0] > 16.1 and m[-1] < 26.75
            assert np.all(q > 0) and np.all(r >= 0) and sum(r) > 0 and np.all(v >= 0)
            # Each saved interval contains exactly n nodes and its full quadrature mass.
            for i, (lo, hi) in enumerate(
                zip(partition[:-1], partition[1:], strict=True)
            ):
                sl = slice(i * order, (i + 1) * order)
                assert np.all((m[sl] > lo) & (m[sl] < hi))
                agree(np.sum(q[sl]), hi - lo)
            agree(np.sum(q), partition[-1] - partition[0])
            if order == 32:
                for actual, expected in (
                    (m, sample["magnitudes"]),
                    (q, sample["quadrature"]),
                    (v, sample["variance"]),
                    (
                        r,
                        np.array(sample["density"])
                        * ((1 + zsource) / 299792.458)
                        * np.array(sample["quadrature"]),
                    ),
                ):
                    discrepancies.append(agree(actual, expected))
            d = b + pixel * v
            assert np.all(np.isfinite(d)) and np.all(d > 0)
            rd, vd, bd = list(map(dec, r)), list(map(dec, v)), dec(b)
            dd = [bd + pd * z for z in vd]
            kd = sum(x / z for x, z in zip(rd, dd, strict=True))
            k = np.sum(r / d)
            neff = length * b * k
            discrepancies.append(agree(neff, ld * bd * kd))
            result["domains"][order] = dict(
                nodes=len(m),
                mass_sum=float(sum(r)),
                zero_masses=int(sum(r == 0)),
                D_min=float(min(d)),
                D_max=float(max(d)),
                K=float(k),
                n_eff=float(neff),
            )
            for label in ("nu", "t3", "t6"):
                w = b / d if label == "nu" else saved[prefix + "weights_" + label[1:]]
                assert w.shape == r.shape and np.all(np.isfinite(w))
                wd = [bd / z for z in dd] if label == "nu" else list(map(dec, w))
                values = moments(r, w, v, length, pixel, b)
                independent, direct = scalar_sums(rd, wd, vd, ld, pd, bd)
                discrepancies.extend(
                    (
                        agree(values, independent),
                        agree(values[-1], direct),
                        agree(independent[-1], direct),
                    )
                )
                if label == "nu":
                    historical = metadata["fixed_coefficients"]
                    discrepancies.extend(
                        (
                            agree(values[-1], b / neff),
                            agree(values[-1], 1 / (length * k)),
                            agree(direct, bd / (ld * bd * kd)),
                        )
                    )
                else:
                    row = next(
                        x for x in metadata["rows"] if x["iterations"] == int(label[1:])
                    )
                    assert row["available"]
                    historical = [row["A"], row["P_pixel"]]
                    assert values[-1] >= float_values[order, "nu"][-1] * (1 - 5e-12)
                    assert independent[-1] >= decimal_values[order, "nu"][-1]
                discrepancies.append(agree(values[3:5], historical))
                float_values[order, label], decimal_values[order, label] = (
                    values,
                    independent,
                )
                result["coefficients"][f"{order}_{label}"] = dict(
                    zip(NAMES, map(float, values), strict=True)
                ) | dict(
                    decimal=list(map(str, independent)), direct_decimal_Q=str(direct)
                )
    result["contrasts"] = {}
    contrast_errors = []

    def contrast(key, numerator, denominator, indices):
        a, e = float_values[numerator], float_values[denominator]
        da, de = decimal_values[numerator], decimal_values[denominator]
        out = {}
        for i in indices:
            value, independent = a[i] / e[i] - 1, da[i] / de[i] - 1
            error = abs(float(value) - float(independent))
            assert error <= 5e-12
            contrast_errors.append(error)
            out[NAMES[i]] = float(value)
        result["contrasts"][key] = out
        return out

    result["stable_over_tested_orders"] = {}
    for label in ("nu", "t3", "t6"):
        changes = [
            contrast(f"{label}_{hi}/{lo}", (hi, label), (lo, label), (3, 4, 6))
            for lo, hi in ((16, 32), (32, 64), (16, 64))
        ]
        result["stable_over_tested_orders"][label] = all(
            abs(x) <= 1e-3 for row in changes for x in row.values()
        )
    for order in (16, 32, 64):
        contrast(f"{order}_t6/t3", (order, "t6"), (order, "t3"), (3, 4, 6))
        for label in ("t3", "t6"):
            contrast(f"{order}_{label}/nu", (order, label), (order, "nu"), (6,))
    result["max_relative_discrepancy"] = max(discrepancies)
    result["max_absolute_contrast_discrepancy"] = max(contrast_errors)
    for p, expected in result["source_sha256"].items():
        assert digest(ROOT / p) == expected
    for p, expected in INPUTS.items():
        assert digest(BASE / p) == expected
    result["stop_reason"] = (
        "Completed prescribed nine coefficient sets and control; stop for review."
    )


def main():
    start = time.monotonic()
    output = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w06-r1-%Y%m%dT%H%M%SZ")
    )
    output.mkdir(parents=True, exist_ok=False)
    result = dict(
        status="incomplete",
        python=platform.python_version(),
        numpy=np.__version__,
        host=platform.node(),
    )
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(30)
    try:
        assert all(
            os.environ.get(k) == "1"
            for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        )
        with np.errstate(all="raise"):
            compare(result)
        result["status"] = "PASS"
    except Exception as exc:
        result["stop_reason"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        signal.alarm(0)
        result["elapsed_seconds"] = time.monotonic() - start
        (output / "summary.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        print(output, result["status"], result["elapsed_seconds"], flush=True)


if __name__ == "__main__":
    main()
