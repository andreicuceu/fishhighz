"""W07: invert four saved QSO-forest Fisher matrices; no forecast imports."""

import hashlib
import json
import os
import platform
import signal
import time
from datetime import datetime, timezone
from decimal import Decimal as D
from decimal import getcontext
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / ".validation/step12-r5-20260914T191855Z"
HASHES = {
    "weight-diagnosis/diagnosis.json": "069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008",
    "weight-diagnosis/bin-0.npz": "d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73",
    "profiles-checked/records-001.report.json": "9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666",
    "profiles-checked/records-001.npz": "e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0",
}
SOURCES = [
    "scripts/compare_forest_weight_bao_errors.py",
    "scripts/diagnose_desi2_weights.py",
    "fishhighz/validation/numerics.py",
    "fishhighz/_information.py",
    "docs/archive/notes/WEIGHTING_DIAGNOSTIC_STEP.md",
    "docs/archive/notes/IMPLEMENTATION_STEP.md",
    "docs/archive/planning/FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md",
    "docs/archive/reviews/step-12-r5.md",
    "docs/archive/reviews/step-12-review-r5.md",
    "docs/archive/reviews/weighting-diagnostics-w06-r1.md",
    "docs/archive/reviews/weighting-diagnostics-w06-review-r1.md",
]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def check(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected, dtype=float)
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    assert actual.shape == expected.shape
    zero = expected == 0
    assert np.all(actual[zero] == 0)
    residual = np.max(abs((actual[~zero] - expected[~zero]) / expected[~zero]))
    assert residual <= 5e-12, residual
    return float(residual)


def invert(saved):
    f = saved.copy()
    assert f.shape == (2, 2) and np.all(np.isfinite(f))
    assert np.all(f.diagonal() > 0)
    scale = np.sqrt(f.diagonal())
    r = f / scale[:, None] / scale[None, :]
    eps = 64 * np.finfo(float).eps * 2
    assert np.all(abs(r - r.T) <= eps * np.maximum(1, np.maximum(abs(r), abs(r.T))))
    asymmetry = float(abs(f[0, 1] - f[1, 0]))
    if asymmetry:
        f[0, 1] += 0.5 * (f[1, 0] - f[0, 1])
        f[1, 0] = f[0, 1]
    r = (r + r.T) / 2
    eig = np.linalg.eigvalsh(r)
    tol = eps * max(1, max(abs(eig)))
    assert np.all(eig > tol)
    cov = np.linalg.solve(f, np.eye(2))
    errors = np.sqrt(cov.diagonal())
    corr = cov / np.outer(errors, errors)
    a, b, d = (D.from_float(float(x)) for x in (f[0, 0], f[0, 1], f[1, 1]))
    det = a * d - b * b
    assert det > 0
    dc = [[d / det, -b / det], [-b / det, a / det]]
    de = [(d / det).sqrt(), (a / det).sqrt()]
    dr = [[D(1), -b / (a * d).sqrt()], [-b / (a * d).sqrt(), D(1)]]
    discrepancy = {
        k: check(x, y)
        for k, x, y in (
            ("covariance", cov, dc),
            ("errors", errors, de),
            ("correlation", corr, dr),
        )
    }
    return dict(
        fisher=saved.tolist(),
        covariance=cov.tolist(),
        errors=errors.tolist(),
        correlation=corr.tolist(),
        normalized_eigenvalues=eig.tolist(),
        rank_tolerance=float(tol),
        rank=2,
        asymmetry=asymmetry,
        decimal_covariance=[[str(x) for x in row] for row in dc],
        decimal_errors=[str(x) for x in de],
        decimal_correlation=str(dr[0][1]),
        solve_decimal_relative=discrepancy,
    ), de


def run(result):
    getcontext().prec = 80
    np.seterr(all="raise")
    threads = {
        k: os.environ.get(k)
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    }
    assert set(threads.values()) == {"1"}
    result.update(
        threads=threads,
        python=platform.python_version(),
        numpy=np.__version__,
        host=platform.node(),
        input_root=str(INPUT),
        input_sha256=HASHES,
        source_sha256={s: digest(ROOT / s) for s in SOURCES},
    )
    assert all(digest(INPUT / p) == h for p, h in HASHES.items())
    diagnosis = json.loads((INPUT / "weight-diagnosis/diagnosis.json").read_text())
    report = json.loads(
        (INPUT / "profiles-checked/records-001.report.json").read_text()
    )
    c, s = report["context"], report["settings"]
    assert c["profile"] == s["profile"] == "accuracy" and c["bin"] == 0
    assert c["bounds"] == s["bounds"] == [2.0, 2.235]
    assert c["parameters"] == s["parameters"] == ["ap_0", "at_0"]
    assert [f["id"] for f in c["fields"]] == s["fields"]
    (fi,) = [i for i, f in enumerate(c["fields"]) if f["id"] == "lya(qso)"]
    assert c["fields"][fi] == dict(
        id="lya(qso)", kind="forest", physical="lya", background="qso"
    )
    (pair,) = [i for i, p in enumerate(c["selected_pairs"]) if p == [fi, fi]]
    assert fi == 0 and len(c["selected_pairs"]) == 15
    assert report["final_controls"]["magnitude_order"] == 32
    assert report["final_controls"]["iterations"] == 6
    bin0 = diagnosis["bins"][0]
    assert bin0["bin"] == 0 and bin0["arrays"] == "bin-0.npz"
    assert bin0["sha256"] == HASHES["weight-diagnosis/bin-0.npz"]
    result.update(
        context=c,
        pair_index=pair,
        final_controls=report["final_controls"],
        source_width_policy=s["source_width_policy"],
    )
    decimals = []
    result["matrices"] = []
    with (
        np.load(INPUT / "weight-diagnosis/bin-0.npz", allow_pickle=False) as data,
        np.load(
            INPUT / "profiles-checked/records-001.npz", allow_pickle=False
        ) as primary,
    ):
        assert np.array_equal(primary["selected_pairs"], c["selected_pairs"])
        for order, count in ((32, 0), (32, 3), (32, 6), (64, 0)):
            (row,) = [
                r
                for r in bin0["rows"]
                if (r["order"], r["iterations"]) == (order, count)
            ]
            assert row["available"] is True
            assert row["fixed_initial_weights"] is (count == 0)
            key = f"order{order}_iterations{count}_pair_fisher"
            assert row["pair_fisher_key"] == key
            assert data[key].shape == (15, 2, 2)
            assert data[key.replace("fisher", "rank")][pair] == 2
            assert np.array_equal(
                data[key.replace("fisher", "constrained")][pair], [1, 1]
            )
            f = data[key][pair]
            derived, dec = invert(f)
            derived.update(
                order=order,
                iterations=count,
                pair_fisher_key=key,
                available=True,
                fixed_initial_weights=count == 0,
            )
            derived["stored_relative"] = {
                name: check(derived[name], data[key.replace("fisher", name)][pair])
                for name in ("covariance", "errors", "correlation")
            }
            if (order, count) == (32, 6):
                assert primary["pair_fisher"].shape == (15, 2, 2)
                derived["primary_relative"] = check(f, primary["pair_fisher"][pair])
            result["matrices"].append(derived)
            decimals.append(dec)
    control, _ = invert(np.array([[4.0, 1.0], [1.0, 9.0]]))
    control["analytic_relative"] = {
        "covariance": check(
            control["covariance"], np.array([[9.0, -1.0], [-1.0, 4.0]]) / 35
        ),
        "errors": check(control["errors"], [3 / np.sqrt(35), 2 / np.sqrt(35)]),
        "correlation": check(control["correlation"], [[1, -1 / 6], [-1 / 6, 1]]),
    }
    assert not np.allclose(control["errors"], [1 / 2, 1 / 3], rtol=5e-12, atol=0)
    result["analytic_control"] = control
    result["contrasts"] = {}
    for label, num, den in (
        ("t3/nu32", 1, 0),
        ("t6/nu32", 2, 0),
        ("t6/t3", 2, 1),
        ("nu64/nu32", 3, 0),
    ):
        x = np.array(result["matrices"][num]["errors"])
        y = np.array(result["matrices"][den]["errors"])
        delta = x / y - 1
        independent = [a / b - 1 for a, b in zip(decimals[num], decimals[den])]
        discrepancy = float(max(abs(delta - np.array(independent, dtype=float))))
        assert discrepancy <= 5e-12
        result["contrasts"][label] = dict(
            fractional_error_change=delta.tolist(),
            decimal_change=[str(x) for x in independent],
            absolute_discrepancy=discrepancy,
        )
    f32, f64 = (np.array(result["matrices"][i]["fisher"]) for i in (0, 3))
    norm = float(np.linalg.norm(f64 - f32) / np.linalg.norm(f32))
    a, b = ([D.from_float(float(v)) for v in f.flat] for f in (f32, f64))
    dn = (sum((y - x) ** 2 for x, y in zip(a, b)) / sum(x * x for x in a)).sqrt()
    assert abs(norm - float(dn)) <= 5e-12
    delta = result["contrasts"]["nu64/nu32"]["fractional_error_change"]
    result["reference_stability"] = dict(
        fisher_relative=norm,
        decimal_fisher_relative=str(dn),
        absolute_discrepancy=abs(norm - float(dn)),
        threshold=1e-3,
        stable=bool(norm <= 1e-3 and max(map(abs, delta)) <= 1e-3),
    )
    assert all(digest(INPUT / p) == h for p, h in HASHES.items())
    assert all(digest(ROOT / p) == h for p, h in result["source_sha256"].items())
    result["status"] = "passed; awaiting independent/user review"


def stop(signum, frame):
    raise TimeoutError("W07 30-second cap reached; no extension")


if __name__ == "__main__":
    output = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w07-r1-%Y%m%dT%H%M%SZ")
    )
    output.mkdir(parents=True, exist_ok=False)
    result = {}
    start = time.monotonic()
    signal.signal(signal.SIGALRM, stop)
    signal.alarm(30)
    try:
        run(result)
    except BaseException as error:
        result["status"] = f"stopped: {type(error).__name__}: {error}"
        raise
    finally:
        signal.alarm(0)
        result["elapsed_seconds"] = time.monotonic() - start
        (output / "summary.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        print(output)
