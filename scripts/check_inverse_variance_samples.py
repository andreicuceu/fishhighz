"""W10 r1: two ordered saved-input coefficient and 2x2 BAO checks."""

import hashlib
import json
import math
import os
import platform
import signal
import time
from datetime import datetime, timezone
from decimal import Decimal as D
from pathlib import Path

import numpy as np

from fishhighz.fields import ObservedField
from fishhighz.geometry import prepare_geometry
from fishhighz.response import InstrumentResponse
from fishhighz.validation._recipes import RECIPES
from fishhighz.weights import prepare_forest_weights

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / ".validation/step12-r5-20260914T191855Z"
HASHES = {
    "weight-diagnosis/diagnosis.json": "069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008",
    "weight-diagnosis/bin-0.npz": "d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73",
    "weight-diagnosis/bin-5.npz": "0246bc1093e98e7ac48e859e1e86a575f10a0bf6afbbd29908c791c35558fe13",
    "profiles-checked/records-001.report.json": "9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666",
    "profiles-checked/records-011.report.json": "f2585acd91a965d90dd52a97887b3365fe0dbb25c64eacb187f181600506f09f",
}
SOURCES = [
    "scripts/check_inverse_variance_samples.py",
    "docs/archive/notes/WEIGHTING_DIAGNOSTIC_STEP.md",
    "fishhighz/weights.py",
    "fishhighz/kernels/weights.py",
    "fishhighz/_arrays.py",
    "fishhighz/geometry.py",
    "fishhighz/response.py",
    "fishhighz/fields.py",
    "fishhighz/validation/_recipes.py",
    "scripts/diagnose_desi2_weights.py",
    "scripts/compare_forest_weight_refinement.py",
    "scripts/compare_forest_weight_bao_errors.py",
    "scripts/replay_inverse_variance_bao.py",
]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def check(out, name, actual, expected):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    assert a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), (
        name
    )
    nz = b != 0
    assert np.all(a[~nz] == 0), name
    residual = float(np.max(abs((a[nz] - b[nz]) / b[nz]), initial=0))
    out.setdefault("residuals", {})[name] = residual
    assert residual <= 5e-12, (name, residual)


def contrast(numerator, denominator):
    values = [
        float(a) / float(b) - 1 for a, b in zip(numerator, denominator, strict=True)
    ]
    independent = [
        float(D.from_float(float(a)) / D.from_float(float(b)) - 1)
        for a, b in zip(numerator, denominator, strict=True)
    ]
    residual = max(abs(a - b) for a, b in zip(values, independent, strict=True))
    assert residual <= 5e-12
    return dict(
        numerator=list(numerator),
        denominator=list(denominator),
        fractional=values,
        independent_fractional=independent,
        absolute_residual=residual,
    )


def invert(data, row, pair, out):
    assert row["available"]
    key = row["pair_fisher_key"]
    assert key == f"order{row['order']}_iterations{row['iterations']}_pair_fisher"
    assert row["fixed_initial_weights"] == (row["iterations"] == 0)
    f = data[key][pair]
    assert f.shape == (2, 2) and np.all(np.isfinite(f)) and np.array_equal(f, f.T)
    a, b, d = map(float, (f[0, 0], f[0, 1], f[1, 1]))
    det = a * d - b * b
    assert a > 0 and d > 0 and det > 0
    eig = 1 - abs(b) / math.sqrt(a * d)
    assert eig > 128 * np.finfo(float).eps * 2
    errors = [math.sqrt(d / det), math.sqrt(a / det)]
    out.update(
        pair_fisher_key=key,
        matrix=f.tolist(),
        determinant=det,
        normalized_min_eigenvalue=eig,
        errors=errors,
    )
    check(
        out, "solve_errors", errors, np.sqrt(np.linalg.solve(f, np.eye(2)).diagonal())
    )
    check(out, "saved_errors", errors, data[key.replace("fisher", "errors")][pair])
    return errors


def sample(
    result, diagnosis, bin_id, label, background, native, record, bounds, intervals
):
    out = dict(
        bin=bin_id,
        field=label,
        native_order=native,
        coefficients={},
        refinement={},
        matrices={},
    )
    result["samples"].append(out)
    result["stage"] = f"{label} bin {bin_id}: identities"
    report = json.loads(
        (BASE / f"profiles-checked/records-{record:03d}.report.json").read_text()
    )
    c, s = report["context"], report["settings"]
    assert c["bin"] == bin_id and c["profile"] == s["profile"] == "accuracy"
    assert c["bounds"] == s["bounds"] == bounds
    assert c["parameters"] == s["parameters"] == [f"ap_{bin_id}", f"at_{bin_id}"]
    assert (
        s["controls"]["magnitude_order"]
        == report["final_controls"]["magnitude_order"]
        == native
    )
    assert [f["id"] for f in c["fields"]] == s["fields"]
    (fi,) = [i for i, f in enumerate(c["fields"]) if f["id"] == label]
    field = dict(id=label, kind="forest", physical="lya", background=background)
    assert c["fields"][fi] == field and fi == (4 if bin_id == 0 else 0)
    (pair,) = [i for i, p in enumerate(c["selected_pairs"]) if p == [fi, fi]]
    (db,) = [b for b in diagnosis["bins"] if b["bin"] == bin_id]
    assert db["arrays"] == f"bin-{bin_id}.npz"
    assert db["sha256"] == HASHES[f"weight-diagnosis/bin-{bin_id}.npz"]
    partition = np.asarray(db["partition"])
    assert np.array_equal(partition, s["partition"]) and np.all(np.diff(partition) > 0)
    assert len(partition) == intervals + 1
    raw = s["samples"][label]
    policies = {
        k: {
            name: value
            for name, value in raw[k + "_diagnostics"]["policies"].items()
            if not isinstance(value, (list, dict))
        }
        for k in ("density", "snr")
    }
    out.update(
        identity=field,
        selected_pair_index=pair,
        field_index=fi,
        bounds=bounds,
        parameters=c["parameters"],
        policies=policies,
        source_width_policy=s["source_width_policy"],
    )
    z, h = s["model"]["z_eval"], s["grid"]["h_fid"]
    av, dd = s["geometry"]["a_v"], s["geometry"]["d_deg"]
    recipe_path = "fishhighz/validation/_recipes.py"
    assert (
        digest(ROOT / recipe_path)
        == report["provenance"]["fishhighz"]["module_hashes"][recipe_path]
    )
    recipe = RECIPES[c["case"]]["survey"]
    resolution = float(recipe["resolution"])
    wavelength = raw["snr_diagnostics"]["original"]["wavelength"]
    pixel = 299792.458 * policies["snr"]["pixel_width_angstrom"] / wavelength
    assert s["resolution"] == "fwhm"
    sigma = 299792.458 / (resolution * 2 * math.sqrt(2 * math.log(2)))
    geometry = prepare_geometry(
        *bounds,
        z_eval=z,
        area_deg2=float(recipe["survey_area"]),
        h_fid=h,
        hubble=lambda x: np.full_like(x, av * h * (1 + z)),
        transverse_distance=lambda x: np.full_like(x, dd / (h * np.pi / 180)),
        z_order=1,
    )
    check(out, "geometry", [geometry.a_v, geometry.d_deg], [av, dd])
    length, zs = raw["length_velocity"], raw["z_source"]
    out["held_scalars"] = dict(
        z_eval=z,
        z_source=zs,
        L=length,
        l_p=pixel,
        h_fid=h,
        a_v=av,
        d_deg=dd,
        wavelength=wavelength,
        resolution=resolution,
        sigma_velocity=sigma,
    )
    out["geometry_policy"] = "W09 constant-background adapter; artificial volume unused"
    with np.load(
        BASE / f"weight-diagnosis/bin-{bin_id}.npz", allow_pickle=False
    ) as data:
        aliases = []
        for order in (16, 32, 64):
            result["stage"] = f"{label} bin {bin_id}: coefficients order {order}"
            rows = [r for r in db["rows"] if r["order"] == order]
            metadata = [
                next(f for f in row["fields"] if f["field"] == label) for row in rows
            ]
            assert metadata and all(f == metadata[0] for f in metadata)
            assert rows[0]["fields"][fi if fi == 0 else 1]["field"] == label
            b = metadata[0]["alias"]
            assert math.isfinite(b) and b > 0
            aliases.append(b)
            assert all(x == b for x in aliases)
            out["held_scalars"]["B_star"] = b
            m, q, r, v = [
                data[f"order{order}_field{fi}_{k}"]
                for k in ("magnitudes", "measure", "masses", "variance")
            ]
            assert all(
                x.shape == (intervals * order,) and np.all(np.isfinite(x))
                for x in (m, q, r, v)
            )
            assert (
                np.all(np.diff(m) > 0)
                and np.all(q > 0)
                and np.all(r >= 0)
                and np.any(r > 0)
                and np.all(v >= 0)
            )
            co = dict(
                nodes=len(m),
                positive_masses=int(np.count_nonzero(r)),
                mass_sum=math.fsum(map(float, r)),
            )
            out["coefficients"][str(order)] = co
            for j, (lo, hi) in enumerate(
                zip(partition[:-1], partition[1:], strict=True)
            ):
                sl = slice(j * order, (j + 1) * order)
                assert np.all((m[sl] > lo) & (m[sl] < hi))
                check(
                    co,
                    f"quadrature_interval_{j}",
                    math.fsum(map(float, q[sl])),
                    hi - lo,
                )
            rho = r / q
            if order == native:
                for key, x in (("magnitudes", m), ("quadrature", q), ("variance", v)):
                    check(co, "native_" + key, x, raw[key])
                rho = np.asarray(raw["density"]) * ((1 + zs) / 299792.458)
            assert np.array_equal(rho > 0, r > 0)
            check(co, "mass_reconstruction", rho * q, r)
            prepared = prepare_forest_weights(
                ObservedField(label, "forest", "lya", background),
                geometry,
                InstrumentResponse(pixel, sigma),
                z_source=zs,
                magnitudes=m,
                quadrature=q,
                rho=rho,
                variance=v,
                length_velocity=length,
                method="inverse_variance",
                alias=b,
            )
            nu = [b / (b + pixel * float(x)) for x in v]
            i1 = math.fsum(float(x) * w for x, w in zip(r, nu, strict=True))
            i2 = math.fsum(float(x) * w * w for x, w in zip(r, nu, strict=True))
            i3 = math.fsum(
                float(x) * w * w * float(y) for x, w, y in zip(r, nu, v, strict=True)
            )
            a, p = i2 / (length * i1**2), pixel * i3 / (length * i1**2)
            public = [prepared.A, prepared.P_pixel, prepared.A * b + prepared.P_pixel]
            scalar = [a, p, a * b + p]
            co.update(
                public=public,
                scalar=scalar,
                scalar_moments=[i1, i2, i3],
                historical=metadata[0]["fixed_coefficients"],
                Q_identity=b / (length * i1),
            )
            check(
                co,
                "moments",
                [prepared.I1[-1], prepared.I2[-1], prepared.I3[-1]],
                [i1, i2, i3],
            )
            check(co, "public_scalar", public, scalar)
            check(co, "historical", public[:2], co["historical"])
            check(co, "Q_identity", scalar[2], co["Q_identity"])
        for lo, hi in ((16, 32), (32, 64), (16, 64)):
            change = contrast(
                out["coefficients"][str(hi)]["public"],
                out["coefficients"][str(lo)]["public"],
            )
            out["refinement"][f"{hi}/{lo}"] = change
            if max(map(abs, change["fractional"])) > 0.001:
                result["status"] = "reference coefficient instability"
                return False
        out["coefficient_stability"] = True
        errors = {}
        for order, count in ((32, 0), (64, 0), (32, 3), (32, 6)):
            result["stage"] = f"{label} bin {bin_id}: historical Fisher {order}/{count}"
            (row,) = [
                r for r in db["rows"] if (r["order"], r["iterations"]) == (order, count)
            ]
            matrix = {}
            out["matrices"][f"{order}/{count}"] = matrix
            if count and not row["available"]:
                matrix["unavailable"] = "historical matrix unavailable; no regeneration"
                continue
            errors[order, count] = invert(data, row, pair, matrix)
            if (order, count) == (64, 0):
                change = contrast(errors[64, 0], errors[32, 0])
                out["reference_error_refinement"] = change
                if max(map(abs, change["fractional"])) > 0.001:
                    result["status"] = "reference BAO error instability"
                    return False
            if count:
                change = contrast(errors[order, count], errors[32, 0])
                change["percent"] = [100 * x for x in change["fractional"]]
                matrix["cumulative_over_reference"] = change
    return True


def stop(_signum, _frame):
    raise TimeoutError("W10 60-second cap; no extension")


if __name__ == "__main__":
    start = time.monotonic()
    output = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w10-r1-%Y%m%dT%H%M%S%fZ")
    )
    output.mkdir(parents=True, exist_ok=False)
    result = dict(
        status="incomplete",
        samples=[],
        input_sha256=HASHES,
        source_sha256={p: digest(ROOT / p) for p in SOURCES},
        python=platform.python_version(),
        numpy=np.__version__,
        host=platform.node(),
    )
    signal.signal(signal.SIGALRM, stop)
    signal.alarm(60)
    try:
        result["threads"] = {
            k: os.environ.get(k)
            for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        }
        assert set(result["threads"].values()) == {"1"}
        assert all(digest(BASE / p) == h for p, h in HASHES.items())
        diagnosis = json.loads((BASE / "weight-diagnosis/diagnosis.json").read_text())
        with np.errstate(all="raise"):
            for args in (
                (0, "lya(lbg)", "lbg", 32, 1, [2.0, 2.235], 162),
                (5, "lya(qso)", "qso", 64, 11, [3.175, 3.410], 165),
            ):
                if not sample(result, diagnosis, *args):
                    break
            else:
                result["status"] = (
                    "numerical and reference-stability checks passed; awaiting review"
                )
    except Exception as error:
        result["status"] = f"stopped: {type(error).__name__}: {error}"
    finally:
        signal.alarm(0)
        result["identities_unchanged"] = all(
            digest(BASE / p) == h for p, h in HASHES.items()
        ) and all(digest(ROOT / p) == h for p, h in result["source_sha256"].items())
        result["elapsed_seconds"] = time.monotonic() - start
        (output / "summary.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        print(output, result["status"], result["elapsed_seconds"], flush=True)
