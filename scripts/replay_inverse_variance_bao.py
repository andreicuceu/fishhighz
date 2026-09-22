"""W09: one saved QSO-forest BAO mode sum, with scalar controls first."""

import hashlib
import json
import math
import os
import platform
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import fisher_matrix
from fishhighz.geometry import SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.models.p1d import default_p1d
from fishhighz.noise import forest_noise
from fishhighz.response import InstrumentResponse, velocity_response
from fishhighz.validation._recipes import RECIPES
from fishhighz.weights import prepare_forest_weights

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / ".validation/step12-r5-20260914T191855Z"
HASHES = {
    "profiles-checked/records-001.report.json": "9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666",
    "profiles-checked/records-001.npz": "e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0",
    "weight-diagnosis/diagnosis.json": "069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008",
    "weight-diagnosis/bin-0.npz": "d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73",
}
SOURCES = [
    "scripts/replay_inverse_variance_bao.py",
    "scripts/diagnose_desi2_weights.py",
    "scripts/compare_forest_weight_bao_errors.py",
    "docs/archive/notes/WEIGHTING_DIAGNOSTIC_STEP.md",
    "fishhighz/weights.py",
    "fishhighz/noise.py",
    "fishhighz/response.py",
    "fishhighz/geometry.py",
    "fishhighz/fields.py",
    "fishhighz/covariance.py",
    "fishhighz/fisher.py",
    "fishhighz/_arrays.py",
    "fishhighz/models/p1d.py",
    "fishhighz/kernels/weights.py",
    "fishhighz/kernels/response.py",
    "fishhighz/kernels/covariance.py",
    "fishhighz/kernels/fisher.py",
    "fishhighz/validation/_recipes.py",
]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def check(result, name, actual, expected):
    a, b = np.asarray(actual), np.asarray(expected)
    assert a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), (
        name
    )
    zero = b == 0
    assert np.all(a[zero] == 0), name
    residual = float(np.max(abs((a[~zero] - b[~zero]) / b[~zero]), initial=0))
    result["residuals"][name] = residual
    assert residual <= 5e-12, (name, residual)


def scalar_fisher(power, noise, modes, jacobian):
    covariance = np.array(
        [
            2 * (float(p) + float(n)) ** 2 / float(m)
            for p, n, m in zip(power, noise, modes)
        ]
    )
    fisher = np.array(
        [
            [
                math.fsum(
                    float(j[a]) * float(j[b]) / float(c)
                    for j, c in zip(jacobian, covariance)
                )
                for b in range(2)
            ]
            for a in range(2)
        ]
    )
    return covariance, fisher


def invert(fisher):
    assert np.array_equal(fisher, fisher.T)
    a, b, d = map(float, (fisher[0, 0], fisher[0, 1], fisher[1, 1]))
    determinant = a * d - b * b
    assert a > 0 and d > 0 and determinant > 0
    covariance = np.array([[d, -b], [-b, a]]) / determinant
    errors = np.array([math.sqrt(d / determinant), math.sqrt(a / determinant)])
    return dict(
        fisher=fisher.tolist(),
        determinant=determinant,
        covariance=covariance.tolist(),
        errors=errors.tolist(),
        normalized_min_eigenvalue=1 - abs(b) / math.sqrt(a * d),
    )


def run(result):
    result["stage"] = "input identities and selection"
    result["input_sha256"] = HASHES
    result["source_sha256"] = {p: digest(ROOT / p) for p in SOURCES}
    assert all(digest(INPUT / p) == h for p, h in HASHES.items())
    report = json.loads(
        (INPUT / "profiles-checked/records-001.report.json").read_text()
    )
    diagnosis = json.loads((INPUT / "weight-diagnosis/diagnosis.json").read_text())
    c, s = report["context"], report["settings"]
    assert c["profile"] == s["profile"] == "accuracy" and c["bin"] == 0
    assert c["bounds"] == s["bounds"] == [2.0, 2.235]
    assert c["parameters"] == s["parameters"] == ["ap_0", "at_0"]
    assert report["final_controls"]["magnitude_order"] == 32
    assert report["final_controls"]["iterations"] == 6
    (fi,) = [i for i, f in enumerate(c["fields"]) if f["id"] == "lya(qso)"]
    assert c["fields"][fi] == dict(
        id="lya(qso)", kind="forest", physical="lya", background="qso"
    )
    (pair,) = [i for i, p in enumerate(c["selected_pairs"]) if p == [fi, fi]]
    (required,) = [i for i, p in enumerate(c["required_pairs"]) if p == [fi, fi]]
    (bin0,) = [b for b in diagnosis["bins"] if b["bin"] == 0]
    assert (
        bin0["arrays"] == "bin-0.npz"
        and bin0["sha256"] == HASHES["weight-diagnosis/bin-0.npz"]
    )
    rows = {}
    for count in (0, 3, 6):
        (row,) = [
            r for r in bin0["rows"] if r["order"] == 32 and r["iterations"] == count
        ]
        assert row["available"] and row["fixed_initial_weights"] == (count == 0)
        assert row["pair_fisher_key"] == f"order32_iterations{count}_pair_fisher"
        rows[count] = row
    (coeff,) = [f for f in rows[6]["fields"] if f["field"] == "lya(qso)"]
    (t6,) = [r for r in coeff["rows"] if r["iterations"] == 6 and r["available"]]
    sample = s["samples"]["lya(qso)"]
    assert sample["density_diagnostics"]["policies"]["negative"] == "floor_negative"
    av, ddeg = s["geometry"]["a_v"], s["geometry"]["d_deg"]
    z, h = s["model"]["z_eval"], s["grid"]["h_fid"]
    # Only these conversion factors are used; the adapter volume is artificial.
    recipe_path = "fishhighz/validation/_recipes.py"
    assert (
        digest(ROOT / recipe_path)
        == report["provenance"]["fishhighz"]["module_hashes"][recipe_path]
    )
    recipe = RECIPES[c["case"]]
    resolution = float(recipe["survey"]["resolution"])
    pixel = (
        SPEED_LIGHT_KMS
        * sample["snr_diagnostics"]["policies"]["pixel_width_angstrom"]
        / sample["snr_diagnostics"]["original"]["wavelength"]
    )
    assert s["resolution"] == "fwhm"
    sigma = SPEED_LIGHT_KMS / (resolution * 2 * math.sqrt(2 * math.log(2)))
    instrument = InstrumentResponse(pixel, sigma)
    geometry = prepare_geometry(
        *c["bounds"],
        z_eval=z,
        area_deg2=float(recipe["survey"]["survey_area"]),
        h_fid=h,
        hubble=lambda x: np.full_like(x, av * h * (1 + z)),
        transverse_distance=lambda x: np.full_like(x, ddeg / (h * np.pi / 180)),
        z_order=1,
    )
    check(result, "geometry_conversion", [geometry.a_v, geometry.d_deg], [av, ddeg])
    result["selection"] = dict(
        bin=0,
        profile="accuracy",
        field="lya(qso)",
        magnitude_order=32,
        field_index=fi,
        selected_pair_index=pair,
        required_pair_index=required,
        remapped_pair=[0, 0],
        parameters=c["parameters"],
        bounds=c["bounds"],
    )
    result["fixed_scalars"] = dict(
        z_eval=z,
        z_source=sample["z_source"],
        h_fid=h,
        a_v=av,
        d_deg=ddeg,
        pixel_velocity=pixel,
        sigma_velocity=sigma,
        resolution=resolution,
        length_velocity=sample["length_velocity"],
        B_star=coeff["alias"],
    )
    result["geometry_adapter"] = (
        "constant H and D_M; artificial integrated volume never used; saved modes only"
    )
    with (
        np.load(INPUT / "profiles-checked/records-001.npz") as p,
        np.load(INPUT / "weight-diagnosis/bin-0.npz") as d,
    ):
        assert np.array_equal(p["selected_pairs"], c["selected_pairs"])
        assert np.array_equal(p["required_pairs"], c["required_pairs"])
        k, mu, modes = p["k"], p["mu"], p["modes"]
        power, jacobian = (
            p["observed_signal"][:, required],
            p["observed_j"][:, required, :],
        )
        assert k.shape == mu.shape == modes.shape == power.shape == (16384,)
        assert jacobian.shape == (16384, 2)
        result["stage"] = "historical control"
        q = k * mu / av
        w = velocity_response(
            q, pixel_width_velocity=pixel, gaussian_sigma_velocity=sigma
        )
        check(result, "historical_response", w, p["response"][:, fi])
        p1d = default_p1d([], z, q)

        def noise(a, pp):
            return (a * p1d * w * w + pp) * ddeg * ddeg / av

        old_noise = noise(t6["A"], t6["P_pixel"])
        check(result, "historical_noise", old_noise, p["noise"][:, required])
        old_cov, old_f = scalar_fisher(power, old_noise, modes, jacobian)
        check(result, "historical_fisher", old_f, d[rows[6]["pair_fisher_key"]][pair])
        check(result, "historical_primary_fisher", old_f, p["pair_fisher"][pair])
        result["historical_control"] = invert(old_f)
        result["historical_control"]["coefficients"] = [t6["A"], t6["P_pixel"]]
        print("Historical response/noise/Fisher control passed.", flush=True)

        result["stage"] = "reference coefficients and noise"
        prefix = f"order32_field{fi}_"
        m, quadrature, variance = (
            np.asarray(sample[key]) for key in ("magnitudes", "quadrature", "variance")
        )
        rho = np.asarray(sample["density"]) * (
            (1 + sample["z_source"]) / SPEED_LIGHT_KMS
        )
        masses = d[prefix + "masses"]
        assert m.shape == (5184,) and np.all(rho >= 0)
        for name, actual, expected in [
            ("masses", rho * quadrature, masses),
            ("magnitudes", m, d[prefix + "magnitudes"]),
            ("quadrature", quadrature, d[prefix + "measure"]),
            ("variance", variance, d[prefix + "variance"]),
        ]:
            check(result, name, actual, expected)
        field = ObservedField("lya(qso)", "forest", "lya", "qso")
        prepared = prepare_forest_weights(
            field,
            geometry,
            instrument,
            z_source=sample["z_source"],
            magnitudes=m,
            quadrature=quadrature,
            rho=rho,
            variance=variance,
            length_velocity=sample["length_velocity"],
            method="inverse_variance",
            alias=coeff["alias"],
        )
        nu = np.array(
            [
                coeff["alias"] / (coeff["alias"] + pixel * float(v)) if r > 0 else 0
                for r, v in zip(masses, variance)
            ]
        )
        moments = [
            math.fsum(float(r) * float(u) for r, u in zip(masses, nu)),
            math.fsum(float(r) * float(u) ** 2 for r, u in zip(masses, nu)),
            math.fsum(
                float(r) * float(u) ** 2 * float(v)
                for r, u, v in zip(masses, nu, variance)
            ),
        ]
        denominator = sample["length_velocity"] * moments[0] ** 2
        independent = [moments[1] / denominator, pixel * moments[2] / denominator]
        public_coeff = [prepared.A, prepared.P_pixel]
        check(result, "reference_weights", prepared.weights, nu)
        check(
            result,
            "reference_moments",
            [prepared.I1[-1], prepared.I2[-1], prepared.I3[-1]],
            moments,
        )
        check(result, "reference_coefficients_scalar", public_coeff, independent)
        check(
            result,
            "reference_coefficients_saved",
            public_coeff,
            coeff["fixed_coefficients"],
        )
        check(
            result,
            "reference_coefficients_W08",
            public_coeff,
            [0.029419845663300446, 0.31435492218062416],
        )
        public_noise = forest_noise(prepared, field, geometry, instrument, k, mu, p1d)
        scalar_noise = noise(*independent)
        check(result, "reference_noise", public_noise.total, scalar_noise)
        check(
            result,
            "pixel_unsmoothed",
            public_noise.pixel,
            np.full(k.shape, independent[1] * ddeg**2 / av),
        )
        result["coefficients"] = dict(
            public=public_coeff, independent=independent, moments=moments
        )

        result["stage"] = "reference covariance and Fisher"
        selection = PairSelection([field], [("lya(qso)", "lya(qso)")])
        total = combine_observed_power(power[:, None], public_noise.total[:, None])
        covariance = gaussian_covariance(total, modes, selection)
        scalar_cov, scalar_f = scalar_fisher(power, scalar_noise, modes, jacobian)
        check(result, "reference_covariance", covariance[:, 0, 0], scalar_cov)
        public_f = fisher_matrix(jacobian[:, None, :], covariance)
        check(result, "reference_fisher_scalar", public_f, scalar_f)
        saved_f = d[rows[0]["pair_fisher_key"]][pair]
        check(result, "reference_fisher_W07", public_f, saved_f)
        new, independent_result, saved = map(invert, (public_f, scalar_f, saved_f))
        check(
            result,
            "reference_errors_scalar",
            new["errors"],
            independent_result["errors"],
        )
        check(result, "reference_errors_W07", new["errors"], saved["errors"])
        check(
            result,
            "reference_errors_stored",
            new["errors"],
            d["order32_iterations0_pair_errors"][pair],
        )
        check(
            result,
            "reference_errors_printed",
            new["errors"],
            [0.029619522112556678, 0.027975323168394890],
        )
        result.update(
            new_reference=new,
            independent_reference=independent_result,
            historical_reference=saved,
        )
        result["stage"] = "historical contrasts"
        result["contrasts"] = {}
        for count, rounded in [
            (3, [7.6305714, 10.3256347]),
            (6, [8.6195623, 11.7086249]),
        ]:
            historical = invert(d[rows[count]["pair_fisher_key"]][pair])
            numerator, denominator = (
                np.array(historical["errors"]),
                np.array(new["errors"]),
            )
            ratio = numerator / denominator
            check(
                result, f"t{count}_ratio", ratio, numerator / np.array(saved["errors"])
            )
            contrast = 100 * (ratio - 1)
            assert np.all(abs(contrast - rounded) <= 5e-8)
            result["contrasts"][f"t{count}"] = dict(
                historical=historical,
                numerator=numerator.tolist(),
                denominator=denominator.tolist(),
                ratio=ratio.tolist(),
                percent=contrast.tolist(),
            )
    assert all(digest(INPUT / p) == h for p, h in HASHES.items())
    assert all(digest(ROOT / p) == h for p, h in result["source_sha256"].items())
    result["status"] = "passed; awaiting user review and acceptance"


def stop(signum, frame):
    raise TimeoutError("W09 60-second cap; no extension")


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, stop)
    signal.alarm(60)
    start = time.monotonic()
    threads = {
        k: os.environ.get(k)
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    }
    assert set(threads.values()) == {"1"}
    assert os.environ.get("FISHHIGHZ_FISHER_BACKEND", "numpy") == "numpy"
    out = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w09-r1-%Y%m%dT%H%M%S%fZ")
    )
    out.mkdir(parents=True, exist_ok=False)
    result = dict(
        status="incomplete",
        residuals={},
        threads=threads,
        backend="numpy",
        python=platform.python_version(),
        numpy=np.__version__,
        host=platform.node(),
        command="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 .venv/bin/python -B scripts/replay_inverse_variance_bao.py",
        git_head=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
    )
    (out / "git-status.txt").write_text(
        subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True)
    )
    try:
        run(result)
    except Exception as error:
        result.update(status="stopped", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic() - start
        (out / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(out, result["status"], flush=True)
