"""W11: bounded sensitivity to the fixed forest reference parallel mode."""

import hashlib
import json
import math
import os
import platform
import signal
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import fisher_matrix
from fishhighz.geometry import SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.models.p1d import default_p1d, p1d_floor
from fishhighz.noise import forest_noise
from fishhighz.response import InstrumentResponse, velocity_response
from fishhighz.validation._recipes import RECIPES
from fishhighz.weights import prepare_forest_weights

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / ".validation/step12-r5-20260914T191855Z"
HASHES = {
    "profiles-checked/records-001.report.json": (
        "9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666"
    ),
    "profiles-checked/records-001.npz": (
        "e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0"
    ),
    "weight-diagnosis/diagnosis.json": (
        "069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008"
    ),
    "weight-diagnosis/bin-0.npz": (
        "d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73"
    ),
}
SOURCES = [
    "WEIGHTING_DIAGNOSTIC_STEP.md",
    "scripts/check_forest_reference_mode.py",
    "fishhighz/weights.py",
    "fishhighz/kernels/weights.py",
    "fishhighz/noise.py",
    "fishhighz/response.py",
    "fishhighz/kernels/response.py",
    "fishhighz/covariance.py",
    "fishhighz/kernels/covariance.py",
    "fishhighz/fisher.py",
    "fishhighz/kernels/fisher.py",
    "fishhighz/models/p1d.py",
    "fishhighz/geometry.py",
    "fishhighz/fields.py",
    "fishhighz/_arrays.py",
    "fishhighz/_information.py",
    "fishhighz/validation/_recipes.py",
]
REFERENCE_MODES = (0.00035, 0.001, 0.003)
RTOL = 5e-12


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def compare(result, name, actual, expected):
    actual = np.asarray(actual, dtype=np.float64)
    expected = np.asarray(expected, dtype=np.float64)
    assert actual.shape == expected.shape, name
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), name
    zero = expected == 0
    assert np.all(actual[zero] == 0), name
    residual = float(
        np.max(np.abs((actual[~zero] - expected[~zero]) / expected[~zero]), initial=0)
    )
    result["residuals"][name] = residual
    assert residual <= RTOL, (name, residual)


def compare_signed(result, name, actual, expected):
    actual = np.asarray(actual, dtype=np.float64)
    expected = np.asarray(expected, dtype=np.float64)
    assert actual.shape == expected.shape, name
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), name
    residual = float(np.max(np.abs(actual - expected), initial=0))
    result["residuals"][name] = residual
    assert residual <= RTOL, (name, residual)


def scalar_p1d_shape(z):
    evolution = math.log1p(z) - math.log(4.0)
    slope = -2.55 - 0.28 * evolution
    floor = 0.009 * math.exp((-0.5 * slope - 1.0) / -0.1)
    return evolution, slope, floor


def scalar_p1d(z, q):
    evolution, slope, floor = scalar_p1d_shape(z)
    u = math.log(max(float(q), floor)) - math.log(0.009)
    return (math.pi * 0.064 / 0.009) * math.exp(
        (2.0 + slope) * u - 0.1 * u * u + 3.55 * evolution
    )


def scalar_response(q, pixel, sigma):
    x = float(q) * pixel / 2.0
    sinc = 1.0 if x == 0.0 else math.sin(x) / x
    return sinc * math.exp(-0.5 * (float(q) * sigma) ** 2)


def scalar_fisher(power, noise, modes, jacobian):
    covariance = np.fromiter(
        (
            2.0 * (float(p) + float(n)) ** 2 / float(count)
            for p, n, count in zip(power, noise, modes)
        ),
        dtype=np.float64,
        count=len(power),
    )
    fisher = np.array(
        [
            [
                math.fsum(
                    float(row[a]) * float(row[b]) / float(covariance_node)
                    for row, covariance_node in zip(jacobian, covariance)
                )
                for b in range(2)
            ]
            for a in range(2)
        ]
    )
    return covariance, fisher


def invert_2x2(fisher):
    assert np.array_equal(fisher, fisher.T)
    a, b, d = map(float, (fisher[0, 0], fisher[0, 1], fisher[1, 1]))
    determinant = a * d - b * b
    assert a > 0 and d > 0 and determinant > 0
    covariance = np.array([[d, -b], [-b, a]]) / determinant
    errors = np.array([math.sqrt(d / determinant), math.sqrt(a / determinant)])
    return {
        "fisher": fisher.tolist(),
        "determinant": determinant,
        "covariance": covariance.tolist(),
        "errors": errors.tolist(),
        "normalized_min_eigenvalue": 1.0 - abs(b) / math.sqrt(a * d),
    }


def prepare_context(result, report):
    context, settings = report["context"], report["settings"]
    assert context["profile"] == settings["profile"] == "accuracy"
    assert context["bin"] == 0 and context["bounds"] == [2.0, 2.235]
    assert settings["bounds"] == context["bounds"]
    assert context["parameters"] == settings["parameters"] == ["ap_0", "at_0"]
    assert report["final_controls"]["magnitude_order"] == 32
    (field_index,) = [
        i for i, field in enumerate(context["fields"]) if field["id"] == "lya(qso)"
    ]
    assert context["fields"][field_index] == {
        "id": "lya(qso)",
        "kind": "forest",
        "physical": "lya",
        "background": "qso",
    }
    (selected_pair,) = [
        i
        for i, pair in enumerate(context["selected_pairs"])
        if pair == [field_index, field_index]
    ]
    (required_pair,) = [
        i
        for i, pair in enumerate(context["required_pairs"])
        if pair == [field_index, field_index]
    ]
    assert selected_pair == required_pair == 0
    sample = settings["samples"]["lya(qso)"]
    assert sample["density_diagnostics"]["policies"]["negative"] == "floor_negative"
    z = float(settings["model"]["z_eval"])
    h_fid = float(settings["grid"]["h_fid"])
    a_v = float(settings["geometry"]["a_v"])
    d_deg = float(settings["geometry"]["d_deg"])
    recipe_path = "fishhighz/validation/_recipes.py"
    assert (
        digest(ROOT / recipe_path)
        == report["provenance"]["fishhighz"]["module_hashes"][recipe_path]
    )
    recipe = RECIPES[context["case"]]
    resolution = float(recipe["survey"]["resolution"])
    pixel = (
        SPEED_LIGHT_KMS
        * sample["snr_diagnostics"]["policies"]["pixel_width_angstrom"]
        / sample["snr_diagnostics"]["original"]["wavelength"]
    )
    assert settings["resolution"] == "fwhm"
    sigma = SPEED_LIGHT_KMS / (resolution * 2.0 * math.sqrt(2.0 * math.log(2.0)))
    instrument = InstrumentResponse(pixel, sigma)
    geometry = prepare_geometry(
        *context["bounds"],
        z_eval=z,
        area_deg2=float(recipe["survey"]["survey_area"]),
        h_fid=h_fid,
        hubble=lambda x: np.full_like(x, a_v * h_fid * (1.0 + z)),
        transverse_distance=lambda x: np.full_like(x, d_deg / (h_fid * np.pi / 180.0)),
        z_order=1,
    )
    compare(result, "geometry_conversion", [geometry.a_v, geometry.d_deg], [a_v, d_deg])
    result["selection"] = {
        "profile": "accuracy",
        "bin": 0,
        "bounds": context["bounds"],
        "field": "lya(qso)",
        "field_index": field_index,
        "selected_pair_index": selected_pair,
        "required_pair_index": required_pair,
        "parameters": context["parameters"],
        "magnitude_order": 32,
        "magnitude_nodes": len(sample["magnitudes"]),
        "fourier_nodes": 16384,
    }
    result["fixed_scalars"] = {
        "z_eval": z,
        "z_source": sample["z_source"],
        "h_fid": h_fid,
        "a_v": a_v,
        "d_deg": d_deg,
        "pixel_velocity": pixel,
        "sigma_velocity": sigma,
        "resolution": resolution,
        "length_velocity": sample["length_velocity"],
    }
    result["geometry_adapter"] = (
        "W09 constant-H and constant-D_M context only; artificial volume discarded"
    )
    return context, settings, sample, field_index, geometry, instrument


def run(result):
    result["stage"] = "identities and saved-input selection"
    result["input_sha256"] = dict(HASHES)
    result["source_sha256"] = {path: digest(ROOT / path) for path in SOURCES}
    assert all(digest(INPUT / path) == expected for path, expected in HASHES.items())
    report = json.loads(
        (INPUT / "profiles-checked/records-001.report.json").read_text()
    )
    diagnosis = json.loads((INPUT / "weight-diagnosis/diagnosis.json").read_text())
    context, settings, sample, field_index, geometry, instrument = prepare_context(
        result, report
    )
    (bin0,) = [entry for entry in diagnosis["bins"] if entry["bin"] == 0]
    assert bin0["arrays"] == "bin-0.npz"
    assert bin0["sha256"] == HASHES["weight-diagnosis/bin-0.npz"]
    (baseline_row,) = [
        row for row in bin0["rows"] if row["order"] == 32 and row["iterations"] == 0
    ]
    assert baseline_row["available"] and baseline_row["fixed_initial_weights"]
    assert baseline_row["pair_fisher_key"] == "order32_iterations0_pair_fisher"

    with (
        np.load(INPUT / "profiles-checked/records-001.npz") as profile_arrays,
        np.load(INPUT / "weight-diagnosis/bin-0.npz") as diagnosis_arrays,
    ):
        assert np.array_equal(
            profile_arrays["selected_pairs"], context["selected_pairs"]
        )
        assert np.array_equal(
            profile_arrays["required_pairs"], context["required_pairs"]
        )
        k = profile_arrays["k"]
        mu = profile_arrays["mu"]
        modes = profile_arrays["modes"]
        power = profile_arrays["observed_signal"][:, 0]
        jacobian = profile_arrays["observed_j"][:, 0, :]
        assert k.shape == mu.shape == modes.shape == power.shape == (16384,)
        assert jacobian.shape == (16384, 2)
        assert np.all(modes > 0)

        prefix = f"order32_field{field_index}_"
        magnitudes, quadrature, variance = (
            np.asarray(sample[name])
            for name in ("magnitudes", "quadrature", "variance")
        )
        rho = np.asarray(sample["density"]) * (
            (1.0 + sample["z_source"]) / SPEED_LIGHT_KMS
        )
        masses = diagnosis_arrays[prefix + "masses"]
        assert (
            magnitudes.shape
            == quadrature.shape
            == variance.shape
            == rho.shape
            == (5184,)
        )
        assert np.all(rho >= 0) and np.all(quadrature > 0)
        for name, actual, expected in (
            ("masses", rho * quadrature, masses),
            ("magnitudes", magnitudes, diagnosis_arrays[prefix + "magnitudes"]),
            ("quadrature", quadrature, diagnosis_arrays[prefix + "measure"]),
            ("variance", variance, diagnosis_arrays[prefix + "variance"]),
        ):
            compare(result, f"saved_{name}", actual, expected)

        z = result["fixed_scalars"]["z_eval"]
        pixel = result["fixed_scalars"]["pixel_velocity"]
        sigma = result["fixed_scalars"]["sigma_velocity"]
        a_v = result["fixed_scalars"]["a_v"]
        d_deg = result["fixed_scalars"]["d_deg"]
        length = result["fixed_scalars"]["length_velocity"]
        q_modes = k * mu / a_v
        physical_p1d = default_p1d([], z, q_modes)
        physical_response = velocity_response(
            q_modes,
            pixel_width_velocity=pixel,
            gaussian_sigma_velocity=sigma,
        )
        compare(
            result,
            "saved_physical_response",
            physical_response,
            profile_arrays["response"][:, field_index],
        )
        scalar_mode_p1d = np.fromiter(
            (scalar_p1d(z, q) for q in q_modes),
            dtype=np.float64,
            count=len(q_modes),
        )
        scalar_mode_response = np.fromiter(
            (scalar_response(q, pixel, sigma) for q in q_modes),
            dtype=np.float64,
            count=len(q_modes),
        )
        compare(result, "physical_p1d_scalar", physical_p1d, scalar_mode_p1d)
        compare(
            result,
            "physical_response_scalar",
            physical_response,
            scalar_mode_response,
        )
        scalar_floor = scalar_p1d_shape(z)[2]
        public_floor = p1d_floor(z)
        compare(result, "p1d_floor_scalar", public_floor, scalar_floor)
        q_min, q_max = map(float, (np.min(q_modes), np.max(q_modes)))
        assert q_min <= REFERENCE_MODES[0] < public_floor
        assert all(public_floor < q_star <= q_max for q_star in REFERENCE_MODES[1:])
        result["reference_domain"] = {
            "saved_q_min": q_min,
            "saved_q_max": q_max,
            "p1d_floor": public_floor,
            "baseline_below_floor": True,
            "alternatives_above_floor": True,
            "all_reference_modes_inside_saved_domain": True,
        }

        field = ObservedField("lya(qso)", "forest", "lya", "qso")
        selection = PairSelection([field], [("lya(qso)", "lya(qso)")])
        baseline_errors = None
        baseline_b = None
        baseline_q0 = None
        runs = []
        omitted = []
        result["stage"] = "ordered public and scalar sensitivity calculations"
        for trial_index, q_star in enumerate(REFERENCE_MODES):
            public_reference_p1d = float(default_p1d([], z, [q_star])[0])
            public_reference_response = float(
                velocity_response(
                    [q_star],
                    pixel_width_velocity=pixel,
                    gaussian_sigma_velocity=sigma,
                )[0]
            )
            public_b = public_reference_p1d * public_reference_response**2
            independent_reference_p1d = scalar_p1d(z, q_star)
            independent_reference_response = scalar_response(q_star, pixel, sigma)
            independent_b = (
                independent_reference_p1d * independent_reference_response**2
            )
            compare(
                result,
                f"q{q_star:g}_reference_p1d_scalar",
                public_reference_p1d,
                independent_reference_p1d,
            )
            compare(
                result,
                f"q{q_star:g}_reference_response_scalar",
                public_reference_response,
                independent_reference_response,
            )
            compare(
                result,
                f"q{q_star:g}_B_star_scalar",
                public_b,
                independent_b,
            )
            if trial_index == 0:
                baseline_b = public_b
                compare(result, "baseline_saved_alias", public_b, 16.356748967852234)

            prepared = prepare_forest_weights(
                field,
                geometry,
                instrument,
                z_source=sample["z_source"],
                magnitudes=magnitudes,
                quadrature=quadrature,
                rho=rho,
                variance=variance,
                length_velocity=length,
                method="inverse_variance",
                alias=public_b,
            )
            independent_weights = np.fromiter(
                (
                    independent_b / (independent_b + pixel * float(v)) if r > 0 else 0.0
                    for r, v in zip(masses, variance)
                ),
                dtype=np.float64,
                count=len(masses),
            )
            moments = [
                math.fsum(
                    float(r) * float(weight)
                    for r, weight in zip(masses, independent_weights)
                ),
                math.fsum(
                    float(r) * float(weight) ** 2
                    for r, weight in zip(masses, independent_weights)
                ),
                math.fsum(
                    float(r) * float(weight) ** 2 * float(v)
                    for r, weight, v in zip(masses, independent_weights, variance)
                ),
            ]
            denominator = length * moments[0] ** 2
            independent_a = moments[1] / denominator
            independent_pixel = pixel * moments[2] / denominator
            compare(
                result,
                f"q{q_star:g}_weights_scalar",
                prepared.weights,
                independent_weights,
            )
            compare(
                result,
                f"q{q_star:g}_moments_scalar",
                [prepared.I1[-1], prepared.I2[-1], prepared.I3[-1]],
                moments,
            )
            compare(
                result,
                f"q{q_star:g}_coefficients_scalar",
                [prepared.A, prepared.P_pixel],
                [independent_a, independent_pixel],
            )

            public_noise = forest_noise(
                prepared,
                field,
                geometry,
                instrument,
                k,
                mu,
                physical_p1d,
            )
            conversion = d_deg * d_deg / a_v
            scalar_noise = np.fromiter(
                (
                    (independent_a * p1d_value * response_value**2 + independent_pixel)
                    * conversion
                    for p1d_value, response_value in zip(
                        scalar_mode_p1d, scalar_mode_response
                    )
                ),
                dtype=np.float64,
                count=len(k),
            )
            compare(
                result,
                f"q{q_star:g}_mode_noise_scalar",
                public_noise.total,
                scalar_noise,
            )
            compare(
                result,
                f"q{q_star:g}_pixel_unsmoothed",
                public_noise.pixel,
                np.full(k.shape, independent_pixel * conversion),
            )
            total = combine_observed_power(power[:, None], public_noise.total[:, None])
            public_covariance = gaussian_covariance(total, modes, selection)
            public_fisher = fisher_matrix(jacobian[:, None, :], public_covariance)
            scalar_covariance, scalar_fisher_matrix = scalar_fisher(
                power, scalar_noise, modes, jacobian
            )
            compare(
                result,
                f"q{q_star:g}_covariance_scalar",
                public_covariance[:, 0, 0],
                scalar_covariance,
            )
            compare(
                result,
                f"q{q_star:g}_fisher_scalar",
                public_fisher,
                scalar_fisher_matrix,
            )
            public_inversion = invert_2x2(public_fisher)
            scalar_inversion = invert_2x2(scalar_fisher_matrix)
            compare(
                result,
                f"q{q_star:g}_errors_scalar",
                public_inversion["errors"],
                scalar_inversion["errors"],
            )

            q_identity = prepared.A * public_b + prepared.P_pixel
            q_identity_rhs = public_b / (length * float(prepared.I1[-1]))
            independent_identity = independent_a * independent_b + independent_pixel
            independent_identity_rhs = independent_b / (length * moments[0])
            compare(
                result,
                f"q{q_star:g}_within_run_Q_identity",
                q_identity,
                q_identity_rhs,
            )
            compare(
                result,
                f"q{q_star:g}_within_run_Q_identity_scalar",
                independent_identity,
                independent_identity_rhs,
            )
            q0 = prepared.A * baseline_b + prepared.P_pixel
            independent_q0 = independent_a * baseline_b + independent_pixel
            compare(result, f"q{q_star:g}_Q0_scalar", q0, independent_q0)
            if trial_index == 0:
                baseline_q0 = q0
                baseline_errors = np.asarray(public_inversion["errors"])
                sensitivity = np.zeros(2)
                saved_fisher = diagnosis_arrays[baseline_row["pair_fisher_key"]][0]
                compare(result, "baseline_saved_fisher", public_fisher, saved_fisher)
                compare(
                    result,
                    "baseline_saved_errors",
                    baseline_errors,
                    diagnosis_arrays["order32_iterations0_pair_errors"][0],
                )
                compare(
                    result,
                    "baseline_printed_errors",
                    baseline_errors,
                    [0.029619522112556678, 0.027975323168394890],
                )
            else:
                assert q0 >= baseline_q0 * (1.0 - RTOL), (q_star, q0, baseline_q0)
                sensitivity = (
                    np.asarray(public_inversion["errors"]) / baseline_errors - 1.0
                )
                independent_sensitivity = (
                    np.asarray(scalar_inversion["errors"]) / baseline_errors - 1.0
                )
                compare_signed(
                    result,
                    f"q{q_star:g}_sensitivity_scalar",
                    sensitivity,
                    independent_sensitivity,
                )

            runs.append(
                {
                    "q_star": q_star,
                    "reference": {
                        "intrinsic_p1d": public_reference_p1d,
                        "response": public_reference_response,
                        "response_squared": public_reference_response**2,
                        "B_star": public_b,
                    },
                    "coefficients": {
                        "A": prepared.A,
                        "P_pixel": prepared.P_pixel,
                        "I1": float(prepared.I1[-1]),
                        "I2": float(prepared.I2[-1]),
                        "I3": float(prepared.I3[-1]),
                    },
                    "Q0_at_baseline_physical_mode": q0,
                    "within_run_Q": q_identity,
                    "public": public_inversion,
                    "independent": scalar_inversion,
                    "sensitivity_fraction_using_baseline_errors": sensitivity.tolist(),
                    "sensitivity_percent_using_baseline_errors": (
                        100.0 * sensitivity
                    ).tolist(),
                }
            )
            print(
                f"q_star={q_star:.8g}: errors={public_inversion['errors']}, "
                f"Delta={sensitivity.tolist()}",
                flush=True,
            )
            if trial_index > 0 and np.any(np.abs(sensitivity) > 0.001):
                omitted = list(REFERENCE_MODES[trial_index + 1 :])
                result["scientific_outcome"] = "reference-mode sensitivity detected"
                break
        else:
            result["scientific_outcome"] = (
                "no sensitivity above 0.1% at the tested points"
            )

        result["runs"] = runs
        result["public_weight_preparations"] = len(runs)
        result["single_spectrum_mode_sums"] = len(runs)
        result["actual_stopping_point"] = runs[-1]["q_star"]
        result["omitted_reference_modes"] = omitted
        result["held_quantities"] = [
            "saved observed signal",
            "saved mean Jacobian",
            "saved k, mu and mode counts",
            "geometry conversions and response parameters",
            "mode-dependent intrinsic P1D and response",
            "magnitude nodes, masses, variance, quadrature, L and pixel width",
        ]
        result["varied_quantity"] = (
            "q_star -> intrinsic P1D(q_star)*W(q_star)^2 -> fixed weights; "
            "noise and covariance recomputed between trials"
        )
        result["worst_residual"] = max(
            result["residuals"].items(), key=lambda item: item[1]
        )

    assert all(digest(INPUT / path) == expected for path, expected in HASHES.items())
    assert all(
        digest(ROOT / path) == expected
        for path, expected in result["source_sha256"].items()
    )
    result["status"] = "completed; awaiting independent and user review"


def stop(signum, frame):
    raise TimeoutError("W11 60-second cap; no extension")


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, stop)
    signal.alarm(60)
    start = time.monotonic()
    threads = {
        name: os.environ.get(name)
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    }
    assert set(threads.values()) == {"1"}
    assert os.environ.get("FISHHIGHZ_FISHER_BACKEND", "numpy") == "numpy"
    output = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w11-r1-%Y%m%dT%H%M%S%fZ")
    )
    output.mkdir(parents=True, exist_ok=False)
    result = {
        "status": "incomplete",
        "stage": "startup",
        "residuals": {},
        "threads": threads,
        "backend": "numpy",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "command": (
            "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
            "timeout 60 .venv/bin/python -B scripts/check_forest_reference_mode.py"
        ),
    }
    try:
        run(result)
    except Exception as error:
        result.update(status="stopped", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        signal.alarm(0)
        result["elapsed_seconds"] = time.monotonic() - start
        (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(output, result["status"], flush=True)
