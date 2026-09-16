"""W12: saved FF/FG/GG fixed-reference calculation at magnitude orders 32/64."""

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
from fishhighz.noise import forest_noise, prepare_noise
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import InstrumentResponse, velocity_response
from fishhighz.validation._recipes import RECIPES
from fishhighz.validation.accuracy import FIXED_REFERENCE, _forest_input
from fishhighz.validation.schema import validate_payload
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
    "scripts/check_accuracy_fixed_reference.py",
    "fishhighz/validation/accuracy.py",
    "fishhighz/validation/study.py",
    "fishhighz/validation/trials.py",
    "fishhighz/validation/schema.py",
    "fishhighz/validation/profiles.py",
    "fishhighz/weights.py",
    "fishhighz/noise.py",
    "fishhighz/covariance.py",
    "fishhighz/fisher.py",
]
UNCHANGED_HISTORICAL_SOURCES = [
    "fishhighz/noise.py",
    "fishhighz/covariance.py",
    "fishhighz/fisher.py",
    "fishhighz/validation/compatibility.py",
    "fishhighz/kernels/weights.py",
    "fishhighz/kernels/covariance.py",
    "fishhighz/kernels/fisher.py",
]
REVIEWED_W09_WEIGHT_SOURCE = (
    "5e9110b63d13eded2cb84dbd0e82b958418ed59602d910bfeb1428a27d67f663"
)


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def check(result, name, actual, expected, *, rtol=5e-12):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    zero = expected == 0
    assert np.all(actual[zero] == 0), name
    residual = float(
        np.max(np.abs((actual[~zero] - expected[~zero]) / expected[~zero]), initial=0)
    )
    result["residuals"][name] = residual
    assert residual <= rtol, (name, residual)


def invert(fisher):
    fisher = np.asarray(fisher)
    assert np.array_equal(fisher, fisher.T)
    a, b, d = map(float, (fisher[0, 0], fisher[0, 1], fisher[1, 1]))
    determinant = a * d - b * b
    assert a > 0 and d > 0 and determinant > 0
    covariance = np.array([[d, -b], [-b, a]]) / determinant
    return dict(
        fisher=fisher.tolist(),
        determinant=determinant,
        covariance=covariance.tolist(),
        errors=[math.sqrt(d / determinant), math.sqrt(a / determinant)],
    )


def independent_covariance(total, modes):
    ff, fg, gg = total.T
    numerator = np.empty((len(modes), 3, 3))
    numerator[:, 0, 0] = 2 * ff**2
    numerator[:, 0, 1] = numerator[:, 1, 0] = 2 * ff * fg
    numerator[:, 0, 2] = numerator[:, 2, 0] = 2 * fg**2
    numerator[:, 1, 1] = ff * gg + fg**2
    numerator[:, 1, 2] = numerator[:, 2, 1] = 2 * gg * fg
    numerator[:, 2, 2] = 2 * gg**2
    return numerator / modes[:, None, None]


def independent_information(covariance, jacobian):
    solved = np.linalg.solve(covariance, jacobian)
    joint = np.einsum("nsa,nsb->ab", jacobian, solved)
    individual = np.array(
        [
            np.einsum(
                "na,nb,n->ab",
                jacobian[:, i],
                jacobian[:, i],
                1 / covariance[:, i, i],
            )
            for i in range(3)
        ]
    )
    return joint, individual


def run(result):
    result["stage"] = "source and input identity"
    assert all(digest(INPUT / name) == sha for name, sha in HASHES.items())
    result["input_sha256"] = dict(HASHES)
    result["source_sha256"] = {name: digest(ROOT / name) for name in SOURCES}
    report = json.loads(
        (INPUT / "profiles-checked/records-001.report.json").read_text()
    )
    diagnosis = json.loads((INPUT / "weight-diagnosis/diagnosis.json").read_text())
    context, settings = report["context"], report["settings"]
    assert context["profile"] == settings["profile"] == "accuracy"
    assert context["bin"] == 0 and context["bounds"] == [2.0, 2.235]
    assert context["parameters"] == ["ap_0", "at_0"]
    assert report["final_controls"]["magnitude_order"] == 32
    assert report["final_controls"]["iterations"] == 6
    field_index = next(
        i for i, field in enumerate(context["fields"]) if field["id"] == "lya(qso)"
    )
    galaxy_index = next(
        i for i, field in enumerate(context["fields"]) if field["id"] == "qso"
    )
    pair_values = [
        [field_index, field_index],
        [field_index, galaxy_index],
        [galaxy_index, galaxy_index],
    ]
    selected_indices = [context["selected_pairs"].index(pair) for pair in pair_values]
    required_indices = [context["required_pairs"].index(pair) for pair in pair_values]
    assert [
        [0 if i == field_index else 1, 0 if j == field_index else 1]
        for i, j in pair_values
    ] == [[0, 0], [0, 1], [1, 1]]
    result["selection"] = dict(
        profile="accuracy",
        bin=0,
        bounds=context["bounds"],
        fields=["lya(qso)", "qso"],
        pairs=["FF", "FG", "GG"],
        selected_indices=selected_indices,
        required_indices=required_indices,
        parameters=context["parameters"],
        magnitude_orders=[32, 64],
    )
    old_modules = report["provenance"]["fishhighz"]["module_hashes"]
    for name in UNCHANGED_HISTORICAL_SOURCES:
        assert digest(ROOT / name) == old_modules[name]
    assert digest(ROOT / "fishhighz/weights.py") == REVIEWED_W09_WEIGHT_SOURCE
    result["historical_source_identity"] = {
        name: old_modules[name] for name in UNCHANGED_HISTORICAL_SOURCES
    }
    result["reviewed_W09_weight_source"] = REVIEWED_W09_WEIGHT_SOURCE

    with np.load(
        INPUT / "profiles-checked/records-001.npz", allow_pickle=False
    ) as primary_file:
        primary = {name: primary_file[name] for name in primary_file.files}
    assert validate_payload(context, primary, report, require_pass=False, schema=3)
    assert report["passed"] is False
    result["historical_record"] = dict(
        validated=True,
        passed=report["passed"],
        trial_contract_version=report["trial_contract"]["version"],
        interpretation="historical legacy cumulative result; unresolved status unchanged",
    )

    (bin_row,) = [row for row in diagnosis["bins"] if row["bin"] == 0]
    assert bin_row["arrays"] == "bin-0.npz"
    assert bin_row["sha256"] == HASHES["weight-diagnosis/bin-0.npz"]
    with np.load(INPUT / "weight-diagnosis/bin-0.npz", allow_pickle=False) as file:
        weight_arrays = {name: file[name] for name in file.files}

    sample = settings["samples"]["lya(qso)"]
    z_eval = settings["model"]["z_eval"]
    a_v, d_deg = settings["geometry"]["a_v"], settings["geometry"]["d_deg"]
    recipe = RECIPES[context["case"]]
    resolution = float(recipe["survey"]["resolution"])
    pixel = (
        SPEED_LIGHT_KMS
        * sample["snr_diagnostics"]["policies"]["pixel_width_angstrom"]
        / sample["snr_diagnostics"]["original"]["wavelength"]
    )
    sigma = SPEED_LIGHT_KMS / (resolution * 2 * math.sqrt(2 * math.log(2)))
    response = InstrumentResponse(pixel, sigma)
    geometry = prepare_geometry(
        *context["bounds"],
        z_eval=z_eval,
        area_deg2=float(recipe["survey"]["survey_area"]),
        h_fid=settings["grid"]["h_fid"],
        hubble=lambda z: np.full_like(
            z, a_v * settings["grid"]["h_fid"] * (1 + z_eval)
        ),
        transverse_distance=lambda z: np.full_like(
            z, d_deg / (settings["grid"]["h_fid"] * np.pi / 180)
        ),
        z_order=1,
    )
    check(result, "geometry", [geometry.a_v, geometry.d_deg], [a_v, d_deg])
    field = ObservedField("lya(qso)", "forest", "lya", "qso")
    galaxy = ObservedField("qso", "galaxy", "qso")
    selection = PairSelection([field, galaxy], [(0, 0), (0, 1), (1, 1)])
    registry = ParameterRegistry(
        [
            Parameter("ap_0", 1.0, "target"),
            Parameter("at_0", 1.0, "target"),
        ]
    )
    k, mu, modes = primary["k"], primary["mu"], primary["modes"]
    signal = primary["observed_signal"][:, required_indices]
    jacobian = primary["observed_j"][:, selected_indices]
    assert signal.shape == (16384, 3) and jacobian.shape == (16384, 3, 2)
    assert np.all(primary["noise"][:, required_indices[1]] == 0)
    galaxy_noise = primary["noise"][:, required_indices[2]]
    q_mode = k * mu / a_v
    physical_p1d = default_p1d([], z_eval, q_mode)
    physical_response = velocity_response(
        q_mode,
        pixel_width_velocity=pixel,
        gaussian_sigma_velocity=sigma,
    )
    check(
        result, "saved_response", physical_response, primary["response"][:, field_index]
    )

    result["orders"] = {}
    calculations = {}
    for order in (32, 64):
        prefix = f"order{order}_field{field_index}_"
        magnitudes = weight_arrays[prefix + "magnitudes"]
        quadrature = weight_arrays[prefix + "measure"]
        masses = weight_arrays[prefix + "masses"]
        variance = weight_arrays[prefix + "variance"]
        rho = masses / quadrature
        density = rho * SPEED_LIGHT_KMS / (1 + sample["z_source"])
        if order == 32:
            check(result, "order32_magnitudes_anchor", magnitudes, sample["magnitudes"])
            check(result, "order32_quadrature_anchor", quadrature, sample["quadrature"])
            check(result, "order32_variance_anchor", variance, sample["variance"])
            check(
                result,
                "order32_masses_anchor",
                masses,
                np.asarray(sample["density"])
                * (1 + sample["z_source"])
                / SPEED_LIGHT_KMS
                * quadrature,
            )
        source, weighting = _forest_input(
            field,
            dict(
                z_source=sample["z_source"],
                density=density,
                variance=variance,
                length_velocity=sample["length_velocity"],
                density_diagnostics={"saved_order": order},
                snr_diagnostics={"saved_order": order},
            ),
            magnitudes,
            quadrature,
            response,
            registry,
            z_eval,
            method="inverse_variance",
        )
        assert source.auxiliary_coordinates is None
        assert source.weight_options["method"] == "inverse_variance"
        assert "iterations" not in source.weight_options
        prepared = prepare_forest_weights(
            field, geometry, response, **dict(source.weight_options)
        )
        check(
            result, f"order{order}_masses", prepared.rho * prepared.quadrature, masses
        )
        b_star = weighting["reference"]["B_star"]
        direct_b = float(
            default_p1d([], z_eval, [FIXED_REFERENCE["q_star"]])[0]
            * velocity_response(
                [FIXED_REFERENCE["q_star"]],
                pixel_width_velocity=pixel,
                gaussian_sigma_velocity=sigma,
            )[0]
            ** 2
        )
        check(result, f"order{order}_B_star", b_star, direct_b)
        reference_weights = np.where(
            masses > 0, b_star / (b_star + pixel * variance), 0.0
        )
        check(result, f"order{order}_weights", prepared.weights, reference_weights)
        i1 = math.fsum(float(r * w) for r, w in zip(masses, reference_weights))
        i2 = math.fsum(float(r * w * w) for r, w in zip(masses, reference_weights))
        i3 = math.fsum(
            float(r * w * w * v) for r, w, v in zip(masses, reference_weights, variance)
        )
        independent_coefficients = np.array(
            [
                i2 / (sample["length_velocity"] * i1**2),
                pixel * i3 / (sample["length_velocity"] * i1**2),
            ]
        )
        check(
            result,
            f"order{order}_coefficients",
            [prepared.A, prepared.P_pixel],
            independent_coefficients,
        )
        noise_f = forest_noise(prepared, field, geometry, response, k, mu, physical_p1d)
        independent_noise = (
            (prepared.A * physical_p1d * physical_response**2 + prepared.P_pixel)
            * d_deg**2
            / a_v
        )
        check(result, f"order{order}_forest_noise", noise_f.total, independent_noise)
        noise = prepare_noise(
            selection,
            len(k),
            diagonal={field.id: noise_f.total, galaxy.id: galaxy_noise},
            independent_sampling=True,
        )
        assert np.all(noise[:, 1] == 0)
        total = combine_observed_power(signal, noise)
        covariance = gaussian_covariance(total, modes, selection)
        direct_covariance = independent_covariance(total, modes)
        check(result, f"order{order}_covariance", covariance, direct_covariance)
        check(
            result,
            f"order{order}_cross_variance",
            covariance[:, 1, 1] * modes,
            total[:, 0] * total[:, 2] + total[:, 1] ** 2,
        )
        public_joint = fisher_matrix(jacobian, covariance)
        public_individual = np.array(
            [
                fisher_matrix(
                    jacobian[:, i : i + 1], covariance[:, i : i + 1, i : i + 1]
                )
                for i in range(3)
            ]
        )
        direct_joint, direct_individual = independent_information(
            direct_covariance, jacobian
        )
        check(result, f"order{order}_joint_fisher", public_joint, direct_joint)
        check(
            result,
            f"order{order}_individual_fisher",
            public_individual,
            direct_individual,
        )
        joint = invert(public_joint)
        individual = [invert(matrix) for matrix in public_individual]
        calculations[order] = dict(
            coefficients=np.array([prepared.A, prepared.P_pixel]),
            joint_errors=np.array(joint["errors"]),
            individual_errors=np.array([row["errors"] for row in individual]),
        )
        result["orders"][str(order)] = dict(
            forest_input=dict(
                method=weighting["method"],
                iterations=weighting["iterations"],
                reference=weighting["reference"],
            ),
            coefficients=dict(A=prepared.A, P_pixel=prepared.P_pixel),
            joint=joint,
            individual=dict(zip(("FF", "FG", "GG"), individual)),
        )

    check(
        result,
        "W09_FF_order32_errors",
        calculations[32]["individual_errors"][0],
        [0.029619522112556678, 0.027975323168394890],
    )
    ratios = {
        "A": calculations[64]["coefficients"][0] / calculations[32]["coefficients"][0]
        - 1,
        "P_pixel": calculations[64]["coefficients"][1]
        / calculations[32]["coefficients"][1]
        - 1,
        "joint_errors": calculations[64]["joint_errors"]
        / calculations[32]["joint_errors"]
        - 1,
        "individual_errors": calculations[64]["individual_errors"]
        / calculations[32]["individual_errors"]
        - 1,
    }
    assert all(np.all(np.abs(value) <= 0.001) for value in ratios.values())
    result["refinement"] = {
        name: dict(
            order32=np.asarray(
                calculations[32]["coefficients" if name in ("A", "P_pixel") else name]
            ).tolist(),
            order64=np.asarray(
                calculations[64]["coefficients" if name in ("A", "P_pixel") else name]
            ).tolist(),
            ratio_minus_one=np.asarray(value).tolist(),
        )
        for name, value in ratios.items()
    }
    result["refinement"]["A"].update(
        order32=calculations[32]["coefficients"][0],
        order64=calculations[64]["coefficients"][0],
    )
    result["refinement"]["P_pixel"].update(
        order32=calculations[32]["coefficients"][1],
        order64=calculations[64]["coefficients"][1],
    )
    assert all(digest(INPUT / name) == sha for name, sha in HASHES.items())
    assert all(
        digest(ROOT / name) == sha for name, sha in result["source_sha256"].items()
    )
    result["limitations"] = [
        "one saved bin and the FF/FG/GG two-field subset only",
        "saved signal, Jacobian, Fourier nodes, modes and galaxy noise held fixed",
        "synthetic/controller and local saved-input checks do not establish full-survey convergence",
        "no P3D weighting query, auxiliary coordinates, real-model controller or full forecast",
    ]
    result["status"] = "passed; ready for independent and user review"


def stop(signum, frame):
    raise TimeoutError("W12 saved-input calculation exceeded its 60-second cap")


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, stop)
    signal.alarm(60)
    started = time.monotonic()
    threads = {
        name: os.environ.get(name)
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    }
    assert set(threads.values()) == {"1"}
    assert os.environ.get("FISHHIGHZ_FISHER_BACKEND", "numpy") == "numpy"
    output = (
        ROOT
        / ".validation/forest-weight-diagnostics"
        / datetime.now(timezone.utc).strftime("w12-r1-%Y%m%dT%H%M%S%fZ")
    )
    output.mkdir(parents=True, exist_ok=False)
    result = dict(
        status="incomplete",
        stage="startup",
        residuals={},
        threads=threads,
        backend="numpy",
        python=platform.python_version(),
        numpy=np.__version__,
        host=platform.node(),
        command=(
            "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
            "timeout 60 .venv/bin/python -B scripts/check_accuracy_fixed_reference.py"
        ),
        git_head=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
    )
    (output / "git-status.txt").write_text(
        subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True)
    )
    try:
        run(result)
    except Exception as error:
        result.update(status="stopped", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        signal.alarm(0)
        result["elapsed_seconds"] = time.monotonic() - started
        (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(output, result["status"], flush=True)
