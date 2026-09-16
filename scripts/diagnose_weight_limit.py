"""Bounded saved-input diagnosis of the cumulative forest-weight limit.

The script reads only the preserved Step 12 arrays.  It has three explicit
phases: initialize an exclusive output directory, run one population/order
batch, and finalize exactly sixty batches into tables and SVG figures.
"""

import argparse
import csv
import hashlib
import json
import platform
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS
from fishhighz.kernels.weights import _integrals
from fishhighz.response import pixel_width_angstrom_to_velocity
from fishhighz.validation.weight_limit import (
    batch_payload,
    relative_log_change,
    signed_log_ratio,
    trajectory,
    validate_batch_payload,
)

CASE = "lya_qso_lbg_lae_15x2pt"
FIELDS = ("lya(qso)", "lya(lbg)")
ORDERS = (4, 8, 16, 32, 64)
CHECKPOINTS = (0, 3, 6, 12, 24, 48, 96, 128, 192, 256, 384, 512, 768, 1024)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_json(path):
    return json.loads(Path(path).read_text())


def _json_default(item):
    if isinstance(item, np.generic):
        return item.item()
    raise TypeError(f"cannot serialize {type(item).__name__}")


def _write_json(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, allow_nan=False, default=_json_default) + "\n"
    )


def _context(source, profiles):
    source = Path(source).resolve()
    profiles = Path(profiles).resolve()
    diagnosis_path = source / "diagnosis.json"
    manifest_path = profiles / "manifest.json"
    diagnosis = _read_json(diagnosis_path)
    manifest = _read_json(manifest_path)
    bins = diagnosis.get("bins", [])
    if [row.get("bin") for row in bins] != list(range(6)):
        raise ValueError("weight diagnosis must contain exactly bins 0 through 5")
    records = [
        row for row in manifest.get("records", []) if row["task"]["case"] == CASE
    ]
    expected = {
        (profile, index)
        for profile in ("accuracy", "compatibility")
        for index in range(6)
    }
    actual = {(row["task"]["profile"], row["task"]["bin"]) for row in records}
    if actual != expected or len(records) != 12:
        raise ValueError("profile bundle must contain six bins in both profiles")
    for row in records:
        arrays = profiles / row["arrays"]
        report = profiles / row["report_file"]
        if digest(arrays) != row["sha256"] or digest(report) != row["report_sha256"]:
            raise ValueError("profile payload hash mismatch")
    for row in bins:
        path = source / row["arrays"]
        if digest(path) != row["sha256"]:
            raise ValueError("weight-diagnosis array hash mismatch")
    return source, profiles, diagnosis, manifest, records


def initialize(source, profiles, output):
    """Create a unique output directory and bind the immutable input context."""
    source, profiles, diagnosis, manifest, records = _context(source, profiles)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    historical = {
        profile: sum(
            row["scientific_passed"]
            for row in records
            if row["task"]["profile"] == profile
        )
        for profile in ("accuracy", "compatibility")
    }
    diagnostic_rows = manifest.get("diagnostics", [])
    expected_batches = []
    profile_matches = []
    for bin_index in range(6):
        for field in FIELDS:
            for order in ORDERS:
                inputs = _field_inputs(
                    source, diagnosis, profiles, records, bin_index, field, order
                )
                expectation = _batch_expectation(
                    inputs, bin_index=bin_index, field=field, order=order
                )
                expected_batches.append(_context_expectation(expectation))
                if inputs["source_match"] is not None:
                    profile_matches.append(
                        dict(bin=bin_index, field=field, order=order)
                    )
    context = dict(
        schema=2,
        kind="step13_weight_limit_context",
        case=CASE,
        source=str(source),
        profiles=str(profiles),
        source_diagnosis_sha256=digest(source / "diagnosis.json"),
        source_manifest_sha256=diagnosis["source_sha256"],
        profiles_manifest_sha256=digest(profiles / "manifest.json"),
        source_arrays={row["arrays"]: row["sha256"] for row in diagnosis["bins"]},
        profile_payloads={row["arrays"]: row["sha256"] for row in records},
        historical_profile_verdicts=dict(
            compatibility_passed=historical["compatibility"],
            accuracy_passed=historical["accuracy"],
            records=len(records),
            diagnostic_requested=len(diagnostic_rows),
            diagnostic_completed=sum(
                row.get("status") == "completed" for row in diagnostic_rows
            ),
            diagnostic_unavailable=sum(
                row.get("status") != "completed" for row in diagnostic_rows
            ),
            interpretation="preserved Step 12 outcomes, not Step 13 results",
        ),
        expected_batches=expected_batches,
        profile_source_matches=profile_matches,
        checkpoints=list(CHECKPOINTS),
        interpreter=sys.executable,
        python=platform.python_version(),
        numpy=np.__version__,
        backend="NumPy logaddexp; one thread required by invocation policy",
    )
    _write_json(output / "context.json", context)
    (output / "batches").mkdir()
    return context


def _verified_context(source, profiles, output, *, verify_expected=True):
    source, profiles, diagnosis, manifest, records = _context(source, profiles)
    output = Path(output).resolve()
    context = _read_json(output / "context.json")
    if (
        context.get("schema") != 2
        or context.get("kind") != "step13_weight_limit_context"
    ):
        raise ValueError("unsupported Step 13 context")
    checks = {
        "source_diagnosis_sha256": digest(source / "diagnosis.json"),
        "source_manifest_sha256": diagnosis["source_sha256"],
        "profiles_manifest_sha256": digest(profiles / "manifest.json"),
        "source_arrays": {row["arrays"]: row["sha256"] for row in diagnosis["bins"]},
        "profile_payloads": {row["arrays"]: row["sha256"] for row in records},
        "checkpoints": list(CHECKPOINTS),
    }
    for name, value in checks.items():
        if context.get(name) != value:
            raise ValueError(f"{name} changed after initialization")
    if context.get("source") != str(source) or context.get("profiles") != str(profiles):
        raise ValueError("input paths differ from initialized context")
    if verify_expected:
        rebuilt = []
        profile_matches = []
        for bin_index in range(6):
            for field in FIELDS:
                for order in ORDERS:
                    inputs = _field_inputs(
                        source, diagnosis, profiles, records, bin_index, field, order
                    )
                    expectation = _batch_expectation(
                        inputs, bin_index=bin_index, field=field, order=order
                    )
                    rebuilt.append(_context_expectation(expectation))
                    if inputs["source_match"] is not None:
                        profile_matches.append(
                            dict(bin=bin_index, field=field, order=order)
                        )
        if context.get("expected_batches") != rebuilt:
            raise ValueError("expected batch context disagrees with immutable sources")
        if context.get("profile_source_matches") != profile_matches:
            raise ValueError("profile/source match inventory disagrees with sources")
    diagnostic_rows = manifest.get("diagnostics", [])
    historical = {
        profile: sum(
            row["scientific_passed"]
            for row in records
            if row["task"]["profile"] == profile
        )
        for profile in ("accuracy", "compatibility")
    }
    expected_historical = dict(
        compatibility_passed=historical["compatibility"],
        accuracy_passed=historical["accuracy"],
        records=len(records),
        diagnostic_requested=len(diagnostic_rows),
        diagnostic_completed=sum(
            row.get("status") == "completed" for row in diagnostic_rows
        ),
        diagnostic_unavailable=sum(
            row.get("status") != "completed" for row in diagnostic_rows
        ),
        interpretation="preserved Step 12 outcomes, not Step 13 results",
    )
    if context.get("historical_profile_verdicts") != expected_historical:
        raise ValueError("historical profile verdicts disagree with immutable sources")
    return source, profiles, diagnosis, manifest, records, output, context


def _accuracy_report(profiles, records, bin_index):
    row = next(
        row
        for row in records
        if row["task"]["profile"] == "accuracy" and row["task"]["bin"] == bin_index
    )
    report = _read_json(profiles / row["report_file"])
    if (
        report["context"]["bin"] != bin_index
        or report["context"]["profile"] != "accuracy"
    ):
        raise ValueError("profile report context mismatch")
    return row, report


def _field_inputs(source, diagnosis, profiles, records, bin_index, field, order):
    if field not in FIELDS or order not in ORDERS or bin_index not in range(6):
        raise ValueError("batch must select one expected bin, forest field and order")
    bin_row = diagnosis["bins"][bin_index]
    profile_row, profile = _accuracy_report(profiles, records, bin_index)
    fields = profile_row["task"]["fields"]
    field_index = next(i for i, item in enumerate(fields) if item["id"] == field)
    prefix = f"order{order}_field{field_index}_"
    arrays_path = source / bin_row["arrays"]
    with np.load(arrays_path, allow_pickle=False) as saved:
        names = ("magnitudes", "measure", "masses", "variance")
        arrays = {name: saved[prefix + name] for name in names}
        historical_weights = {
            count: saved[prefix + f"weights_{count}"]
            for count in (3, 6, 12, 24)
            if prefix + f"weights_{count}" in saved.files
        }
    candidates = [row for row in bin_row["rows"] if row["order"] == order]
    field_metadata = next(
        item for item in candidates[0]["fields"] if item["field"] == field
    )
    sample = profile["settings"]["samples"][field]
    policies = sample["snr_diagnostics"]["policies"]
    z_eval = profile["settings"]["model"]["z_eval"]
    pixel = pixel_width_angstrom_to_velocity(
        policies["pixel_width_angstrom"],
        lambda_obs_angstrom=LYA_REST_ANGSTROM * (1 + z_eval),
    )
    length = sample["length_velocity"]
    final_order = profile["settings"]["controls"]["magnitude_order"]
    source_match = None
    if order == final_order:
        expected_mass = (
            np.asarray(sample["density"])
            * ((1 + sample["z_source"]) / SPEED_LIGHT_KMS)
            * np.asarray(sample["quadrature"])
        )
        comparisons = {
            "magnitudes": np.asarray(sample["magnitudes"]),
            "measure": np.asarray(sample["quadrature"]),
            "masses": expected_mass,
            "variance": np.asarray(sample["variance"]),
        }
        for name, expected_values in comparisons.items():
            np.testing.assert_allclose(
                arrays[name], expected_values, rtol=5e-12, atol=0
            )
        source_match = dict(
            profile_arrays=profile_row["arrays"],
            profile_sha256=profile_row["sha256"],
            magnitude_order=final_order,
            matched=True,
            tolerance=5e-12,
        )
    return dict(
        **arrays,
        historical_weights=historical_weights,
        historical_rows=field_metadata["rows"],
        field_index=field_index,
        signal=field_metadata["signal"],
        alias=field_metadata["alias"],
        length=length,
        pixel=pixel,
        source_match=source_match,
        arrays_path=arrays_path,
        arrays_sha256=bin_row["sha256"],
        profile_report=profile_row["report_file"],
        profile_report_sha256=profile_row["report_sha256"],
    )


def _input_identity(inputs):
    digest_value = hashlib.sha256()
    for name in (
        "magnitudes",
        "masses",
        "variance",
        "length",
        "pixel",
        "signal",
        "alias",
    ):
        array = np.ascontiguousarray(np.atleast_1d(inputs[name]))
        digest_value.update(name.encode())
        digest_value.update(array.dtype.str.encode())
        digest_value.update(str(array.shape).encode())
        digest_value.update(array.tobytes())
    return digest_value.hexdigest()


def _source_record(inputs):
    return dict(
        arrays=inputs["arrays_path"].name,
        sha256=inputs["arrays_sha256"],
        profile_report=inputs["profile_report"],
        profile_report_sha256=inputs["profile_report_sha256"],
        source_match=inputs["source_match"],
    )


def _batch_expectation(inputs, *, bin_index, field, order):
    source_inputs = {
        name: np.asarray(np.atleast_1d(inputs[name]), dtype=np.float64)
        for name in (
            "magnitudes",
            "masses",
            "variance",
            "length",
            "pixel",
            "signal",
            "alias",
        )
    }
    return dict(
        bin=bin_index,
        field=field,
        order=order,
        inputs=source_inputs,
        input_sha256=_input_identity(source_inputs),
        source=_source_record(inputs),
        attempts=[dict(iterations=count, outcome="completed") for count in CHECKPOINTS],
        capped_at=CHECKPOINTS[-1],
    )


def _context_expectation(expectation):
    return dict(
        bin=expectation["bin"],
        field=expectation["field"],
        order=expectation["order"],
        input_sha256=expectation["input_sha256"],
        source=expectation["source"],
        attempts=expectation["attempts"],
        capped_at=expectation["capped_at"],
    )


def _historical_comparison(inputs, result):
    states = {state.iterations: state for state in result.snapshots}
    rows = []
    for historical in inputs["historical_rows"]:
        count = historical["iterations"]
        state = states[count]
        row = dict(iterations=count, historical_available=historical["available"])
        if historical["available"]:
            row.update(
                A_relative=abs(state.coefficients.A / historical["A"] - 1),
                P_pixel_relative=abs(
                    state.coefficients.P_pixel / historical["P_pixel"] - 1
                ),
            )
        else:
            row["historical_failure"] = {
                key: value
                for key, value in historical.items()
                if key not in ("iterations", "available")
            }
        if count in inputs["historical_weights"]:
            old = inputs["historical_weights"][count]
            support = inputs["masses"] > 0
            positive = support & (old > 0)
            row["saved_positive_support"] = int(np.count_nonzero(positive))
            if np.all(positive == support):
                row["log_weight_max_absolute"] = float(
                    np.max(
                        np.abs(
                            np.log(old[positive])
                            - state.log_amplitude
                            - state.log_shape[positive]
                        )
                    )
                )
            try:
                _integrals(
                    inputs["masses"],
                    old,
                    inputs["variance"],
                    inputs["length"],
                    inputs["pixel"],
                )
                row["production_integrals"] = "available"
            except (FloatingPointError, ValueError) as error:
                row["production_integrals"] = "rejected"
                row["production_error"] = str(error)
        rows.append(row)
    return rows


def _finite_trajectory_comparison(result, report):
    state = result.snapshots[-1]
    limit_state = report["linearized"]["fixed_grid_limit"]
    return dict(
        iterations=state.iterations,
        A_relative=relative_log_change(state.coefficients.log_A, limit_state["log_A"]),
        P_pixel_relative=relative_log_change(
            state.coefficients.log_P_pixel, limit_state["log_P_pixel"]
        ),
        interpretation=(
            "finite nonlinear checkpoint versus the exact linearized fixed-grid "
            "shape; this is not an asymptotic error estimate"
        ),
    )


def _semantic_close(actual, expected, name):
    """Compare derived report content quantity by quantity."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f"{name}: inconsistent derived schema")
        for key, value in expected.items():
            _semantic_close(actual[key], value, f"{name} {key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(f"{name}: inconsistent derived inventory")
        for index, value in enumerate(expected):
            _semantic_close(actual[index], value, f"{name} {index}")
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if actual != expected:
            raise ValueError(f"{name}: inconsistent derived content")
    elif not np.isclose(actual, expected, rtol=5e-12, atol=0, equal_nan=False):
        raise ValueError(f"{name}: inconsistent derived numerical content")


def run_batch(source, profiles, output, bin_index, field, order):
    """Run and save exactly one population/order batch without overwriting."""
    source, profiles, diagnosis, _, records, output, context = _verified_context(
        source, profiles, output, verify_expected=False
    )
    inputs = _field_inputs(
        source, diagnosis, profiles, records, bin_index, field, order
    )
    expected = _batch_expectation(inputs, bin_index=bin_index, field=field, order=order)
    context_expected = next(
        row
        for row in context["expected_batches"]
        if (row["bin"], row["field"], row["order"]) == (bin_index, field, order)
    )
    if context_expected != _context_expectation(expected):
        raise ValueError("batch expectation disagrees with immutable sources")
    started = time.perf_counter()
    result = trajectory(
        inputs["masses"],
        inputs["variance"],
        inputs["length"],
        inputs["pixel"],
        inputs["signal"],
        inputs["alias"],
        checkpoints=CHECKPOINTS,
        coordinates=inputs["magnitudes"],
    )
    elapsed = time.perf_counter() - started
    arrays, report = batch_payload(
        result,
        magnitudes=inputs["magnitudes"],
        masses=inputs["masses"],
        variance=inputs["variance"],
        length=inputs["length"],
        pixel=inputs["pixel"],
        signal=inputs["signal"],
        alias=inputs["alias"],
        bin_index=bin_index,
        field=field,
        order=order,
        source=expected["source"],
    )
    report["historical_comparison"] = _historical_comparison(inputs, result)
    report["trajectory_to_fixed_grid_limit"] = _finite_trajectory_comparison(
        result, report
    )
    report["seconds"] = elapsed
    report["within_batch_wall_target"] = elapsed <= 30
    slug = f"bin-{bin_index}-{field[4:-1]}-order-{order}"
    array_path = output / "batches" / f"{slug}.npz"
    report_path = output / "batches" / f"{slug}.json"
    if array_path.exists() or report_path.exists():
        raise FileExistsError(f"refusing to overwrite batch {slug}")
    np.savez_compressed(array_path, **arrays)
    report["arrays"] = array_path.name
    report["arrays_sha256"] = digest(array_path)
    _write_json(report_path, report)
    with np.load(array_path, allow_pickle=False) as saved:
        validate_batch_payload(
            {name: saved[name] for name in saved.files},
            report,
            expected=expected,
        )
    return report


def _grid_screen(reports, bin_index, field):
    selected = {
        report["order"]: report
        for report in reports
        if report["bin"] == bin_index and report["field"] == field
    }
    finite_cap_comparisons = []
    asymptotic_comparisons = []
    for lower, upper in ((16, 32), (32, 64)):
        first = selected[lower]["checkpoints"][-1]
        second = selected[upper]["checkpoints"][-1]
        finite_cap_comparisons.append(
            dict(
                orders=[lower, upper],
                A_signed_ratio=signed_log_ratio(first["log_A"], second["log_A"]),
                A_relative=relative_log_change(first["log_A"], second["log_A"]),
                P_pixel_signed_ratio=signed_log_ratio(
                    first["log_P_pixel"], second["log_P_pixel"]
                ),
                P_pixel_relative=relative_log_change(
                    first["log_P_pixel"], second["log_P_pixel"]
                ),
            )
        )
        first_limit = selected[lower]["linearized"]["fixed_grid_limit"]
        second_limit = selected[upper]["linearized"]["fixed_grid_limit"]
        asymptotic_comparisons.append(
            dict(
                orders=[lower, upper],
                A_signed_ratio=signed_log_ratio(
                    first_limit["log_A"], second_limit["log_A"]
                ),
                A_relative=relative_log_change(
                    first_limit["log_A"], second_limit["log_A"]
                ),
                P_pixel_signed_ratio=signed_log_ratio(
                    first_limit["log_P_pixel"], second_limit["log_P_pixel"]
                ),
                P_pixel_relative=relative_log_change(
                    first_limit["log_P_pixel"], second_limit["log_P_pixel"]
                ),
                effective_measure_ratio=np.exp(
                    second_limit["log_effective_measure"]
                    - first_limit["log_effective_measure"]
                ),
                dominant_coordinate_shift=abs(
                    selected[upper]["linearized"]["dominant_coordinate"]
                    - selected[lower]["linearized"]["dominant_coordinate"]
                ),
            )
        )
    finite_cap_grid_passed = all(
        row["A_relative"] <= 1e-4 and row["P_pixel_relative"] <= 1e-4
        for row in finite_cap_comparisons
    )
    asymptotic_grid_passed = all(
        row["A_relative"] <= 1e-4 and row["P_pixel_relative"] <= 1e-4
        for row in asymptotic_comparisons
    )
    iteration_passed = all(
        selected[order]["fixed_mesh_screen"]["passed"] for order in (16, 32, 64)
    )
    rising_pixel_scaling = all(
        row["P_pixel_signed_ratio"] > 0.8 and row["effective_measure_ratio"] < 1
        for row in asymptotic_comparisons
    )
    lower_limit = selected[16]["linearized"]["fixed_grid_limit"]
    upper_limit = selected[64]["linearized"]["fixed_grid_limit"]
    return dict(
        bin=bin_index,
        field=field,
        orders=[16, 32, 64],
        finite_cap_comparisons=finite_cap_comparisons,
        fixed_grid_asymptotic_limits=[
            dict(
                order=order,
                dominant_coordinate=selected[order]["linearized"][
                    "dominant_coordinate"
                ],
                **selected[order]["linearized"]["fixed_grid_limit"],
            )
            for order in (16, 32, 64)
        ],
        asymptotic_comparisons=asymptotic_comparisons,
        iteration_passed=iteration_passed,
        finite_cap_grid_passed=finite_cap_grid_passed,
        asymptotic_grid_passed=asymptotic_grid_passed,
        rising_pixel_scaling=rising_pixel_scaling,
        P_pixel_order_exponent=(
            (upper_limit["log_P_pixel"] - lower_limit["log_P_pixel"]) / np.log(64 / 16)
        ),
        effective_measure_order_exponent=(
            (
                upper_limit["log_effective_measure"]
                - lower_limit["log_effective_measure"]
            )
            / np.log(64 / 16)
        ),
        candidate_finite_limit=False,
        outcome="unresolved_within_stated_bounds",
        continuum_assessment=dict(
            measure="dR = rho(m) dm in deg^-2 (km/s)^-1",
            support="bounded increasing physical magnitude interval",
            variance="strictly positive at every saved supported node",
            limit_ordering=(
                "iteration to infinity at each fixed grid, followed by a "
                "continuum refinement; no joint-limit interchange is assumed"
            ),
            fixed_grid_iteration_limit_established=True,
            nested_refinement_family_established=False,
            non_atomic_measure_limit_established=False,
            uniform_nonlinear_remainder_control_established=False,
            tail_diagnostics=[
                dict(
                    order=order,
                    mean_coordinate=selected[order]["linearized"]["fixed_grid_limit"][
                        "mean_coordinate"
                    ],
                    sigma_coordinate=selected[order]["linearized"]["fixed_grid_limit"][
                        "sigma_coordinate"
                    ],
                    **selected[order]["fixed_grid_tail"],
                )
                for order in (16, 32, 64)
            ],
            conclusion=(
                "the sufficient continuum concentration assumptions cannot be "
                "verified from the five finite, non-nested saved quadratures"
            ),
        ),
        basis=(
            "strong signed empirical refinement trend over orders 16/32/64; "
            "finite saved refinements do not establish either a continuum "
            "divergence theorem or a finite limit"
        ),
    )


def _polyline(values, x0, y0, width, height, x_limits, y_limits):
    x_min, x_max = x_limits
    y_min, y_max = y_limits
    points = []
    for x, y in values:
        px = x0 + width * (x - x_min) / (x_max - x_min)
        py = y0 + height * (1 - (y - y_min) / (y_max - y_min))
        points.append(f"{px:.2f},{py:.2f}")
    return " ".join(points)


def _write_svg(path, reports, quantity, title):
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    values = []
    for field_index, field in enumerate(FIELDS):
        for report in reports:
            if report["field"] != field or report["order"] != 64:
                continue
            curve = [
                (np.log2(max(row["iterations"], 1)), row[quantity] / np.log(10))
                for row in report["checkpoints"]
                if row["iterations"] > 0 and row[quantity] is not None
            ]
            values.append((field_index, report["bin"], curve))
    all_y = [y for _, _, curve in values for _, y in curve]
    y_min, y_max = min(all_y), max(all_y)
    if y_min == y_max:
        y_min, y_max = y_min - 0.5, y_max + 0.5
    width, height = 900, 650
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="450" y="28" text-anchor="middle" font-size="18">{title}</text>',
    ]
    for panel, field in enumerate(FIELDS):
        x0, y0, panel_width, panel_height = 70 + panel * 420, 70, 350, 500
        elements.extend(
            [
                f'<rect x="{x0}" y="{y0}" width="{panel_width}" height="{panel_height}" fill="none" stroke="black"/>',
                f'<text x="{x0 + panel_width / 2}" y="{y0 - 12}" text-anchor="middle">{field}, order 64</text>',
                f'<text x="{x0 + panel_width / 2}" y="610" text-anchor="middle">log2(iterations)</text>',
            ]
        )
        for _, bin_index, curve in (row for row in values if row[0] == panel):
            points = _polyline(
                curve, x0, y0, panel_width, panel_height, (0, 10), (y_min, y_max)
            )
            elements.append(
                f'<polyline points="{points}" fill="none" stroke="{colors[bin_index]}" stroke-width="1.5"/>'
            )
    elements.append(
        '<text x="18" y="320" transform="rotate(-90 18 320)" '
        'text-anchor="middle">log10 quantity</text>'
    )
    for index, color in enumerate(colors):
        elements.append(
            f'<text x="{675 + (index % 2) * 70}" y="{590 + (index // 2) * 16}" fill="{color}">bin {index}</text>'
        )
    elements.append("</svg>")
    Path(path).write_text("\n".join(elements) + "\n")


def finalize(source, profiles, output):
    """Validate all sixty batches and write the bounded scientific summary."""
    source, profiles, diagnosis, _, records, output, context = _verified_context(
        source, profiles, output
    )
    final_paths = (
        output / "summary.json",
        output / "decision-table.csv",
        output / "A-concentration.svg",
        output / "P-pixel-concentration.svg",
        output / "effective-measure.svg",
    )
    if any(path.exists() for path in final_paths):
        raise FileExistsError("refusing to overwrite finalized evidence")
    expected = context["expected_batches"]
    json_paths = sorted((output / "batches").glob("*.json"))
    npz_paths = sorted((output / "batches").glob("*.npz"))
    if len(json_paths) != 60 or len(npz_paths) != 60:
        raise ValueError("finalization requires exactly 60 JSON and 60 NPZ batches")
    reports = []
    identities = set()
    for path in json_paths:
        report = _read_json(path)
        identity = (report["bin"], report["field"], report["order"])
        if identity in identities:
            raise ValueError("duplicate batch identity")
        identities.add(identity)
        slug = f"bin-{identity[0]}-{identity[1][4:-1]}-order-{identity[2]}"
        if path.name != f"{slug}.json" or report.get("arrays") != f"{slug}.npz":
            raise ValueError("batch filenames disagree with source identity")
        arrays_path = path.with_name(report["arrays"])
        if digest(arrays_path) != report["arrays_sha256"]:
            raise ValueError("batch array file hash mismatch")
        inputs = _field_inputs(
            source, diagnosis, profiles, records, identity[0], identity[1], identity[2]
        )
        expectation = _batch_expectation(
            inputs, bin_index=identity[0], field=identity[1], order=identity[2]
        )
        with np.load(arrays_path, allow_pickle=False) as saved:
            validate_batch_payload(
                {name: saved[name] for name in saved.files},
                report,
                expected=expectation,
            )
        recomputed = trajectory(
            inputs["masses"],
            inputs["variance"],
            inputs["length"],
            inputs["pixel"],
            inputs["signal"],
            inputs["alias"],
            checkpoints=CHECKPOINTS,
            coordinates=inputs["magnitudes"],
        )
        _semantic_close(
            report.get("historical_comparison"),
            _historical_comparison(inputs, recomputed),
            "historical comparison",
        )
        _semantic_close(
            report.get("trajectory_to_fixed_grid_limit"),
            _finite_trajectory_comparison(recomputed, report),
            "finite trajectory comparison",
        )
        seconds = report.get("seconds")
        if (
            not isinstance(seconds, (int, float))
            or not np.isfinite(seconds)
            or seconds < 0
        ):
            raise ValueError("batch timing is invalid")
        if report.get("within_batch_wall_target") is not (seconds <= 30):
            raise ValueError("batch wall-target verdict is inconsistent")
        reports.append(report)
    if identities != {(row["bin"], row["field"], row["order"]) for row in expected}:
        raise ValueError("batch identities do not cover the expected scope")
    decisions = [
        _grid_screen(reports, bin_index, field)
        for bin_index in range(6)
        for field in FIELDS
    ]
    summary = dict(
        schema=2,
        kind="step13_weight_limit_summary",
        context_sha256=digest(output / "context.json"),
        batches=60,
        population_diagnoses=12,
        checkpoints=list(CHECKPOINTS),
        max_batch_seconds=max(report["seconds"] for report in reports),
        all_batches_within_wall_target=all(
            report["within_batch_wall_target"] for report in reports
        ),
        maximum_historical_A_relative=max(
            row.get("A_relative", 0)
            for report in reports
            for row in report["historical_comparison"]
        ),
        maximum_historical_P_pixel_relative=max(
            row.get("P_pixel_relative", 0)
            for report in reports
            for row in report["historical_comparison"]
        ),
        decisions=decisions,
        finite_limit_candidates=sum(row["candidate_finite_limit"] for row in decisions),
        conclusion=(
            "all twelve saved forest populations remain unresolved in the "
            "continuum limit within the stated bounds"
        ),
        scope="saved forest-weight coefficients only; no forecast reassembly or production adoption",
    )
    _write_json(output / "summary.json", summary)
    with (output / "decision-table.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "bin",
                "field",
                "iteration_passed",
                "finite_cap_grid_passed",
                "asymptotic_grid_passed",
                "rising_pixel_scaling",
                "candidate_finite_limit",
                "outcome",
            ),
        )
        writer.writeheader()
        writer.writerows(
            {name: row[name] for name in writer.fieldnames} for row in decisions
        )
    _write_svg(output / "A-concentration.svg", reports, "log_A", "Forest coefficient A")
    _write_svg(
        output / "P-pixel-concentration.svg",
        reports,
        "log_P_pixel",
        "Forest pixel-noise coefficient P_pixel",
    )
    _write_svg(
        output / "effective-measure.svg",
        reports,
        "log_effective_measure",
        "Effective sightline measure",
    )
    return summary


def check_finalized(source, profiles, output):
    """Reconstruct a finalized bundle and compare every final artifact byte."""
    output = Path(output).resolve()
    final_names = (
        "summary.json",
        "decision-table.csv",
        "A-concentration.svg",
        "P-pixel-concentration.svg",
        "effective-measure.svg",
    )
    if any(not (output / name).is_file() for name in final_names):
        raise ValueError("finalized bundle is incomplete")
    with tempfile.TemporaryDirectory(
        prefix="step13-r2-check-", dir=output.parent
    ) as name:
        rebuilt = Path(name)
        shutil.copy2(output / "context.json", rebuilt / "context.json")
        shutil.copytree(output / "batches", rebuilt / "batches")
        summary = finalize(source, profiles, rebuilt)
        for final_name in final_names:
            if (output / final_name).read_bytes() != (
                rebuilt / final_name
            ).read_bytes():
                raise ValueError(f"finalized {final_name} is inconsistent")
    return dict(
        checked=True,
        batches=summary["batches"],
        population_diagnoses=summary["population_diagnoses"],
        conclusion=summary["conclusion"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--profiles", required=True)
    parser.add_argument("--output", required=True)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--initialize", action="store_true")
    modes.add_argument("--run-batch", action="store_true")
    modes.add_argument("--finalize", action="store_true")
    modes.add_argument("--check-finalized", action="store_true")
    parser.add_argument("--bin", type=int, dest="bin_index")
    parser.add_argument("--field", choices=FIELDS)
    parser.add_argument("--order", type=int, choices=ORDERS)
    args = parser.parse_args()
    if args.initialize:
        result = initialize(args.source, args.profiles, args.output)
    elif args.run_batch:
        if args.bin_index is None or args.field is None or args.order is None:
            parser.error("--run-batch requires --bin, --field and --order")
        result = run_batch(
            args.source,
            args.profiles,
            args.output,
            args.bin_index,
            args.field,
            args.order,
        )
    elif args.finalize:
        result = finalize(args.source, args.profiles, args.output)
    else:
        result = check_finalized(args.source, args.profiles, args.output)
    print(json.dumps(result, indent=2, allow_nan=False, default=_json_default))


if __name__ == "__main__":
    main()
