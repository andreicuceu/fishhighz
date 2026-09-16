#!/usr/bin/env python3
"""Offline W01 replay of the signed legacy forest-weight recurrence.

This diagnostic is intentionally separate from the production weighting API.
It consumes one saved compatibility record, preserves the literal legacy
operation order, and stops after the first independently confirmed coefficient
change.  No survey reader or model controller is called.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import sys
import time
from dataclasses import dataclass
from decimal import Decimal, DecimalException, DivisionByZero, localcontext
from pathlib import Path

import numpy as np

CASE = "lya_qso_lbg_lae_15x2pt"
PROFILE = "compatibility"
POPULATION = "lya(qso)"
PAIR_LABEL = "lya(qso)_lya(qso)"
BIN_INDEX = 0
BIN_BOUNDS = (2.0, 2.235)
CHECKPOINTS = (0, 3, 6, 12, 24)
BASELINE_COUNT = 3
COEFFICIENT_TRIGGER = 1.0e-3
RELATIVE_TOLERANCE = 5.0e-12
EVIDENCE_SCHEMA = "fishhighz-weighting-diagnostic-w01-r2-v1"
MANIFEST_SCHEMA = "fishhighz-weighting-diagnostic-w01-r2-manifest-v1"


class DiagnosticArithmeticError(ValueError):
    """A classified arithmetic failure in the literal signed recurrence."""

    def __init__(self, kind, operation, detail, *, index=None):
        super().__init__(f"{kind}: {operation}: {detail}")
        self.kind = kind
        self.operation = operation
        self.detail = detail
        self.index = index

    def record(self):
        """Return a JSON-compatible failure description."""
        return {
            "kind": self.kind,
            "operation": self.operation,
            "detail": self.detail,
            "index": self.index,
        }


@dataclass(frozen=True)
class LegacyInputs:
    """Exact saved inputs to the finite signed recurrence."""

    magnitudes: np.ndarray
    density: np.ndarray
    variance: np.ndarray
    magnitude_step: float
    forest_length: float
    pixel_width: float
    auxiliary_signal: float
    auxiliary_p1d: float


@dataclass(frozen=True)
class ProjectionInputs:
    """Saved arrays and scalars for one auto-spectrum Fisher projection."""

    k: np.ndarray
    mu: np.ndarray
    modes: np.ndarray
    baseline_total: np.ndarray
    observed_j: np.ndarray
    saved_reference_fisher: np.ndarray
    saved_pair_fisher: np.ndarray
    distance_to_velocity: float
    angle_to_distance: float
    covariance_redshift: float
    resolution_width: float
    selected_pairs: np.ndarray
    required_pairs: np.ndarray
    field_ids: tuple[str, ...]
    selected_column: int
    required_column: int


@dataclass(frozen=True)
class SourceInputs:
    """Verified W01 inputs plus their immutable provenance."""

    legacy: LegacyInputs
    projection: ProjectionInputs
    captured_weights: np.ndarray
    captured_a: float
    captured_p_pixel: float
    provenance: dict


def sha256_file(path):
    """Return the SHA-256 digest of one bounded file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def _require_hash(path, expected, context):
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{context} SHA-256 mismatch: {actual} != {expected}")
    return actual


def _task_matches(task):
    return (
        task.get("case") == CASE
        and task.get("profile") == PROFILE
        and task.get("bin") == BIN_INDEX
        and tuple(task.get("bounds", ())) == BIN_BOUNDS
        and task.get("kind") == "real_bao"
        and task.get("diagnostic_id") is None
    )


def _resolve_auto_columns(field_ids, selected_pairs, required_pairs):
    """Resolve selected-J and required-power columns independently by identity."""
    if tuple(field_ids).count(POPULATION) != 1:
        raise ValueError("expected exactly one lya(qso) field")
    field_index = tuple(field_ids).index(POPULATION)
    target = np.array([field_index, field_index])

    def resolve(pairs, name):
        pairs = np.asarray(pairs)
        match = np.flatnonzero(np.all(pairs == target, axis=1))
        if len(match) != 1:
            raise ValueError(f"expected one {PAIR_LABEL} {name} column")
        return int(match[0])

    return resolve(selected_pairs, "selected"), resolve(required_pairs, "required")


def _verify_live_files(mapping, *, relevant_names=None):
    verified = {}
    for raw_path, expected in mapping.items():
        path = Path(raw_path)
        if relevant_names is not None and path.name not in relevant_names:
            continue
        if not path.is_file():
            raise ValueError(f"recorded source is unavailable: {path}")
        verified[str(path)] = _require_hash(path, expected, "recorded source")
    return verified


def load_source(input_directory):
    """Resolve and verify the unique saved W01 record and its producer metadata."""
    root = Path(input_directory).resolve()
    manifest_path = root / "manifest.json"
    manifest = _read_json(manifest_path)
    if (
        manifest.get("schema") != 3
        or manifest.get("kind") != "fishhighz-validation"
        or manifest.get("payload_kind") != "real_bao"
    ):
        raise ValueError("input is not the required schema-3 real-BAO bundle")
    records = [row for row in manifest.get("records", ()) if _task_matches(row["task"])]
    if len(records) != 1:
        raise ValueError("expected exactly one compatibility bin-0 manifest record")
    record = records[0]
    if record.get("status") != "completed" or not record.get("scientific_passed"):
        raise ValueError("selected compatibility record is not complete and passed")

    report_path = root / record["report_file"]
    arrays_path = root / record["arrays"]
    _require_hash(report_path, record["report_sha256"], "record report")
    _require_hash(arrays_path, record["sha256"], "record arrays")
    report = _read_json(report_path)
    if report.get("context") != record["task"] or not _task_matches(report["context"]):
        raise ValueError("record task/report context mismatch")
    if report["settings"].get("profile") != PROFILE:
        raise ValueError("record settings profile mismatch")

    field_ids = tuple(report["settings"]["fields"])
    task_fields = tuple(row["id"] for row in record["task"]["fields"])
    if field_ids != task_fields:
        raise ValueError("settings/task field identity mismatch")
    pair_row = report["settings"]["pair_inputs"].get(PAIR_LABEL)
    if pair_row is None:
        raise ValueError(f"missing {PAIR_LABEL} saved inputs")

    reference_manifests = [
        (Path(path), digest)
        for path, digest in manifest.get("inputs", {}).items()
        if path.endswith("/reference-01/manifest.json")
    ]
    if len(reference_manifests) != 1:
        raise ValueError("expected one recorded reference manifest")
    reference_manifest_path, reference_manifest_hash = reference_manifests[0]
    _require_hash(
        reference_manifest_path, reference_manifest_hash, "reference manifest"
    )
    reference_manifest = _read_json(reference_manifest_path)
    references = [
        row for row in reference_manifest.get("records", ()) if row.get("case") == CASE
    ]
    if len(references) != 1 or references[0].get("status") != "completed":
        raise ValueError("required reference case is unavailable")
    reference = references[0]
    metadata_path = reference_manifest_path.parent / reference["metadata"]
    metadata_hash = _require_hash(
        metadata_path, reference["sha256"], "reference metadata"
    )
    if metadata_hash != report.get("reference_metadata_sha256"):
        raise ValueError("report does not bind the selected reference metadata")
    metadata = _read_json(metadata_path)
    source_records = [
        row
        for row in metadata.get("records", ())
        if row.get("case") == CASE
        and row.get("bin") == BIN_INDEX
        and tuple(row.get("bounds", ())) == BIN_BOUNDS
    ]
    if len(source_records) != 1:
        raise ValueError("expected one producer bin-0 record")
    source_record = source_records[0]
    source_array = metadata_path.parent / source_record["array"]
    _require_hash(source_array, source_record["sha256"], "producer bin arrays")
    if source_record["fields"] != list(field_ids):
        raise ValueError("producer/report field identity mismatch")
    if source_record["pair_inputs"][PAIR_LABEL] != pair_row:
        raise ValueError("producer/report weighting inputs differ")

    reference_provenance = report["provenance"]["reference"]
    if metadata["identity"]["sources"] != reference_provenance["sources"]:
        raise ValueError("producer/report source-hash inventories differ")
    live_reference = _verify_live_files(
        reference_provenance["sources"],
        relevant_names={
            "weights.py",
            "covariance.py",
            "power_spectrum.py",
            "fisher.py",
        },
    )
    live_resources = _verify_live_files(report["provenance"]["resources"])

    with np.load(arrays_path, allow_pickle=False) as saved:
        arrays = {name: saved[name].copy() for name in saved.files}
    expected_inventory = record["inventory"]
    for name in (
        "k",
        "mu",
        "modes",
        "total",
        "observed_j",
        "selected_pairs",
        "required_pairs",
        "pair_fisher",
        "reference_pair_fisher",
    ):
        if list(arrays[name].shape) != expected_inventory[name]:
            raise ValueError(f"saved {name} shape does not match manifest inventory")
    if len(arrays["k"]) != 5000:
        raise ValueError("W01 requires the saved 5000-node Fourier grid")
    if not np.array_equal(arrays["selected_pairs"], record["task"]["selected_pairs"]):
        raise ValueError("selected-pair indices differ between task and arrays")
    if not np.array_equal(arrays["required_pairs"], record["task"]["required_pairs"]):
        raise ValueError("required-pair indices differ between task and arrays")
    selected_column, required_column = _resolve_auto_columns(
        field_ids, arrays["selected_pairs"], arrays["required_pairs"]
    )

    magnitudes = np.asarray(pair_row["magnitudes"], dtype=np.float64)
    density = np.asarray(pair_row["density"], dtype=np.float64)
    variance = np.asarray(pair_row["variance"], dtype=np.float64)
    captured_weights = np.asarray(pair_row["_w_lya"], dtype=np.float64)
    if not (
        magnitudes.shape
        == density.shape
        == variance.shape
        == captured_weights.shape
        == (107,)
    ):
        raise ValueError("W01 requires four aligned 107-sample weighting arrays")
    if not np.array_equal(magnitudes, np.linspace(16.1, 26.75, 107)):
        raise ValueError("saved magnitude grid is not the literal uniform grid")
    if (magnitudes[0], magnitudes[-1]) != (16.1, 26.75):
        raise ValueError("unexpected saved magnitude support")
    if np.count_nonzero(density < 0) != 9 or np.count_nonzero(density == 0) != 0:
        raise ValueError("saved signed-density inventory differs from W01")

    legacy = LegacyInputs(
        magnitudes=magnitudes,
        density=density,
        variance=variance,
        magnitude_step=float(magnitudes[1] - magnitudes[0]),
        forest_length=float(pair_row["forest_length"]),
        pixel_width=float(pair_row["_pix_kms"]),
        auxiliary_signal=float(pair_row["auxiliary_signal"]),
        auxiliary_p1d=float(pair_row["auxiliary_p1d"]),
    )
    projection = ProjectionInputs(
        k=arrays["k"],
        mu=arrays["mu"],
        modes=arrays["modes"],
        baseline_total=arrays["total"][:, required_column],
        observed_j=arrays["observed_j"][:, selected_column, :],
        saved_reference_fisher=arrays["reference_pair_fisher"][selected_column],
        saved_pair_fisher=arrays["pair_fisher"][selected_column],
        distance_to_velocity=float(pair_row["_distance_to_velocity"]),
        angle_to_distance=float(pair_row["_angle_to_distance"]),
        covariance_redshift=float(pair_row["_z_mean"]),
        resolution_width=float(pair_row["_res_kms"]),
        selected_pairs=arrays["selected_pairs"],
        required_pairs=arrays["required_pairs"],
        field_ids=field_ids,
        selected_column=selected_column,
        required_column=required_column,
    )
    provenance = {
        "input_directory": str(root),
        "manifest": {"path": str(manifest_path), "sha256": sha256_file(manifest_path)},
        "report": {"path": str(report_path), "sha256": record["report_sha256"]},
        "arrays": {"path": str(arrays_path), "sha256": record["sha256"]},
        "reference_manifest": {
            "path": str(reference_manifest_path),
            "sha256": reference_manifest_hash,
        },
        "reference_metadata": {"path": str(metadata_path), "sha256": metadata_hash},
        "producer_arrays": {
            "path": str(source_array),
            "sha256": source_record["sha256"],
        },
        "record_effective_hash": record["effective_hash"],
        "original_recipe_hash": record["original_recipe_hash"],
        "reference_versions": reference_provenance["versions"],
        "verified_live_reference_sources": live_reference,
        "verified_live_resources": live_resources,
    }
    return SourceInputs(
        legacy=legacy,
        projection=projection,
        captured_weights=captured_weights,
        captured_a=float(pair_row["_aliasing_weights"][-1]),
        captured_p_pixel=float(pair_row["_effective_noise_power"][-1]),
        provenance=provenance,
    )


def _floating_failure(error, operation, values=None):
    index = None
    if values is not None:
        bad = np.flatnonzero(~np.isfinite(values))
        if len(bad):
            index = int(bad[0])
    return DiagnosticArithmeticError(
        "arithmetic_failure", operation, str(error), index=index
    )


def _literal_terms(inputs, weights):
    dm = inputs.magnitude_step
    try:
        with np.errstate(all="raise"):
            i1_terms = np.multiply(np.multiply(inputs.density, weights), dm)
            squared = np.square(weights)
            i2_terms = np.multiply(np.multiply(inputs.density, squared), dm)
            i3_terms = np.multiply(
                np.multiply(np.multiply(inputs.density, squared), inputs.variance), dm
            )
            i1_prefix = np.cumsum(i1_terms)
            i2_prefix = np.cumsum(i2_terms)
            i3_prefix = np.cumsum(i3_terms)
    except FloatingPointError as error:
        raise _floating_failure(error, "literal moment products") from error
    return i1_prefix, i2_prefix, i3_prefix


def float_snapshot(inputs, count):
    """Evaluate one absolute update count from the shared literal initialization."""
    if not isinstance(count, int) or count < 0:
        raise ValueError("count must be a nonnegative integer")
    try:
        with np.errstate(all="raise"):
            b = inputs.auxiliary_p1d / inputs.pixel_width
            denominator = b + inputs.variance
            if np.any(denominator == 0):
                index = int(np.flatnonzero(denominator == 0)[0])
                raise DiagnosticArithmeticError(
                    "singular_signed_recurrence",
                    "initialization",
                    "zero denominator",
                    index=index,
                )
            weights = b / denominator
    except FloatingPointError as error:
        raise _floating_failure(error, "initialization") from error

    for iteration in range(1, count + 1):
        i1_prefix, _, _ = _literal_terms(inputs, weights)
        try:
            with np.errstate(all="raise"):
                scale = inputs.forest_length / inputs.pixel_width
                n_eff = i1_prefix * scale
                if np.any(n_eff == 0):
                    index = int(np.flatnonzero(n_eff == 0)[0])
                    raise DiagnosticArithmeticError(
                        "singular_signed_recurrence",
                        f"update {iteration}",
                        "zero cumulative effective density",
                        index=index,
                    )
                noise_ratio = inputs.variance / n_eff
                denominator = inputs.auxiliary_signal + noise_ratio
                if np.any(denominator == 0):
                    index = int(np.flatnonzero(denominator == 0)[0])
                    raise DiagnosticArithmeticError(
                        "singular_signed_recurrence",
                        f"update {iteration}",
                        "zero weight denominator",
                        index=index,
                    )
                weights = inputs.auxiliary_signal / denominator
        except FloatingPointError as error:
            raise _floating_failure(error, f"update {iteration}") from error
        if not np.all(np.isfinite(weights)):
            raise DiagnosticArithmeticError(
                "arithmetic_failure", f"update {iteration}", "nonfinite weights"
            )

    i1_prefix, i2_prefix, i3_prefix = _literal_terms(inputs, weights)
    i1, i2, i3 = (float(array[-1]) for array in (i1_prefix, i2_prefix, i3_prefix))
    try:
        with np.errstate(all="raise"):
            i1_squared = np.square(np.float64(i1))
            denominator = inputs.forest_length * i1_squared
            if denominator == 0:
                raise DiagnosticArithmeticError(
                    "singular_signed_recurrence",
                    "coefficients",
                    "zero I1 normalization",
                )
            coefficient_a = i2 / denominator
            p_pixel = inputs.pixel_width * i3 / denominator
    except FloatingPointError as error:
        raise _floating_failure(error, "coefficients") from error
    scalars = np.array([i1, i2, i3, coefficient_a, p_pixel], dtype=np.float64)
    if not np.all(np.isfinite(scalars)):
        raise DiagnosticArithmeticError(
            "arithmetic_failure", "coefficients", "nonfinite moment or coefficient"
        )
    scale = inputs.forest_length / inputs.pixel_width
    return {
        "count": count,
        "weights": weights,
        "i1_prefix": i1_prefix,
        "i2_prefix": i2_prefix,
        "i3_prefix": i3_prefix,
        "n_eff": i1_prefix * scale,
        "i1": i1,
        "i2": i2,
        "i3": i3,
        "a": float(coefficient_a),
        "p_pixel": float(p_pixel),
    }


def decimal_snapshot(inputs, count, precision):
    """Independent scalar Decimal recurrence from the exact saved float64 values."""
    try:
        with localcontext() as context:
            context.prec = precision
            density = [Decimal.from_float(float(value)) for value in inputs.density]
            variance = [Decimal.from_float(float(value)) for value in inputs.variance]
            dm = Decimal.from_float(inputs.magnitude_step)
            length = Decimal.from_float(inputs.forest_length)
            pixel = Decimal.from_float(inputs.pixel_width)
            signal = Decimal.from_float(inputs.auxiliary_signal)
            alias = Decimal.from_float(inputs.auxiliary_p1d)
            b = alias / pixel
            weights = [b / (b + noise) for noise in variance]
            for _ in range(count):
                cumulative = Decimal(0)
                updated = []
                for rho, weight, noise in zip(density, weights, variance):
                    cumulative += rho * weight * dm
                    n_eff = cumulative * (length / pixel)
                    updated.append(signal / (signal + noise / n_eff))
                weights = updated
            i1 = sum(rho * weight * dm for rho, weight in zip(density, weights))
            i2 = sum(rho * weight**2 * dm for rho, weight in zip(density, weights))
            i3 = sum(
                rho * weight**2 * noise * dm
                for rho, weight, noise in zip(density, weights, variance)
            )
            denominator = length * i1**2
            coefficient_a = i2 / denominator
            p_pixel = pixel * i3 / denominator
    except (DecimalException, ZeroDivisionError) as error:
        kind = (
            "singular_signed_recurrence"
            if isinstance(error, (DivisionByZero, ZeroDivisionError))
            else "arithmetic_failure"
        )
        raise DiagnosticArithmeticError(
            kind, f"Decimal count {count}", str(error)
        ) from error
    return {
        "count": count,
        "precision": precision,
        "weights": weights,
        "i1": i1,
        "i2": i2,
        "i3": i3,
        "a": coefficient_a,
        "p_pixel": p_pixel,
    }


def _relative_change(value, reference):
    if reference == 0:
        return None if value != 0 else 0.0
    return float(value / reference - 1.0)


def _zero_safe_relative_discrepancy(value, reference):
    """Return a relative discrepancy, requiring exact equality at zero."""
    if reference == 0:
        return 0.0 if value == 0 else float("inf")
    return abs(float(value / reference - 1.0))


def _scalar_baseline_check(value, reference):
    discrepancy = _zero_safe_relative_discrepancy(value, reference)
    return {
        "reference_zero": bool(reference == 0),
        "exact_zero_match": bool(reference == 0 and value == 0),
        "relative": discrepancy if np.isfinite(discrepancy) else None,
        "passed": bool(discrepancy <= RELATIVE_TOLERANCE),
    }


def coefficient_comparison(current, reference):
    """Compare A and P_pixel separately, without a combined coefficient norm."""
    a_change = _relative_change(current["a"], reference["a"])
    p_change = _relative_change(current["p_pixel"], reference["p_pixel"])
    trigger = any(
        value is None or abs(value) > COEFFICIENT_TRIGGER
        for value in (a_change, p_change)
    )
    return {
        "a_relative": a_change,
        "p_pixel_relative": p_change,
        "a_absolute": float(current["a"] - reference["a"]),
        "p_pixel_absolute": float(current["p_pixel"] - reference["p_pixel"]),
        "weight_max_absolute": float(
            np.max(np.abs(current["weights"] - reference["weights"]))
        ),
        "triggered": bool(trigger),
    }


def _decimal_confirmation(float_state, decimal_states):
    rows = []
    previous = None
    for state in decimal_states:
        a = float(state["a"])
        p_pixel = float(state["p_pixel"])
        rows.append(
            {
                "precision": state["precision"],
                "a": str(state["a"]),
                "p_pixel": str(state["p_pixel"]),
                "float_a_relative": _relative_change(float_state["a"], a),
                "float_p_pixel_relative": _relative_change(
                    float_state["p_pixel"], p_pixel
                ),
            }
        )
        if previous is not None:
            for name in ("a", "p_pixel"):
                if previous[name] == 0:
                    refinement = "0" if state[name] == 0 else None
                else:
                    refinement = str(state[name] / previous[name] - Decimal(1))
                rows[-1][f"precision_refinement_{name}_relative"] = refinement
        previous = state
    return rows


def select_bounded_attempts(evaluator, confirmer, baseline_validator=None):
    """Evaluate absolute counts with the W01 first-discrepancy stopping rule."""
    states = {}
    confirmations = {}
    attempts = []
    terminal = None
    for count in CHECKPOINTS:
        if terminal is not None:
            suffix = (
                "identification" if terminal == "coefficient_sensitive" else "failure"
            )
            attempts.append({"count": count, "status": f"not_attempted_after_{suffix}"})
            continue
        try:
            state = evaluator(count)
        except DiagnosticArithmeticError as error:
            attempts.append(
                {"count": count, "status": error.kind, "failure": error.record()}
            )
            confirmations[count] = confirmer(count, error)
            terminal = error.kind
            continue
        states[count] = state
        confirmations[count] = confirmer(count, None)
        if count == BASELINE_COUNT and baseline_validator is not None:
            baseline_validator(state)
        attempt = {"count": count, "status": "completed"}
        if count not in (0, BASELINE_COUNT):
            comparison = coefficient_comparison(state, states[BASELINE_COUNT])
            attempt["relative_to_t3"] = comparison
            if comparison["triggered"]:
                terminal = "coefficient_sensitive"
        attempts.append(attempt)
    if terminal is None:
        terminal = "no_sensitivity_detected_within_cap"
    return states, confirmations, attempts, terminal


def _max_elementwise_relative(value, reference):
    value = np.asarray(value)
    reference = np.asarray(reference)
    if value.shape != reference.shape or not np.all(np.isfinite(value)):
        return float("inf")
    zero = reference == 0
    if np.any(value[zero] != 0):
        return float("inf")
    nonzero = ~zero
    if not np.any(nonzero):
        return 0.0
    return float(np.max(np.abs(value[nonzero] / reference[nonzero] - 1.0)))


def baseline_gate(source, state):
    """Check each captured three-update quantity at its own scale."""
    weights_relative = _max_elementwise_relative(
        state["weights"], source.captured_weights
    )
    result = {
        "weights_relative": (
            weights_relative if np.isfinite(weights_relative) else None
        ),
        "weights_passed": bool(weights_relative <= RELATIVE_TOLERANCE),
        "a": _scalar_baseline_check(state["a"], source.captured_a),
        "p_pixel": _scalar_baseline_check(state["p_pixel"], source.captured_p_pixel),
    }
    result["passed"] = bool(
        result["weights_passed"]
        and result["a"]["passed"]
        and result["p_pixel"]["passed"]
    )
    return result


def legacy_p1d(redshift, k_parallel):
    """Literal PD2013 P1D with the lyaforecast stationary low-k floor."""
    k_parallel = np.asarray(k_parallel, dtype=np.float64)
    n_f_z = -2.55 - 0.28 * np.log((1 + redshift) / 4.0)
    floor = 0.009 * np.exp((-0.5 * n_f_z - 1) / -0.1)
    q = np.fmax(k_parallel, floor)
    exponent = 3 + n_f_z - 0.1 * np.log(q / 0.009)
    return (
        np.pi
        * 0.064
        / 0.009
        * (q / 0.009) ** (exponent - 1)
        * ((1 + redshift) / 4.0) ** 3.55
    )


def legacy_response(k_parallel, pixel_width, resolution_width):
    """Literal sinc-times-Gaussian field response, including the exact q=0 limit."""
    q = np.asarray(k_parallel, dtype=np.float64)
    x = 0.5 * q * pixel_width
    pixel = np.ones_like(x)
    nonzero = x != 0
    pixel[nonzero] = np.sin(x[nonzero]) / x[nonzero]
    return pixel * np.exp(-0.5 * q**2 * resolution_width**2)


def _rank_aware_summary(fisher):
    """Apply the established normalized rank rule without finite null error bars."""
    fisher = np.asarray(fisher, dtype=np.float64)
    diagonal = np.diag(fisher)
    if (
        fisher.shape != (2, 2)
        or np.any(diagonal < 0)
        or not np.all(np.isfinite(fisher))
    ):
        raise ValueError("invalid one-spectrum Fisher matrix")
    scale = np.sqrt(diagonal)
    safe = np.where(scale == 0, 1.0, scale)
    normalized = fisher / safe[:, None] / safe[None, :]
    eigenvalues, eigenvectors = np.linalg.eigh(normalized)
    tolerance = 64 * np.finfo(float).eps * 2 * max(1.0, np.max(np.abs(eigenvalues)))
    if eigenvalues[0] < -tolerance:
        raise ValueError("one-spectrum Fisher matrix is indefinite")
    good = eigenvalues > tolerance
    null = eigenvectors[:, ~good]
    available = np.sum(null**2, axis=1) <= 64 * np.finfo(float).eps * 2
    errors = np.full(2, np.nan)
    covariance = np.zeros((2, 2))
    if np.any(good):
        basis = eigenvectors[:, good] / safe[:, None]
        partial = (basis / eigenvalues[good]) @ basis.T
        mask = available[:, None] & available[None, :]
        covariance[mask] = partial[mask]
        errors[available] = np.sqrt(np.diag(covariance)[available])
    rank = int(np.count_nonzero(good))
    if rank == 2:
        direct_errors = np.sqrt(np.diag(np.linalg.inv(fisher)))
        np.testing.assert_allclose(errors, direct_errors, rtol=5e-13, atol=0)
    return {
        "fisher": fisher,
        "covariance": covariance,
        "errors": errors,
        "constrained": available.astype(np.int64),
        "rank": rank,
        "normalized_eigenvalues": eigenvalues,
        "rank_tolerance": float(tolerance),
    }


def scalar_fisher_oracle(observed_j, covariance):
    """Independent scalar accumulation and analytic 2x2 inversion."""
    f00 = f01 = f11 = 0.0
    for derivative, variance in zip(observed_j, covariance):
        inverse = 1.0 / float(variance)
        f00 += float(derivative[0]) ** 2 * inverse
        f01 += float(derivative[0]) * float(derivative[1]) * inverse
        f11 += float(derivative[1]) ** 2 * inverse
    determinant = f00 * f11 - f01 * f01
    if determinant <= 0 or not np.isfinite(determinant):
        return np.array([[f00, f01], [f01, f11]]), None
    inverse = np.array([[f11, -f01], [-f01, f00]]) / determinant
    return np.array([[f00, f01], [f01, f11]]), np.sqrt(np.diag(inverse))


def project_one_spectrum(source, state, baseline_state):
    """Project one coefficient state with saved mean Jacobian and Fourier nodes."""
    projection = source.projection
    q = projection.k * projection.mu / projection.distance_to_velocity
    p1d = legacy_p1d(projection.covariance_redshift, q)
    response = legacy_response(
        q, source.legacy.pixel_width, projection.resolution_width
    )
    conversion = projection.angle_to_distance**2 / projection.distance_to_velocity

    def noise(item):
        return (item["a"] * p1d * response**2 + item["p_pixel"]) * conversion

    total = projection.baseline_total + noise(state) - noise(baseline_state)
    if not np.all(np.isfinite(total)) or np.any(total <= 0):
        raise ValueError(
            "one-spectrum physical total auto-power is not positive finite"
        )
    covariance = 2 * total**2 / projection.modes
    if not np.all(np.isfinite(covariance)) or np.any(covariance <= 0):
        raise ValueError("one-spectrum covariance is not positive finite")
    fisher = np.einsum(
        "ni,nj,n->ij", projection.observed_j, projection.observed_j, 1 / covariance
    )
    summary = _rank_aware_summary(fisher)
    scalar_fisher, scalar_errors = scalar_fisher_oracle(
        projection.observed_j, covariance
    )
    fisher_discrepancy = _max_elementwise_relative(fisher, scalar_fisher)
    if fisher_discrepancy > 5e-12:
        raise ValueError("vectorized/scalar Fisher calculations disagree")
    if summary["rank"] == 2:
        error_discrepancy = _max_elementwise_relative(summary["errors"], scalar_errors)
        if error_discrepancy > 5e-12:
            raise ValueError("matrix/analytic 2x2 inversions disagree")
    else:
        error_discrepancy = None
    summary.update(
        {
            "total": total,
            "covariance_node": covariance,
            "noise": noise(state),
            "p1d": p1d,
            "response": response,
            "scalar_fisher_discrepancy": fisher_discrepancy,
            "analytic_error_discrepancy": error_discrepancy,
        }
    )
    return summary


def _array_digest_update(digest, name, value):
    array = np.ascontiguousarray(np.asarray(value))
    digest.update(name.encode())
    digest.update(array.dtype.str.encode())
    digest.update(json.dumps(array.shape).encode())
    digest.update(array.tobytes())


def source_fingerprint(source):
    """Bind all selected recurrence and one-spectrum inputs to their identities."""
    digest = hashlib.sha256()
    identity = {
        "case": CASE,
        "profile": PROFILE,
        "population": POPULATION,
        "pair": PAIR_LABEL,
        "bin": BIN_INDEX,
        "bounds": BIN_BOUNDS,
        "field_ids": source.projection.field_ids,
        "selected_column": source.projection.selected_column,
        "required_column": source.projection.required_column,
    }
    digest.update(json.dumps(identity, sort_keys=True).encode())
    for name, value in (
        ("magnitude_step", source.legacy.magnitude_step),
        ("forest_length", source.legacy.forest_length),
        ("pixel_width", source.legacy.pixel_width),
        ("auxiliary_signal", source.legacy.auxiliary_signal),
        ("auxiliary_p1d", source.legacy.auxiliary_p1d),
        ("distance_to_velocity", source.projection.distance_to_velocity),
        ("angle_to_distance", source.projection.angle_to_distance),
        ("covariance_redshift", source.projection.covariance_redshift),
        ("resolution_width", source.projection.resolution_width),
        ("captured_a", source.captured_a),
        ("captured_p_pixel", source.captured_p_pixel),
    ):
        _array_digest_update(digest, name, np.float64(value))
    for name, value in (
        ("magnitudes", source.legacy.magnitudes),
        ("density", source.legacy.density),
        ("variance", source.legacy.variance),
        ("captured_weights", source.captured_weights),
        ("k", source.projection.k),
        ("mu", source.projection.mu),
        ("modes", source.projection.modes),
        ("baseline_total", source.projection.baseline_total),
        ("observed_j", source.projection.observed_j),
        ("saved_reference_fisher", source.projection.saved_reference_fisher),
        ("saved_pair_fisher", source.projection.saved_pair_fisher),
        ("selected_pairs", source.projection.selected_pairs),
        ("required_pairs", source.projection.required_pairs),
    ):
        _array_digest_update(digest, name, value)
    return digest.hexdigest()


def _state_signs(state):
    def counts(value):
        value = np.asarray(value)
        return {
            "negative": int(np.count_nonzero(value < 0)),
            "zero": int(np.count_nonzero(value == 0)),
            "positive": int(np.count_nonzero(value > 0)),
        }

    return {
        "weights": counts(state["weights"]),
        "i1_prefix": counts(state["i1_prefix"]),
        "i2_prefix": counts(state["i2_prefix"]),
        "i3_prefix": counts(state["i3_prefix"]),
        "moments": {
            name: int(np.sign(state[name]))
            for name in ("i1", "i2", "i3", "a", "p_pixel")
        },
    }


def _decimal_control(inputs, count, float_state=None, float_error=None):
    """Evaluate both Decimal precisions and classify float/Decimal agreement."""
    decimal_states = []
    failures = []
    for precision in (80, 160):
        try:
            decimal_states.append(decimal_snapshot(inputs, count, precision))
            failures.append(None)
        except DiagnosticArithmeticError as error:
            failures.append(error.record())

    rows = []
    if len(decimal_states) == 2:
        if float_state is None:
            rows = [
                {
                    "precision": state["precision"],
                    "a": str(state["a"]),
                    "p_pixel": str(state["p_pixel"]),
                }
                for state in decimal_states
            ]
        else:
            rows = _decimal_confirmation(float_state, decimal_states)

    discrepancies = (
        [
            row[name]
            for row in rows
            for name in ("float_a_relative", "float_p_pixel_relative")
        ]
        if float_state is not None
        else []
    )
    agreement = bool(
        float_state is not None
        and len(rows) == 2
        and all(
            value is not None and abs(value) <= RELATIVE_TOLERANCE
            for value in discrepancies
        )
    )
    classification = None
    if float_error is not None:
        decimal_kinds = {failure["kind"] for failure in failures if failure is not None}
        classification = (
            "singular_signed_recurrence"
            if decimal_kinds == {"singular_signed_recurrence"}
            and all(failure is not None for failure in failures)
            else "arithmetic_failure"
        )
    elif not agreement:
        classification = "arithmetic_failure"

    return {
        "results": rows,
        "failures": failures,
        "float64_failure": None if float_error is None else float_error.record(),
        "agreement": agreement,
        "classification": classification,
    }


def run_diagnostic(source):
    """Run W01, retaining bounded evidence for every terminal outcome."""
    start = time.monotonic()
    states = {}
    confirmations = {}
    attempts = []
    terminal = None
    terminal_event = None
    baseline = None

    for count in CHECKPOINTS:
        if terminal is not None:
            suffix = (
                "identification" if terminal == "coefficient_sensitive" else "failure"
            )
            attempts.append({"count": count, "status": f"not_attempted_after_{suffix}"})
            continue
        try:
            state = float_snapshot(source.legacy, count)
        except DiagnosticArithmeticError as error:
            confirmation = _decimal_control(source.legacy, count, float_error=error)
            confirmations[count] = confirmation
            terminal = confirmation["classification"]
            terminal_event = {
                "count": count,
                "kind": terminal,
                "reason": "float64_evaluation_failure",
                "float64_failure": error.record(),
                "decimal_failures": confirmation["failures"],
            }
            attempts.append(
                {
                    "count": count,
                    "status": terminal,
                    "failure": terminal_event,
                }
            )
            continue

        states[count] = state
        confirmation = _decimal_control(source.legacy, count, float_state=state)
        confirmations[count] = confirmation
        if confirmation["classification"] is not None:
            terminal = confirmation["classification"]
            terminal_event = {
                "count": count,
                "kind": terminal,
                "reason": "decimal_float64_disagreement",
                "decimal_failures": confirmation["failures"],
            }
            attempts.append(
                {
                    "count": count,
                    "status": terminal,
                    "failure": terminal_event,
                    "float_state_available": True,
                }
            )
            continue

        attempt = {"count": count, "status": "completed"}
        if count == BASELINE_COUNT:
            baseline = baseline_gate(source, state)
            if not baseline["passed"]:
                terminal = "source_discrepancy"
                terminal_event = {
                    "count": count,
                    "kind": terminal,
                    "reason": "captured_three_update_baseline_mismatch",
                    "baseline_gate": baseline,
                }
                attempt = {
                    "count": count,
                    "status": terminal,
                    "failure": terminal_event,
                    "float_state_available": True,
                }
        elif count != 0:
            comparison = coefficient_comparison(state, states[BASELINE_COUNT])
            attempt["relative_to_t3"] = comparison
            if comparison["triggered"]:
                terminal = "coefficient_sensitive"
                terminal_event = {
                    "count": count,
                    "kind": terminal,
                    "reason": "first_coefficient_trigger",
                    "comparison": comparison,
                }
        attempts.append(attempt)

    if terminal is None:
        terminal = "no_sensitivity_detected_within_cap"
        terminal_event = {
            "count": CHECKPOINTS[-1],
            "kind": terminal,
            "reason": "finite_cap_completed_without_trigger",
        }

    projected_counts = []
    projections = {}
    if terminal == "coefficient_sensitive":
        identified = next(
            row["count"]
            for row in attempts
            if row.get("relative_to_t3", {}).get("triggered")
        )
        projected_counts = [BASELINE_COUNT, identified]
    elif terminal == "no_sensitivity_detected_within_cap":
        projected_counts = [BASELINE_COUNT, CHECKPOINTS[-1]]
    for count in projected_counts:
        projections[count] = project_one_spectrum(
            source, states[count], states[BASELINE_COUNT]
        )
    if projected_counts:
        baseline_fisher = projections[BASELINE_COUNT]["fisher"]
        reference_discrepancy = _max_elementwise_relative(
            baseline_fisher, source.projection.saved_reference_fisher
        )
        pair_discrepancy = _max_elementwise_relative(
            baseline_fisher, source.projection.saved_pair_fisher
        )
        if reference_discrepancy > RELATIVE_TOLERANCE:
            raise ValueError("t=3 Fisher does not reproduce reference_pair_fisher")
        if pair_discrepancy > RELATIVE_TOLERANCE:
            raise ValueError("t=3 Fisher does not reproduce pair_fisher")
    else:
        reference_discrepancy = pair_discrepancy = None

    forecast_sensitive = False
    forecast_comparison = None
    if len(projected_counts) == 2:
        reference_count, comparison_count = projected_counts
        first = projections[reference_count]
        second = projections[comparison_count]
        fisher_relative = float(
            np.linalg.norm(second["fisher"] - first["fisher"])
            / np.linalg.norm(first["fisher"])
        )
        error_relative = [
            _relative_change(new, old)
            if np.isfinite(new) and np.isfinite(old)
            else None
            for new, old in zip(second["errors"], first["errors"])
        ]
        forecast_sensitive = fisher_relative > COEFFICIENT_TRIGGER or any(
            value is not None and abs(value) > COEFFICIENT_TRIGGER
            for value in error_relative
        )
        forecast_comparison = {
            "reference_count": reference_count,
            "comparison_count": comparison_count,
            "fisher_relative": fisher_relative,
            "error_relative": error_relative,
        }

    verdicts = []
    if baseline is not None and baseline["passed"]:
        verdicts.append("reproduced_legacy")
    if terminal == "coefficient_sensitive":
        verdicts.append("coefficient_sensitive")
    elif terminal in ("arithmetic_failure", "singular_signed_recurrence"):
        verdicts.append(terminal)
    elif terminal == "source_discrepancy":
        verdicts.append("source_discrepancy")
    else:
        verdicts.append("no_sensitivity_detected_within_cap")
    if forecast_sensitive:
        verdicts.append("forecast_sensitive")

    table = []
    previous = None
    for attempt in attempts:
        count = attempt["count"]
        state = states.get(count)
        projection = projections.get(count)
        row = {"count": count, "status": attempt["status"]}
        if state is not None:
            baseline_state = states.get(BASELINE_COUNT)
            row.update(
                {
                    "a": state["a"],
                    "p_pixel": state["p_pixel"],
                    "a_relative_to_t3": None
                    if baseline_state is None
                    else _relative_change(state["a"], baseline_state["a"]),
                    "p_pixel_relative_to_t3": None
                    if baseline_state is None
                    else _relative_change(state["p_pixel"], baseline_state["p_pixel"]),
                    "a_relative_to_previous": None
                    if previous is None
                    else _relative_change(state["a"], previous["a"]),
                    "p_pixel_relative_to_previous": None
                    if previous is None
                    else _relative_change(state["p_pixel"], previous["p_pixel"]),
                    "signs": _state_signs(state),
                }
            )
            previous = state
        if projection is not None:
            row.update(
                {
                    "errors": [
                        float(value) if np.isfinite(value) else None
                        for value in projection["errors"]
                    ],
                    "rank": projection["rank"],
                }
            )
        else:
            row.update({"errors": [None, None], "rank": None})
        table.append(row)

    return {
        "states": states,
        "confirmations": confirmations,
        "attempts": attempts,
        "terminal": terminal,
        "terminal_event": terminal_event,
        "baseline_gate": baseline,
        "projections": projections,
        "projected_counts": projected_counts,
        "reference_fisher_discrepancy": reference_discrepancy,
        "saved_pair_fisher_discrepancy": pair_discrepancy,
        "forecast_comparison": forecast_comparison,
        "verdicts": verdicts,
        "table": table,
        "seconds": time.monotonic() - start,
    }


def _source_from_payload(summary, arrays):
    scalars = summary["inputs"]["scalars"]
    legacy = LegacyInputs(
        magnitudes=arrays["source_magnitudes"],
        density=arrays["source_density"],
        variance=arrays["source_variance"],
        magnitude_step=scalars["magnitude_step"],
        forest_length=scalars["forest_length"],
        pixel_width=scalars["pixel_width"],
        auxiliary_signal=scalars["auxiliary_signal"],
        auxiliary_p1d=scalars["auxiliary_p1d"],
    )
    identity = summary["identity"]
    projection = ProjectionInputs(
        k=arrays["k"],
        mu=arrays["mu"],
        modes=arrays["modes"],
        baseline_total=arrays["baseline_total"],
        observed_j=arrays["observed_j"],
        saved_reference_fisher=arrays["saved_reference_fisher"],
        saved_pair_fisher=arrays["saved_pair_fisher"],
        distance_to_velocity=scalars["distance_to_velocity"],
        angle_to_distance=scalars["angle_to_distance"],
        covariance_redshift=scalars["covariance_redshift"],
        resolution_width=scalars["resolution_width"],
        selected_pairs=arrays["selected_pairs"],
        required_pairs=arrays["required_pairs"],
        field_ids=tuple(identity["field_ids"]),
        selected_column=identity["selected_column"],
        required_column=identity["required_column"],
    )
    return SourceInputs(
        legacy=legacy,
        projection=projection,
        captured_weights=arrays["captured_weights"],
        captured_a=scalars["captured_a"],
        captured_p_pixel=scalars["captured_p_pixel"],
        provenance=summary["provenance"],
    )


def validate_attempt_inventory(attempts, terminal):
    """Require every prescribed count and an honest first-terminal stopping record."""
    if [row.get("count") for row in attempts] != list(CHECKPOINTS):
        raise ValueError(
            "attempt inventory must contain every prescribed count in order"
        )
    terminal_seen = None
    for row in attempts:
        status = row.get("status")
        if terminal_seen is not None:
            suffix = (
                "identification"
                if terminal_seen == "coefficient_sensitive"
                else "failure"
            )
            if status != f"not_attempted_after_{suffix}":
                raise ValueError("later count was not marked as stopped")
            continue
        if status in (
            "arithmetic_failure",
            "singular_signed_recurrence",
            "source_discrepancy",
        ):
            terminal_seen = status
        elif status == "completed" and row.get("relative_to_t3", {}).get("triggered"):
            terminal_seen = "coefficient_sensitive"
        elif status != "completed":
            raise ValueError("attempt status is inconsistent with the W01 sequence")
    if terminal == "no_sensitivity_detected_within_cap":
        if terminal_seen is not None or attempts[-1]["status"] != "completed":
            raise ValueError("no-sensitivity verdict does not cover the finite cap")
    elif terminal_seen != terminal:
        raise ValueError("terminal outcome is not the first actual terminal event")
    if terminal not in (
        "coefficient_sensitive",
        "no_sensitivity_detected_within_cap",
        "arithmetic_failure",
        "singular_signed_recurrence",
        "source_discrepancy",
    ):
        raise ValueError("unknown terminal outcome")
    return True


def _assert_semantic_equal(actual, expected, path="value"):
    """Compare reconstructed JSON content with zero-safe numerical scaling."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise ValueError(f"{path} mapping inventory mismatch")
        for name in expected:
            _assert_semantic_equal(actual[name], expected[name], f"{path}.{name}")
        return
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(f"{path} sequence inventory mismatch")
        for index, (item, reference) in enumerate(zip(actual, expected)):
            _assert_semantic_equal(item, reference, f"{path}[{index}]")
        return
    if isinstance(expected, bool):
        if not isinstance(actual, bool) or actual != expected:
            raise ValueError(f"{path} differs from reconstructed value")
        return
    if expected is None or isinstance(expected, str):
        if actual != expected:
            raise ValueError(f"{path} differs from reconstructed value")
        return
    if isinstance(expected, int):
        if (
            not isinstance(actual, int)
            or isinstance(actual, bool)
            or actual != expected
        ):
            raise ValueError(f"{path} differs from reconstructed value")
        return
    if isinstance(expected, float):
        if not isinstance(actual, (int, float)):
            raise ValueError(f"{path} is not numerical")
        discrepancy = _zero_safe_relative_discrepancy(actual, expected)
        if discrepancy > RELATIVE_TOLERANCE:
            raise ValueError(f"{path} differs from reconstructed value")
        return
    if actual != expected:
        raise ValueError(f"{path} differs from reconstructed value")


def _assert_array_equal(actual, expected, name):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if actual.shape != expected.shape or actual.dtype != expected.dtype:
        raise ValueError(f"saved {name} shape/dtype mismatch")
    if np.issubdtype(expected.dtype, np.inexact):
        np.testing.assert_allclose(
            actual,
            expected,
            rtol=RELATIVE_TOLERANCE,
            atol=0,
            equal_nan=True,
            err_msg=f"saved {name} differs from reconstruction",
        )
    elif not np.array_equal(actual, expected):
        raise ValueError(f"saved {name} differs from reconstruction")


def validate_payload(summary, arrays, *, authoritative_source=None):
    """Reconstruct every reported W01 result from authenticated source inputs."""
    if summary.get("schema") != EVIDENCE_SCHEMA:
        raise ValueError("unsupported W01 evidence schema")
    identity = summary.get("identity", {})
    expected = {
        "case": CASE,
        "profile": PROFILE,
        "population": POPULATION,
        "pair": PAIR_LABEL,
        "bin": BIN_INDEX,
        "bounds": list(BIN_BOUNDS),
    }
    for name, value in expected.items():
        if identity.get(name) != value:
            raise ValueError(f"incorrect W01 {name} identity")
    source = _source_from_payload(summary, arrays)
    selected, required = _resolve_auto_columns(
        source.projection.field_ids,
        source.projection.selected_pairs,
        source.projection.required_pairs,
    )
    if (selected, required) != (
        source.projection.selected_column,
        source.projection.required_column,
    ):
        raise ValueError("recorded pair columns do not match pair identity")
    fingerprint = source_fingerprint(source)
    if fingerprint != summary.get("source_fingerprint"):
        raise ValueError("selected scientific inputs do not match their fingerprint")
    if authoritative_source is not None:
        authoritative_fingerprint = source_fingerprint(authoritative_source)
        _assert_semantic_equal(
            summary["provenance"],
            authoritative_source.provenance,
            "provenance",
        )
    else:
        authoritative_fingerprint = None
    if (
        authoritative_fingerprint is not None
        and fingerprint != authoritative_fingerprint
    ):
        raise ValueError(
            "selected scientific inputs do not match the authoritative source"
        )
    validate_attempt_inventory(summary["attempts"], summary["terminal"])
    expected_result = run_diagnostic(source)
    expected_summary = _summary(source, expected_result, "validation reconstruction")
    for name in ("identity", "scope", "inputs", "source_fingerprint", "provenance"):
        _assert_semantic_equal(summary[name], expected_summary[name], name)
    for name in (
        "attempts",
        "terminal",
        "terminal_event",
        "verdicts",
        "baseline_gate",
        "projected_counts",
        "reference_fisher_discrepancy",
        "saved_pair_fisher_discrepancy",
        "forecast_comparison",
        "projection_checks",
        "table",
    ):
        _assert_semantic_equal(summary[name], expected_summary[name], name)
    expected_confirmations = {
        str(count): value for count, value in expected_result["confirmations"].items()
    }
    _assert_semantic_equal(
        summary["decimal_confirmations"],
        expected_confirmations,
        "decimal_confirmations",
    )
    expected_arrays = _payload_arrays(source, expected_result)
    if set(arrays) != set(expected_arrays):
        raise ValueError("saved NPZ array inventory mismatch")
    for name, expected_array in expected_arrays.items():
        _assert_array_equal(arrays[name], expected_array, name)
    return True


def _payload_arrays(source, result):
    state_counts = [count for count in CHECKPOINTS if count in result["states"]]
    states = [result["states"][count] for count in state_counts]
    projected = result["projected_counts"]
    magnitude_nodes = len(source.legacy.magnitudes)
    fourier_nodes = len(source.projection.k)

    def state_stack(name):
        if states:
            return np.stack([state[name] for state in states])
        return np.empty((0, magnitude_nodes), dtype=np.float64)

    def projection_stack(name, shape):
        if projected:
            return np.stack([result["projections"][count][name] for count in projected])
        return np.empty((0, *shape), dtype=np.float64)

    return {
        "source_magnitudes": source.legacy.magnitudes,
        "source_density": source.legacy.density,
        "source_variance": source.legacy.variance,
        "captured_weights": source.captured_weights,
        "k": source.projection.k,
        "mu": source.projection.mu,
        "modes": source.projection.modes,
        "baseline_total": source.projection.baseline_total,
        "observed_j": source.projection.observed_j,
        "saved_reference_fisher": source.projection.saved_reference_fisher,
        "saved_pair_fisher": source.projection.saved_pair_fisher,
        "selected_pairs": source.projection.selected_pairs,
        "required_pairs": source.projection.required_pairs,
        "state_counts": np.array(state_counts, dtype=np.int64),
        "weights": state_stack("weights"),
        "i1_prefix": state_stack("i1_prefix"),
        "i2_prefix": state_stack("i2_prefix"),
        "i3_prefix": state_stack("i3_prefix"),
        "n_eff": state_stack("n_eff"),
        "moments": np.asarray(
            [
                [state[name] for name in ("i1", "i2", "i3", "a", "p_pixel")]
                for state in states
            ],
            dtype=np.float64,
        ).reshape(len(states), 5),
        "projected_counts": np.array(projected, dtype=np.int64),
        "projected_fisher": projection_stack("fisher", (2, 2)),
        "projected_errors": projection_stack("errors", (2,)),
        "projected_total": projection_stack("total", (fourier_nodes,)),
        "projected_covariance": projection_stack("covariance_node", (fourier_nodes,)),
        "p1d": (
            result["projections"][projected[0]]["p1d"]
            if projected
            else np.empty(0, dtype=np.float64)
        ),
        "response": (
            result["projections"][projected[0]]["response"]
            if projected
            else np.empty(0, dtype=np.float64)
        ),
    }


def _summary(source, result, command):
    projection = source.projection
    return {
        "schema": EVIDENCE_SCHEMA,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "identity": {
            "case": CASE,
            "profile": PROFILE,
            "population": POPULATION,
            "pair": PAIR_LABEL,
            "bin": BIN_INDEX,
            "bounds": list(BIN_BOUNDS),
            "field_ids": list(projection.field_ids),
            "selected_column": projection.selected_column,
            "required_column": projection.required_column,
        },
        "scope": {
            "magnitude_nodes": len(source.legacy.magnitudes),
            "fourier_nodes": len(projection.k),
            "requested_counts": list(CHECKPOINTS),
            "coefficient_trigger": COEFFICIENT_TRIGGER,
            "relative_tolerance": RELATIVE_TOLERANCE,
            "new_forecast_run": False,
            "populations_evaluated": [POPULATION],
            "bins_evaluated": [BIN_INDEX],
            "spectra_projected": [PAIR_LABEL],
        },
        "inputs": {
            "scalars": {
                "forest_length": source.legacy.forest_length,
                "pixel_width": source.legacy.pixel_width,
                "magnitude_step": source.legacy.magnitude_step,
                "auxiliary_signal": source.legacy.auxiliary_signal,
                "auxiliary_p1d": source.legacy.auxiliary_p1d,
                "distance_to_velocity": projection.distance_to_velocity,
                "angle_to_distance": projection.angle_to_distance,
                "covariance_redshift": projection.covariance_redshift,
                "resolution_width": projection.resolution_width,
                "captured_a": source.captured_a,
                "captured_p_pixel": source.captured_p_pixel,
            },
            "negative_density_nodes": int(np.count_nonzero(source.legacy.density < 0)),
            "zero_density_nodes": int(np.count_nonzero(source.legacy.density == 0)),
        },
        "source_fingerprint": source_fingerprint(source),
        "provenance": source.provenance,
        "attempts": result["attempts"],
        "terminal": result["terminal"],
        "terminal_event": result["terminal_event"],
        "verdicts": result["verdicts"],
        "baseline_gate": result["baseline_gate"],
        "decimal_confirmations": {
            str(count): value for count, value in result["confirmations"].items()
        },
        "projected_counts": result["projected_counts"],
        "reference_fisher_discrepancy": result["reference_fisher_discrepancy"],
        "saved_pair_fisher_discrepancy": result["saved_pair_fisher_discrepancy"],
        "forecast_comparison": result["forecast_comparison"],
        "projection_checks": {
            str(count): {
                "rank": result["projections"][count]["rank"],
                "constrained": result["projections"][count]["constrained"].tolist(),
                "normalized_eigenvalues": result["projections"][count][
                    "normalized_eigenvalues"
                ].tolist(),
                "rank_tolerance": result["projections"][count]["rank_tolerance"],
                "scalar_fisher_discrepancy": result["projections"][count][
                    "scalar_fisher_discrepancy"
                ],
                "analytic_error_discrepancy": result["projections"][count][
                    "analytic_error_discrepancy"
                ],
            }
            for count in result["projected_counts"]
        },
        "table": result["table"],
        "execution": {
            "command": command,
            "seconds": result["seconds"],
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "thread_environment": {
                name: __import__("os").environ.get(name)
                for name in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                )
            },
        },
    }


def write_evidence(output_directory, instructions_path, source, result, command):
    """Create the numerical portion of a unique W01 evidence directory."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=False)
    instructions = Path(instructions_path).resolve()
    if not instructions.is_file():
        raise ValueError("instructions path is not a file")
    arrays = _payload_arrays(source, result)
    arrays_path = output / "arrays.npz"
    np.savez_compressed(arrays_path, **arrays)
    snapshot_path = output / "WEIGHTING_DIAGNOSTIC_STEP.md"
    shutil.copyfile(instructions, snapshot_path)
    script_snapshot = output / "diagnose_legacy_forest_iterations.py"
    shutil.copyfile(Path(__file__).resolve(), script_snapshot)
    summary = _summary(source, result, command)
    summary["artifacts"] = {
        "arrays.npz": sha256_file(arrays_path),
        "WEIGHTING_DIAGNOSTIC_STEP.md": sha256_file(snapshot_path),
        "diagnose_legacy_forest_iterations.py": sha256_file(script_snapshot),
        "script_source_path": str(Path(__file__).resolve()),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    validate_payload(summary, arrays, authoritative_source=source)
    return summary_path


def _write_json(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def finalize_evidence(output_directory, input_directory, handoff_path, checks_path):
    """Bind numerical evidence, checks, script, instructions, and handoff."""
    output = Path(output_directory)
    manifest_path = output / "manifest.json"
    if manifest_path.exists():
        raise FileExistsError(f"evidence is already finalized: {manifest_path}")
    copies = {
        "checks.json": Path(checks_path).resolve(),
        "weighting-diagnostics-w01-r2.md": Path(handoff_path).resolve(),
    }
    for name, source in copies.items():
        if not source.is_file():
            raise ValueError(f"final evidence source is unavailable: {source}")
        destination = output / name
        if destination.exists():
            raise FileExistsError(
                f"final evidence artifact already exists: {destination}"
            )
        shutil.copyfile(source, destination)
    names = (
        "arrays.npz",
        "summary.json",
        "WEIGHTING_DIAGNOSTIC_STEP.md",
        "diagnose_legacy_forest_iterations.py",
        "checks.json",
        "weighting-diagnostics-w01-r2.md",
    )
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "files": {name: sha256_file(output / name) for name in names},
    }
    _write_json(manifest_path, manifest)
    validate_saved_evidence(
        output,
        input_directory,
        instructions_path=output / "WEIGHTING_DIAGNOSTIC_STEP.md",
        handoff_path=handoff_path,
    )
    return manifest_path


def validate_saved_evidence(
    output_directory,
    input_directory,
    *,
    instructions_path=None,
    handoff_path=None,
):
    """Authenticate and reconstruct a finalized W01 evidence directory."""
    output = Path(output_directory)
    manifest = _read_json(output / "manifest.json")
    expected_names = {
        "arrays.npz",
        "summary.json",
        "WEIGHTING_DIAGNOSTIC_STEP.md",
        "diagnose_legacy_forest_iterations.py",
        "checks.json",
        "weighting-diagnostics-w01-r2.md",
    }
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unsupported W01 final manifest schema")
    if set(manifest.get("files", {})) != expected_names:
        raise ValueError("final evidence artifact inventory mismatch")
    for name, digest in manifest["files"].items():
        _require_hash(output / name, digest, f"final evidence {name}")

    checks = _read_json(output / "checks.json")
    if checks.get("schema") != "fishhighz-weighting-diagnostic-checks-v1":
        raise ValueError("unsupported final check inventory schema")
    if not checks.get("checks") or not all(
        row.get("exit_code") == 0 and row.get("command") for row in checks["checks"]
    ):
        raise ValueError("final check inventory contains an unsuccessful check")

    summary = _read_json(output / "summary.json")
    for name in (
        "arrays.npz",
        "WEIGHTING_DIAGNOSTIC_STEP.md",
        "diagnose_legacy_forest_iterations.py",
    ):
        _require_hash(output / name, summary["artifacts"][name], f"saved {name}")
    _require_hash(
        output / "diagnose_legacy_forest_iterations.py",
        sha256_file(__file__),
        "executing diagnostic script",
    )
    if instructions_path is not None:
        _require_hash(
            instructions_path,
            manifest["files"]["WEIGHTING_DIAGNOSTIC_STEP.md"],
            "current instructions",
        )
    if handoff_path is not None:
        _require_hash(
            handoff_path,
            manifest["files"]["weighting-diagnostics-w01-r2.md"],
            "current handoff",
        )
    with np.load(output / "arrays.npz", allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in saved.files}
    authoritative = load_source(input_directory)
    return validate_payload(summary, arrays, authoritative_source=authoritative)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, help="saved profiles-checked directory"
    )
    parser.add_argument(
        "--output", required=True, help="new exclusive W01 evidence directory"
    )
    parser.add_argument(
        "--instructions", required=True, help="approved WEIGHTING_DIAGNOSTIC_STEP.md"
    )
    parser.add_argument(
        "--check-only", action="store_true", help="validate an existing --output"
    )
    parser.add_argument(
        "--finalize", action="store_true", help="bind checks and handoff to --output"
    )
    parser.add_argument("--handoff", help="revision-2 handoff to authenticate")
    parser.add_argument("--checks", help="final check inventory to authenticate")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.check_only and args.finalize:
        raise ValueError("--check-only and --finalize are mutually exclusive")
    if args.check_only:
        validate_saved_evidence(
            args.output,
            args.input,
            instructions_path=args.instructions,
            handoff_path=args.handoff,
        )
        print(f"validated {Path(args.output) / 'summary.json'}")
        return 0
    if args.finalize:
        if args.handoff is None or args.checks is None:
            raise ValueError("--finalize requires --handoff and --checks")
        manifest = finalize_evidence(args.output, args.input, args.handoff, args.checks)
        print(manifest)
        return 0
    source = load_source(args.input)
    result = run_diagnostic(source)
    command = " ".join(
        [sys.executable, str(Path(__file__).resolve()), *(argv or sys.argv[1:])]
    )
    summary = write_evidence(args.output, args.instructions, source, result, command)
    print(summary)
    print(", ".join(result["verdicts"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
