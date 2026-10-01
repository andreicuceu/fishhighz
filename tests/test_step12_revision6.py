"""Per-spectrum trial binding across sixteen decades of information."""

import copy
import json
import runpy
from pathlib import Path

import numpy as np
import pytest

from fishhighz.validation import schema, trials
from fishhighz.validation.cases import selection
from fishhighz.validation.evidence import check, execute
from fishhighz.validation.numerics import relative
from fishhighz.validation.study import DEFAULT, study

CASE = "lya_qso_lbg_lae_15x2pt"


def provenance():
    """Construct deterministic synthetic source and distribution identities.

    Returns
    -------
    provenance : dict
        Fixed source, wheel, reference, and resource hashes.
    """
    return dict(
        fishhighz=dict(origin="fixture", module_hashes={"a": "a" * 64}),
        wheel=dict(origin="fixture", modules={"a": "a" * 64}, sha256="b" * 64),
        reference=dict(
            reference_origin="/fixture",
            sources={"/fixture/a": "c" * 64},
            versions={"fixture": "1"},
        ),
        resources={"fixture": "d" * 64},
    )


@pytest.fixture(autouse=True)
def tiny_nodes(monkeypatch):
    # Only quadrature construction is replaced, with four deterministic nodes.
    # The trial controller retains every real bounded control and its identity.
    """Replace production quadrature with four synthetic Fourier nodes.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture restoring patched attributes and environment variables after the
        test.
    """
    monkeypatch.setattr(
        schema,
        "grid_nodes",
        lambda _: (
            np.array([0.02, 0.1, 0.2, 0.4]),
            np.array([0.1, 0.3, 0.6, 0.9]),
            np.ones(4),
        ),
    )


def payload(ratio=1e-16, pair=0, *, converged=False, task=None):
    """Construct convergence evidence with one weak selected spectrum.

    Parameters
    ----------
    ratio : float, optional
        Dimensionless relative Fisher strength assigned to one spectrum. Default
        is 1e-16.
    pair : int, optional
        Selected-spectrum index whose information is weakened. Default is 0.
    converged : bool, optional
        Whether the synthetic weighting response is independent of iteration
        count. Default is False.
    task : dict, optional
        Synthetic forecast request including case, redshift-bin index, profile,
        and pair selection. Default is None.

    Returns
    -------
    fixture : tuple
        Task, numerical evidence arrays, and convergence report.
    """
    task = task or schema.request(CASE, 0, "accuracy")
    settings = dict(
        profile=task["profile"],
        parameters=task["parameters"],
        bounds=task["bounds"],
        fields=[f["id"] for f in task["fields"]],
        controls=dict(DEFAULT),
        grid=dict(
            kind="gauss_legendre",
            volume=1000.0,
            h_fid=0.7,
            k_intervals=128,
            k_order=4,
            mu_order=32,
        ),
    )
    k_grid, mu_grid, _ = schema.grid_nodes(settings)
    pairs = np.array(task["required_pairs"])
    total = np.tile(np.eye(5)[pairs[:, 0], pairs[:, 1]], (len(k_grid), 1))
    observed_jacobian = np.ones((len(k_grid), 15, 2))
    observed_jacobian[:, :, 1] = mu_grid[:, None] ** 2
    observed_jacobian[:, pair] *= np.sqrt(ratio)
    base, report = schema.assemble(task, total, observed_jacobian, settings)

    class WeakStudy:
        selection = selection(CASE)
        _prepared = {}

        def evaluate(self, task, controls):
            """Evaluate the synthetic Fisher trial for the supplied numerical controls.

            Parameters
            ----------
            task : dict
                Synthetic forecast request including case, redshift-bin index, profile,
                and pair selection.
            controls : dict
                Numerical quadrature and weighting controls for this synthetic trial.

            Returns
            -------
            payload : tuple
                Synthetic numerical arrays and validation report, including the
                requested test modification.

            Notes
            -----
            Uses small analytic matrices to exercise validation control flow; it does not run a survey forecast.
            """
            evidence_arrays, evidence_report = copy.deepcopy((base, report))
            factor = 1 if converged else (1 + 1 / controls["iterations"]) / (1 + 1 / 24)
            evidence_arrays["pair_fisher"][pair] *= factor
            evidence_arrays["fisher"] += (
                evidence_arrays["pair_fisher"][pair] - base["pair_fisher"][pair]
            )
            return evidence_arrays, evidence_report

    evidence_arrays, evidence_report = study(WeakStudy(), task)
    evidence_report["settings"]["controls"] = evidence_report["final_controls"].copy()
    evidence_report["effective_hash"] = schema.canonical(evidence_report["settings"])
    evidence_arrays["effective_token"] = schema.token(evidence_report["settings"])
    evidence_report["provenance"] = provenance()
    return task, evidence_arrays, evidence_report


def rejected_everywhere(tmp_path, task, arrays, report, reason, record_property):
    """Require matching semantic rejection by validation, writing, and offline checking.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest for generated inputs and results.
    task : dict
        Synthetic forecast request including case, redshift-bin index, profile,
        and pair selection.
    arrays : dict of ndarray
        Named numerical arrays in the synthetic evidence payload.
    report : dict
        Synthetic validation report with settings, provenance, and convergence
        diagnostics.
    reason : str
        Expected diagnostic pattern identifying the rejected evidence.
    record_property : callable
        Pytest callback recording diagnostic evidence in the test report.
    """
    with pytest.raises(ValueError, match=reason) as error:
        schema.validate_payload(task, arrays, report, require_pass=False)
    record_property("validator_rejection", str(error.value))
    out = tmp_path / "writer"
    manifest = execute(
        out,
        suite="full",
        cases=[CASE],
        profiles=("accuracy",),
        bin_indices=[0],
        worker=lambda _: (arrays, report),
    )
    row = manifest["records"][0]
    assert row["status"] == "failed" and reason in row["error"]
    record_property("writer_rejection", row["error"])
    # Writer saved the mutated arrays/report and their actual hashes. Promote
    # only outer flags to reach the offline reader's semantic validation.
    row.update(status="completed", scientific_passed=True)
    manifest["complete"] = True
    (out / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match=reason) as error:
        check(out, verify_sources=False)
    record_property("offline_rejection", str(error.value))


@pytest.mark.parametrize("ratio", [1.0, 1e-8, 1e-16])
@pytest.mark.parametrize("pair", [0, 1, 8, 14])
@pytest.mark.parametrize("family", ["weights", "combined", "all"])
def test_weak_lower_operand(tmp_path, ratio, pair, family, record_property):
    """Check weak lower operand.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    ratio : float
        Dimensionless ratio tested across numerical scales, supplied by pytest
        parametrization.
    pair : int
        One selected tracer pair, supplied by pytest parametrization.
    family : str
        Scientific recipe family, supplied by pytest parametrization.
    record_property : callable
        Pytest callback used to attach validation measurements to the test
        report.
    """
    task, evidence_arrays, evidence_report = payload(ratio, pair)
    assert schema.validate_payload(
        task, evidence_arrays, evidence_report, require_pass=False
    )
    i = evidence_report["metric_names"].index("weights" if family == "all" else family)
    # Analytic inverse of s*F gives sigma/sqrt(s), independent of amplitude.
    expected = 1 - np.sqrt((1 + 1 / 24) / (1 + 1 / 12))
    assert evidence_report["metrics"][i]["pair_error_relative"] == pytest.approx(
        expected, abs=2e-15
    )
    original = evidence_arrays["study_pair_fisher"].copy()
    indices = range(len(evidence_report["metric_names"])) if family == "all" else [i]
    for index in indices:
        evidence_arrays["metric_pair_fisher"][index, 0, pair] = evidence_arrays[
            "metric_pair_fisher"
        ][index, 1, pair]
    evidence_report["metrics"] = schema.metric_values(
        evidence_arrays, evidence_report["metric_names"]
    )
    evidence_report["passed"] = True
    assert np.array_equal(original, evidence_arrays["study_pair_fisher"])
    rejected_everywhere(
        tmp_path,
        task,
        evidence_arrays,
        evidence_report,
        "metric/trial operand",
        record_property,
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "primary",
        "final",
        "replay",
        "metric",
        "verdict",
        "false_pass",
        "pair_fisher",
        "pair_errors",
        "pair_covariance",
        "pair_correlation",
    ],
)
def test_weak_other_bindings(tmp_path, mutation, record_property):
    """Check weak other bindings.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    mutation : str
        Modification applied to the otherwise valid fixture, supplied by pytest
        parametrization.
    record_property : callable
        Pytest callback used to attach validation measurements to the test
        report.
    """
    task, evidence_arrays, evidence_report = payload(converged=True)
    if mutation == "primary":
        evidence_arrays["pair_fisher"][0] *= 1.04
        with pytest.raises(ValueError, match="primary/final pair trial"):
            trials.validate(evidence_arrays, evidence_report)
        reason = "pair_fisher block"
    elif mutation == "final":
        i = evidence_report["trial_contract"]["successful_ids"].index(
            evidence_report["trial_contract"]["final_id"]
        )
        evidence_arrays["study_pair_fisher"][i, 0] *= 1.04
        reason = "primary/final pair trial"
    elif mutation == "replay":
        # Failed combined trial has no successful operand binding. Its saved
        # placeholder must nevertheless agree with the controller replay.
        refs = evidence_report["trial_contract"]["metric_trials"][-1]
        identity = refs[0]
        i = evidence_report["trial_contract"]["successful_ids"].index(identity)
        row = next(
            o
            for o in evidence_report["trial_contract"]["outcomes"]
            if o["id"] == identity
        )
        row.update(outcome="failed", error="unrepresentable arithmetic")
        del row["array_index"]
        evidence_report["trial_contract"]["successful_ids"].pop(i)
        evidence_report["study_controls"].pop(i)
        for key in ("study_fisher", "study_pair_fisher", "study_volume"):
            evidence_arrays[key] = np.delete(evidence_arrays[key], i, axis=0)
        for row in evidence_report["trial_contract"]["outcomes"]:
            if row.get("array_index", -1) > i:
                row["array_index"] -= 1
        evidence_report["unresolved_controls"] = [
            dict(control="combined", error="unrepresentable arithmetic")
        ]
        evidence_report["passed"] = False
        assert schema.validate_payload(
            task, evidence_arrays, evidence_report, require_pass=False
        )
        evidence_arrays["metric_pair_fisher"][-1, 0, 0] *= 1.04
        evidence_report["metrics"] = schema.metric_values(
            evidence_arrays, evidence_report["metric_names"]
        )
        reason = "replayed trial operands metric_pair_fisher"
    elif mutation == "metric":
        evidence_report["metrics"][0]["pair_error_relative"] = 0.0001
        reason = "replayed trial metrics differ"
    elif mutation == "false_pass":
        task, evidence_arrays, evidence_report = payload()
        evidence_report["passed"] = True
        reason = "replayed trial convergence verdict differs"
    elif mutation == "verdict":
        evidence_report["passed"] = False
        reason = "replayed trial convergence verdict differs"
    else:
        evidence_arrays[mutation][0] *= 1.04
        reason = mutation + " block"
    rejected_everywhere(
        tmp_path, task, evidence_arrays, evidence_report, reason, record_property
    )


@pytest.mark.parametrize("ratio", [1.0, 1e-8, 1e-16])
def test_valid_weak_converged_roundoff_and_unconverged(tmp_path, ratio):
    """Check valid weak converged roundoff and unconverged.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    ratio : float
        Dimensionless ratio tested across numerical scales, supplied by pytest
        parametrization.
    """
    task, evidence_arrays, evidence_report = payload(ratio, converged=True)
    evidence_arrays["study_pair_fisher"][:, 0] *= 1 + np.finfo(float).eps
    assert schema.validate_payload(task, evidence_arrays, evidence_report)
    manifest = execute(
        tmp_path / "valid",
        suite="full",
        cases=[CASE],
        profiles=("accuracy",),
        bin_indices=[0],
        worker=lambda _: (evidence_arrays, evidence_report),
    )
    assert (
        manifest["complete"]
        and check(tmp_path / "valid", verify_sources=False)["complete"]
    )
    task, evidence_arrays, evidence_report = payload(ratio)
    assert not evidence_report["passed"] and schema.validate_payload(
        task, evidence_arrays, evidence_report, require_pass=False
    )
    with pytest.raises(ValueError, match="scientific validation failed"):
        schema.validate_payload(task, evidence_arrays, evidence_report)


@pytest.mark.parametrize("scale", [1e-300, 1e-150, 1.0, 1e150, 1e300])
def test_range_safe_block_comparison(scale):
    """Check range safe block comparison.

    Parameters
    ----------
    scale : float
        Multiplicative scale used to test invariance or numerical range,
        supplied by pytest parametrization.
    """
    fisher_matrix = np.diag([scale, scale / 2])
    schema._close_blocks(fisher_matrix * (1 + 1e-13), fisher_matrix, "range")
    with pytest.raises(ValueError, match="range block"):
        schema._close_blocks(fisher_matrix * 1.04, fisher_matrix, "range")
    assert relative(fisher_matrix, np.zeros_like(fisher_matrix)) == float("inf")
    assert relative(np.zeros_like(fisher_matrix), np.zeros_like(fisher_matrix)) == 0


@pytest.mark.parametrize("value", [np.inf, np.nan])
def test_nonfinite_block_rejected(value):
    """Check nonfinite block rejected.

    Parameters
    ----------
    value : float
        Value at the tested validation boundary, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError, match="finite block"):
        schema._close_blocks(np.array([[value]]), np.ones((1, 1)), "finite")


def test_null_and_partially_constrained_trials():
    """Check null and partially constrained trials."""
    for null in (np.zeros((2, 2)), np.diag([1e-250, 0.0])):
        task, evidence_arrays, evidence_report = payload(converged=True)
        evidence_arrays["pair_fisher"][0] = null
        evidence_arrays["study_pair_fisher"][:, 0] = null
        evidence_arrays["metric_pair_fisher"][:, :, 0] = null
        evidence_report["metrics"] = schema.metric_values(
            evidence_arrays, evidence_report["metric_names"]
        )
        assert trials.validate(evidence_arrays, evidence_report)


def test_full_scoped_inventory(tmp_path, monkeypatch):
    """Check full scoped inventory.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    module = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts/check_desi2_15x2pt.py")
    )
    diagnostics = [
        schema.request(CASE, i, "accuracy", kind="diagnostic", diagnostic_id=p)
        for i in range(6)
        for p in module["POLICIES"]
    ]

    def worker(t):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        t : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.
        """
        _, evidence_arrays, evidence_report = payload(converged=True, task=t)
        if t["profile"] == "compatibility":
            from fishhighz.validation.numerics import change

            evidence_arrays["reference_fisher"] = evidence_arrays["fisher"].copy()
            evidence_arrays["reference_pair_fisher"] = evidence_arrays[
                "pair_fisher"
            ].copy()
            evidence_report["comparison"] = change(
                evidence_arrays["fisher"],
                evidence_arrays["fisher"],
                evidence_arrays["pair_fisher"],
                evidence_arrays["pair_fisher"],
                1,
                1,
            )
        return evidence_arrays, evidence_report

    out = tmp_path / "scope"
    manifest = execute(
        out,
        suite="full",
        cases=[CASE],
        profiles=("compatibility", "accuracy"),
        diagnostic_requests=diagnostics,
        worker=worker,
    )
    assert (
        manifest["complete"]
        and len(manifest["records"]) == 12
        and len(manifest["diagnostics"]) == 72
    )
    assert module["scoped_gate"](out)["complete"]
    for kind in ("omission", "extra", "duplicate", "kind", "diagnostic"):
        bad = copy.deepcopy(manifest)
        if kind == "omission":
            bad["records"].pop()
        elif kind == "extra":
            bad["records"].append(bad["records"][0])
        elif kind == "duplicate":
            bad["records"][1] = bad["records"][0]
        elif kind == "kind":
            bad["requested"][0]["kind"] = "synthetic_bao"
        else:
            bad["diagnostics"].pop()
        (out / "manifest.json").write_text(json.dumps(bad))
        with pytest.raises(ValueError, match="12 primary|missing/extra"):
            module["scoped_gate"](out)


@pytest.mark.parametrize("offset", [1e-16, 1e-4])
def test_compatibility_metric_roundoff(offset):
    """Check compatibility metric roundoff.

    Parameters
    ----------
    offset : float
        Perturbation or index offset, supplied by pytest parametrization.
    """
    task = schema.request(CASE, 0, "compatibility")
    _, evidence_arrays, evidence_report = payload(converged=True, task=task)
    evidence_arrays["reference_fisher"] = evidence_arrays["fisher"].copy()
    evidence_arrays["reference_pair_fisher"] = evidence_arrays["pair_fisher"].copy()
    evidence_report["comparison"] = dict(
        fisher_relative=offset,
        error_relative=0.0,
        pair_error_relative=0.0,
        volume_relative=0.0,
    )
    if offset < 1e-15:
        assert schema.validate_payload(task, evidence_arrays, evidence_report)
    else:
        with pytest.raises(ValueError, match="incorrect reference comparison metrics"):
            schema.validate_payload(
                task, evidence_arrays, evidence_report, require_pass=False
            )
