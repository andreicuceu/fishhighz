"""Actual-control and failed-trial mutations reach schema-3 semantics."""

import copy
import json

import numpy as np
import pytest

from fishhighz.validation.evidence import (
    check,
    execute,
)
from fishhighz.validation.schema import metric_values, request, validate_payload
from fishhighz.validation.synthetic import convergence_payload
from fishhighz.validation.trials import validate


@pytest.fixture(scope="module")
def original():
    """Construct a valid synthetic accuracy payload with producer provenance.

    Returns
    -------
    fixture : tuple
        Task, numerical evidence arrays, and convergence report.
    """
    task = request("lbg_lae_3x2pt", 0, "accuracy")
    arrays, report = convergence_payload(task)
    report["provenance"] = dict(
        fishhighz=dict(origin="fixture", module_hashes={"a": "a" * 64}),
        wheel=dict(origin="fixture", modules={"a": "a" * 64}, sha256="b" * 64),
        reference=dict(
            reference_origin="/fixture",
            sources={"/fixture/a": "c" * 64},
            versions={"fixture": "1"},
        ),
        resources={"fixture": "d" * 64},
    )
    return task, arrays, report


def test_distinct_controls_identical_information(original):
    """Check distinct controls identical information.

    Parameters
    ----------
    original : tuple of dict
        Validated task, numerical evidence arrays, and convergence report
        supplied by the original fixture.
    """
    task, evidence_arrays, evidence_report = original
    assert validate_payload(task, evidence_arrays, evidence_report)
    assert all(m["fisher_relative"] == 0 for m in evidence_report["metrics"])


@pytest.mark.parametrize(
    "mutation",
    [
        "duplicate_operands",
        "wrong_ref",
        "repeated_ref",
        "missing_trial",
        "reordered_trials",
        "controls",
        "primary",
        "pair",
        "volume",
        "lower",
        "missing_contract",
        "missing_outcome",
    ],
)
def test_trial_mutations_writer_and_checker(tmp_path, original, mutation):
    """Check trial mutations writer and checker.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    original : tuple of dict
        Validated task, numerical evidence arrays, and convergence report
        supplied by the original fixture.
    mutation : str
        Modification applied to the otherwise valid fixture, supplied by pytest
        parametrization.
    """
    task, evidence_arrays, evidence_report = copy.deepcopy(original)
    # Begin with a failed metric and retain the contradictory actual trial.
    if mutation == "duplicate_operands":
        index = evidence_report["trial_contract"]["successful_ids"].index(
            evidence_report["trial_contract"]["metric_trials"][0][0]
        )
        evidence_arrays["study_fisher"][index] *= 2
        evidence_arrays["study_pair_fisher"][index] *= 2
    elif mutation == "wrong_ref":
        evidence_report["trial_contract"]["metric_trials"][0][0] = evidence_report[
            "trial_contract"
        ]["metric_trials"][1][0]
    elif mutation == "repeated_ref":
        evidence_report["trial_contract"]["metric_trials"][0][0] = evidence_report[
            "trial_contract"
        ]["metric_trials"][0][1]
    elif mutation == "missing_trial":
        evidence_arrays["study_fisher"] = evidence_arrays["study_fisher"][:-1]
    elif mutation == "reordered_trials":
        evidence_report["study_controls"].reverse()
    elif mutation == "controls":
        evidence_report["study_controls"][0]["iterations"] = 24
    elif mutation == "primary":
        evidence_arrays["study_fisher"][
            evidence_report["trial_contract"]["successful_ids"].index(
                evidence_report["trial_contract"]["final_id"]
            )
        ] *= 2
    elif mutation == "pair":
        evidence_arrays["study_pair_fisher"][
            evidence_report["trial_contract"]["successful_ids"].index(
                evidence_report["trial_contract"]["metric_trials"][0][0]
            )
        ] *= 2
    elif mutation == "volume":
        evidence_arrays["study_volume"][
            evidence_report["trial_contract"]["successful_ids"].index(
                evidence_report["trial_contract"]["metric_trials"][0][0]
            )
        ] *= 2
    elif mutation == "lower":
        evidence_report["combined_lower_controls"]["step"] *= 2
    elif mutation == "missing_contract":
        del evidence_report["trial_contract"]
    elif mutation == "missing_outcome":
        evidence_report["trial_contract"]["outcomes"].pop()
    evidence_report["metrics"] = metric_values(
        evidence_arrays, evidence_report["metric_names"]
    )
    with pytest.raises((ValueError, KeyError)):
        validate_payload(task, evidence_arrays, evidence_report)
    out = tmp_path / "writer"
    manifest = execute(
        out,
        suite="full",
        cases=[task["case"]],
        bin_indices=[0],
        worker=lambda _: (evidence_arrays, evidence_report),
    )
    assert not manifest["complete"] and manifest["records"][0]["status"] == "failed"
    # Consistently rehash and force all manifest flags to reach offline semantics.
    row = manifest["records"][0]
    row.update(status="completed", scientific_passed=True)
    manifest["complete"] = True
    (out / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises((ValueError, KeyError)):
        check(out, verify_sources=False)


def test_failed_attempt_cannot_be_suppressed():
    """Check failed attempt cannot be suppressed."""
    from types import SimpleNamespace

    from fishhighz.validation.study import study

    class FailedStudy:
        selection = SimpleNamespace(
            fields=[SimpleNamespace(kind="forest")], selected_pairs=np.array([[0, 0]])
        )
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

            Raises
            ------
            ValueError
                Deliberately raised to exercise the rejection path in the enclosing
                test.

            Notes
            -----
            Uses small analytic matrices to exercise validation control flow; it does not run a survey forecast.
            """
            if controls["iterations"] == 24:
                raise ValueError("first unrepresentable operation")
            fisher = np.eye(2) * controls["iterations"]
            return dict(fisher=fisher, pair_fisher=fisher[None]), dict(
                settings=dict(grid=dict(volume=1.0))
            )

    evidence_arrays, evidence_report = study(FailedStudy(), {})
    # Replay itself must require the attempted 24-update failure even if both
    # the outcomes and unbound legacy failure list were maliciously cleared.
    assert any(
        o["outcome"] == "failed" for o in evidence_report["trial_contract"]["outcomes"]
    )
    evidence_report["trial_contract"]["outcomes"] = [
        o
        for o in evidence_report["trial_contract"]["outcomes"]
        if o["outcome"] == "success"
    ]
    evidence_report["unresolved_controls"] = []
    evidence_report["settings"]["controls"] = evidence_report["final_controls"]
    evidence_report["settings"]["grid"].update(k_intervals=128, mu_order=32, k_order=4)
    with pytest.raises(ValueError, match="replay"):
        validate(evidence_arrays, evidence_report)


def test_exact_scoped_dispatch(tmp_path):
    """Check exact scoped dispatch.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    from fishhighz.validation.synthetic import evidence_payload

    seen = []

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

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        seen.append(t)
        return evidence_payload(t)

    manifest = execute(
        tmp_path / "b",
        suite="full",
        cases=["lya_qso_lbg_lae_15x2pt"],
        profiles=("compatibility", "accuracy"),
        kind="synthetic_bao",
        worker=worker,
    )
    assert len(seen) == 12 and sum(len(t["selected_pairs"]) for t in seen) == 180
    assert {t["case"] for t in seen} == {"lya_qso_lbg_lae_15x2pt"}
    assert manifest["complete"]
    assert len(manifest["not_run_cases"]) == 6


def test_cumulative_first_cell_fixed_point_and_refinement():
    """Check cumulative first cell fixed point and refinement."""
    from fishhighz.kernels.weights import _iterate

    # At the first cell the map is w -> c*w/(c*w+variance/S).
    # Subdivision reduces c; a positive fixed point disappears below threshold.
    for mass in (0.25, 0.5, 2.0):
        source_weights, _ = _iterate(np.array([mass]), np.ones(1), 1, 1, 1, 1, 100)
        np.testing.assert_allclose(source_weights, [max(0, 1 - 1 / mass)], atol=1e-25)


def test_guard_diagnosis_identifies_mixed_product():
    """Check guard diagnosis identifies mixed product."""
    from fishhighz.validation.weight_diagnosis import first_underflow

    result = first_underflow(np.array([1e-100]), np.array([1e-110]), np.array([1e100]))
    assert result["operation"] == "r*w*w"
    assert result["index"] == 0 and "underflow" in result["error"]
    assert (
        first_underflow(np.array([0.0, 1.0]), np.array([1.0, 1.0]), np.ones(2)) is None
    )


@pytest.mark.parametrize("null", [False, True])
def test_rank_and_roundoff_trial_controls(null):
    """Check rank and roundoff trial controls.

    Parameters
    ----------
    null : bool
        Choice of exact null direction or zero-response case, supplied by pytest
        parametrization.
    """
    task = request("lbg_lae_3x2pt", 0, "accuracy", kind="synthetic_bao")
    evidence_arrays, evidence_report = convergence_payload(task, null=null)
    evidence_arrays["study_fisher"][:, 0, 0] *= 1 + np.finfo(float).eps
    assert validate(evidence_arrays, evidence_report)


def test_scoped_diagnostic_inventory_and_gate(tmp_path):
    """Check scoped diagnostic inventory and gate.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    import runpy
    from pathlib import Path

    from fishhighz.validation.schema import request
    from fishhighz.validation.synthetic import evidence_payload

    script = Path(__file__).resolve().parents[1] / "scripts/check_desi2_15x2pt.py"
    module = runpy.run_path(str(script))
    diagnostics = [
        request(
            "lya_qso_lbg_lae_15x2pt", i, "accuracy", kind="diagnostic", diagnostic_id=p
        )
        for i in range(6)
        for p in module["POLICIES"]
    ]
    manifest = execute(
        tmp_path / "b",
        suite="full",
        cases=["lya_qso_lbg_lae_15x2pt"],
        profiles=("compatibility", "accuracy"),
        kind="synthetic_bao",
        diagnostic_requests=diagnostics,
        worker=evidence_payload,
    )
    assert (
        len(manifest["records"]) == 12
        and len(manifest["diagnostics"]) == 72
        and manifest["complete"]
    )
    # Synthetic coverage is never relabeled real scientific acceptance.
    with pytest.raises(ValueError, match="12 primary"):
        module["scoped_gate"](tmp_path / "b")
    manifest["schema"] = 2
    (tmp_path / "b/manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="schema 3"):
        module["scoped_gate"](tmp_path / "b")
