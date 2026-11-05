"""Rank-aware report errors and working-directory-independent example loading."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose


@pytest.fixture(scope="module")
def report_module():
    """Import the optional report without reading or generating study outputs."""
    pytest.importorskip("eftkernels")
    pytest.importorskip("matplotlib")
    path = Path(__file__).resolve().parents[1] / "examples/eft_3x2_report.py"
    spec = importlib.util.spec_from_file_location("eft_report", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fixed_higher_singular_errors(report_module):
    """Retain identifiable errors and expose an exact f/eta-bias degeneracy."""
    design = np.eye(7)
    design[:, 6] = design[:, 4]
    row = dict(
        data="FF",
        kmax=0.08,
        names=["alpha_w", "phi_w", "alpha_s", "phi_s", "f", "bias_F", "bias_eta_F"],
    )
    record = report_module.fixed_higher_forecast(
        row, design[:, None, :], np.full(7, 0.05), 0.3
    )
    assert record["rank"] == 6
    assert record["errors"][:4] == [1.0] * 4
    assert record["errors"][4] is None
    assert record["errors"][5] == 1.0
    assert record["errors"][6] is None
    assert record["fsigma8_error"] is None
    assert len(record["null_directions"]) == 1
    json.dumps(record, allow_nan=False)


def test_fixed_higher_full_rank_and_cut(report_module):
    """Agree with an ordinary inverse when the selected design is well posed."""
    design = np.arange(1, 36, dtype=float).reshape(7, 5) + np.eye(7, 5) * 10
    row = dict(data="GG", kmax=0.08, names=list(report_module.study.TARGETS))
    k_grid = np.array([0.05] * 6 + [0.2])
    record = report_module.fixed_higher_forecast(row, design[:, None, :], k_grid, 0.3)
    expected = np.sqrt(np.diag(np.linalg.inv(design[:6].T @ design[:6])))
    assert record["rank"] == 5
    assert_allclose(record["errors"], expected, rtol=1e-13)
    assert_allclose(record["fsigma8_error"], expected[4] * 0.3, rtol=1e-13)


def test_identifiable_errors_underdetermined(report_module):
    """Complete the null basis when the design has fewer rows than parameters."""
    errors, rank, null, _ = report_module.study.identifiable_errors(np.eye(2, 5))
    assert rank == 2
    assert null.shape == (3, 5)
    assert_allclose(errors[:2], 1.0)
    assert np.isinf(errors[2:]).all()


def test_report_missing_inputs(report_module, tmp_path):
    """List all absent prerequisites before creating a partial report."""
    with pytest.raises(FileNotFoundError, match="reproduction commands") as error:
        report_module.build_report(tmp_path / "results", tmp_path / "report")
    assert "halo_loop_sums.json" in str(error.value)
    assert "step_half" in str(error.value)
    assert not (tmp_path / "report").exists()


def test_report_cli_outside_checkout(report_module, tmp_path):
    """Load the sibling study module even from an unrelated directory."""
    result = subprocess.run(
        [sys.executable, report_module.__file__, "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "--results" in result.stdout
