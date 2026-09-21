"""Focused checks for the DESI Run-2 compatibility examples."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from fishhighz.adapters import desi2_compatibility as adapter
from fishhighz.validation.cases import selection
from fishhighz.validation.profile_definitions import forecast_selection

ROOT = Path(__file__).parents[1]
EXAMPLES = ROOT / "examples"


def _load_result_helper():
    spec = importlib.util.spec_from_file_location(
        "desi2_result_helper", EXAMPLES / "_desi2_results.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _captured_bins():
    original = selection(adapter.CASE)
    rows = []
    for index in range(6):
        selected = forecast_selection(adapter.CASE, index)
        required = [tuple(pair) for pair in selected.required_pairs.tolist()]
        pairs = [tuple(pair) for pair in selected.selected_pairs.tolist()]
        total = np.zeros((4, len(required)))
        for column, (i, j) in enumerate(required):
            if i == j:
                total[:, column] = 3.0 + i
        jacobian = np.empty((4, len(pairs), 2))
        for pair_index in range(len(pairs)):
            jacobian[:, pair_index, 0] = [1, 2, 1, 2]
            jacobian[:, pair_index, 1] = [2, 1, 2, 1]
        rows.append(
            {
                "bin": index,
                "mean_z": 2.1 + 0.2 * index,
                "covariance_z": 2.09 + 0.2 * index,
                "k": np.ones(4),
                "mu": np.zeros(4),
                "modes": np.ones(4),
                "total": total,
                "observed_j": jacobian,
                "required_pairs": required,
                "selected_pairs": pairs,
                "pair_inputs": {},
            }
        )
    fields = [field.id for field in original.fields]
    return rows, fields, original, {"versions": {}}, Path("15x2pt.ini")


def test_full_profile_selection_and_result_counts(monkeypatch, tmp_path):
    monkeypatch.setattr(adapter, "_capture", lambda *args: _captured_bins())
    settings, records = adapter.run_desi2_compatibility(
        reference=tmp_path,
        profile="full-compatibility",
        work_directory=tmp_path / "work",
    )
    assert settings["individual_result_count"] == 78
    assert settings["joint_result_count"] == 6
    assert len(settings["bin1_excluded_pairs"]) == 12
    assert sum(row["kind"] == "individual" for row in records[:4]) == 3
    assert all(row["status"] == "available" for row in records)


def test_fixed_profile_changes_only_forest_auto_noise(monkeypatch):
    fields = selection(adapter.CASE).fields
    required = [(0, 0), (0, 1), (1, 1)]
    row = {
        "required_pairs": required,
        "total": np.arange(18.0).reshape(6, 3),
        "k": np.ones(6),
        "mu": np.zeros(6),
        "pair_inputs": {
            "lya(qso)_lya(qso)": {
                "_aliasing_weights": [2.0],
                "_effective_noise_power": [3.0],
            }
        },
    }
    monkeypatch.setattr(
        adapter.WeightInputs,
        "from_pair",
        classmethod(
            lambda cls, context: type("Inputs", (), {"magnitudes": range(107)})()
        ),
    )
    monkeypatch.setattr(
        adapter,
        "solve",
        lambda *args, **kwargs: {
            "status": "converged",
            "updates": 8,
            "candidate": 4,
            "coefficients": np.array([0, 0, 0, 7.0, 11.0]),
        },
    )
    calls = iter((np.full(6, 2.0), np.full(6, 5.0)))
    monkeypatch.setattr(adapter, "forest_noise", lambda *args, **kwargs: next(calls))
    changed, states = adapter._fixed_total(row, required, fields)
    np.testing.assert_array_equal(changed[:, 0], row["total"][:, 0] + 3)
    np.testing.assert_array_equal(changed[:, 1:], row["total"][:, 1:])
    assert set(states) == {"lya(qso)"}


def test_result_output_round_trip_and_new_directory(tmp_path):
    helper = _load_result_helper()
    output = helper.create_output(tmp_path / "result")
    with pytest.raises(FileExistsError):
        helper.create_output(output)
    records = [
        {
            "bin_index": 0,
            "bounds": (2.0, 2.2),
            "redshift": 2.1,
            "kind": "joint",
            "pair": None,
            "pair_indices": None,
            "status": "unavailable_rank_deficient",
            "fisher": np.eye(2),
            "sigma_ap": None,
            "sigma_at": None,
            "correlation": None,
        }
    ]
    helper.write_results(output, settings={"array": np.array([1, 2])}, records=records)
    document = json.loads((output / "settings.json").read_text())
    assert document["array"] == [1, 2]
    assert document["results"][0]["sigma_ap"] is None
    with np.load(output / "results.npz", allow_pickle=False) as arrays:
        np.testing.assert_array_equal(arrays["fisher"], [np.eye(2)])
        assert arrays["available"].tolist() == [0]
        assert np.isnan(arrays["sigma_ap"][0])


@pytest.mark.parametrize(
    "script", ["desi2_full_compatibility.py", "desi2_fixed_compatibility.py"]
)
def test_short_script_help(script):
    result = subprocess.run(
        [sys.executable, str(EXAMPLES / script), "--help"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "--reference" in result.stdout
    assert "--output" in result.stdout
