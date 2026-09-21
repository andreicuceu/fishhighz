"""Scientific contract for the standalone research-baseline API example."""

import runpy
from pathlib import Path

import numpy as np


def test_research_bao_example_is_joint_adaptive_and_synthetic():
    namespace = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "examples/research_bao_forecast.py")
    )
    report = namespace["run"]()

    assert report["interpretation"] == "synthetic API example; not a DESI-2 forecast"
    assert report["selected_spectra"] == [
        "forest x forest",
        "forest x galaxy",
        "galaxy x galaxy",
    ]
    assert report["weighting"]["method"] == "early_lyaforecast"
    assert report["weighting"]["status"] == "converged"
    assert report["weighting"]["updates"] >= 6
    assert report["weighting"]["stopping_controls"] == {
        "rtol": 1e-4,
        "min_updates": 3,
        "stable_steps": 3,
        "max_updates": 96,
    }
    assert report["fixed_preparation"]["selected_pair_count"] == 3
    assert report["fixed_preparation"]["required_pair_count"] == 3
    assert report["fixed_preparation"]["noise_convention"] == "independent sampling"
    assert report["joint_rank"] == 2
    assert np.all(np.isfinite(report["joint_ap_at_errors"]))
    assert all(
        np.all(np.isfinite(errors))
        for errors in report["individual_ap_at_errors"].values()
    )
    assert not any(
        np.allclose(report["joint_ap_at_errors"], errors)
        for errors in report["individual_ap_at_errors"].values()
    )
    assert report["joint_vs_sum_max_relative_difference"] > 1e-3
    assert not np.allclose(
        report["joint_data_fisher"], report["summed_individual_data_fisher"]
    )
