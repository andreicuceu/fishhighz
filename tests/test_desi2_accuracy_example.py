"""Focused contract tests for the documented DESI-2 accuracy example."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.forecast import prepare_bin
from fishhighz.geometry import prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import InstrumentResponse
from fishhighz.results import FisherResult
from fishhighz.survey import BinSpec


def _example():
    examples = Path(__file__).resolve().parents[1] / "examples"
    sys.path.insert(0, str(examples))
    try:
        spec = importlib.util.spec_from_file_location(
            "desi2_accuracy_example", examples / "desi2_accuracy.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def test_accuracy_example_controls_and_selection():
    example = _example()
    assert example.CONTROLS == {
        "k_intervals": 128,
        "k_order": 4,
        "k_min": 0.01,
        "k_max": 0.5,
        "mu_order": 32,
        "z_order": 32,
        "magnitude_order": 16,
        "ap_step": 2.5e-4,
        "weight_rtol": 1e-5,
    }

    inputs = object.__new__(example.SurveyInputs)
    inputs.fields = (
        ObservedField("lya(qso)", "forest", "lya", background="qso"),
        ObservedField("qso", "galaxy", "qso"),
        ObservedField("lbg", "galaxy", "lbg"),
        ObservedField("lae", "galaxy", "lae"),
        ObservedField("lya(lbg)", "forest", "lya", background="lbg"),
    )
    assert len(inputs.selection(0).selected_pairs) == 3
    assert len(inputs.selection(0).required_pairs) == 3
    assert len(inputs.selection(1).selected_pairs) == 15
    assert len(inputs.selection(1).required_pairs) == 15


def test_accuracy_example_constraints_report_rank():
    example = _example()
    registry = ParameterRegistry(
        [Parameter("ap", 1, "target"), Parameter("at", 1, "target")]
    )
    available = FisherResult(registry, [[4.0, 1.0], [1.0, 9.0]])
    status, sigma_ap, sigma_at, correlation = example._constraints(available)
    assert status == "available"
    assert np.all(np.isfinite([sigma_ap, sigma_at, correlation]))

    unavailable = FisherResult(registry, [[1.0, 1.0], [1.0, 1.0]])
    assert example._constraints(unavailable) == (
        "unavailable",
        None,
        None,
        None,
    )


def test_accuracy_example_refuses_existing_output_before_forecast(
    tmp_path, monkeypatch
):
    example = _example()
    reference = tmp_path / "reference"
    reference.mkdir()
    template = tmp_path / "template.fits"
    template.touch()
    output = tmp_path / "existing"
    output.mkdir()

    def unexpected(*args, **kwargs):
        raise AssertionError("forecast preparation must not start")

    monkeypatch.setattr(example, "_prepare_bins", unexpected)
    with pytest.raises(FileExistsError):
        example.run(reference=reference, template=template, output=output)


def test_accuracy_example_forecast_records_two_bins_end_to_end():
    example = _example()
    registry = ParameterRegistry(
        [
            Parameter(f"{name}_{index}", 1.0, "target", step=0.001)
            for index in range(2)
            for name in ("ap", "at")
        ]
    )
    field = ObservedField("galaxy", "galaxy", "galaxy")
    selection = PairSelection([field])
    grid = gauss_legendre_grid([0.03, 0.1, 0.2], k_order=2, mu_order=3, h_fid=0.7)
    prepared = []
    for index, z_min in enumerate((2.0, 3.0)):
        geometry = prepare_geometry(
            z_min,
            z_min + 0.5,
            z_eval=z_min + 0.25,
            area_deg2=100,
            h_fid=0.7,
            z_order=3,
            hubble=lambda z: np.full_like(z, 200.0),
            transverse_distance=lambda z: 1000.0 * (1 + z),
        )

        def model(theta, redshift, k, mu, pairs):
            value = 5.0 + theta[0] * (1 + k) + theta[1] * (1 + mu)
            return np.repeat(value[:, None], len(pairs), axis=1)

        binding = BoundParameters(
            registry,
            ("ap", "at"),
            {"ap": f"ap_{index}", "at": f"at_{index}"},
        )
        p3d = PreparedP3D(
            registry,
            selection,
            [P3DProvider(f"bin-{index}", model, binding, selection.required_pairs)],
        )
        prepared.append(
            prepare_bin(
                BinSpec(
                    f"bin-{index}",
                    geometry,
                    grid,
                    p3d,
                    {"galaxy": InstrumentResponse(0, 0)},
                    full_noise=np.ones((len(grid.k_flat), 1)),
                )
            )
        )

    records, combined = example._forecast_records(
        SimpleNamespace(registry=registry), prepared
    )
    assert combined.shape == (4, 4)
    assert len(records) == 4
    assert [(row["bin_index"], row["kind"]) for row in records] == [
        (0, "joint"),
        (0, "individual"),
        (1, "joint"),
        (1, "individual"),
    ]
    assert all(row["status"] == "available" for row in records)
    assert all(np.asarray(row["fisher"]).shape == (2, 2) for row in records)
    assert records[0]["pair_indices"] is None
    assert records[1]["pair_indices"] == [0, 0]
