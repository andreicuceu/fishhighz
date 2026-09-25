"""Synthetic fixed-growth BAO parameter and result checks."""

import json
from pathlib import Path

import numpy as np
import pytest
from test_public_forecast import Background, inputs

from fishhighz import Forecast, parse_survey_ini
from fishhighz.bao_marginalized import active_names, make_registry, reported_constraint
from fishhighz.kernels.kaiser import _scales
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult


def test_dilation_basis_and_covariance_jacobian():
    alpha, phi = 1.07, 0.92
    ap, at, q = _scales((alpha, phi), 3)
    assert ap == pytest.approx(alpha * phi ** (-2 / 3))
    assert at == pytest.approx(alpha * phi ** (1 / 3))
    assert q == pytest.approx(1 / alpha**3)
    jacobian = np.array(
        [
            [phi ** (-2 / 3), -2 * alpha * phi ** (-5 / 3) / 3],
            [phi ** (1 / 3), alpha * phi ** (-2 / 3) / 3],
        ]
    )
    covariance = np.array([[0.02, -0.003], [-0.003, 0.01]])
    np.testing.assert_allclose(
        jacobian @ covariance @ jacobian.T,
        np.array(
            [
                [
                    np.sum(jacobian[0, :, None] * covariance * jacobian[0, None, :]),
                    jacobian[0] @ covariance @ jacobian[1],
                ],
                [
                    jacobian[1] @ covariance @ jacobian[0],
                    jacobian[1] @ covariance @ jacobian[1],
                ],
            ]
        ),
    )


def test_standard_registry_has_five_nuisances_and_fixed_growth(tmp_path):
    source = tmp_path / "bao.ini"
    source.write_text(
        Path("fishhighz/data/desi2_accuracy.ini").read_text()
        + "\n[model]\nmode = bao_marginalized\n"
    )
    config = parse_survey_ini(source)
    registry = make_registry(config, Background())
    assert len(active_names(config, 0)) == 5
    assert len(active_names(config, 1)) == 7
    assert active_names(config, 1) == (
        "alpha_iso",
        "phi",
        "b_lya",
        "b_qso",
        "b_lbg",
        "b_lae",
        "beta_lya",
    )
    assert not any(parameter.id.startswith("f_") for parameter in registry.parameters)
    assert len(registry.ids) == 40


def test_empty_bin_and_subgroup_nuisances(tmp_path):
    source = tmp_path / "subgroups.ini"
    content = Path("fishhighz/data/desi2_accuracy.ini").read_text()
    content = (
        content.replace(
            "selected = lya(qso)xlya(qso), lya(qso)xqso, qsoxqso", "selected = "
        )
        .replace(
            "[pairs bin 2]\nselected = all",
            "[pairs bin 2]\nselected = lbgxlbg, lbgxlae, laexlae",
        )
        .replace(
            "[pairs bin 3]\nselected = all",
            "[pairs bin 3]\nselected = lya(qso)xlya(qso), lya(qso)xlya(lbg), lya(lbg)xlya(lbg)",
        )
    )
    source.write_text(content + "\n[model]\nmode = bao_marginalized\n")
    config = parse_survey_ini(source)
    assert active_names(config, 0) == ("alpha_iso", "phi")
    assert active_names(config, 1) == ("alpha_iso", "phi", "b_lbg", "b_lae")
    assert active_names(config, 2) == ("alpha_iso", "phi", "b_lya", "beta_lya")
    registry = make_registry(config, Background())
    assert not any(name.endswith("_0") for name in registry.ids)


def test_dense_marginalization_and_singular_status():
    registry = ParameterRegistry(
        [
            Parameter(name, 1, "target" if index < 2 else "nuisance")
            for index, name in enumerate(
                ("alpha_iso", "phi", "b_qso", "b_lya", "beta_lya")
            )
        ]
    )
    matrix = np.random.default_rng(53).normal(size=(12, 5))
    fisher = matrix.T @ matrix
    status, covariance, errors, correlation = reported_constraint(
        FisherResult(registry, fisher), ("alpha_iso", "phi")
    )
    assert status == "available"
    np.testing.assert_allclose(covariance, np.linalg.inv(fisher)[:2, :2], rtol=1e-13)
    np.testing.assert_allclose(errors**2, np.diag(covariance))
    np.testing.assert_allclose(np.diag(correlation), 1)
    fisher[-1, :] = 0
    fisher[:, -1] = 0
    assert reported_constraint(
        FisherResult(registry, fisher), ("alpha_iso", "phi")
    ) == ("unavailable", None, None, None)


def test_joint_only_save_and_fixed_kaiser_growth(tmp_path):
    source, template, readers = inputs(tmp_path)
    source.write_text(
        source.read_text()
        .replace("parameterization = ap_at", "parameterization = alpha_iso_phi")
        .replace("parameter_names = ap, at", "parameter_names = alpha_iso, phi")
        .replace("wiggle_scaling = ap_at", "wiggle_scaling = alpha_iso_phi")
        .replace("[model]", "[model]\nmode = bao_marginalized")
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    spec = forecast.prepare().survey.bins[0]
    model = spec.p3d.routes[0].provider.model
    assert "f" not in model.local_names
    assert model.bases == ("ap_at", "alpha_iso_phi")
    assert model.fixed[2 * len(model.fields)] == pytest.approx(0.8)
    result = forecast.run(individuals=False)
    assert len(result.joint) == 1 and not result.individual
    assert result.joint[0].target_ids == ("alpha_iso_0", "phi_0")
    assert result.joint[0].parameter_ids == (
        "alpha_iso_0",
        "phi_0",
        "b_qso_0",
        "b_lbg_0",
    )
    output = result.save(tmp_path / "result")
    settings = json.loads((output / "settings.json").read_text())
    assert settings["schema"] == {
        "name": "fishhighz-bao-marginalized-result",
        "version": 1,
    }
    assert settings["resolved_settings"]["fixed_growth_by_bin"]["0"] == 0.8
    with np.load(output / "results.npz", allow_pickle=False) as arrays:
        assert arrays["record_0000_fisher"].shape == (4, 4)
