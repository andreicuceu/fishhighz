"""Synthetic full-shape inference and Fourier-cut controls."""

import json
from pathlib import Path

import numpy as np
import pytest
from test_public_forecast import Background, inputs

from fishhighz import Forecast, parse_survey_ini
from fishhighz.full_shape import (
    ISO_TARGETS,
    TARGETS,
    active_names,
    make_registry,
    reported_constraint,
)
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult


def full_ini(source):
    text = source.read_text().replace("[model]", "[model]\nmode = full_shape")
    for before, after in [
        ("parameterization = ap_at", "parameterization = alpha_phi"),
        ("parameter_names = ap, at", "parameter_names = alpha_w,phi_w,alpha_s,phi_s,f"),
        ("smooth_scaling = identity", "smooth_scaling = alpha_phi"),
        ("wiggle_scaling = ap_at", "wiggle_scaling = alpha_phi"),
    ]:
        text = text.replace(before, after)
    source.write_text(text)


def test_isotropic_basis_native_fisher_jacobian(tmp_path):
    old_source, template, readers = inputs(tmp_path)
    full_ini(old_source)
    old = Forecast(
        old_source, background=Background(), template=template, readers=readers
    )
    old_result = old.run(individuals=False)
    iso_source = tmp_path / "iso.ini"
    iso_source.write_text(
        old_source.read_text()
        .replace("alpha_phi", "alpha_iso_phi")
        .replace("alpha_w", "alpha_iso_w")
        .replace("alpha_s", "alpha_iso_s")
    )
    iso = Forecast(
        iso_source, background=Background(), template=template, readers=readers
    )
    new_result = iso.run(individuals=False)
    old_record, new_record = old_result.joint[0], new_result.joint[0]
    assert new_record.target_ids == tuple(f"{name}_0" for name in ISO_TARGETS[:4]) + (
        "fsigma8_0",
    )
    assert new_record.status == old_record.status
    old_ids = old_record.fisher.registry.ids
    new_ids = new_record.fisher.registry.ids
    assert len(old_ids) == len(new_ids)
    jacobian = np.eye(len(old_ids))
    for component in ("w", "s"):
        old_alpha = old_ids.index(f"alpha_{component}_0")
        new_phi = new_ids.index(f"phi_{component}_0")
        jacobian[old_alpha, new_phi] = -1 / 6
    expected = jacobian.T @ old_record.fisher.data_fisher @ jacobian
    np.testing.assert_allclose(
        new_record.fisher.data_fisher, expected, rtol=3e-5, atol=1e-7
    )
    if new_record.status == "available":
        inverse = np.linalg.inv(jacobian)
        old_cov = np.linalg.inv(old_record.fisher.data_fisher)
        target_indices = [new_ids.index(f"{name}_0") for name in ISO_TARGETS]
        expected_cov = inverse @ old_cov @ inverse.T
        sigma8 = new_record.sigma8_fid
        scaling = np.diag([1, 1, 1, 1, sigma8])
        np.testing.assert_allclose(
            new_record.target_covariance,
            scaling @ expected_cov[np.ix_(target_indices, target_indices)] @ scaling,
            rtol=3e-5,
        )


def test_forest_only_isotropic_fisher_jacobian_with_shared_nuisances(tmp_path):
    """Four dilation targets retain the complete forest nuisance covariance."""
    from fishhighz.full_shape import make_model
    from fishhighz.models.templates import prepare_template

    text = Path("fishhighz/data/desi2_accuracy.ini").read_text()
    forest_pairs = "lya(qso)xlya(qso), lya(qso)xlya(lbg), lya(lbg)xlya(lbg)"
    text = text.replace("selected = all", f"selected = {forest_pairs}")
    text = text.replace(
        "selected = lya(qso)xlya(qso), lya(qso)xqso, qsoxqso",
        "selected = lya(qso)xlya(qso)",
    )
    path = tmp_path / "forest_old.ini"
    path.write_text(text + "\n[model]\nmode = full_shape\ntarget_set = dilation_only\n")
    old = parse_survey_ini(path)
    path.write_text(path.read_text() + "parameterization = alpha_iso_phi\n")
    new = parse_survey_ini(path)
    k = np.linspace(0.001, 1.0, 600)
    template = prepare_template(
        k,
        100 + np.sin(100 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    pairs = [(0, 0), (0, 4), (4, 4)]
    kval = np.array([0.04, 0.08, 0.12])
    mu = np.array([0.25, 0.55, 0.85])
    covariance = np.eye(3) * 0.7 + np.ones((3, 3)) * 0.1
    fishers = []
    for config in (old, new):
        registry = make_registry(config, Background())
        fields = config.observed_fields
        model, binding = make_model(
            config,
            config.bins[1],
            registry,
            template,
            fields,
            {field.id: 1.0 for field in fields},
            {field.id: 1.45 for field in fields if field.kind == "forest"},
            {field.id: (0.0, 0.0) for field in fields},
            0.8,
            1.0,
        )
        indices = binding.binding.local_to_global
        active = tuple(f"{name}_1" for name in active_names(config, 1))
        assert len(active) == 6
        gradients = []
        for name in active:
            slot = registry.ids.index(name)
            step = registry.parameters[slot].step
            plus = registry.fiducials.copy()
            minus = plus.copy()
            plus[slot] += step
            minus[slot] -= step
            gradients.append(
                (
                    model(plus[indices], config.bins[1].z_eval, kval, mu, pairs)
                    - model(minus[indices], config.bins[1].z_eval, kval, mu, pairs)
                )
                / (2 * step)
            )
        gradients = np.asarray(gradients).transpose(1, 2, 0)
        fisher = sum(g.T @ np.linalg.solve(covariance, g) for g in gradients)
        fishers.append(fisher)
    jacobian = np.eye(6)
    jacobian[0, 1] = jacobian[2, 3] = -1 / 6
    np.testing.assert_allclose(
        fishers[1], jacobian.T @ fishers[0] @ jacobian, rtol=3e-5, atol=1e-7
    )


def test_registry_selection(tmp_path):
    source = tmp_path / "full.ini"
    source.write_text(
        Path("fishhighz/data/desi2_accuracy.ini").read_text()
        + "\n[model]\nmode = full_shape\n"
    )
    config = parse_survey_ini(source)
    registry = make_registry(config, Background())
    assert [len(active_names(config, i)) for i in range(6)] == [8, 10, 10, 10, 10, 10]
    assert len(registry.ids) == 58
    assert active_names(config, 1, [("lya(qso)", "lya(lbg)")]) == TARGETS + (
        "b_lya",
        "beta_lya",
    )


def test_dense_marginalization_and_growth_cross_covariance():
    registry = ParameterRegistry(
        [
            Parameter(name, 1.0, "target" if i < 5 else "nuisance")
            for i, name in enumerate((*TARGETS, "b"))
        ]
    )
    a = np.random.default_rng(42).normal(size=(16, 6))
    fisher = a.T @ a
    status, covariance, errors, correlations = reported_constraint(
        FisherResult(registry, fisher), TARGETS, 0.3
    )
    j = np.diag([1.0, 1.0, 1.0, 1.0, 0.3])
    np.testing.assert_allclose(
        covariance, j @ np.linalg.inv(fisher)[:5, :5] @ j, rtol=1e-13, atol=1e-15
    )
    assert status == "available"
    np.testing.assert_allclose(errors**2, np.diag(covariance))
    np.testing.assert_allclose(np.diag(correlations), 1.0)
    fisher[4] = 0
    fisher[:, 4] = 0
    assert reported_constraint(FisherResult(registry, fisher), TARGETS, 0.3) == (
        "unavailable",
        None,
        None,
        None,
    )


@pytest.mark.parametrize("lo,hi", [(0.02, 0.1), (0.01, 0.15), (0.025, 0.2)])
def test_public_cuts_and_numeric_serialization(tmp_path, lo, hi):
    source, template, readers = inputs(tmp_path)
    full_ini(source)
    source.write_text(
        source.read_text()
        .replace("k_min = 0.01", f"k_min = {lo}")
        .replace("k_max = 0.1", f"k_max = {hi}")
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    grid = forecast.prepare().survey.bins[0].grid
    assert grid.k.min() > lo and grid.k.max() < hi
    assert grid.q_mode.sum() == pytest.approx(
        (hi**3 - lo**3) / (6 * np.pi**2), rel=1e-14
    )
    result = forecast.run()
    grid_metadata = result.resolved_settings["grid_by_bin"]["desi2_accuracy-1"]
    assert len(grid_metadata["k_nodes"]) == grid.n_k
    assert len(grid_metadata["mu_nodes"]) == grid.n_mu
    assert (
        sum(item["radial_nodes"] for item in grid_metadata["category_intervals"])
        == grid.n_k
    )
    assert len(result.joint[0].parameter_ids) == 7
    assert len(result.individual[0].parameter_ids) == 6
    assert result.joint[0].target_ids[-1] == "fsigma8_0"
    output = result.save(tmp_path / "out")
    assert (
        json.loads((output / "settings.json").read_text())["schema"]["name"]
        == "fishhighz-full-shape-result"
    )
    with np.load(output / "results.npz", allow_pickle=False) as arrays:
        assert arrays["record_0000_fisher"].shape == (6, 6)
        assert arrays["combined_fisher"].shape == (7, 7)
        assert all(arrays[key].dtype != object for key in arrays)


@pytest.mark.parametrize(
    "lo,hi", [(0, 0.1), (0.1, 0.1), (0.2, 0.1), ("nan", 0.1), (0.01, "inf")]
)
def test_invalid_limits(tmp_path, lo, hi):
    source, _, _ = inputs(tmp_path)
    full_ini(source)
    source.write_text(
        source.read_text()
        .replace("k_min = 0.01", f"k_min = {lo}")
        .replace("k_max = 0.1", f"k_max = {hi}")
    )
    with pytest.raises(ValueError, match="invalid domains"):
        Forecast(source)


def test_bao_rejects_full_shape_controls(tmp_path):
    source, _, _ = inputs(tmp_path)
    source.write_text(
        source.read_text().replace("[model]", "[model]\nbiases = marginalized")
    )
    with pytest.raises(ValueError, match="full_shape mode required"):
        Forecast(source)


def test_coverage_is_rejected_at_derivative_stencil(tmp_path):
    source, template, readers = inputs(tmp_path)
    full_ini(source)
    source.write_text(
        source.read_text()
        .replace("k_max = 0.1", "k_max = 0.99")
        .replace("k_order = 2", "k_order = 20")
        .replace("[numerical]", "[numerical]\nscale_step = 0.1")
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    with pytest.raises(ValueError, match="outside template domain"):
        forecast.run()


def test_shared_forest_binding_and_singular_spectra(tmp_path):
    from fishhighz.full_shape import make_model
    from fishhighz.models.templates import prepare_template

    source = tmp_path / "full.ini"
    source.write_text(
        Path("fishhighz/data/desi2_accuracy.ini").read_text()
        + "\n[model]\nmode = full_shape\n"
    )
    config = parse_survey_ini(source)
    registry = make_registry(config, Background())
    k = np.linspace(0.001, 1.0, 1000)
    template = prepare_template(
        k,
        100 + np.sin(100 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    fields = config.observed_fields
    biases = {field.id: 1.0 for field in fields}
    betas = {field.id: 1.45 for field in fields if field.kind == "forest"}
    model, binding = make_model(
        config,
        config.bins[1],
        registry,
        template,
        fields,
        biases,
        betas,
        {field.id: (0.0, 0.0) for field in fields},
        0.8,
        1.0,
    )
    global_indices = binding.binding.local_to_global
    assert (
        registry.ids[global_indices[model.local_names.index("bias_field_0")]]
        == "b_lya_1"
    )
    assert (
        registry.ids[global_indices[model.local_names.index("bias_field_4")]]
        == "b_lya_1"
    )
    assert (
        registry.ids[global_indices[model.local_names.index("beta_field_4")]]
        == "beta_lya_1"
    )
    kval = np.repeat(np.linspace(0.02, 0.2, 30), 12)
    mu = np.tile(np.linspace(0.02, 0.98, 12), 30)
    theta = registry.fiducials.copy()
    for pair, names in [
        ((0, 4), ("lya(qso)", "lya(lbg)")),
        ((0, 1), ("lya(qso)", "qso")),
    ]:
        ids = tuple(f"{name}_1" for name in active_names(config, 1, [names]))
        columns = []
        for name in ids:
            idx = registry.ids.index(name)
            step = registry.parameters[idx].step
            plus = theta.copy()
            minus = theta.copy()
            plus[idx] += step
            minus[idx] -= step
            columns.append(
                (
                    (
                        model(
                            plus[global_indices],
                            config.bins[1].z_eval,
                            kval,
                            mu,
                            [pair],
                        )
                        - model(
                            minus[global_indices],
                            config.bins[1].z_eval,
                            kval,
                            mu,
                            [pair],
                        )
                    )
                    / (2 * step)
                ).ravel()
            )
        jac = np.asarray(columns).T
        active = ParameterRegistry(
            [registry.parameters[registry.ids.index(name)] for name in ids]
        )
        result = FisherResult(active, jac.T @ jac)
        assert (
            reported_constraint(result, tuple(f"{name}_1" for name in TARGETS), 0.8)[0]
            == "unavailable"
        )


def test_default_steps_halved_on_synthetic_full_shape(tmp_path):
    from fishhighz.derivatives import check_convergence
    from fishhighz.models.templates import prepare_template

    source, _, readers = inputs(tmp_path)
    full_ini(source)
    k = np.linspace(0.001, 1.0, 1000)
    template = prepare_template(
        k,
        100 + np.sin(100 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    prepared = Forecast(
        source, background=Background(), template=template, readers=readers
    ).prepare()
    spec = prepared.survey.bins[0]
    result = check_convergence(
        spec.p3d,
        spec.p3d.registry.fiducials,
        prepared.config.bins[0].z_eval,
        spec.grid.k_flat,
        spec.grid.mu_flat,
        rtol=1e-3,
        atol=1e-9,
    )
    assert result.passed


def test_template_covering_all_nodes_but_not_cuts_is_rejected(tmp_path):
    from fishhighz.models.templates import prepare_template

    source, _, readers = inputs(tmp_path)
    full_ini(source)
    k = np.linspace(0.019, 0.091, 200)
    template = prepare_template(
        k,
        100 + np.sin(30 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    with pytest.raises(ValueError, match="smooth, fiducial: fixed observed cuts"):
        Forecast(
            source, background=Background(), template=template, readers=readers
        ).prepare()


@pytest.mark.parametrize("component,suffix", [("smooth", "s"), ("wiggle", "w")])
@pytest.mark.parametrize("one_sided", [False, True])
def test_perturbed_cut_coverage_uses_actual_schedule(
    tmp_path, component, suffix, one_sided
):
    from dataclasses import replace
    from types import SimpleNamespace

    from fishhighz.full_shape import validate_template_coverage
    from fishhighz.models.templates import prepare_template

    source, _, readers = inputs(tmp_path)
    full_ini(source)
    k = np.linspace(0.009999, 0.10001, 200)
    template = prepare_template(
        k,
        100 + np.sin(30 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    spec = forecast.prepare().survey.bins[0]
    # Keep all inactive dilation stencils safely inside support; test each
    # component independently, including a forward stencil at a lower bound.
    registry = ParameterRegistry(
        [
            replace(
                p,
                step=(0.00025 if p.id == f"alpha_{suffix}_0" else 1e-8),
                bounds=(1.0, None)
                if one_sided and p.id == f"alpha_{suffix}_0"
                else None,
            )
            for p in spec.p3d.registry.parameters
        ]
    )
    probe = SimpleNamespace(
        grid=spec.grid, p3d=SimpleNamespace(registry=registry, routes=spec.p3d.routes)
    )
    validate_template_coverage(probe, step_scale=0.1)
    with pytest.raises(
        ValueError,
        match=f"{component}, parameter alpha_{suffix}_0, "
        + ("forward" if one_sided else "central"),
    ):
        validate_template_coverage(probe, step_scale=1.0)


def test_public_custom_step_scale_changes_coverage(tmp_path):
    from fishhighz.models.templates import prepare_template

    source, _, readers = inputs(tmp_path)
    full_ini(source)
    k = np.linspace(0.009999, 0.10001, 200)
    template = prepare_template(
        k,
        100 + np.sin(30 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    forecast.run(step_scale=0.1)
    with pytest.raises(ValueError, match="parameter alpha_w_0"):
        forecast.run(step_scale=1.0)


def test_new_controls_are_optional_and_empty_bin_is_explicit(tmp_path):
    source, template, readers = inputs(tmp_path)
    full_ini(source)
    source.write_text(
        source.read_text()
        .replace("z_edges = 2.0, 2.2", "z_edges = 2.0, 2.2, 2.4")
        .replace("num_z_bins = 1", "num_z_bins = 2")
        .replace("selected = qsoxqso, qsoxlbg, lbgxlbg", "selected = ")
        + "\n[pairs bin 2]\nselected = qsoxqso, qsoxlbg, lbgxlbg\n"
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    result = forecast.run(individuals=False)
    assert result.joint[0].status == "excluded"
    assert result.joint[0].fisher is None
    assert result.joint[0].parameter_ids == ()
    assert result.joint[1].bin_index == 1
    assert len(result.individual) == 0
    assert len(result.excluded) == 3
    assert not any(name.endswith("_0") for name in result.combined.registry.ids)
    output = result.save(tmp_path / "empty_out")
    settings = json.loads((output / "settings.json").read_text())
    assert [item["bin_index"] for item in settings["records"][:2]] == [0, 1]
    assert settings["resolved_settings"]["selected_pairs_by_bin"]["0"] == []


def test_native_category_cut_sets_observed_domain_and_is_saved(tmp_path):
    source, template, readers = inputs(tmp_path)
    full_ini(source)
    source.write_text(
        source.read_text().replace(
            "[numerical]", "[numerical]\nk_max_galaxy_galaxy = 0.15"
        )
    )
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    prepared = forecast.prepare()
    assert prepared.survey.bins[0].grid.k_max == pytest.approx(0.15)
    result = forecast.run(individuals=False)
    assert len(result.joint) == 1 and len(result.individual) == 0
    output = result.save(tmp_path / "category_out")
    settings = json.loads((output / "settings.json").read_text())
    assert settings["resolved_settings"]["effective_k_max"]["galaxy_galaxy"] == 0.15
    assert settings["resolved_settings"]["run"]["individuals"] is False


def _three_category_bin():
    from fishhighz.fields import ObservedField, PairSelection
    from fishhighz.forecast import prepare_bin
    from fishhighz.geometry import prepare_geometry
    from fishhighz.grids import gauss_legendre_grid
    from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
    from fishhighz.response import InstrumentResponse
    from fishhighz.survey import BinSpec

    fields = (
        ObservedField("forest", "forest", "forest", background="qso"),
        ObservedField("galaxy", "galaxy", "galaxy"),
    )
    selection = PairSelection(fields)
    registry = ParameterRegistry(
        [
            Parameter("amplitude", 1.0, "target", step=1e-4),
            Parameter("forest_bias", 1.0, "nuisance", step=1e-4),
        ]
    )
    binding = BoundParameters(
        registry,
        ("amplitude", "forest_bias"),
        {"amplitude": "amplitude", "forest_bias": "forest_bias"},
    )
    amplitudes = {(0, 0): 1.2, (0, 1): 0.35, (1, 1): 0.9}
    forest_responses = {(0, 0): 0.7, (0, 1): 0.2, (1, 1): 0.0}

    def power(theta, redshift, k, mu, pairs):
        return np.column_stack(
            [
                np.full(
                    len(k),
                    theta[0] * amplitudes[tuple(pair)]
                    + theta[1] * forest_responses[tuple(pair)],
                )
                for pair in pairs
            ]
        )

    p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("linear", power, binding, selection.required_pairs)],
    )
    grid = gauss_legendre_grid(
        [0.01, 0.10, 0.15, 0.20], k_order=2, mu_order=2, h_fid=0.7
    )
    geometry = prepare_geometry(
        2.0,
        2.2,
        z_eval=2.1,
        area_deg2=100.0,
        h_fid=0.7,
        z_order=2,
        hubble=lambda z: np.full_like(z, 100.0),
        transverse_distance=lambda z: np.full_like(z, 1000.0),
    )
    spec = BinSpec(
        "categories",
        geometry,
        grid,
        p3d,
        {field.id: InstrumentResponse(0, 0) for field in fields},
        full_noise=np.tile([2.0, 0.1, 3.0], (len(grid.k_flat), 1)),
    )
    return prepare_bin(spec), amplitudes, forest_responses


@pytest.mark.parametrize("forest_cut,cross_cut", [(0.20, 0.15), (0.15, 0.20)])
def test_category_intervals_match_direct_covariance(forest_cut, cross_cut):
    from types import SimpleNamespace

    from fishhighz.covariance import gaussian_covariance
    from fishhighz.fields import PairSelection
    from fishhighz.forecast import run_bin
    from fishhighz.full_shape import _interval_result

    prepared, amplitudes, forest_responses = _three_category_bin()
    config = SimpleNamespace(
        numerical={
            "k_min": 0.01,
            "k_max": 0.20,
            "k_max_galaxy_galaxy": 0.10,
            "k_max_galaxy_forest": cross_cut,
            "k_max_forest_forest": forest_cut,
        }
    )
    actual = _interval_result(
        prepared, config, batch_size=4, step_scale=1.0, numerical=False
    )
    fields = prepared.p3d.selection.fields
    expected = np.zeros((2, 2))
    old_pairs = [tuple(pair) for pair in prepared.p3d.selection.required_pairs]
    for node, k in enumerate(prepared.k):
        selected_pairs = [
            pair
            for pair, cut in (((0, 0), forest_cut), ((0, 1), cross_cut), ((1, 1), 0.10))
            if k < cut
        ]
        selection = PairSelection(fields, selected_pairs)
        columns = [old_pairs.index(tuple(pair)) for pair in selection.required_pairs]
        covariance = gaussian_covariance(
            prepared.total[node, columns][None, :],
            prepared.modes[node : node + 1],
            selection,
        )[0]
        derivative = np.array(
            [
                [amplitudes[tuple(pair)], forest_responses[tuple(pair)]]
                for pair in selected_pairs
            ]
        )
        expected += derivative.T @ np.linalg.solve(covariance, derivative)
    np.testing.assert_allclose(actual.data_fisher, expected, rtol=2e-11)
    assert actual.marginalized_covariance(("amplitude",))[0, 0] == pytest.approx(
        np.linalg.inv(expected)[0, 0], rel=2e-11
    )
    own = _interval_result(
        prepared,
        config,
        pairs=[(0, 1)],
        batch_size=4,
        step_scale=1.0,
        numerical=False,
    ).data_fisher[0, 0]
    assert own > 0 and own != pytest.approx(actual)

    equal = SimpleNamespace(numerical={"k_min": 0.01, "k_max": 0.20})
    np.testing.assert_allclose(
        _interval_result(
            prepared, equal, batch_size=4, step_scale=1.0, numerical=False
        ).data_fisher,
        run_bin(prepared, batch_size=4).result.data_fisher,
        rtol=1e-12,
    )


@pytest.mark.parametrize("first_bin", [True, False])
def test_five_field_category_covariance_closure_and_first_bin(first_bin):
    """Retain covariance autos above their observable cuts in five-field fits."""
    from types import SimpleNamespace

    from fishhighz.covariance import gaussian_covariance
    from fishhighz.fields import ObservedField, PairSelection
    from fishhighz.forecast import prepare_bin
    from fishhighz.full_shape import _interval_result
    from fishhighz.geometry import prepare_geometry
    from fishhighz.grids import gauss_legendre_grid
    from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
    from fishhighz.response import InstrumentResponse
    from fishhighz.survey import BinSpec

    names = ("lya(qso)", "qso", "lbg", "lae", "lya(lbg)")
    fields = tuple(
        ObservedField(
            name,
            "forest" if name.startswith("lya(") else "galaxy",
            name,
            background="qso"
            if name == "lya(qso)"
            else "lbg"
            if name == "lya(lbg)"
            else None,
        )
        for name in names
    )
    all_pairs = [tuple(pair) for pair in PairSelection(fields).selected_pairs.tolist()]
    first_pairs = [(0, 0), (0, 1), (1, 1)]
    selection = PairSelection(fields, first_pairs if first_bin else all_pairs)
    registry = ParameterRegistry(
        [
            Parameter("dilation", 1.0, "target", step=1e-4),
            Parameter("forest_bias", 1.0, "nuisance", step=1e-4),
        ]
    )
    binding = BoundParameters(
        registry,
        ("dilation", "forest_bias"),
        {"dilation": "dilation", "forest_bias": "forest_bias"},
    )

    def coefficients(pair):
        i, j = pair
        return np.array([0.01 * (1 + i + j), 0.01 * (1 + i * j)])

    def power(theta, redshift, k, mu, pairs):
        return np.column_stack(
            [
                np.full(
                    len(k),
                    (3.0 if pair[0] == pair[1] else 0.1) + coefficients(pair) @ theta,
                )
                for pair in pairs
            ]
        )

    provider = P3DProvider("linear", power, binding, selection.required_pairs)
    p3d = PreparedP3D(registry, selection, [provider])
    grid = gauss_legendre_grid(
        [0.01, 0.10, 0.15, 0.20], k_order=2, mu_order=2, h_fid=0.7
    )
    geometry = prepare_geometry(
        2.0,
        2.2,
        z_eval=2.1,
        area_deg2=100.0,
        h_fid=0.7,
        z_order=2,
        hubble=lambda z: np.full_like(z, 100.0),
        transverse_distance=lambda z: np.full_like(z, 1000.0),
    )
    required = [tuple(pair) for pair in selection.required_pairs.tolist()]
    noise = np.tile(
        [0.5 if i == j else 0.0 for i, j in required], (len(grid.k_flat), 1)
    )
    spec = BinSpec(
        "five-fields",
        geometry,
        grid,
        p3d,
        {field.id: InstrumentResponse(0, 0) for field in fields},
        full_noise=noise,
    )
    prepared = prepare_bin(spec)
    config = SimpleNamespace(
        numerical={
            "k_min": 0.01,
            "k_max": 0.20,
            "k_max_galaxy_galaxy": 0.10,
            "k_max_galaxy_forest": 0.20,
            "k_max_forest_forest": 0.10,
        }
    )
    actual = _interval_result(
        prepared, config, batch_size=4, step_scale=1.0, numerical=False
    ).data_fisher
    expected = np.zeros((2, 2))
    for node, k in enumerate(prepared.k):
        active = [
            pair
            for pair in selection.selected_pairs.tolist()
            if k < (0.10 if (fields[pair[0]].kind == fields[pair[1]].kind) else 0.20)
        ]
        active = [tuple(pair) for pair in active]
        selected = PairSelection(fields, active)
        columns = [
            required.index(tuple(pair)) for pair in selected.required_pairs.tolist()
        ]
        covariance = gaussian_covariance(
            prepared.total[node, columns][None, :],
            prepared.modes[node : node + 1],
            selected,
        )[0]
        gradient = np.array([coefficients(pair) for pair in active])
        expected += gradient.T @ np.linalg.solve(covariance, gradient)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=1e-13)
    if first_bin:
        assert len(selection.selected_pairs) == 3
        assert any(k > 0.10 for k in prepared.k)
        assert (0, 0) in required and (1, 1) in required
    else:
        assert len(selection.selected_pairs) == 15
        assert any(k > 0.10 for k in prepared.k)


def test_forest_dilation_only_keeps_fiducial_growth_fixed(tmp_path):
    from fishhighz.full_shape import make_model
    from fishhighz.models.templates import prepare_template

    text = Path("fishhighz/data/desi2_accuracy.ini").read_text()
    text = text.replace(
        "selected = all",
        "selected = lya(qso)xlya(qso), lya(qso)xlya(lbg), lya(lbg)xlya(lbg)",
    )
    text = text.replace(
        "selected = lya(qso)xlya(qso), lya(qso)xqso, qsoxqso",
        "selected = lya(qso)xlya(qso)",
    )
    full_source = tmp_path / "forest_full.ini"
    full_source.write_text(text + "\n[model]\nmode = full_shape\n")
    reduced_source = tmp_path / "forest_dilation.ini"
    reduced_source.write_text(
        text + "\n[model]\nmode = full_shape\ntarget_set = dilation_only\n"
    )
    full_config = parse_survey_ini(full_source)
    reduced_config = parse_survey_ini(reduced_source)
    full_registry = make_registry(full_config, Background())
    reduced_registry = make_registry(reduced_config, Background())
    assert "f_0" in full_registry.ids and "f_0" not in reduced_registry.ids
    assert (
        sum(parameter.role == "target" for parameter in reduced_registry.parameters)
        == 24
    )
    k = np.linspace(0.001, 1.0, 300)
    template = prepare_template(
        k,
        100 + np.sin(100 * k),
        100 * np.ones_like(k),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    fields = reduced_config.observed_fields
    biases = {field.id: 1.0 for field in fields}
    betas = {field.id: 1.45 for field in fields if field.kind == "forest"}
    widths = {field.id: (0.0, 0.0) for field in fields}
    models = []
    for config, registry in (
        (full_config, full_registry),
        (reduced_config, reduced_registry),
    ):
        model, binding = make_model(
            config,
            config.bins[0],
            registry,
            template,
            fields,
            biases,
            betas,
            widths,
            0.8,
            1.0,
        )
        models.append(
            model(
                registry.fiducials[binding.binding.local_to_global],
                config.bins[0].z_eval,
                np.array([0.05]),
                np.array([0.5]),
                [(0, 0)],
            )
        )
    np.testing.assert_allclose(models[0], models[1], rtol=1e-13)
    invalid_source = tmp_path / "forest_invalid.ini"
    invalid_source.write_text(
        Path("fishhighz/data/desi2_accuracy.ini").read_text()
        + "\n[model]\nmode = full_shape\ntarget_set = dilation_only\n"
    )
    with pytest.raises(ValueError, match="forest-only"):
        parse_survey_ini(invalid_source)


@pytest.mark.parametrize("amplitude", ["alpha", "alpha_iso"])
def test_four_target_forest_result_round_trip(tmp_path, amplitude):
    from types import SimpleNamespace

    from fishhighz.public import SpectrumConstraint, SurveyResult

    names = (f"{amplitude}_w", "phi_w", f"{amplitude}_s", "phi_s", "b_lya", "beta_lya")
    registry = ParameterRegistry(
        [
            Parameter(name, 1.0, "target" if index < 4 else "nuisance")
            for index, name in enumerate(names)
        ]
    )
    fisher = FisherResult(registry, np.eye(6) * 3.0 + np.ones((6, 6)) * 0.1)
    status, covariance, errors, correlations = reported_constraint(
        fisher, names[:4], 0.8
    )
    record = SpectrumConstraint(
        1,
        "forest-2",
        (2.2, 2.4),
        2.3,
        "joint",
        None,
        names,
        status,
        fisher,
        None,
        None,
        None,
        names[:4],
        covariance,
        errors,
        correlations,
        (1.0, 1.0, 1.0, 1.0),
        0.8,
    )
    result = SurveyResult(
        SimpleNamespace(model={"mode": "full_shape"}),
        None,
        (),
        (record,),
        (),
        {},
        {"target_names": names[:4]},
        fisher,
    )
    output = result.save(tmp_path / "four_target")
    settings = json.loads((output / "settings.json").read_text())
    assert settings["records"][0]["target_ids"] == list(names[:4])
    assert len(settings["records"][0]["target_fiducials"]) == 4
    with np.load(output / "results.npz", allow_pickle=False) as arrays:
        np.testing.assert_allclose(arrays["record_0000_target_covariance"], covariance)
        assert arrays["record_0000_target_covariance"].shape == (4, 4)
