"""Synthetic checks for the public INI-driven forecast facade."""

import json
from types import SimpleNamespace

import numpy as np

from fishhighz import Forecast
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.forecast import prepare_bin
from fishhighz.geometry import prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.models.templates import prepare_template
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.public import SurveyResult, _pair_spec
from fishhighz.response import InstrumentResponse
from fishhighz.survey import BinSpec

INI = """\
[schema]
name = fishhighz-native-survey
version = 1
[cosmology]
camb_ini = package:camb_configs/Planck18.ini
template = package:templates/Planck18_z_2.406.fits
template_redshift = 2.406
damping_reference_redshift = 2.3
[survey]
area_deg2 = 10
band = r
min_band_mag = 16.1
max_band_mag = 26.75
z_edges = 2.0, 2.2
num_z_bins = 1
resolution = 2500
reconstruction_factor = 2
lya_rest_angstrom = 1215.67
evaluation_redshift = geometric_1plusz
[fields]
ids = qso, lbg
[field qso]
kind = galaxy
physical_model = qso
tracer = qso
density = package:not-used-qso
target_density = 90
min_band_mag = 16.1
max_band_mag = 26.75
density_magnitude_bounds = none
[field lbg]
kind = galaxy
physical_model = lbg
tracer = lbg
density = package:not-used-lbg
target_density = 380
min_band_mag = 16.1
max_band_mag = 26.75
density_magnitude_bounds = none
[model]
parameterization = ap_at
parameter_names = ap, at
smooth_scaling = identity
wiggle_scaling = ap_at
damping = mixed_squared_width
growth = camb_sigma8_ratio
reconstruction = per_field
forest_beta = 1.45
damping_amplitude = 3.26
forest_reconstruction_factor = 1.0
[input policies]
density_semantics = cell_count_per_deg2
density_width_policy = legacy_first_spacing
density_redshift_normalization = target_density
density_negative_policy = floor_negative
density_interpolation = RectBivariateSpline_kx2_ky2_s0
snr_smoothing = legacy
snr_interpolation = linear_RegularGridInterpolator
snr_bright_policy = clamp_to_brightest_tabulated_magnitude
snr_clamp = 1e-10
snr_sentinel = 1e20
magnitude_partition = density_knots_support_snr_nodes_negative_roots
weighting_method = early_lyaforecast
weighting_reference = fiducial_auto_p3d_and_p1d_times_response_squared
weighting_reference_k_t_deg = 2.4
weighting_reference_k_p_velocity = 0.00035
weighting_rtol = 1e-5
weighting_min_updates = 3
weighting_stable_steps = 3
weighting_max_updates = 96
sampling_noise = independent_diagonal
[numerical]
k_intervals = 2
k_order = 2
k_min = 0.01
k_max = 0.1
mu_order = 3
z_order = 2
magnitude_order = 4
ap_step = 2.5e-4
weight_rtol = 1e-5
[pairs bin 1]
selected = qsoxqso, qsoxlbg, lbgxlbg
"""


class Background:
    h_fid = 0.7
    template_growth_redshift = 2.406
    damping_reference_redshift = 2.3
    sigma8_damping_reference = 0.8

    def sigma8_at(self, redshift):
        return 0.8

    def growth_rate_at(self, redshift):
        return 0.8

    def hubble_parameter(self, redshift):
        return np.full(np.asarray(redshift).shape, 100.0)

    def transverse_comoving_distance(self, redshift):
        return np.full(np.asarray(redshift).shape, 1000.0)


class Density:
    magnitudes = np.array([16.1, 20.0, 26.75])

    def sample(self, redshift, magnitudes):
        return {"values": np.ones(len(magnitudes)), "provenance": {"synthetic": True}}


def inputs(tmp_path):
    source = tmp_path / "survey.ini"
    source.write_text(INI)
    k = np.linspace(0.001, 1.0, 20)
    template = prepare_template(
        k,
        np.ones_like(k),
        np.full_like(k, 0.5),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    readers = {name: {"density": Density()} for name in ("qso", "lbg")}
    return source, template, readers


def test_forecast_parsing_is_lazy_and_prepare_accepts_injected_factories(
    tmp_path, monkeypatch
):
    source, template, readers = inputs(tmp_path)
    calls = []

    def factory(config):
        calls.append(config.schema_name)
        return Background()

    forecast = Forecast(
        source,
        background_factory=factory,
        template=template,
        readers=readers,
    )
    assert forecast._prepared is None
    assert calls == []
    prepared = forecast.prepare()
    assert calls == ["fishhighz-native-survey"]
    assert prepared.prepared_bins[0].p3d.registry is prepared.survey.registry
    monkeypatch.setattr(
        "fishhighz.public.prepare_camb",
        lambda **_: (_ for _ in ()).throw(AssertionError()),
    )


def test_bundled_recipe_name_is_resolved_without_cwd_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    forecast = Forecast("desi2_accuracy.ini")
    assert forecast.config.path is None
    (tmp_path / "desi2_accuracy.ini").write_text("[cosmo]\n")
    try:
        Forecast("desi2_accuracy.ini")
    except ValueError as error:
        assert "lyaforecast INI" in str(error)
    else:
        raise AssertionError(
            "an existing relative INI must not be replaced by the bundle"
        )


def test_forecast_run_retains_fisher_results_and_own_covariance(tmp_path, monkeypatch):
    source, template, readers = inputs(tmp_path)
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    original = __import__("fishhighz.public", fromlist=["prepare_bin"]).prepare_bin
    selected_lengths = []

    def checked(spec):
        selected_lengths.append(len(spec.p3d.selection.selected_pairs))
        return original(spec)

    monkeypatch.setattr("fishhighz.public.prepare_bin", checked)
    result = forecast.run(batch_size=32)
    assert isinstance(result, SurveyResult)
    assert len(result.individual) == 3
    assert len(result.joint) == 1
    assert len(result.excluded) == 0
    assert all(item.fisher is not None for item in result.individual + result.joint)
    # One three-spectrum joint preparation is followed by three one-spectrum
    # preparations; the latter own their covariance factors.
    assert selected_lengths[0] == 3
    assert selected_lengths[1:] == [1, 1, 1]
    assert result.joint[0].fisher.data_fisher.shape == (2, 2)
    output = result.save(tmp_path / "result")
    assert (output / "settings.json").is_file()
    assert (output / "results.npz").is_file()
    settings = json.loads((output / "settings.json").read_text())
    assert settings["resolved_settings"]["run"] == {
        "batch_size": 32,
        "step_scale": 0.25,
        "numerical": False,
    }


def test_pair_spec_remaps_packed_full_noise_to_new_required_pairs():
    registry = ParameterRegistry([Parameter("A", 1.0, "target", step=0.01)])
    fields = tuple(ObservedField(name, "galaxy", name) for name in ("a", "b", "c"))
    selection = PairSelection(fields)

    def model(theta, redshift, k, mu, pairs):
        return np.full((len(k), len(pairs)), theta[0])

    binding = BoundParameters(registry, ("A",), {"A": "A"})
    p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("synthetic", model, binding, selection.required_pairs)],
    )
    grid = gauss_legendre_grid([0.02, 0.1], k_order=2, mu_order=2, h_fid=0.7)
    geometry = prepare_geometry(
        2.0,
        3.0,
        z_eval=2.5,
        area_deg2=10.0,
        h_fid=0.7,
        z_order=2,
        hubble=lambda z: np.full_like(z, 100.0),
        transverse_distance=lambda z: np.full_like(z, 1000.0),
    )
    spec = BinSpec(
        "packed",
        geometry,
        grid,
        p3d,
        {field.id: InstrumentResponse(0, 0) for field in fields},
        full_noise=np.tile([1.0, 0.5, 2.0, 4.0, 0.5, 6.0], (len(grid.k_flat), 1)),
    )
    one = _pair_spec(spec, (0, 2))
    np.testing.assert_array_equal(one.full_noise, spec.full_noise[:, [0, 2, 5]])
    prepared = prepare_bin(one)
    np.testing.assert_array_equal(prepared.noise, one.full_noise)


def test_convergence_serialization_retains_solver_diagnostics():
    from fishhighz.public import _convergence

    diagnostics = {
        "status": "capped",
        "reason": "no confirmed convergence",
        "updates": 96,
        "state_updates": 96,
        "candidate": None,
        "last_step": {"amplitude": 1e-6},
        "confirmation": None,
        "forward_residual": None,
    }
    value = _convergence(
        [
            SimpleNamespace(
                id="bin",
                weights={
                    "forest": SimpleNamespace(convergence=diagnostics),
                },
            )
        ]
    )
    assert dict(value["bin"]["forest"]) == diagnostics


def test_public_cli_uses_forecast_result_save(tmp_path, monkeypatch):
    from fishhighz import cli

    saved = []

    class FakeResult:
        def save(self, path):
            saved.append(path)

    class FakeForecast:
        def __init__(self, source):
            assert source == "input.ini"

        def run(self):
            return FakeResult()

    monkeypatch.setattr(cli, "Forecast", FakeForecast)
    assert cli.main(["input.ini", "--output", str(tmp_path / "out")]) == 0
    assert saved == [str(tmp_path / "out")]


def test_saved_inputs_survive_deleted_ini_and_distinguish_injection(tmp_path):
    import hashlib

    source, template, readers = inputs(tmp_path)
    source.write_text(INI.replace("target_density = 90", "target_density = 91"))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    source.unlink()
    result = forecast.run()
    output = result.save(tmp_path / "captured")
    settings = json.loads((output / "settings.json").read_text())["resolved_settings"]
    assert settings["config"]["fields"][0]["target_density"] == 91
    assert settings["inputs"]["ini"]["sha256"] == digest
    assert settings["inputs"]["background"]["source"] == "injected"
    assert settings["inputs"]["template"]["source"] == "injected"
    assert settings["inputs"]["fields"]["qso"]["density"] == {
        "source": "injected",
        "type": "Density",
        "sha256": None,
    }


def test_configured_readers_capture_all_source_hashes():
    import hashlib
    from contextlib import ExitStack

    from fishhighz.public import _reader_identity
    from fishhighz.resources import bundled_resource
    from fishhighz.survey_config import _normalise_readers, parse_survey_ini

    config = parse_survey_ini()
    with ExitStack() as stack:
        readers = _normalise_readers(config, config.fields, None, stack)
        identities = _reader_identity(config, readers, "configured")
    for field in config.fields:
        density = identities[field.observed.id]["density"]
        assert density["identifier"] == field.density
        assert (
            density["sha256"]
            == hashlib.sha256(
                bundled_resource(field.density.removeprefix("package:")).read_bytes()
            ).hexdigest()
        )
        if field.observed.kind == "forest":
            snr = identities[field.observed.id]["snr"]
            assert len(snr["identifiers"]) == len(
                readers[field.observed.id]["snr"].reader.provenance["paths"]
            )
            for identifier, digest in zip(
                snr["identifiers"], snr["sha256"], strict=True
            ):
                assert identifier.startswith(field.snr + "/")
                assert (
                    digest
                    == hashlib.sha256(
                        bundled_resource(
                            identifier.removeprefix("package:")
                        ).read_bytes()
                    ).hexdigest()
                )


def test_default_cosmology_and_template_input_hashes(tmp_path, monkeypatch):
    import hashlib
    from contextlib import ExitStack

    from fishhighz.public import (
        _default_background,
        _default_template,
        _object_identity,
    )
    from fishhighz.resources import bundled_resource
    from fishhighz.survey_config import parse_survey_ini

    config = parse_survey_ini()
    monkeypatch.setattr("fishhighz.public.prepare_camb", lambda **kwargs: Background())
    identity = {}
    with ExitStack() as stack:
        background = _default_background(config, stack, identity)
        template = _default_template(config, background, stack)
    assert identity["identifier"] == config.cosmology["camb_ini"]
    assert (
        identity["sha256"]
        == hashlib.sha256(
            bundled_resource(
                config.cosmology["camb_ini"].removeprefix("package:")
            ).read_bytes()
        ).hexdigest()
    )
    assert (
        _object_identity(template, "configured")["sha256"]
        == hashlib.sha256(
            bundled_resource(
                config.cosmology["template"].removeprefix("package:")
            ).read_bytes()
        ).hexdigest()
    )


def test_supplied_instances_take_precedence_over_factories(tmp_path):
    source, template, readers = inputs(tmp_path)

    def unused(*args):
        raise AssertionError("factory must not replace an explicit input")

    prepared = Forecast(
        source,
        background=Background(),
        template=template,
        readers=readers,
        background_factory=unused,
        template_factory=unused,
        readers_factory=unused,
    ).prepare()
    assert prepared.input_identity["background"]["source"] == "injected"
    assert prepared.input_identity["template"]["source"] == "injected"
    assert prepared.input_identity["fields"]["qso"]["density"]["source"] == "injected"


def test_saved_exposure_and_bias_settings(tmp_path):
    source, template, readers = inputs(tmp_path)
    text = INI.replace(
        "[field qso]\nkind = galaxy\nphysical_model = qso\ntracer = qso",
        "[field qso]\nkind = forest\nphysical_model = lya\ntracer = lya\n"
        "background = qso\nsnr = package:not-used-snr\nnum_exposures = 7\n"
        "pix_width_angstrom = 0.8\nmin_rest_frame_lya = 1040\nmax_rest_frame_lya = 1205",
    ).replace(
        "target_density = 380",
        "target_density = 380\nbias_z = 2, 4\nbias_values = 3.7, 3.8",
    )
    source.write_text(text)

    class SNR:
        magnitudes = Density.magnitudes

        def sample(self, **kwargs):
            return {"values": np.ones(len(kwargs["magnitudes"])), "provenance": {}}

    readers["qso"]["snr"] = SNR()
    result = Forecast(
        source, background=Background(), template=template, readers=readers
    ).run()
    output = result.save(tmp_path / "survey-settings")
    fields = json.loads((output / "settings.json").read_text())["resolved_settings"][
        "config"
    ]["fields"]
    assert fields[0]["num_exposures"] == 7
    assert fields[1]["bias_values"] == [3.7, 3.8]


def test_reader_digest_is_not_recomputed_when_saving(tmp_path):
    import hashlib

    source, template, readers = inputs(tmp_path)
    table = tmp_path / "density.dat"
    table.write_text("original input")
    digest = hashlib.sha256(table.read_bytes()).hexdigest()
    readers["qso"]["density"].provenance = {"path": str(table), "sha256": digest}
    forecast = Forecast(
        source, background=Background(), template=template, readers=readers
    )
    forecast.prepare()
    table.write_text("modified input")
    readers["qso"]["density"].provenance["sha256"] = "mutated reader"
    table.unlink()
    result = forecast.run()
    output = result.save(tmp_path / "reader-identity")
    identity = json.loads((output / "settings.json").read_text())["resolved_settings"][
        "inputs"
    ]["fields"]["qso"]["density"]
    assert identity["sha256"] == digest
    assert identity["identifier"] == str(table)
