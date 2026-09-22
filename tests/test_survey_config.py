"""Deterministic controls for the native Step-3 survey boundary."""

from contextlib import ExitStack
from pathlib import Path

import numpy as np
import pytest

import fishhighz
from fishhighz.magnitude import composite
from fishhighz.models.templates import prepare_template
from fishhighz.survey_config import _build_readers, _resolve_resource


def test_native_recipe_identity_and_selection():
    config = fishhighz.parse_survey_ini()
    assert config.schema_version == 1
    assert [field.observed.id for field in config.fields] == [
        "lya(qso)",
        "qso",
        "lbg",
        "lae",
        "lya(lbg)",
    ]
    assert [len(item.selected_pairs) for item in config.bins] == [3, 15, 15, 15, 15, 15]
    expected = [(2.0 + 0.235 * index, 2.0 + 0.235 * (index + 1)) for index in range(6)]
    assert np.allclose(
        [(item.z_min, item.z_max) for item in config.bins], expected, rtol=0, atol=1e-14
    )
    for item in config.bins:
        assert item.z_eval == pytest.approx(
            np.sqrt((1 + item.z_min) * (1 + item.z_max)) - 1
        )
    assert config.numerical["k_intervals"] == "128"
    assert config.numerical["k_order"] == "4"
    assert config.numerical["mu_order"] == "32"
    assert config.numerical["z_order"] == "32"
    assert config.numerical["magnitude_order"] == "16"
    assert config.numerical["ap_step"] == "2.5e-4"
    assert config.numerical["weight_rtol"] == "1e-5"
    assert config.input_policies["density_negative_policy"] == "floor_negative"
    assert config.input_policies["snr_smoothing"] == "legacy"
    assert config.input_policies["weighting_reference_k_t_deg"] == "2.4"
    assert config.input_policies["weighting_reference_k_p_velocity"] == "0.00035"
    assert config.fields[0].density_magnitude_bounds == "survey"
    assert config.fields[1].density_magnitude_bounds == "none"
    assert config.fields[2].z_norm_min is None


def test_bundled_reader_policies_match_accuracy_example():
    config = fishhighz.parse_survey_ini()
    with ExitStack() as stack:
        readers = _build_readers(config, config.fields, stack)
        qso = readers["qso"]["density"].reader.provenance
        lbg = readers["lbg"]["density"].reader.provenance
        lae = readers["lae"]["density"].reader.provenance
        lya_qso = readers["lya(qso)"]["density"].reader.provenance
        lya_lbg = readers["lya(lbg)"]["density"].reader.provenance
    assert qso["z_norm_min"] == pytest.approx(2.15)
    assert qso["magnitude_bounds"] is None
    assert lbg["z_norm_min"] is None
    assert lbg["magnitude_bounds"] is None
    assert lae["z_norm_min"] is None
    assert lae["magnitude_bounds"] is None
    assert lya_qso["z_norm_min"] == pytest.approx(2.15)
    assert lya_qso["magnitude_bounds"] == pytest.approx((16.1, 26.75))
    assert lya_lbg["z_norm_min"] == pytest.approx(2.15)
    assert lya_lbg["magnitude_bounds"] == pytest.approx((16.1, 26.75))


def test_legacy_ini_is_not_translated(tmp_path):
    legacy = tmp_path / "legacy.ini"
    legacy.write_text("[cosmo]\nfilename = Planck18.ini\n[tracer 1]\ntracer = qso\n")
    with pytest.raises(fishhighz.UnsupportedSchemaError, match="lyaforecast INI"):
        fishhighz.parse_survey_ini(legacy)


def test_schema_rejects_missing_sections_garbage_lists_and_fixed_labels(tmp_path):
    source = Path("tests/data/desi2_accuracy_expanded.ini").read_text()
    missing = tmp_path / "missing.ini"
    missing.write_text(source.replace("\n[numerical]\n", "\n[removed]\n", 1))
    with pytest.raises(fishhighz.UnsupportedSchemaError, match="missing required"):
        fishhighz.parse_survey_ini(missing)
    garbage = tmp_path / "garbage.ini"
    garbage.write_text(
        source.replace("z_edges = 2.0, 2.235", "z_edges = 2.0, garbage", 1)
    )
    with pytest.raises(ValueError, match="without garbage"):
        fishhighz.parse_survey_ini(garbage)
    altered = tmp_path / "altered.ini"
    altered.write_text(
        source.replace(
            "snr_interpolation = linear_RegularGridInterpolator",
            "snr_interpolation = cubic",
            1,
        )
    )
    with pytest.raises(ValueError, match="snr_interpolation"):
        fishhighz.parse_survey_ini(altered)
    altered = tmp_path / "altered-snr-policy.ini"
    altered.write_text(source.replace("snr_clamp = 1e-10", "snr_clamp = 1e-9", 1))
    with pytest.raises(ValueError, match="snr_clamp"):
        fishhighz.parse_survey_ini(altered)


def test_package_and_ini_relative_resource_resolution(tmp_path):
    config = fishhighz.parse_survey_ini()
    assert _resolve_resource(config.cosmology["camb_ini"], config) == (
        "package",
        "camb_configs/Planck18.ini",
    )
    external = tmp_path / "survey.ini"
    external.write_text(Path("fishhighz/data/desi2_accuracy.ini").read_text())
    external_config = fishhighz.parse_survey_ini(external)
    # Native package references remain package references after relocation;
    # an ordinary relative path is anchored at the INI directory.
    assert _resolve_resource("tables/density.dat", external_config) == (
        "external",
        tmp_path / "tables/density.dat",
    )


def test_composite_magnitude_quadrature_is_ordered_and_exact():
    nodes, weights = composite((0.0, 1.0, 3.0), 4)
    assert np.all(np.diff(nodes) > 0)
    assert np.sum(weights) == pytest.approx(3.0)
    assert np.sum(nodes * weights) == pytest.approx(4.5)


class _Background:
    h_fid = 0.7
    sigma8_damping_reference = 0.8
    template_growth_redshift = 2.406
    damping_reference_redshift = 2.3

    def sigma8_at(self, redshift):
        assert redshift > 0
        return 0.8

    def growth_rate_at(self, redshift):
        assert redshift > 0
        return 0.8

    def hubble_parameter(self, redshift):
        return np.full(np.asarray(redshift).shape, 100.0)

    def transverse_comoving_distance(self, redshift):
        return np.full(np.asarray(redshift).shape, 1000.0)


class _Density:
    magnitudes = np.array([16.1, 20.0, 26.75])

    def sample(self, redshift, magnitudes):
        return {"values": np.ones(len(magnitudes)), "provenance": {"source": "fake"}}


class _SNR:
    magnitudes = np.array([16.1, 20.0, 26.75])

    def sample(self, **kwargs):
        return {
            "values": np.ones(len(kwargs["magnitudes"])),
            "provenance": {"source": "fake"},
        }


class _BackgroundMissingTemplateGrowth(_Background):
    @property
    def template_growth_redshift(self):
        raise AttributeError("template-growth metadata is absent")


class _BackgroundMissingDampingReference(_Background):
    @property
    def damping_reference_redshift(self):
        raise AttributeError("damping-reference metadata is absent")


class _BackgroundMissingDampingSigma(_Background):
    @property
    def sigma8_damping_reference(self):
        raise AttributeError("damping sigma8 metadata is absent")


def _synthetic_template_and_readers(config):
    k = np.linspace(0.001, 1.0, 20)
    template = prepare_template(
        k,
        np.ones_like(k),
        np.full_like(k, 0.5),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    readers = {
        field.observed.id: {
            "density": _Density(),
            "snr": _SNR() if field.observed.kind == "forest" else None,
        }
        for field in config.fields
    }
    return template, readers


def test_small_injected_preparation_preserves_covariance_closure_and_registry():
    config = fishhighz.parse_survey_ini()
    k = np.linspace(0.001, 1.0, 20)
    template = prepare_template(
        k,
        np.ones_like(k),
        np.full_like(k, 0.5),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    readers = {
        field.observed.id: {
            "density": _Density(),
            "snr": _SNR() if field.observed.kind == "forest" else None,
        }
        for field in config.fields
    }
    prepared = fishhighz.prepare_survey(
        config, background=_Background(), template=template, readers=readers
    )
    assert len(prepared.bins) == 6
    assert prepared.registry.ids == (
        "ap_0",
        "at_0",
        "ap_1",
        "at_1",
        "ap_2",
        "at_2",
        "ap_3",
        "at_3",
        "ap_4",
        "at_4",
        "ap_5",
        "at_5",
    )
    for item in prepared.bins:
        assert item.p3d.registry is prepared.registry
        assert len(item.p3d.selection.required_pairs) >= len(
            item.p3d.selection.selected_pairs
        )
        assert set(item.responses) == {field.id for field in prepared.fields}
    assert len(prepared.bins[0].p3d.selection.required_pairs) == 3
    assert len(prepared.bins[1].p3d.selection.required_pairs) == 15


def test_injected_background_and_template_must_match_configured_redshifts():
    config = fishhighz.parse_survey_ini()
    k = np.linspace(0.001, 1.0, 20)
    template = prepare_template(
        k,
        np.ones_like(k),
        np.full_like(k, 0.5),
        z_ref=2.407,
        h_template=0.7,
        h_fid=0.7,
    )
    readers = {
        field.observed.id: {
            "density": _Density(),
            "snr": _SNR() if field.observed.kind == "forest" else None,
        }
        for field in config.fields
    }
    with pytest.raises(ValueError, match="template z_ref/configured"):
        fishhighz.prepare_survey(
            config, background=_Background(), template=template, readers=readers
        )


def test_injected_background_reference_redshifts_are_checked():
    config = fishhighz.parse_survey_ini()
    k = np.linspace(0.001, 1.0, 20)
    template = prepare_template(
        k,
        np.ones_like(k),
        np.full_like(k, 0.5),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )
    readers = {
        field.observed.id: {
            "density": _Density(),
            "snr": _SNR() if field.observed.kind == "forest" else None,
        }
        for field in config.fields
    }
    background = _Background()
    background.damping_reference_redshift = 2.31
    with pytest.raises(ValueError, match="damping_reference_redshift"):
        fishhighz.prepare_survey(
            config, background=background, template=template, readers=readers
        )


@pytest.mark.parametrize(
    ("background_type", "missing"),
    [
        (_BackgroundMissingTemplateGrowth, "template_growth_redshift"),
        (_BackgroundMissingDampingReference, "damping_reference_redshift"),
        (_BackgroundMissingDampingSigma, "sigma8_damping_reference"),
    ],
)
def test_background_damping_metadata_is_required(background_type, missing):
    config = fishhighz.parse_survey_ini()
    template, readers = _synthetic_template_and_readers(config)
    with pytest.raises(ValueError, match=missing):
        fishhighz.prepare_survey(
            config,
            background=background_type(),
            template=template,
            readers=readers,
        )


@pytest.mark.parametrize(
    ("attribute", "value", "message"),
    [
        ("template_growth_redshift", np.nan, "template_growth_redshift"),
        ("damping_reference_redshift", np.inf, "damping_reference_redshift"),
        ("sigma8_damping_reference", 0.0, "sigma8_damping_reference"),
    ],
)
def test_background_damping_metadata_must_be_finite_and_positive(
    attribute, value, message
):
    config = fishhighz.parse_survey_ini()
    template, readers = _synthetic_template_and_readers(config)
    background = _Background()
    setattr(background, attribute, value)
    with pytest.raises(ValueError, match=message):
        fishhighz.prepare_survey(
            config, background=background, template=template, readers=readers
        )


def test_native_modules_do_not_import_validation_or_external_forecast_packages():
    source = "\n".join(
        Path(path).read_text()
        for path in (
            "fishhighz/accuracy.py",
            "fishhighz/magnitude.py",
            "fishhighz/survey_config.py",
        )
    )
    assert "import fishhighz.validation" not in source
    assert "from fishhighz.validation" not in source
    assert "import lyaforecast" not in source
    assert "import vega" not in source


def _write_modified_ini(tmp_path, changes, *, expanded=False):
    import configparser

    source = (
        "tests/data/desi2_accuracy_expanded.ini"
        if expanded
        else "fishhighz/data/desi2_accuracy.ini"
    )
    parser = configparser.ConfigParser(interpolation=None)
    parser.read(source)
    for section, values in changes.items():
        if section not in parser:
            parser.add_section(section)
        for key, value in values.items():
            if value is None:
                parser.remove_option(section, key)
            else:
                parser[section][key] = value
    path = tmp_path / "modified.ini"
    with path.open("w") as stream:
        parser.write(stream)
    return fishhighz.parse_survey_ini(path)


def test_named_prescription_matches_original_expanded_ini():
    from fishhighz.accuracy import REVISION

    compact = fishhighz.parse_survey_ini()
    original = fishhighz.parse_survey_ini("tests/data/desi2_accuracy_expanded.ini")
    for name in (
        "cosmology",
        "survey",
        "model",
        "input_policies",
        "numerical",
        "fields",
        "bins",
    ):
        assert getattr(compact, name) == getattr(original, name)
    assert compact.provenance["prescription"] == {
        "name": "accuracy",
        "revision": REVISION,
    }
    assert original.prescription is None


@pytest.mark.parametrize(
    "section,key", [("numerical", "weight_rtol"), ("input policies", "weighting_rtol")]
)
def test_single_weight_tolerance_override(tmp_path, section, key):
    config = _write_modified_ini(tmp_path, {section: {key: "2e-5"}})
    assert config.numerical["weight_rtol"] == "2e-5"
    assert config.input_policies["weighting_rtol"] == "2e-5"


def test_conflicting_weight_tolerances_rejected(tmp_path):
    with pytest.raises(ValueError, match="must match"):
        _write_modified_ini(
            tmp_path,
            {
                "numerical": {"weight_rtol": "2e-5"},
                "input policies": {"weighting_rtol": "3e-5"},
            },
        )


@pytest.mark.parametrize(
    "changes",
    [
        {"numerical": {"mu_orders": "4"}},
        {"prescription": {"name": "unknown"}},
        {"prescription": {"revision": "future"}},
        {"model": {"damping": "unknown"}},
    ],
)
def test_named_prescription_rejects_unknown_choices(tmp_path, changes):
    with pytest.raises(ValueError):
        _write_modified_ini(tmp_path, changes)


@pytest.mark.parametrize("field", ["qso", "lya(qso)"])
def test_per_field_magnitude_selection_not_silently_ignored(tmp_path, field):
    with pytest.raises(ValueError, match="unsupported per-field magnitude selection"):
        _write_modified_ini(
            tmp_path, {f"field {field}": {"min_band_mag": "22", "max_band_mag": "23"}}
        )
    config = _write_modified_ini(
        tmp_path, {f"field {field}": {"min_band_mag": "16.1", "max_band_mag": "26.75"}}
    )
    assert config.fields == fishhighz.parse_survey_ini().fields


def test_field_provenance_records_all_scientific_attributes(tmp_path):
    from dataclasses import fields

    config = _write_modified_ini(
        tmp_path,
        {
            "field qso": {"target_density": "80"},
            "field lya(qso)": {"num_exposures": "2"},
            "field lbg": {"bias_values": "3.5, 3.6"},
        },
    )
    for field, recorded in zip(config.fields, config.provenance["fields"]):
        for attribute in fields(field):
            if attribute.name != "observed":
                assert recorded[attribute.name] == getattr(field, attribute.name)
    assert config.provenance["fields"][0]["num_exposures"] == 2
    assert config.provenance["fields"][1]["target_density"] == 80
    assert config.provenance["fields"][2]["bias_values"] == (3.5, 3.6)


def test_compact_and_expanded_prepare_identical_small_numerical_inputs(tmp_path):
    changes = {
        "numerical": {
            "k_intervals": "2",
            "k_order": "2",
            "mu_order": "2",
            "z_order": "2",
            "magnitude_order": "2",
        }
    }
    compact = _write_modified_ini(tmp_path, changes)
    expanded = _write_modified_ini(tmp_path, changes, expanded=True)
    preparations = []
    for config in (compact, expanded):
        template, readers = _synthetic_template_and_readers(config)
        preparations.append(
            fishhighz.prepare_survey(
                config, background=_Background(), template=template, readers=readers
            )
        )
    for left, right in zip(preparations[0].bins, preparations[1].bins):
        for name, value in vars(left.geometry).items():
            np.testing.assert_equal(value, getattr(right.geometry, name))
        assert left.galaxies == right.galaxies
        for name in ("k", "mu", "weights", "q_mode"):
            np.testing.assert_array_equal(
                getattr(left.grid, name), getattr(right.grid, name)
            )
        for field in left.forests:
            for name, value in left.forests[field].weight_options.items():
                np.testing.assert_equal(
                    value, right.forests[field].weight_options[name]
                )
    assert preparations[0].bins[0].grid.k.size == 4


@pytest.mark.parametrize(
    "omitted_section,omitted_key,override_section,override_key",
    [
        ("numerical", "weight_rtol", "input policies", "weighting_rtol"),
        ("input policies", "weighting_rtol", "numerical", "weight_rtol"),
    ],
)
def test_expanded_ini_can_specify_weight_tolerance_once(
    tmp_path, omitted_section, omitted_key, override_section, override_key
):
    config = _write_modified_ini(
        tmp_path,
        {
            omitted_section: {omitted_key: None},
            override_section: {override_key: "2e-5"},
        },
        expanded=True,
    )
    assert config.numerical["weight_rtol"] == "2e-5"
    assert config.input_policies["weighting_rtol"] == "2e-5"


def test_omitted_field_magnitude_limits_inherit_survey_normalization(tmp_path):
    config = _write_modified_ini(
        tmp_path, {"survey": {"min_band_mag": "17", "max_band_mag": "25"}}
    )
    assert config.fields[0].magnitude_bounds == (17.0, 25.0)
    assert config.fields[1].magnitude_bounds is None
    assert config.provenance["survey"]["min_band_mag"] == "17"
