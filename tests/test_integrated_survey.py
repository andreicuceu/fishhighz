"""Integrated forest-source mode: INI resolution, survey preparation and limits."""

import configparser

import numpy as np
import pytest
from test_survey_config import _Background
from test_weights import FIELD, RESPONSE, geometry

import fishhighz
from fishhighz.accuracy import (
    INTEGRATED_DEFAULTS,
    INTEGRATED_REVISION,
    NATIVE_REVISION,
    REVISION,
    REVISION_DENSITY_INTERPOLATION,
)
from fishhighz.forecast import prepare_bin
from fishhighz.forest_integration import integrated_forest_source, integration_nodes
from fishhighz.geometry import SPEED_LIGHT_KMS
from fishhighz.models.templates import prepare_template
from fishhighz.public import PreparedForecast
from fishhighz.weights import (
    density_per_velocity,
    prepare_forest_weights,
    prepare_integrated_forest_weights,
)

BUNDLED = "fishhighz/data/desi2_accuracy.ini"
EXPANDED = "tests/data/desi2_accuracy_expanded.ini"
LYA_REST = 1215.67
SMALL_NUMERICAL = {
    "k_intervals": "2",
    "k_order": "2",
    "mu_order": "2",
    "z_order": "2",
    "magnitude_order": "4",
}


def _write_ini(tmp_path, changes, *, source=BUNDLED, name="modified.ini"):
    """Apply option edits to an INI and write the result.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory receiving the INI.
    changes : mapping
        Section-to-option mappings; a None value removes the option and a None
        section removes the section.
    source : str, optional
        Starting INI, default the bundled recipe.
    name : str, optional
        File name of the written INI.

    Returns
    -------
    path : pathlib.Path
        Written INI.
    """
    parser = configparser.ConfigParser(interpolation=None)
    parser.read(source)
    for section, values in changes.items():
        if values is None:
            parser.remove_section(section)
            continue
        if section not in parser:
            parser.add_section(section)
        for key, value in values.items():
            if value is None:
                parser.remove_option(section, key)
            else:
                parser[section][key] = value
    path = tmp_path / name
    with path.open("w") as stream:
        parser.write(stream)
    return path


def _parse(tmp_path, changes, **kwargs):
    """Write and parse a modified INI.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory receiving the INI.
    changes : mapping
        Edits, see ``_write_ini``.
    **kwargs : dict
        Forwarded to ``_write_ini``.

    Returns
    -------
    config : SurveyConfig
        Parsed configuration.
    """
    return fishhighz.parse_survey_ini(_write_ini(tmp_path, changes, **kwargs))


INTEGRATED = {"input policies": {"forest_source_integration": "integrated"}}


# ---------------------------------------------------------------------------
# Defaults: absent keys leave the parsed configuration and provenance unchanged.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("expanded", [False, True])
def test_central_configuration_has_no_integration_state(tmp_path, expanded):
    """No new key, field or provenance entry appears without the integrated mode."""
    config = _parse(tmp_path, {}, source=EXPANDED if expanded else BUNDLED)
    assert config.forest_integration is None
    assert "forest_source_integration" not in config.provenance
    assert "forest_source_integration" not in config.input_policies
    assert not set(INTEGRATED_DEFAULTS) & set(config.numerical)
    assert config.prescription is None or config.prescription["revision"] == (
        NATIVE_REVISION
    )
    for field in config.provenance["fields"]:
        assert "min_zq_forest" not in field and "max_zq_forest" not in field
    # The integrated defaults are not INI_DEFAULTS, so compact == expanded.
    compact = fishhighz.parse_survey_ini(BUNDLED)
    full = fishhighz.parse_survey_ini(EXPANDED)
    assert dict(compact.numerical) == dict(full.numerical)


@pytest.mark.parametrize(
    "revision", [None, NATIVE_REVISION, REVISION], ids=["default", "native", "hist"]
)
def test_existing_revisions_unchanged(tmp_path, revision):
    """Existing revisions keep their density treatment and an unchanged identity."""
    changes = {} if revision is None else {"prescription": {"revision": revision}}
    config = _parse(tmp_path, changes)
    expected = revision or NATIVE_REVISION
    assert config.prescription["revision"] == expected
    assert config.forest_integration is None
    assert (
        config.input_policies["density_interpolation"]
        == REVISION_DENSITY_INTERPOLATION[expected]
    )


def test_integrated_revision_density_map():
    """The integrated revision uses piecewise-constant cells."""
    assert (
        REVISION_DENSITY_INTERPOLATION[INTEGRATED_REVISION]
        == (REVISION_DENSITY_INTERPOLATION[NATIVE_REVISION])
    )
    assert INTEGRATED_REVISION == "early-lyaforecast-integrated-2026-10-02"


# ---------------------------------------------------------------------------
# Resolution and conflict matrix.
# ---------------------------------------------------------------------------


def test_integrated_key_implies_revision_and_records_block(tmp_path):
    """The key alone selects the revision, the default orders and the block."""
    config = _parse(tmp_path, INTEGRATED)
    assert config.prescription["revision"] == INTEGRATED_REVISION
    block = config.forest_integration
    assert block["mode"] == "integrated"
    assert block["revision"] == INTEGRATED_REVISION
    for key, value in INTEGRATED_DEFAULTS.items():
        assert block[key] == value
    assert set(block["fields"]) == {"lya(qso)", "lya(lbg)"}
    assert all(
        limits["min_zq_forest"] is None and limits["max_zq_forest"] is None
        for limits in block["fields"].values()
    )
    assert config.provenance["forest_source_integration"]["revision"] == (
        INTEGRATED_REVISION
    )
    assert config.provenance["prescription"]["revision"] == INTEGRATED_REVISION
    # The expanded policy and numerical mappings never contain the new keys.
    assert "forest_source_integration" not in config.input_policies
    assert not set(INTEGRATED_DEFAULTS) & set(config.numerical)
    assert config.input_policies["density_interpolation"] == "piecewise_constant_cells"


def test_integrated_revision_implies_integrated_mode(tmp_path):
    """The integrated revision alone, or with the key, selects the mode."""
    only_revision = _parse(
        tmp_path, {"prescription": {"revision": INTEGRATED_REVISION}}
    )
    both = _parse(
        tmp_path,
        {
            "prescription": {"revision": INTEGRATED_REVISION},
            **INTEGRATED,
        },
    )
    assert only_revision.forest_integration is not None
    assert dict(only_revision.provenance["forest_source_integration"]) == dict(
        both.provenance["forest_source_integration"]
    )


@pytest.mark.parametrize(
    "changes",
    [
        {
            "prescription": {"revision": INTEGRATED_REVISION},
            "input policies": {"forest_source_integration": "central"},
        },
        {
            "prescription": {"revision": NATIVE_REVISION},
            **INTEGRATED,
        },
        {
            "prescription": {"revision": REVISION},
            **INTEGRATED,
        },
    ],
    ids=["central+integrated revision", "integrated+native", "integrated+historical"],
)
def test_mode_revision_conflicts_rejected(tmp_path, changes):
    """Conflicting mode and revision raise ValueError."""
    with pytest.raises(ValueError):
        _parse(tmp_path, changes)


def test_unknown_mode_and_missing_prescription_rejected(tmp_path):
    """Unknown mode values and a missing [prescription] are rejected."""
    with pytest.raises(ValueError, match="unsupported"):
        _parse(tmp_path, {"input policies": {"forest_source_integration": "both"}})
    with pytest.raises(ValueError, match="requires a \\[prescription\\]"):
        _parse(tmp_path, {"prescription": None, **INTEGRATED}, source=EXPANDED)


def test_central_key_equals_absent_key(tmp_path):
    """Explicit central mode is the same as no key, with no new provenance."""
    explicit = _parse(
        tmp_path, {"input policies": {"forest_source_integration": "central"}}
    )
    default = _parse(tmp_path, {})
    assert explicit.forest_integration is None
    assert explicit.provenance["numerical"] == default.provenance["numerical"]
    assert explicit.provenance["input_policies"] == default.provenance["input_policies"]
    assert "forest_source_integration" not in explicit.provenance


def test_spline_density_with_integrated_mode_rejected(tmp_path):
    """Integrated mode needs piecewise-constant cells."""
    spline = {
        "input policies": {
            "forest_source_integration": "integrated",
            "density_interpolation": "RectBivariateSpline_kx2_ky2_s0",
            "magnitude_partition": "density_knots_support_snr_nodes_negative_roots",
        }
    }
    with pytest.raises(ValueError, match="piecewise-constant"):
        _parse(tmp_path, spline)


@pytest.mark.parametrize(
    "changes",
    [
        {"numerical": {"forest_zq_order": "8"}},
        {"numerical": {"forest_lambda_order": "8"}},
        {"numerical": {"forest_lambda_panels": "2"}},
        {"field lya(qso)": {"min_zq_forest": "2.1"}},
        {"field lya(lbg)": {"max_zq_forest": "3.9"}},
    ],
)
def test_integrated_only_keys_rejected_in_central_mode(tmp_path, changes):
    """Integrated-only keys require the integrated mode."""
    with pytest.raises(ValueError, match="integrated-only"):
        _parse(tmp_path, changes)
    with pytest.raises(ValueError, match="integrated-only"):
        _parse(
            tmp_path,
            {**changes, "input policies": {"forest_source_integration": "central"}},
        )


@pytest.mark.parametrize("field", ["qso", "lbg", "lae"])
def test_zq_keys_rejected_on_galaxy_fields(tmp_path, field):
    """Source-redshift limits are only defined for forest fields."""
    changes = {**INTEGRATED, f"field {field}": {"min_zq_forest": "2.1"}}
    with pytest.raises(ValueError, match="forest fields only"):
        _parse(tmp_path, changes)


@pytest.mark.parametrize(
    "changes,match",
    [
        ({"numerical": {"forest_zq_order": "0"}}, "positive integer"),
        ({"numerical": {"forest_lambda_order": "3.5"}}, "positive integer"),
        ({"numerical": {"forest_lambda_panels": "x"}}, "positive integer"),
        ({"field lya(qso)": {"min_zq_forest": "nan"}}, "finite"),
        ({"field lya(qso)": {"min_zq_forest": "-1"}}, "nonnegative"),
        (
            {"field lya(qso)": {"min_zq_forest": "3", "max_zq_forest": "2.5"}},
            "must exceed",
        ),
    ],
)
def test_integrated_values_validated(tmp_path, changes, match):
    """Orders and limits must be valid."""
    with pytest.raises(ValueError, match=match):
        _parse(tmp_path, {**INTEGRATED, **changes})


def test_explicit_orders_and_limits_are_recorded(tmp_path):
    """Explicit settings reach the config block and its provenance."""
    config = _parse(
        tmp_path,
        {
            **INTEGRATED,
            "numerical": {
                "forest_zq_order": "6",
                "forest_lambda_order": "7",
                "forest_lambda_panels": "3",
            },
            "field lya(qso)": {"min_zq_forest": "2.1", "max_zq_forest": "3.9"},
        },
    )
    block = config.provenance["forest_source_integration"]
    assert (
        block["forest_zq_order"],
        block["forest_lambda_order"],
        block["forest_lambda_panels"],
    ) == (6, 7, 3)
    assert block["fields"]["lya(qso)"]["min_zq_forest"] == 2.1
    assert block["fields"]["lya(qso)"]["max_zq_forest"] == 3.9
    assert block["fields"]["lya(lbg)"]["min_zq_forest"] is None


# ---------------------------------------------------------------------------
# Synthetic readers with an exact dn/dv_q = const limit.
# ---------------------------------------------------------------------------

MAGNITUDE_AXIS = np.array([16.1, 26.75])


def _amplitude(magnitudes):
    """Magnitude dependence of the synthetic dn/dz (deg^-2 z^-1 mag^-1)."""
    return 0.5 * 10 ** (0.15 * (np.asarray(magnitudes) - 20.0))


def _variance(magnitudes):
    """Magnitude-only dimensionless pixel variance."""
    return 0.2 * 10 ** (0.3 * (np.asarray(magnitudes) - 20.0))


class _ConstantVelocityDensity:
    """dn/dz = K(m)/(1+z): constant dn/dv_q, supporting scalar and grid queries."""

    magnitudes = MAGNITUDE_AXIS
    z_edges = np.array([0.0, 10.0])

    def sample(self, z, magnitudes):
        """Return the scalar-redshift density."""
        values = _amplitude(magnitudes) / (1 + z)
        return {"values": values, "provenance": {"source": "exact"}}

    def sample_grid(self, z, magnitudes):
        """Return the density on a (z, magnitude) grid."""
        values = _amplitude(magnitudes)[None, :] / (1 + np.asarray(z))[:, None]
        return {"values": values, "provenance": {"counts": {}}}


class _MagnitudeOnlySNR:
    """Variance independent of redshift, wavelength and pixel width."""

    magnitudes = MAGNITUDE_AXIS
    z = np.array([0.0, 10.0])

    def sample(self, *, magnitudes, **_):
        """Return the scalar-query variance."""
        return {"values": _variance(magnitudes), "provenance": {}}

    def variance_grid(self, *, z_source, magnitudes, **_):
        """Return the variance on paired points times magnitudes."""
        values = np.tile(_variance(magnitudes), (len(z_source), 1))
        return {"values": values, "provenance": {"counts": {}}}


def _template():
    """Return a constant synthetic template."""
    k_grid = np.linspace(0.001, 1.0, 20)
    return prepare_template(
        k_grid,
        np.ones_like(k_grid),
        np.full_like(k_grid, 0.5),
        z_ref=2.406,
        h_template=0.7,
        h_fid=0.7,
    )


def _small_config(tmp_path, *, integrated, extra=None):
    """Parse a one-bin forest-by-forest and forest-by-QSO configuration.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory receiving the INI.
    integrated : bool
        Whether to request the integrated mode (with no source-redshift limit,
        so that the window is untruncated for the synthetic density).
    extra : mapping, optional
        Further edits applied after the defaults of this helper.

    Returns
    -------
    config : SurveyConfig
        One-bin (2.0-2.26) configuration with lya(qso) and qso.
    """
    changes = {
        "survey": {"z_edges": "2.0, 2.26", "num_z_bins": "1"},
        "fields": {"ids": "lya(qso), qso"},
        "numerical": dict(SMALL_NUMERICAL),
        "pairs bin 1": {"selected": "lya(qso)xlya(qso), lya(qso)xqso"},
    }
    for section in ("lbg", "lae", "lya(lbg)"):
        changes[f"field {section}"] = None
    for index in range(2, 7):
        changes[f"pairs bin {index}"] = None
    changes["field lya(qso)"] = {"z_norm_min": None}
    if integrated:
        changes["input policies"] = {"forest_source_integration": "integrated"}
        changes["field lya(qso)"]["min_zq_forest"] = "0"
    for section, values in (extra or {}).items():
        changes.setdefault(section, {}).update(values)
    return _parse(tmp_path, changes)


def _readers(config):
    """Return synthetic readers for every configured field."""
    return {
        field.observed.id: {
            "density": _ConstantVelocityDensity(),
            "snr": _MagnitudeOnlySNR() if field.observed.kind == "forest" else None,
        }
        for field in config.fields
    }


def _prepare(config):
    """Prepare the survey and its first bin with the synthetic readers."""
    survey = fishhighz.prepare_survey(
        config,
        background=_Background(),
        template=_template(),
        readers=_readers(config),
    )
    return survey, prepare_bin(survey.bins[0])


def test_integrated_equals_central_for_constant_velocity_density(tmp_path):
    """Constant dn/dv_q and magnitude-only variance: integrated == central.

    With dn/dz = K/(1+z) the integrand (1+z) dn/dz is constant, so the
    integrals over y = ln(1+z_q) and ln(lambda) are exact for Gauss-Legendre
    with breakpoints at the overlap kinks, and the integrated mode must
    reproduce L I_1, A, P_pixel and the magnitude weights of the central mode.
    The residual is floating-point summation order; the expected agreement is
    6e-14 (relative, measured for A and P_pixel; 2e-16 for the weights) and the
    asserted tolerance 1e-11.
    """
    central_survey, central_bin = _prepare(_small_config(tmp_path, integrated=False))
    integrated_survey, integrated_bin = _prepare(
        _small_config(tmp_path, integrated=True)
    )
    reference = central_bin.weights["lya(qso)"]
    prepared = integrated_bin.weights["lya(qso)"]
    tolerance = dict(rtol=1e-11, atol=0)
    np.testing.assert_allclose(prepared.A, reference.A, **tolerance)
    np.testing.assert_allclose(prepared.P_pixel, reference.P_pixel, **tolerance)
    length = SPEED_LIGHT_KMS * np.log(1205.0 / 1040.0)
    np.testing.assert_allclose(prepared.N1, length * reference.I1[-1], **tolerance)
    for row in prepared.weights:
        np.testing.assert_allclose(row, reference.weights, **tolerance)

    # The noise power entering the covariance is therefore the same.
    np.testing.assert_allclose(
        integrated_bin.noise, central_bin.noise, rtol=1e-10, atol=0
    )
    assert integrated_survey.bins[0].forests["lya(qso)"].integrated is not None
    assert central_survey.bins[0].forests["lya(qso)"].integrated is None


def test_integrated_survey_provenance_and_z_eff(tmp_path):
    """Per-bin integration records, resolved limits and z_eff are in provenance."""
    config = _small_config(tmp_path, integrated=True)
    survey, prepared_bin = _prepare(config)
    provenance = survey.provenance
    assert provenance["forest_source_integration"]["revision"] == INTEGRATED_REVISION
    assert provenance["forest_source_integration"]["resolved_zq_limits"][
        "lya(qso)"
    ] == (0.0, 10.0)
    block = provenance["bins"][0]["forest_source_integration"]
    assert block["spec_id"] == survey.bins[0].id
    info = survey.bins[0].forests["lya(qso)"].integrated.info
    assert block["n_pixel"]["lya(qso)"] == info["n_pixel"]
    assert block["n_zq_nodes"]["lya(qso)"] == info["n_zq_nodes"]
    assert "z_eff" not in block

    forecast = PreparedForecast(
        config, survey, _Background(), _template(), (prepared_bin,)
    )
    completed = forecast.provenance["bins"][0]["forest_source_integration"]
    weights = prepared_bin.weights["lya(qso)"]
    assert completed["z_eff"]["lya(qso)"] == weights.z_eff
    assert completed["N1"]["lya(qso)"] == weights.N1
    # z_eff lies inside the bin redshift range of the absorbing gas.
    assert 2.0 < weights.z_eff < 2.26
    # Central-mode provenance is returned unchanged.
    central_survey, central_bin = _prepare(_small_config(tmp_path, integrated=False))
    central_forecast = PreparedForecast(
        central_survey.config,
        central_survey,
        _Background(),
        _template(),
        (central_bin,),
    )
    assert central_forecast.provenance is central_survey.provenance
    assert "forest_source_integration" not in central_survey.provenance


def test_central_survey_provenance_has_no_integration_keys(tmp_path):
    """Central preparation records no integrated-mode keys anywhere."""
    survey, _ = _prepare(_small_config(tmp_path, integrated=False))
    assert "forest_source_integration" not in survey.provenance
    for record in survey.provenance["bins"]:
        assert "forest_source_integration" not in record


def test_integrated_mode_requires_grid_readers(tmp_path):
    """Readers without grid queries are rejected with a field and bin label."""

    class _ScalarOnly:
        magnitudes = MAGNITUDE_AXIS
        z_edges = np.array([0.0, 10.0])

        def sample(self, z, magnitudes):
            """Return the scalar-redshift density."""
            return {"values": _amplitude(magnitudes), "provenance": {}}

    config = _small_config(tmp_path, integrated=True)
    readers = _readers(config)
    readers["lya(qso)"]["density"] = _ScalarOnly()
    with pytest.raises(ValueError, match=r"field lya\(qso\), bin 1.*sample_grid"):
        fishhighz.prepare_survey(
            config, background=_Background(), template=_template(), readers=readers
        )


def test_no_forest_coverage_names_field_and_bin(tmp_path):
    """A source-redshift window outside the tabulated density raises."""
    config = _small_config(
        tmp_path,
        integrated=True,
        extra={"field lya(qso)": {"min_zq_forest": "9", "max_zq_forest": "9.5"}},
    )
    with pytest.raises(
        ValueError, match=r"field lya\(qso\), bin 1.*no forest coverage"
    ):
        fishhighz.prepare_survey(
            config,
            background=_Background(),
            template=_template(),
            readers=_readers(config),
        )


# ---------------------------------------------------------------------------
# Physical normalisation at the WP1/WP2 boundary (independent of survey_config).
# ---------------------------------------------------------------------------


def test_real_integrated_source_matches_central_normalisation():
    """(1+z_q), c and 1/L_bin factors of the measure, end to end.

    A real IntegratedForestSource (untruncated window, dn/dz_q = K(m)/(1+z_q)
    evaluated analytically at the y nodes, variance depending on magnitude
    only) must give A, P_pixel and the converged magnitude weights of the
    central weights with rho = dn/dv_q = K/c and L = c ln(rest_max/rest_min).
    A missing (1+z_q), c or 1/L_bin factor would rescale A and P_pixel.
    Agreement is at floating-point summation level (measured 6e-14 for A and
    P_pixel, 3e-14 for N1); asserted at rtol = 1e-12.
    """
    rest_min, rest_max = 1040.0, 1205.0
    lambda_min, lambda_max = LYA_REST * 3.0, LYA_REST * 3.26
    nodes = integration_nodes(
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        rest_min=rest_min,
        rest_max=rest_max,
        zq_min=0.0,
        zq_max=10.0,
        zq_breaks=[],
        zq_order=6,
        lambda_order=5,
        lambda_panels=2,
        lya_rest_angstrom=LYA_REST,
    )
    magnitudes = np.array([19.0, 21.0, 23.0, 24.5])
    quadrature = np.array([0.4, 0.9, 1.1, 0.6])
    amplitude, variance = _amplitude(magnitudes), _variance(magnitudes)
    density = amplitude[None, :] / (1 + nodes.zq_nodes)[:, None]
    pixel_width = RESPONSE.pixel_width_velocity
    source = integrated_forest_source(
        nodes,
        density,
        magnitudes,
        quadrature,
        np.tile(variance, (len(nodes.lam_obs), 1)),
        pixel_width,
    )
    options = dict(method="early_lyaforecast", signal=3.0, alias=2.0)
    prepared = prepare_integrated_forest_weights(
        FIELD, geometry(), RESPONSE, source, **options
    )

    z_source = 2.9
    reference = prepare_forest_weights(
        FIELD,
        geometry(),
        RESPONSE,
        z_source=z_source,
        magnitudes=magnitudes,
        quadrature=quadrature,
        rho=density_per_velocity(amplitude / (1 + z_source), z_source=z_source),
        variance=variance,
        length_velocity=SPEED_LIGHT_KMS * np.log(rest_max / rest_min),
        **options,
    )
    tolerance = dict(rtol=1e-12, atol=0)
    np.testing.assert_allclose(prepared.A, reference.A, **tolerance)
    np.testing.assert_allclose(prepared.P_pixel, reference.P_pixel, **tolerance)
    length = SPEED_LIGHT_KMS * np.log(rest_max / rest_min)
    np.testing.assert_allclose(prepared.N1, length * reference.I1[-1], **tolerance)
    for row in prepared.weights:
        np.testing.assert_allclose(row, reference.weights, **tolerance)
    # The window is untruncated: z_eff is the mean absorption redshift of the bin.
    assert 2.0 < prepared.z_eff < 2.26


# ---------------------------------------------------------------------------
# Breakpoints of the survey wiring and refinement ladder.
# ---------------------------------------------------------------------------


class _TableDensity(_ConstantVelocityDensity):
    """Cell-edge density with a nonuniform redshift tiling."""

    z_edges = np.array([1.9, 2.15, 2.35, 2.5, 2.9, 3.4])


class _TableSNR(_MagnitudeOnlySNR):
    """S/N source-redshift nodes distinct from the density cell edges."""

    z = np.array([2.0, 2.3, 2.45, 2.7, 3.1, 3.3])

    def variance_grid(self, *, z_source, wavelength, magnitudes, **_):
        """Return a variance with smooth pixel and source-redshift dependence."""
        factor = np.exp(
            10 * (np.asarray(wavelength) / (LYA_REST * 3.13) - 1)
            + 4 * (np.asarray(z_source) - 2.6)
        )
        values = _variance(magnitudes)[None, :] * factor[:, None]
        return {"values": values, "provenance": {"counts": {}}}


def test_survey_breakpoints_include_density_edges_and_snr_nodes(tmp_path):
    """Both table sets of an SRD-like bin are panel boundaries of the y nodes."""
    config = _small_config(
        tmp_path,
        integrated=True,
        extra={"field lya(qso)": {"min_zq_forest": "1.9"}},
    )
    readers = _readers(config)
    for field in ("lya(qso)", "qso"):
        readers[field]["density"] = _TableDensity()
    readers["lya(qso)"]["snr"] = _TableSNR()
    survey = fishhighz.prepare_survey(
        config, background=_Background(), template=_template(), readers=readers
    )
    info = survey.bins[0].forests["lya(qso)"].integrated.info["nodes"]
    y_breakpoints = np.asarray(info["y_breakpoints"])
    window = info["window_zq"]
    wanted = np.r_[_TableDensity.z_edges, _TableSNR.z]
    inside = wanted[(wanted > window[0]) & (wanted < window[1])]
    assert len(inside) >= 4
    for z_break in inside:
        assert np.min(np.abs(y_breakpoints - np.log1p(z_break))) < 1e-12


def _bundled_density_sources(order):
    """Integrated N1..N3 of a smooth synthetic source for one order triple.

    Parameters
    ----------
    order : tuple of int
        (zq_order, lambda_order, lambda_panels).

    Returns
    -------
    moments : ndarray of shape (3,)
        N1, N2 and N3 of the early-lyaforecast weights.
    """
    zq_order, lambda_order, lambda_panels = order
    # Density cells (z_q) with piecewise-constant steps and a smooth
    # pixel-dependent variance: only the quadrature, not the integrand, changes.
    z_edges = np.array([1.9, 2.15, 2.35, 2.5, 2.9, 3.4])
    steps = np.array([0.6, 1.0, 1.7, 1.1, 0.5])
    nodes = integration_nodes(
        lambda_min=LYA_REST * 3.0,
        lambda_max=LYA_REST * 3.26,
        rest_min=1040.0,
        rest_max=1205.0,
        zq_min=1.9,
        zq_max=3.4,
        zq_breaks=z_edges,
        zq_order=zq_order,
        lambda_order=lambda_order,
        lambda_panels=lambda_panels,
        lya_rest_angstrom=LYA_REST,
    )
    magnitudes = np.array([19.0, 21.0, 23.0, 24.5])
    quadrature = np.array([0.4, 0.9, 1.1, 0.6])
    cell = np.clip(np.searchsorted(z_edges, nodes.zq_nodes, side="right") - 1, 0, 4)
    density = steps[cell][:, None] * _amplitude(magnitudes)[None, :]
    # Steep but smooth pixel and source-redshift dependence of the variance.
    pixel_factor = np.exp(
        10 * (nodes.lam_obs / (LYA_REST * 3.13) - 1) + 4 * (nodes.z_q - 2.6)
    )
    variance = _variance(magnitudes)[None, :] * pixel_factor[:, None]
    source = integrated_forest_source(
        nodes,
        density,
        magnitudes,
        quadrature,
        variance,
        RESPONSE.pixel_width_velocity,
    )
    prepared = prepare_integrated_forest_weights(
        FIELD,
        geometry(),
        RESPONSE,
        source,
        method="early_lyaforecast",
        signal=3.0,
        alias=2.0,
        iterations=4,
    )
    return np.array([prepared.N1, prepared.N2, prepared.N3])


def test_synthetic_refinement_ladder_converges():
    """N1-N3 converge along (4,8,1) -> (8,16,1) -> (16,32,2) at fixed weights.

    Breakpoints at the density cell edges make the integrand smooth in every
    panel, so the quadrature converges geometrically; the last step must change
    the moments by less than 1e-9 relative, far below the 0.1 per cent target.
    """
    ladder = [(4, 8, 1), (8, 16, 1), (16, 32, 2)]
    moments = np.array([_bundled_density_sources(order) for order in ladder])
    changes = np.abs(moments[1:] / moments[:-1] - 1)
    assert np.all(changes[0] > changes[1])
    assert np.all(changes[1] < 1e-9)


def _ladder_sigma(tmp_path, order):
    """BAO errors of the one-bin synthetic survey at one order triple.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory receiving the INI.
    order : tuple of int
        (zq_order, lambda_order, lambda_panels).

    Returns
    -------
    sigma : ndarray of shape (2,)
        Joint sigma(ap), sigma(at) of the bin.
    """
    from fishhighz.forecast import run_bin

    config = _small_config(
        tmp_path,
        integrated=True,
        extra={
            "numerical": {
                "forest_zq_order": str(order[0]),
                "forest_lambda_order": str(order[1]),
                "forest_lambda_panels": str(order[2]),
            },
            "field lya(qso)": {"min_zq_forest": "0"},
        },
    )
    readers = _readers(config)
    for field in ("lya(qso)", "qso"):
        readers[field]["density"] = _TableDensity()
    readers["lya(qso)"]["snr"] = _TableSNR()
    survey = fishhighz.prepare_survey(
        config, background=_Background(), template=_template(), readers=readers
    )
    run = run_bin(prepare_bin(survey.bins[0]), step_scale=1.0)
    return np.asarray(run.result.fix_except(("ap_0", "at_0")).marginalized_errors())


def test_survey_refinement_ladder_bao_errors_converge(tmp_path):
    """BAO errors change by < 0.1 per cent between the last two orders."""
    ladder = [(4, 8, 1), (8, 16, 1), (16, 32, 2)]
    sigma = np.array([_ladder_sigma(tmp_path, order) for order in ladder])
    relative = np.abs(sigma[2] / sigma[1] - 1)
    assert np.all(relative < 1e-3)
    assert np.all(np.isfinite(sigma)) and np.all(sigma > 0)
