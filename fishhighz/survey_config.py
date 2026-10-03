"""Strict native INI parsing and survey preparation.

This module is the small configuration boundary for FishHighz.  It owns no
forecast result or command-line state: parsing validates the versioned native
schema, while preparation turns it into the existing immutable survey records.
Neither operation imports ``fishhighz.validation``, lyaforecast, or Vega.
"""

from __future__ import annotations

import configparser
import hashlib
from contextlib import ExitStack
from dataclasses import dataclass
from dataclasses import fields as dataclass_fields
from pathlib import Path
from types import MappingProxyType

import numpy as np

from .accuracy import (
    CONTROLS,
    DENSITY_INTERPOLATION_POLICIES,
    INI_DEFAULTS,
    INTEGRATED_DEFAULTS,
    INTEGRATED_REVISION,
    NATIVE_REVISION,
    REVISION_DENSITY_INTERPOLATION,
)
from .fields import ObservedField, PairSelection
from .forest_integration import (
    integrated_forest_source,
    integration_nodes,
    pixel_width_angstrom,
)
from .geometry import SPEED_LIGHT_KMS, prepare_geometry
from .grids import gauss_legendre_grid
from .magnitude import breakpoints, composite
from .models.biases import linear_tabulated_bias
from .models.external import BoundParameters, P3DProvider, PreparedP3D
from .models.kaiser import KaiserModel, Scaling
from .models.p1d import default_p1d
from .noise import local_galaxy_density
from .parameters import Parameter, ParameterRegistry
from .resources import bundled_path, bundled_paths, resolve_input_path
from .response import (
    InstrumentResponse,
    pixel_width_angstrom_to_velocity,
    resolving_power_fwhm_to_sigma,
)
from .survey import BinSpec, ForestInput, freeze
from .weights import density_per_velocity

SCHEMA_NAME = "fishhighz-native-survey"
SCHEMA_VERSION = 1


class UnsupportedSchemaError(ValueError):
    """Raised when an INI is not the versioned native FishHighz schema."""


def _plain_mapping(value):
    """Copy a mapping into a read-only view.

    Parameters
    ----------
    value : mapping or iterable of pairs
        Configuration values to snapshot.

    Returns
    -------
    mapping : mappingproxy
        Read-only shallow copy of the input.
    """
    return MappingProxyType(dict(value))


def _float_list(value, name, minimum=1):
    """Parse a finite comma- or whitespace-separated numeric list.

    Parameters
    ----------
    value : str
        Numeric tokens separated by commas or whitespace.
    name : str
        Configuration label for errors.
    minimum : int, optional
        Minimum number of entries; default 1.

    Returns
    -------
    values : tuple of float
        Finite values in input order, with units set by the configuration
        option.

    Raises
    ------
    ValueError
        If parsing fails or too few finite values are supplied.
    """
    try:
        tokens = value.replace(",", " ").split()
        result = np.asarray([float(token) for token in tokens], dtype=float)
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError(
            f"{name}: require a finite numeric list without garbage"
        ) from error
    if result.ndim != 1 or len(result) < minimum or not np.all(np.isfinite(result)):
        raise ValueError(
            f"{name}: require a finite list with at least {minimum} values"
        )
    return tuple(float(x) for x in result)


def _tokens(value, name):
    """Parse a nonempty list of configuration tokens.

    Parameters
    ----------
    value : str
        Tokens separated by commas or whitespace.
    name : str
        Configuration label for errors.

    Returns
    -------
    tokens : tuple of str
        Nonempty tokens in input order.

    Raises
    ------
    ValueError
        If the token list is empty.
    """
    values = tuple(x.strip() for x in value.replace(",", " ").split() if x.strip())
    if not values:
        raise ValueError(f"{name}: require a nonempty list")
    return values


def _section_options(parser, section, allowed, required=()):
    """Validate required and allowed options in one INI section.

    Parameters
    ----------
    parser : configparser.ConfigParser
        Parsed INI configuration.
    section : str
        Existing section to validate.
    allowed : iterable of str
        Supported option names.
    required : iterable of str, optional
        Mandatory option names; default empty.

    Returns
    -------
    None
        No value is returned.

    Raises
    ------
    ValueError
        If options are unsupported or required options are missing.
    """
    actual = set(parser[section])
    unknown = actual - set(allowed)
    if unknown:
        raise ValueError(
            f"[{section}]: unsupported option(s): {', '.join(sorted(unknown))}"
        )
    missing = set(required) - actual
    if missing:
        raise ValueError(
            f"[{section}]: missing required option(s): {', '.join(sorted(missing))}"
        )


@dataclass(frozen=True)
class FieldConfig:
    """One observed field and its explicit raw-input conventions."""

    observed: ObservedField
    tracer: str
    density: str
    target_density: float
    z_norm_min: float | None
    magnitude_bounds: tuple[float, float] | None
    density_magnitude_bounds: str
    snr: str | None
    num_exposures: float | None
    pixel_width_angstrom: float | None
    min_rest_frame_lya: float | None
    max_rest_frame_lya: float | None
    bias_redshifts: tuple[float, ...] | None
    bias_values: tuple[float, ...] | None


@dataclass(frozen=True)
class BinConfig:
    index: int
    z_min: float
    z_max: float
    z_eval: float
    selected_pairs: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class SurveyConfig:
    """Validated native INI values, independent of any prepared background."""

    path: Path | None
    schema_name: str
    schema_version: int
    cosmology: object
    survey: object
    model: object
    input_policies: object
    numerical: object
    fields: tuple[FieldConfig, ...]
    bins: tuple[BinConfig, ...]
    prescription: object = None
    source_identity: object = None
    forest_integration: object = None

    @property
    def observed_fields(self):
        """Return observed fields in configured order.

        Returns
        -------
        fields : tuple of ObservedField
            Ordered observed populations, including separate forest backgrounds.
        """
        return tuple(field.observed for field in self.fields)

    @property
    def provenance(self):
        """Snapshot all resolved native configuration values.

        Returns
        -------
        provenance : mappingproxy
            Immutable schema, input identity, model, field and bin metadata.
            A ``forest_source_integration`` block (mode, revision, quadrature
            orders and requested per-field source-redshift limits) is present
            only for the integrated forest-source mode; central-mode provenance
            has no such key.
        """
        record = {
            "schema": {"name": self.schema_name, "version": self.schema_version},
            "prescription": (
                None if self.prescription is None else dict(self.prescription)
            ),
            "ini_input": self.source_identity,
            "ini_path": None if self.path is None else str(self.path),
            "cosmology": dict(self.cosmology),
            "survey": dict(self.survey),
            "model": dict(self.model),
            "input_policies": dict(self.input_policies),
            "numerical": dict(self.numerical),
            "fields": [
                {
                    "id": field.observed.id,
                    "kind": field.observed.kind,
                    "physical_model": field.observed.physical_model,
                    "background": field.observed.background,
                    **{
                        attribute.name: getattr(field, attribute.name)
                        for attribute in dataclass_fields(field)
                        if attribute.name != "observed"
                    },
                }
                for field in self.fields
            ],
            "bins": [
                {
                    "index": item.index,
                    "bounds": [item.z_min, item.z_max],
                    "z_eval": item.z_eval,
                    "selected_pairs": [list(pair) for pair in item.selected_pairs],
                }
                for item in self.bins
            ],
        }
        if self.forest_integration is not None:
            record["forest_source_integration"] = self.forest_integration
        return freeze(record)


@dataclass(frozen=True)
class PreparedSurvey:
    """Native survey records ready for the Step-4 forecast facade."""

    config: SurveyConfig
    fields: tuple[ObservedField, ...]
    registry: ParameterRegistry
    bins: tuple[BinSpec, ...]
    provenance: object

    @property
    def bin_specs(self):
        """Return the prepared survey bin specifications.

        Returns
        -------
        bins : tuple of BinSpec
            Nonempty configured bins and their prepared model/noise inputs.
        """

        return self.bins


def _read_parser(source):
    """Read an INI and retain its exact input digest.

    Parameters
    ----------
    source : path-like or None
        Native INI path, or None for bundled desi2_accuracy.ini.

    Returns
    -------
    parser : configparser.ConfigParser
        Strict parser with interpolation disabled.
    path : pathlib.Path or None
        Resolved external path, or None for the package resource.
    identity : mappingproxy
        Stable input identifier and SHA256 digest.

    Raises
    ------
    FileNotFoundError
        If the external input does not exist.
    """
    if source is None:
        with bundled_path("desi2_accuracy.ini") as path:
            parser = configparser.ConfigParser(interpolation=None, strict=True)
            content = path.read_bytes()
            parser.read_string(content.decode("utf-8"))
            return (
                parser,
                None,
                _plain_mapping(
                    {
                        "identifier": "package:desi2_accuracy.ini",
                        "sha256": hashlib.sha256(content).hexdigest(),
                    }
                ),
            )
    path = Path(source).expanduser().resolve(strict=True)
    parser = configparser.ConfigParser(interpolation=None, strict=True)
    content = path.read_bytes()
    parser.read_string(content.decode("utf-8"))
    return (
        parser,
        path,
        _plain_mapping(
            {"identifier": str(path), "sha256": hashlib.sha256(content).hexdigest()}
        ),
    )


def _density_partition(density_interpolation):
    """Return the magnitude partition paired with a density interpolation.

    Parameters
    ----------
    density_interpolation : str
        Value of ``[input policies] density_interpolation``.

    Returns
    -------
    partition : str
        Magnitude-partition policy that must accompany the interpolation.

    Raises
    ------
    ValueError
        If the interpolation is not one of ``DENSITY_INTERPOLATION_POLICIES``.
    """
    if density_interpolation not in DENSITY_INTERPOLATION_POLICIES:
        raise ValueError(
            f"[input policies] density_interpolation={density_interpolation!r} is "
            f"unsupported; choose from {tuple(DENSITY_INTERPOLATION_POLICIES)}"
        )
    return DENSITY_INTERPOLATION_POLICIES[density_interpolation][1]


_INTEGRATION_MODES = ("central", "integrated")
_INTEGRATION_NUMERICAL_KEYS = tuple(INTEGRATED_DEFAULTS)
_INTEGRATION_FIELD_KEYS = ("min_zq_forest", "max_zq_forest")


def _pop_integration_keys(parser):
    """Remove the integrated forest-source keys from a parsed INI.

    Parameters
    ----------
    parser : configparser.ConfigParser
        Parsed INI; modified in place.

    Returns
    -------
    request : dict
        ``mode`` (raw ``[input policies] forest_source_integration`` string or
        None), ``numerical`` (raw ``[numerical] forest_*`` strings present) and
        ``fields`` (raw ``min_zq_forest``/``max_zq_forest`` strings present,
        keyed by field section name).

    Notes
    -----
    The keys are removed so that the central-mode option validation, the
    expanded policy mappings and all existing provenance are exactly those of an
    INI that never mentioned the integrated mode.
    """
    request = {"mode": None, "numerical": {}, "fields": {}}
    if parser.has_option("input policies", "forest_source_integration"):
        request["mode"] = parser.get("input policies", "forest_source_integration")
        parser.remove_option("input policies", "forest_source_integration")
    for key in _INTEGRATION_NUMERICAL_KEYS:
        if parser.has_option("numerical", key):
            request["numerical"][key] = parser.get("numerical", key)
            parser.remove_option("numerical", key)
    for section in parser.sections():
        if not section.startswith("field "):
            continue
        for key in _INTEGRATION_FIELD_KEYS:
            if parser.has_option(section, key):
                request["fields"].setdefault(section, {})[key] = parser.get(
                    section, key
                )
                parser.remove_option(section, key)
    return request


def _resolve_integration_mode(request, explicit_revision, has_prescription):
    """Resolve the forest-source mode and prescription revision.

    Parameters
    ----------
    request : dict
        Output of ``_pop_integration_keys``.
    explicit_revision : str or None
        ``[prescription] revision`` as written, or None when absent.
    has_prescription : bool
        Whether the INI has a ``[prescription]`` section.

    Returns
    -------
    integrated : bool
        True for the integrated forest-source mode.
    revision : str or None
        Revision to use: ``INTEGRATED_REVISION`` in integrated mode, else the
        explicit revision (None when absent).

    Raises
    ------
    ValueError
        If the mode key is unknown, the key and the revision disagree
        (central with the integrated revision, integrated with any other
        explicit revision), integrated mode lacks a ``[prescription]``, or an
        integrated-only key appears in central mode.

    Notes
    -----
    An absent key selects the central mode unless the integrated revision is
    named; the integrated key without a revision implies the integrated
    revision.
    """
    mode = request["mode"]
    if mode is not None and mode not in _INTEGRATION_MODES:
        raise ValueError(
            f"[input policies] forest_source_integration={mode!r} is unsupported; "
            f"choose from {_INTEGRATION_MODES}"
        )
    if explicit_revision == INTEGRATED_REVISION:
        if mode == "central":
            raise ValueError(
                f"[prescription] revision={INTEGRATED_REVISION!r} conflicts with "
                "forest_source_integration=central"
            )
        integrated = True
    elif mode == "integrated":
        if explicit_revision is not None:
            raise ValueError(
                "forest_source_integration=integrated requires the revision "
                f"{INTEGRATED_REVISION!r} or none, not {explicit_revision!r}"
            )
        integrated = True
    else:
        integrated = False
    if integrated:
        if not has_prescription:
            raise ValueError(
                "the integrated forest-source mode requires a [prescription] section"
            )
        return True, INTEGRATED_REVISION
    if request["numerical"] or request["fields"]:
        keys = sorted(
            [f"[numerical] {key}" for key in request["numerical"]]
            + [
                f"[{section}] {key}"
                for section, values in request["fields"].items()
                for key in values
            ]
        )
        raise ValueError(
            "integrated-only option(s) require forest_source_integration=integrated: "
            + ", ".join(keys)
        )
    return False, explicit_revision


def _integration_settings(request, fields):
    """Parse the integrated-mode numerical orders and source-redshift limits.

    Parameters
    ----------
    request : dict
        Output of ``_pop_integration_keys``.
    fields : sequence of FieldConfig
        Parsed fields, in configured order.

    Returns
    -------
    settings : mappingproxy
        ``mode``, ``revision``, the three quadrature controls (defaults from
        ``INTEGRATED_DEFAULTS``) and ``fields``: requested ``min_zq_forest`` and
        ``max_zq_forest`` of every forest field (None selects the density-table
        default at preparation).

    Raises
    ------
    ValueError
        If an order is not a positive integer, a limit is not finite and
        nonnegative, the limits are not ordered, or a limit is set on a
        galaxy field.
    """
    orders = {}
    for key, default in INTEGRATED_DEFAULTS.items():
        text = request["numerical"].get(key)
        try:
            value = default if text is None else int(text)
        except ValueError as error:
            raise ValueError(f"[numerical] {key} must be a positive integer") from error
        if value < 1:
            raise ValueError(f"[numerical] {key} must be a positive integer")
        orders[key] = value

    by_section = request["fields"]
    forest_limits = {}
    for field in fields:
        section = f"field {field.observed.id}"
        values = by_section.get(section, {})
        if field.observed.kind != "forest":
            if values:
                raise ValueError(
                    f"[{section}] {', '.join(sorted(values))} apply to forest "
                    "fields only"
                )
            continue
        limits = {}
        for key in _INTEGRATION_FIELD_KEYS:
            if key not in values:
                limits[key] = None
                continue
            try:
                value = float(values[key])
            except ValueError as error:
                raise ValueError(f"[{section}] {key} must be numeric") from error
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"[{section}] {key} must be finite and nonnegative")
            limits[key] = value
        if (
            limits["min_zq_forest"] is not None
            and limits["max_zq_forest"] is not None
            and limits["max_zq_forest"] <= limits["min_zq_forest"]
        ):
            raise ValueError(f"[{section}] max_zq_forest must exceed min_zq_forest")
        forest_limits[field.observed.id] = _plain_mapping(limits)
    return _plain_mapping(
        {
            "mode": "integrated",
            "revision": INTEGRATED_REVISION,
            **orders,
            "fields": _plain_mapping(forest_limits),
        }
    )


def parse_survey_ini(source=None):
    """Parse and strictly validate a native FishHighz survey INI.

    Parameters
    ----------
    source : path-like or None, optional
        Native INI path; None selects the bundled DESI-2 accuracy recipe.

    Returns
    -------
    config : SurveyConfig
        Immutable validated configuration, preserving explicit spectrum
        selections.

    Raises
    ------
    UnsupportedSchemaError
        If the input is not a supported native schema.
    ValueError
        If options or physical/numerical domains are invalid.
    FileNotFoundError
        If the external INI is absent.

    Notes
    -----
    The parser expands the selected accuracy prescription and checks supported model, input and numerical conventions. The ``[prescription] revision`` (default ``early-lyaforecast-2026-10-01``; the historical ``early-lyaforecast-2026-09-18`` is also accepted) sets the default source-density interpolation: piecewise-constant cells for the native revision, the quadratic ``RectBivariateSpline_kx2_ky2_s0`` for the historical one. An explicit ``[input policies] density_interpolation`` overrides it, and the revision and effective policies are recorded in the provenance. It reads no survey tables and invokes no CAMB calculation.

    The opt-in integrated forest-source mode (``[input policies]
    forest_source_integration = integrated``, or ``[prescription] revision =
    early-lyaforecast-integrated-2026-10-02``) integrates every forest field
    over source redshift and forest pixel; its keys are parsed here, kept out of
    the expanded policy and numerical mappings, and recorded only in the
    ``forest_source_integration`` provenance block. Without them the parsed
    configuration and provenance are those of the central mode.
    """

    parser, path, source_identity = _read_parser(source)
    if "schema" not in parser:
        legacy = "cosmo" in parser or any(s.startswith("tracer") for s in parser)
        detail = (
            "lyaforecast INI schema is unsupported"
            if legacy
            else "missing [schema] section"
        )
        raise UnsupportedSchemaError(
            f"unsupported FishHighz INI schema: {detail}; expected "
            f"[{SCHEMA_NAME}] version {SCHEMA_VERSION}"
        )
    _section_options(parser, "schema", ("name", "version"), ("name", "version"))
    if parser["schema"]["name"] != SCHEMA_NAME:
        raise UnsupportedSchemaError(
            f"unsupported FishHighz INI schema name {parser['schema']['name']!r}; "
            f"expected {SCHEMA_NAME!r} version {SCHEMA_VERSION}"
        )
    try:
        version = parser["schema"].getint("version")
    except ValueError as error:
        raise UnsupportedSchemaError(
            "FishHighz schema version must be an integer"
        ) from error
    if version != SCHEMA_VERSION:
        raise UnsupportedSchemaError(
            f"unsupported FishHighz INI schema version {version}; supported version is {SCHEMA_VERSION}"
        )

    # Expand defaults only within the explicitly selected scientific model.
    mode = parser.get("model", "mode", fallback="bao")
    full_shape = mode == "full_shape"
    bao_marginalized = mode == "bao_marginalized"
    if full_shape:
        dilation_only = parser["model"].get("target_set", "full") == "dilation_only"
        full_basis = parser["model"].get("parameterization", "alpha_phi")
        amplitude = "alpha_iso" if full_basis == "alpha_iso_phi" else "alpha"
        for key, value in {
            "parameterization": full_basis,
            "parameter_names": (
                f"{amplitude}_w,phi_w,{amplitude}_s,phi_s"
                if dilation_only
                else f"{amplitude}_w,phi_w,{amplitude}_s,phi_s,f"
            ),
            "smooth_scaling": full_basis,
            "wiggle_scaling": full_basis,
            "growth_rate": "fixed_fiducial_f" if dilation_only else "free_f",
            "biases": "marginalized",
            "forest_bindings": "shared",
            "reported_growth": "f_sigma8_fid",
        }.items():
            if key not in parser["model"]:
                parser["model"][key] = value
    if bao_marginalized:
        for key, value in {
            "parameterization": "alpha_iso_phi",
            "parameter_names": "alpha_iso,phi",
            "smooth_scaling": "identity",
            "wiggle_scaling": "alpha_iso_phi",
            "growth_rate": "fixed_fiducial_f",
            "biases": "marginalized",
            "forest_bindings": "shared",
        }.items():
            if key not in parser["model"]:
                parser["model"][key] = value
    prescription = None
    if parser.defaults():
        raise ValueError("[DEFAULT] options are unsupported; use explicit sections")

    # Integrated forest-source keys are resolved first and removed, so the
    # remaining validation and the recorded mappings are those of the central
    # mode unless the integrated mode was requested.
    integration_request = _pop_integration_keys(parser)
    integrated, resolved_revision = _resolve_integration_mode(
        integration_request,
        parser["prescription"].get("revision") if "prescription" in parser else None,
        "prescription" in parser,
    )
    explicit_rtol = {
        section: parser.get(section, key)
        for section, key in (
            ("numerical", "weight_rtol"),
            ("input policies", "weighting_rtol"),
        )
        if parser.has_option(section, key)
    }
    if "prescription" in parser:
        _section_options(parser, "prescription", ("name", "revision"), ("name",))
        if parser["prescription"]["name"] != "accuracy":
            raise ValueError("[prescription] name must be accuracy")
        revision = (
            resolved_revision if resolved_revision is not None else NATIVE_REVISION
        )
        if revision not in REVISION_DENSITY_INTERPOLATION:
            raise ValueError(
                f"[prescription] unsupported revision {revision!r}; "
                f"supported revisions: {tuple(REVISION_DENSITY_INTERPOLATION)}"
            )
        prescription = _plain_mapping({"name": "accuracy", "revision": revision})
        for section, defaults in INI_DEFAULTS.items():
            if section not in parser:
                parser.add_section(section)
            for key, value in defaults.items():
                if key not in parser[section]:
                    parser[section][key] = value
        policies = parser["input policies"]

        # The revision fixes the default source-density treatment: the native
        # revision uses piecewise-constant cells, the historical revision the
        # quadratic spline. An explicit [input policies] entry overrides it and
        # is recorded with the other effective policies. Unknown values are
        # rejected here and mismatched interpolation/partition pairs below,
        # which also covers INIs without a [prescription] section.
        if "density_interpolation" not in policies:
            policies["density_interpolation"] = REVISION_DENSITY_INTERPOLATION[revision]
        partition = _density_partition(policies["density_interpolation"])
        if "magnitude_partition" not in policies:
            policies["magnitude_partition"] = partition
        if "num_z_bins" not in parser["survey"] and "z_edges" in parser["survey"]:
            parser["survey"]["num_z_bins"] = str(
                len(_float_list(parser["survey"]["z_edges"], "[survey] z_edges", 2)) - 1
            )
    if len(explicit_rtol) == 1:
        value = next(iter(explicit_rtol.values()))
        for section, key in (
            ("numerical", "weight_rtol"),
            ("input policies", "weighting_rtol"),
        ):
            if section in parser:
                parser[section][key] = value

    # Validate section structure before reading physical and numerical controls.
    required_sections = {
        "schema",
        "cosmology",
        "survey",
        "fields",
        "model",
        "input policies",
        "numerical",
    }
    missing_sections = required_sections - set(parser.sections())
    if missing_sections:
        raise UnsupportedSchemaError(
            "native FishHighz INI is missing required section(s): "
            + ", ".join(sorted(missing_sections))
        )
    unknown_sections = set(parser.sections()) - required_sections - {"prescription"}
    field_names = (
        _tokens(parser["fields"].get("ids", ""), "[fields] ids")
        if "fields" in parser
        else ()
    )
    expected_fields = {f"field {name}" for name in field_names}
    try:
        declared_bins = parser["survey"].getint("num_z_bins", fallback=0)
    except ValueError as error:
        raise ValueError("[survey] num_z_bins must be an integer") from error
    expected_bins = {f"pairs bin {index}" for index in range(1, declared_bins + 1)}
    unknown_sections -= expected_fields | expected_bins
    if unknown_sections:
        raise ValueError(
            f"unsupported native INI section(s): {', '.join(sorted(unknown_sections))}"
        )

    _section_options(
        parser,
        "cosmology",
        ("camb_ini", "template", "template_redshift", "damping_reference_redshift"),
        ("camb_ini", "template", "template_redshift", "damping_reference_redshift"),
    )
    _section_options(
        parser,
        "survey",
        (
            "area_deg2",
            "band",
            "min_band_mag",
            "max_band_mag",
            "z_edges",
            "num_z_bins",
            "resolution",
            "reconstruction_factor",
            "lya_rest_angstrom",
            "evaluation_redshift",
        ),
        (
            "area_deg2",
            "band",
            "min_band_mag",
            "max_band_mag",
            "z_edges",
            "num_z_bins",
            "resolution",
            "reconstruction_factor",
            "lya_rest_angstrom",
            "evaluation_redshift",
        ),
    )
    _section_options(parser, "fields", ("ids",), ("ids",))
    _section_options(
        parser,
        "model",
        (
            "mode",
            "biases",
            "forest_bindings",
            "reported_growth",
            "target_set",
            "growth_rate",
            "parameterization",
            "parameter_names",
            "smooth_scaling",
            "wiggle_scaling",
            "damping",
            "growth",
            "reconstruction",
            "forest_beta",
            "damping_amplitude",
            "forest_reconstruction_factor",
        ),
        (
            "parameterization",
            "parameter_names",
            "smooth_scaling",
            "wiggle_scaling",
            "damping",
            "growth",
            "reconstruction",
            "forest_beta",
            "damping_amplitude",
            "forest_reconstruction_factor",
        ),
    )
    _section_options(
        parser,
        "input policies",
        (
            "density_semantics",
            "density_width_policy",
            "density_redshift_normalization",
            "density_negative_policy",
            "density_interpolation",
            "snr_smoothing",
            "snr_interpolation",
            "snr_bright_policy",
            "snr_clamp",
            "snr_sentinel",
            "magnitude_partition",
            "weighting_method",
            "weighting_reference",
            "weighting_reference_k_t_deg",
            "weighting_reference_k_p_velocity",
            "weighting_rtol",
            "weighting_min_updates",
            "weighting_stable_steps",
            "weighting_max_updates",
            "sampling_noise",
        ),
        (
            "density_semantics",
            "density_width_policy",
            "density_redshift_normalization",
            "density_negative_policy",
            "density_interpolation",
            "snr_smoothing",
            "snr_interpolation",
            "snr_bright_policy",
            "snr_clamp",
            "snr_sentinel",
            "magnitude_partition",
            "weighting_method",
            "weighting_reference",
            "weighting_reference_k_t_deg",
            "weighting_reference_k_p_velocity",
            "weighting_rtol",
            "weighting_min_updates",
            "weighting_stable_steps",
            "weighting_max_updates",
            "sampling_noise",
        ),
    )
    _section_options(
        parser,
        "numerical",
        tuple(CONTROLS)
        + (
            "weight_rtol",
            "scale_step",
            "relative_step",
            "k_max_galaxy_galaxy",
            "k_max_galaxy_forest",
            "k_max_forest_forest",
        ),
        tuple(CONTROLS) + ("weight_rtol",),
    )

    if mode not in ("bao", "full_shape", "bao_marginalized"):
        raise ValueError("[model] mode must be bao, full_shape or bao_marginalized")
    if full_shape and parser["model"]["parameterization"] not in (
        "alpha_phi",
        "alpha_iso_phi",
    ):
        raise ValueError(
            "[model] full_shape parameterization must be alpha_phi or alpha_iso_phi"
        )
    full_basis = parser["model"]["parameterization"] if full_shape else None
    fixed_labels = {
        "parameterization": full_basis
        if full_shape
        else "alpha_iso_phi"
        if bao_marginalized
        else "ap_at",
        "smooth_scaling": full_basis if full_shape else "identity",
        "wiggle_scaling": full_basis
        if full_shape
        else "alpha_iso_phi"
        if bao_marginalized
        else "ap_at",
        "damping": "mixed_squared_width",
        "growth": "camb_sigma8_ratio",
        "reconstruction": "per_field",
    }
    if mode == "bao":
        inappropriate = set(parser["model"]) & {
            "biases",
            "forest_bindings",
            "reported_growth",
            "growth_rate",
            "target_set",
        }
        inappropriate |= set(parser["numerical"]) & {
            "scale_step",
            "relative_step",
            "k_max_galaxy_galaxy",
            "k_max_galaxy_forest",
            "k_max_forest_forest",
        }
        if inappropriate:
            raise ValueError(
                "full_shape mode required for options: "
                + ", ".join(sorted(inappropriate))
            )
    expected_names = (
        (
            "alpha_iso_w,phi_w,alpha_iso_s,phi_s"
            if full_basis == "alpha_iso_phi"
            else "alpha_w,phi_w,alpha_s,phi_s"
        )
        if full_shape and parser["model"].get("target_set") == "dilation_only"
        else (
            "alpha_iso_w,phi_w,alpha_iso_s,phi_s,f"
            if full_basis == "alpha_iso_phi"
            else "alpha_w,phi_w,alpha_s,phi_s,f"
        )
        if full_shape
        else "alpha_iso,phi"
        if bao_marginalized
        else "ap,at"
    )
    if parser["model"]["parameter_names"].replace(" ", "") != expected_names:
        raise ValueError(f"[model] parameter_names must be {expected_names}")
    if full_shape:
        target_set = parser["model"].get("target_set", "full")
        if target_set not in ("full", "dilation_only"):
            raise ValueError("[model] target_set must be full or dilation_only")
        parser["model"]["target_set"] = target_set
        fixed_labels.update(
            growth_rate=(
                "fixed_fiducial_f" if target_set == "dilation_only" else "free_f"
            ),
            biases="marginalized",
            forest_bindings="shared",
            reported_growth="f_sigma8_fid",
        )
    if bao_marginalized:
        if "target_set" in parser["model"] or "reported_growth" in parser["model"]:
            raise ValueError(
                "[model] target_set/reported_growth require full_shape mode"
            )
        fixed_labels.update(
            growth_rate="fixed_fiducial_f",
            biases="marginalized",
            forest_bindings="shared",
        )
    for key, expected in fixed_labels.items():
        if parser["model"][key] != expected:
            raise ValueError(f"[model] {key}={parser['model'][key]!r} is unsupported")
    for key, default in (("scale_step", "0.00025"), ("relative_step", "0.001")):
        value = parser["numerical"].get(key, default)
        if not np.isfinite(float(value)) or float(value) <= 0:
            raise ValueError(f"[numerical] {key} must be positive finite")
        if full_shape or bao_marginalized:
            parser["numerical"][key] = value
    fixed_policies = {
        "density_semantics": "cell_count_per_deg2",
        "density_width_policy": "legacy_first_spacing",
        "density_redshift_normalization": "target_density",
        "density_negative_policy": "floor_negative",
        "snr_smoothing": "legacy",
        "snr_interpolation": "linear_RegularGridInterpolator",
        "snr_bright_policy": "clamp_to_brightest_tabulated_magnitude",
        "weighting_method": "early_lyaforecast",
        "weighting_reference": "fiducial_auto_p3d_and_p1d_times_response_squared",
        "sampling_noise": "independent_diagonal",
    }
    for key, expected in fixed_policies.items():
        if parser["input policies"][key] != expected:
            raise ValueError(
                f"[input policies] {key}={parser['input policies'][key]!r} is unsupported"
            )
    density_interpolation = parser["input policies"]["density_interpolation"]
    partition = _density_partition(density_interpolation)
    if (
        integrated
        and DENSITY_INTERPOLATION_POLICIES[density_interpolation][0]
        != "piecewise_constant"
    ):
        raise ValueError(
            "the integrated forest-source mode requires piecewise-constant density "
            f"cells, not density_interpolation={density_interpolation!r}"
        )
    if parser["input policies"]["magnitude_partition"] != partition:
        raise ValueError(
            f"[input policies] magnitude_partition="
            f"{parser['input policies']['magnitude_partition']!r} is inconsistent "
            f"with density_interpolation={density_interpolation!r}; expected "
            f"{partition!r}"
        )
    try:
        model_beta = float(parser["model"]["forest_beta"])
        damping_amplitude = float(parser["model"]["damping_amplitude"])
        forest_reconstruction = float(parser["model"]["forest_reconstruction_factor"])
        reference_k_t = float(parser["input policies"]["weighting_reference_k_t_deg"])
        reference_k_p = float(
            parser["input policies"]["weighting_reference_k_p_velocity"]
        )
        snr_clamp = float(parser["input policies"]["snr_clamp"])
        snr_sentinel = float(parser["input policies"]["snr_sentinel"])
        weighting_rtol = float(parser["input policies"]["weighting_rtol"])
    except ValueError as error:
        raise ValueError(
            "native INI contains a malformed numeric policy/model value"
        ) from error
    for name, value in (
        ("forest_beta", model_beta),
        ("damping_amplitude", damping_amplitude),
        ("forest_reconstruction_factor", forest_reconstruction),
        ("weighting_reference_k_t_deg", reference_k_t),
        ("weighting_reference_k_p_velocity", reference_k_p),
        ("snr_clamp", snr_clamp),
        ("snr_sentinel", snr_sentinel),
        ("weighting_rtol", weighting_rtol),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"[model/input policies] {name} must be positive finite")
    if snr_clamp != 1e-10 or snr_sentinel != 1e20:
        raise ValueError(
            "[input policies] snr_clamp and snr_sentinel must match the adopted "
            "legacy SNR policy"
        )
    for key in ("template_redshift", "damping_reference_redshift"):
        try:
            redshift = float(parser["cosmology"][key])
        except ValueError as error:
            raise ValueError(f"[cosmology] {key} must be numeric") from error
        if not np.isfinite(redshift) or redshift < 0:
            raise ValueError(f"[cosmology] {key} must be finite and nonnegative")
    if float(parser["numerical"]["weight_rtol"]) != weighting_rtol:
        raise ValueError(
            "[numerical] weight_rtol must match [input policies] weighting_rtol"
        )

    # Survey limits and quadrature controls are validated in their native units.
    edges = _float_list(parser["survey"]["z_edges"], "[survey] z_edges", 2)
    try:
        area = float(parser["survey"]["area_deg2"])
        min_mag = float(parser["survey"]["min_band_mag"])
        max_mag = float(parser["survey"]["max_band_mag"])
        n_bins = parser["survey"].getint("num_z_bins")
        resolution = float(parser["survey"]["resolution"])
        reconstruction_factor = float(parser["survey"]["reconstruction_factor"])
        lya_rest = float(parser["survey"]["lya_rest_angstrom"])
        numerical = {
            key: float(parser["numerical"][key])
            for key in ("k_min", "k_max", "ap_step", "weight_rtol")
        }
        if full_shape or bao_marginalized:
            for key in (
                "k_max_galaxy_galaxy",
                "k_max_galaxy_forest",
                "k_max_forest_forest",
            ):
                if key in parser["numerical"]:
                    numerical[key] = float(parser["numerical"][key])
        integer_controls = {
            key: parser["numerical"].getint(key)
            for key in (
                "k_intervals",
                "k_order",
                "mu_order",
                "z_order",
                "magnitude_order",
            )
        }
        stopping = {
            key: parser["input policies"].getint(key)
            for key in (
                "weighting_min_updates",
                "weighting_stable_steps",
                "weighting_max_updates",
            )
        }
    except ValueError as error:
        raise ValueError(
            "native INI contains a malformed survey or numerical value"
        ) from error
    if not np.isfinite(area) or area <= 0 or area > 41252.96124941927:
        raise ValueError(
            "[survey] area_deg2 must be positive and no larger than the full sky"
        )
    if not np.isfinite(min_mag) or not np.isfinite(max_mag) or min_mag >= max_mag:
        raise ValueError(
            "[survey] magnitude limits must be finite and strictly ordered"
        )
    if n_bins < 1 or not np.isfinite(resolution) or resolution <= 0:
        raise ValueError(
            "[survey] num_z_bins must be positive and resolution must be positive finite"
        )
    if not np.isfinite(reconstruction_factor) or reconstruction_factor <= 0:
        raise ValueError("[survey] reconstruction_factor must be positive finite")
    if not np.isfinite(lya_rest) or lya_rest <= 0:
        raise ValueError("[survey] lya_rest_angstrom must be positive finite")
    if not np.all(np.isfinite(tuple(numerical.values()))) or (
        numerical["k_min"] <= 0
        or numerical["k_max"] <= numerical["k_min"]
        or numerical["ap_step"] <= 0
        or numerical["weight_rtol"] <= 0
    ):
        raise ValueError("[numerical] k/ap/weight controls have invalid domains")
    if any(
        value <= numerical["k_min"]
        for key, value in numerical.items()
        if key.startswith("k_max_")
    ):
        raise ValueError("[numerical] category k_max must exceed k_min")
    if any(value < 1 for value in integer_controls.values()):
        raise ValueError(
            "[numerical] quadrature and interval controls must be positive integers"
        )
    if (
        stopping["weighting_min_updates"] < 1
        or stopping["weighting_stable_steps"] < 1
        or stopping["weighting_max_updates"] < stopping["weighting_min_updates"]
    ):
        raise ValueError(
            "[input policies] weighting stopping controls have invalid domains"
        )
    if len(edges) != n_bins + 1 or any(b <= a for a, b in zip(edges[:-1], edges[1:])):
        raise ValueError(
            "[survey] z_edges must contain num_z_bins+1 strictly increasing edges"
        )
    if parser["survey"]["evaluation_redshift"] != "geometric_1plusz":
        raise ValueError("[survey] evaluation_redshift must be geometric_1plusz")
    if parser["survey"]["band"] != "r":
        raise ValueError("[survey] band must be the adopted r-band input")
    if len(set(field_names)) != len(field_names):
        raise ValueError("[fields] ids must be unique")

    # Field declarations distinguish source populations from physical tracers.
    allowed_field = {
        "kind",
        "physical_model",
        "background",
        "tracer",
        "density",
        "target_density",
        "z_norm_min",
        "density_magnitude_bounds",
        "min_band_mag",
        "max_band_mag",
        "snr",
        "num_exposures",
        "pix_width_angstrom",
        "min_rest_frame_lya",
        "max_rest_frame_lya",
        "bias_z",
        "bias_values",
    }
    fields = []
    for name in field_names:
        section = f"field {name}"
        if section not in parser:
            raise ValueError(f"missing [{section}] section")
        _section_options(
            parser,
            section,
            allowed_field,
            (
                "kind",
                "physical_model",
                "tracer",
                "density",
                "target_density",
                "density_magnitude_bounds",
            ),
        )
        values = parser[section]
        kind = values["kind"]
        if kind not in ("forest", "galaxy"):
            raise ValueError(f"[{section}] kind must be forest or galaxy")
        background = values.get("background")
        if kind == "forest":
            required = (
                "background",
                "snr",
                "num_exposures",
                "pix_width_angstrom",
                "min_rest_frame_lya",
                "max_rest_frame_lya",
            )
            if any(key not in values for key in required):
                raise ValueError(
                    f"[{section}] forest fields require {', '.join(required)}"
                )
        elif background is not None:
            raise ValueError(f"[{section}] galaxy fields cannot define background")
        target_density = values.getfloat("target_density")
        if not np.isfinite(target_density) or target_density <= 0:
            raise ValueError(f"[{section}] target_density must be positive finite")
        z_norm_min = values.getfloat("z_norm_min") if "z_norm_min" in values else None
        if z_norm_min is not None and (not np.isfinite(z_norm_min) or z_norm_min < 0):
            raise ValueError(f"[{section}] z_norm_min must be finite and nonnegative")
        density_bounds_policy = values["density_magnitude_bounds"]
        if density_bounds_policy not in ("survey", "none"):
            raise ValueError(
                f"[{section}] density_magnitude_bounds must be survey or none"
            )
        field_min_mag = values.getfloat("min_band_mag", fallback=min_mag)
        field_max_mag = values.getfloat("max_band_mag", fallback=max_mag)
        if field_min_mag != min_mag or field_max_mag != max_mag:
            raise ValueError(
                f"[{section}] unsupported per-field magnitude selection; "
                "min_band_mag and max_band_mag must match [survey] limits or be omitted"
            )
        bounds = (min_mag, max_mag) if density_bounds_policy == "survey" else None
        if kind == "forest":
            forest_values = {
                key: values.getfloat(key)
                for key in (
                    "num_exposures",
                    "pix_width_angstrom",
                    "min_rest_frame_lya",
                    "max_rest_frame_lya",
                )
            }
            if (
                not np.all(np.isfinite(tuple(forest_values.values())))
                or any(value <= 0 for value in forest_values.values())
                or forest_values["min_rest_frame_lya"]
                >= forest_values["max_rest_frame_lya"]
            ):
                raise ValueError(
                    f"[{section}] forest exposure, pixel and wavelength controls have invalid domains"
                )
        bias_z = (
            _float_list(values["bias_z"], f"[{section}] bias_z", 2)
            if "bias_z" in values
            else None
        )
        bias_values = (
            _float_list(values["bias_values"], f"[{section}] bias_values", 2)
            if "bias_values" in values
            else None
        )
        if (bias_z is None) != (bias_values is None) or (
            bias_z is not None and len(bias_z) != len(bias_values)
        ):
            raise ValueError(
                f"[{section}] bias_z and bias_values must be paired with equal length"
            )
        if bias_z is not None and (
            bias_z[0] < 0 or any(b <= a for a, b in zip(bias_z[:-1], bias_z[1:]))
        ):
            raise ValueError(f"[{section}] bias_z must be nonnegative and increasing")
        observed = ObservedField(
            name, kind, values["physical_model"], background=background
        )
        fields.append(
            FieldConfig(
                observed,
                values["tracer"],
                values["density"],
                target_density,
                z_norm_min,
                bounds,
                density_bounds_policy,
                values.get("snr"),
                values.getfloat("num_exposures") if "num_exposures" in values else None,
                values.getfloat("pix_width_angstrom")
                if "pix_width_angstrom" in values
                else None,
                values.getfloat("min_rest_frame_lya")
                if "min_rest_frame_lya" in values
                else None,
                values.getfloat("max_rest_frame_lya")
                if "max_rest_frame_lya" in values
                else None,
                bias_z,
                bias_values,
            )
        )

    forest_integration = (
        _integration_settings(integration_request, fields) if integrated else None
    )

    # Preserve each bin selection in canonical observed-field order.
    bins = []
    ids = {field.observed.id for field in fields}
    field_order = [field.observed.id for field in fields]
    for index, (z_min, z_max) in enumerate(zip(edges[:-1], edges[1:]), 1):
        section = f"pairs bin {index}"
        if section not in parser:
            raise ValueError(f"missing [{section}] section")
        _section_options(parser, section, ("selected",), ("selected",))
        pairs = []
        selected = parser[section]["selected"]
        if selected.strip() == "all":
            selected = ",".join(
                f"{left}x{right}"
                for position, left in enumerate(field_order)
                for right in field_order[position:]
            )
        if not selected.strip() and not (full_shape or bao_marginalized):
            raise ValueError(
                f"[{section}] empty selected pairs require full_shape or bao_marginalized mode"
            )
        for token in () if not selected.strip() else selected.split(","):
            pair = tuple(x.strip() for x in token.split("x"))
            if len(pair) != 2 or any(x not in ids for x in pair):
                raise ValueError(f"[{section}] invalid selected pair {token!r}")
            pair = tuple(sorted(pair, key=field_order.index))
            if pair in pairs:
                raise ValueError(f"[{section}] duplicate selected pair {pair}")
            pairs.append(pair)
        z_eval = float(np.sqrt((1 + z_min) * (1 + z_max)) - 1)
        bins.append(BinConfig(index, z_min, z_max, z_eval, tuple(pairs)))
    if (full_shape or bao_marginalized) and not any(
        item.selected_pairs for item in bins
    ):
        raise ValueError(f"{mode} requires at least one nonempty redshift bin")
    if full_shape and parser["model"]["target_set"] == "dilation_only":
        if any(
            any(
                field.observed.kind != "forest"
                for field in fields
                if field.observed.id in pair
            )
            for item in bins
            for pair in item.selected_pairs
        ):
            raise ValueError(
                "dilation_only target_set requires forest-only selected spectra"
            )

    return SurveyConfig(
        path,
        SCHEMA_NAME,
        version,
        _plain_mapping(dict(parser["cosmology"])),
        _plain_mapping(dict(parser["survey"])),
        _plain_mapping(dict(parser["model"])),
        _plain_mapping(dict(parser["input policies"])),
        _plain_mapping(dict(parser["numerical"])),
        tuple(fields),
        tuple(bins),
        prescription,
        source_identity,
        forest_integration,
    )


def parse_ini(source=None):
    """Parse and strictly validate a native FishHighz survey INI.

    Parameters
    ----------
    source : path-like or None, optional
        Native INI path; None selects the bundled DESI-2 accuracy recipe.

    Returns
    -------
    config : SurveyConfig
        Immutable validated configuration, preserving explicit spectrum
        selections.

    Raises
    ------
    UnsupportedSchemaError
        If the input is not a supported native schema.
    ValueError
        If options or physical/numerical domains are invalid.
    FileNotFoundError
        If the external INI is absent.

    Notes
    -----
    The parser expands the selected accuracy prescription and checks supported model, input and numerical conventions. The ``[prescription] revision`` (default ``early-lyaforecast-2026-10-01``; the historical ``early-lyaforecast-2026-09-18`` is also accepted) sets the default source-density interpolation: piecewise-constant cells for the native revision, the quadratic ``RectBivariateSpline_kx2_ky2_s0`` for the historical one. An explicit ``[input policies] density_interpolation`` overrides it, and the revision and effective policies are recorded in the provenance. It reads no survey tables and invokes no CAMB calculation.
    """

    return parse_survey_ini(source)


def _resolve_resource(value, config):
    """Classify package resources and resolve external input paths.

    Parameters
    ----------
    value : str or path-like
        Configured input reference.
    config : SurveyConfig
        Parsed native survey configuration.

    Returns
    -------
    kind : str
        package or external resource category.
    path : str or pathlib.Path
        Package-relative name or resolved external path.
    """
    value = str(value)
    if value.startswith("package:"):
        return "package", value[len("package:") :]
    if config.path is None:
        return "package", value
    return "external", resolve_input_path(value, config.path)


def _reader_samples(field, density, snr, z_source, magnitudes, wavelength):
    """Sample density and forest noise at the explicit source coordinates.

    Parameters
    ----------
    field : FieldConfig
        Observed population and instrument settings.
    density : object
        Density adapter exposing sample or query.
    snr : object or None
        Forest SNR adapter exposing sample or variance; unused for galaxies.
    z_source : float
        Dimensionless source redshift.
    magnitudes : array_like
        Magnitude nodes, shape (n_magnitude,).
    wavelength : float
        Observed Ly-alpha wavelength in angstrom.

    Returns
    -------
    sample : dict
        Density in deg^-2 redshift^-1 mag^-1 and, for forests, dimensionless
        pixel variance; arrays have shape (n_magnitude,). Reader provenance is
        retained.

    Raises
    ------
    ValueError
        If a required reader method is unavailable or reader validation fails.
    """
    if hasattr(density, "sample"):
        density_sample = density.sample(z_source, magnitudes)
        density_values, density_provenance = (
            density_sample["values"],
            density_sample["provenance"],
        )
    elif hasattr(density, "query"):
        density_values = density.query(z_source, magnitudes)
        density_provenance = getattr(density, "provenance", {})
    else:
        raise ValueError(
            f"{field.observed.id}: density reader has no sample/query method"
        )
    result = {
        "density": np.asarray(density_values),
        "density_diagnostics": density_provenance,
    }
    if field.observed.kind == "forest":
        if hasattr(snr, "sample"):
            snr_sample = snr.sample(
                z_source=z_source,
                magnitudes=magnitudes,
                wavelength=wavelength,
                pixel_width_angstrom=field.pixel_width_angstrom,
                exposure_count=field.num_exposures,
            )
            result["variance"] = np.asarray(snr_sample["values"])
            result["snr_diagnostics"] = snr_sample["provenance"]
        elif hasattr(snr, "variance"):
            result["variance"] = np.asarray(
                snr.variance(
                    z_source=z_source,
                    magnitudes=magnitudes,
                    wavelength=wavelength,
                    pixel_width_angstrom=field.pixel_width_angstrom,
                    exposure_count=field.num_exposures,
                )
            )
            result["snr_diagnostics"] = getattr(snr, "provenance", {})
        else:
            raise ValueError(
                f"{field.observed.id}: SNR reader has no sample/variance method"
            )
    return result


def _grid_result(result):
    """Split a grid-query result into values and fallback counts.

    Parameters
    ----------
    result : mapping or array_like
        Legacy adapters return a mapping with ``values`` and ``provenance``;
        strict readers return the array itself.

    Returns
    -------
    values : ndarray
        Query values.
    provenance : dict
        Adapter provenance (empty for strict readers).
    """
    if hasattr(result, "items"):
        return np.asarray(result["values"]), dict(result["provenance"])
    return np.asarray(result), {}


def _integrated_zq_limits(item, density, limits):
    """Resolve the source-redshift limits of one integrated forest field.

    Parameters
    ----------
    item : FieldConfig
        Forest field.
    density : object
        Density adapter exposing the cell edges ``z_edges`` (directly or through
        ``reader``).
    limits : mapping
        Requested ``min_zq_forest`` and ``max_zq_forest`` (None for the default).

    Returns
    -------
    zq_min, zq_max : float
        Source-redshift limits. The defaults are ``z_norm_min`` (else the lowest
        density-table redshift edge) and the highest edge; every limit is
        clipped to the table, outside which the density is not tabulated.

    Raises
    ------
    ValueError
        If the density has no redshift cell edges (spline density) or the
        limits leave an empty interval.
    """
    edges = getattr(getattr(density, "reader", density), "z_edges", None)
    if edges is None:
        raise ValueError(
            f"{item.observed.id}: the integrated forest-source mode requires a "
            "piecewise-constant density exposing z_edges"
        )
    table_min, table_max = float(edges[0]), float(edges[-1])
    requested_min, requested_max = limits["min_zq_forest"], limits["max_zq_forest"]
    if requested_min is None:
        requested_min = table_min if item.z_norm_min is None else item.z_norm_min
    if requested_max is None:
        requested_max = table_max
    zq_min, zq_max = max(requested_min, table_min), min(requested_max, table_max)
    if not zq_max > zq_min:
        raise ValueError(
            f"{item.observed.id}: source-redshift limits [{requested_min}, "
            f"{requested_max}] do not overlap the density table [{table_min}, "
            f"{table_max}]"
        )
    return zq_min, zq_max


def _integrated_source(
    item,
    density,
    snr,
    response,
    bin_config,
    *,
    integration,
    zq_limits,
    magnitudes,
    quadrature,
    lya_rest,
):
    """Build the integrated forest source of one field in one redshift bin.

    Parameters
    ----------
    item : FieldConfig
        Forest field (rest-frame forest limits, exposures).
    density : object
        Density adapter with ``sample_grid`` (legacy) or ``query_grid`` (strict).
    snr : object
        SNR adapter with ``variance_grid``.
    response : InstrumentResponse
        Bin response; its pixel width is the fixed velocity pixel width l_pix.
    bin_config : BinConfig
        Redshift bin; its actual bounds set the observed-wavelength slice.
    integration : mapping
        Resolved ``SurveyConfig.forest_integration`` (quadrature orders).
    zq_limits : tuple of float
        Resolved (zq_min, zq_max) of the field.
    magnitudes, quadrature : ndarray of shape (n_magnitude,)
        Shared magnitude nodes and weights.
    lya_rest : float
        Lyman-alpha rest wavelength in angstrom.

    Returns
    -------
    source : IntegratedForestSource
        Nodes, density (n_y, n_m), variance (n_p, n_m) and provenance info.
    density_provenance, snr_provenance : dict
        Adapter provenance of the grid queries (empty for strict readers).

    Raises
    ------
    ValueError
        If the bin has no forest coverage for this field, a reader lacks the
        vectorised grid queries, or a query leaves its table domain under a
        strict reader.

    Notes
    -----
    The slice is [lya_rest (1+z_min), lya_rest (1+z_max)] of the bin bounds, not
    the evaluation redshift. Breakpoints in ln(1+z_q) are the density-table
    redshift cell edges and the S/N-table source redshifts. The pixel variance
    is queried at (z_q, lambda, pixel width l_pix lambda/c in angstrom) with
    the field's exposure count, the arguments of the central path except for
    the coordinates.
    """
    label = f"field {item.observed.id}, bin {bin_config.index}"
    density_axis = getattr(getattr(density, "reader", density), "z_edges", None)
    snr_axis = getattr(getattr(snr, "reader", snr), "z", None)
    breaks = np.concatenate(
        [
            np.asarray(axis, dtype=float).ravel()
            for axis in (density_axis, snr_axis)
            if axis is not None
        ]
    )
    nodes = integration_nodes(
        lambda_min=lya_rest * (1 + bin_config.z_min),
        lambda_max=lya_rest * (1 + bin_config.z_max),
        rest_min=item.min_rest_frame_lya,
        rest_max=item.max_rest_frame_lya,
        zq_min=zq_limits[0],
        zq_max=zq_limits[1],
        zq_breaks=breaks,
        zq_order=integration["forest_zq_order"],
        lambda_order=integration["forest_lambda_order"],
        lambda_panels=integration["forest_lambda_panels"],
        lya_rest_angstrom=lya_rest,
        label=label,
    )
    for reader, method in (
        (density, ("sample_grid", "query_grid")),
        (snr, ("variance_grid",)),
    ):
        if not any(hasattr(reader, name) for name in method):
            raise ValueError(
                f"{label}: the integrated forest-source mode needs a reader with "
                f"{' or '.join(method)}"
            )
    grid_density = (
        density.sample_grid if hasattr(density, "sample_grid") else density.query_grid
    )
    density_values, density_provenance = _grid_result(
        grid_density(nodes.zq_nodes, magnitudes)
    )

    # Fixed velocity pixel width of the bin, converted to angstrom at each pixel.
    variance_values, snr_provenance = _grid_result(
        snr.variance_grid(
            z_source=nodes.z_q,
            wavelength=nodes.lam_obs,
            magnitudes=magnitudes,
            pixel_width_angstrom=pixel_width_angstrom(
                nodes, response.pixel_width_velocity
            ),
            exposure_count=item.num_exposures,
        )
    )
    info = {
        "zq_limits": [float(zq_limits[0]), float(zq_limits[1])],
        "n_zq_nodes": int(len(nodes.zq_nodes)),
        "n_pixel": int(len(nodes.lam_obs)),
        "nodes": dict(nodes.info),
        "fallback_counts": {
            "density": dict(density_provenance.get("counts", {})),
            "snr": dict(snr_provenance.get("counts", {})),
        },
    }
    source = integrated_forest_source(
        nodes,
        density_values,
        magnitudes,
        quadrature,
        variance_values,
        response.pixel_width_velocity,
        label=label,
        info=info,
    )
    return source, density_provenance, snr_provenance


def _integration_bin_record(spec_id, integration_record):
    """Collect per-field integrated-source metadata of one bin.

    Parameters
    ----------
    spec_id : str
        Identifier of the bin specification (the prepared bin carries the same).
    integration_record : mapping
        Source ``info`` of each active integrated forest field.

    Returns
    -------
    record : dict
        ``spec_id`` and, per quantity, a mapping from field ID to value:
        ``n_zq_nodes``, ``n_pixel``, ``zq_limits``, ``zq_window`` and
        ``fallback_counts`` of the legacy grid queries. The effective redshift,
        ``N1``, ``A`` and ``P_pixel`` are added by
        ``PreparedForecast.provenance`` once the weights exist.
    """
    return {
        "spec_id": spec_id,
        "n_zq_nodes": {k: v["n_zq_nodes"] for k, v in integration_record.items()},
        "n_pixel": {k: v["n_pixel"] for k, v in integration_record.items()},
        "zq_limits": {k: v["zq_limits"] for k, v in integration_record.items()},
        "zq_window": {
            k: v["nodes"]["window_zq"] for k, v in integration_record.items()
        },
        "fallback_counts": {
            k: v["fallback_counts"] for k, v in integration_record.items()
        },
    }


def _background_values(background, z):
    """Read sigma8 and logarithmic growth rate at an exact redshift.

    Parameters
    ----------
    background : object
        Background with scalar growth methods or exact-redshift arrays.
    z : float
        Dimensionless evaluation redshift.

    Returns
    -------
    sigma8 : float
        Dimensionless rms density fluctuation amplitude.
    growth_rate : float
        Dimensionless logarithmic growth rate f.

    Raises
    ------
    ValueError
        If exact redshift growth values cannot be obtained.
    """
    for sigma_name, growth_name in (
        ("sigma8_at", "growth_rate_at"),
        ("sigma8", "growth_rate"),
    ):
        sigma = getattr(background, sigma_name, None)
        growth = getattr(background, growth_name, None)
        if callable(sigma) and callable(growth):
            return float(sigma(z)), float(growth(z))
    z_bins = getattr(background, "z_bins", None)
    sigma_bins = getattr(background, "sigma8_zbins", None)
    growth_bins = getattr(background, "growth_rate_zbins", None)
    if z_bins is not None and sigma_bins is not None and growth_bins is not None:
        indices = np.flatnonzero(np.asarray(z_bins) == z)
        if len(indices) == 1:
            i = indices[0]
            return float(sigma_bins[i]), float(growth_bins[i])
    raise ValueError(f"background must provide exact sigma8/growth at z={z}")


def _require_scalar_match(actual, expected, name):
    """Require agreement within a float64 roundoff allowance.

    Parameters
    ----------
    actual : float
        Prepared value.
    expected : float
        Configured reference, in the same units.
    name : str
        Quantity label for errors.

    Returns
    -------
    None
        No value is returned.

    Raises
    ------
    ValueError
        If inputs are not finite or differ beyond the allowance.

    Notes
    -----
    The absolute allowance is 64*eps64 times the larger of unity and both magnitudes.
    """
    try:
        actual, expected = float(actual), float(expected)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite numeric values") from error
    if not np.isfinite(actual) or not np.isfinite(expected):
        raise ValueError(f"{name} must be finite")
    tolerance = 64 * np.finfo(float).eps * max(1.0, abs(actual), abs(expected))
    if abs(actual - expected) > tolerance:
        raise ValueError(f"{name} mismatch: prepared={actual}, configured={expected}")


def _require_finite_scalar(value, name, *, positive=False):
    """Validate a finite scalar configuration value.

    Parameters
    ----------
    value : float
        Candidate value; units depend on the named quantity.
    name : str
        Quantity label for errors.
    positive : bool, optional
        Also require strict positivity; default False.

    Returns
    -------
    result : float
        Validated finite value.

    Raises
    ------
    ValueError
        If conversion, finiteness or requested positivity fails.
    """
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite numeric") from error
    if not np.isfinite(result) or (positive and result <= 0):
        qualifier = "positive finite" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _build_readers(config, fields, stack):
    """Build density and SNR adapters from materialized resources.

    Parameters
    ----------
    config : SurveyConfig
        Parsed native survey configuration.
    fields : sequence of FieldConfig
        Ordered field configurations.
    stack : contextlib.ExitStack
        Context retaining materialized package resources.

    Returns
    -------
    readers : dict
        Density/SNR adapters keyed by observed-field ID.

    Raises
    ------
    ValueError
        If reader data or the density policy is unsupported.

    Notes
    -----
    Reads configured survey tables and constructs interpolation objects with the adopted explicit compatibility policies.
    """

    from .adapters.legacy_compat import LegacyDensity, LegacySNR
    from .adapters.legacy_inputs import DensityReader, SNRReader

    result = {}
    policy = config.input_policies
    for field in fields:
        kind, density_path = _resolve_resource(field.density, config)
        if kind == "package":
            density_path = stack.enter_context(bundled_path(density_path))
        density_reader = DensityReader(
            density_path,
            semantics=policy["density_semantics"],
            target_density=field.target_density,
            z_norm_min=field.z_norm_min,
            magnitude_bounds=field.magnitude_bounds,
            width_policy=policy["density_width_policy"],
            label=field.observed.id,
            interpolation=DENSITY_INTERPOLATION_POLICIES[
                policy["density_interpolation"]
            ][0],
        )
        if policy["density_negative_policy"] == "floor_negative":
            density = LegacyDensity(density_reader, "floor_negative")
        elif policy["density_negative_policy"] == "reject":
            density = LegacyDensity(density_reader, "reject")
        else:
            raise ValueError("unsupported density_negative_policy")
        snr = None
        if field.observed.kind == "forest":
            kind, snr_path = _resolve_resource(field.snr, config)
            paths = sorted(Path(snr_path).glob("*.dat")) if kind == "external" else None
            if kind == "package":
                paths = stack.enter_context(bundled_paths(snr_path))
            snr = LegacySNR(
                SNRReader(
                    paths, smoothing=policy["snr_smoothing"], label=field.observed.id
                )
            )
        result[field.observed.id] = {"density": density, "snr": snr}
    return result


def _normalise_readers(config, fields, readers, stack):
    """Validate supplied reader coverage or construct configured readers.

    Parameters
    ----------
    config : SurveyConfig
        Parsed native survey configuration.
    fields : sequence of FieldConfig
        Fields that must be covered exactly.
    readers : mapping or None
        Injected readers, or None to load configured resources.
    stack : contextlib.ExitStack
        Context retaining materialized package resources.

    Returns
    -------
    normalized : dict
        Per-field dictionaries containing density and snr adapters.

    Raises
    ------
    ValueError
        If injected readers do not cover exactly the configured fields.
    """
    if readers is None:
        return _build_readers(config, fields, stack)
    if set(readers) != {field.observed.id for field in fields}:
        raise ValueError("readers must cover exactly all configured field IDs")
    result = {}
    for field in fields:
        value = readers[field.observed.id]
        if hasattr(value, "sample") or hasattr(value, "query"):
            density, snr = value, None
        else:
            if not hasattr(value, "keys") or "density" not in value:
                raise ValueError(f"{field.observed.id}: reader requires density")
            density, snr = value["density"], value.get("snr")
        result[field.observed.id] = {"density": density, "snr": snr}
    return result


def prepare_survey(config, *, background, template, readers=None, prepare=False):
    """Assemble native bin geometry, spectra and noise inputs.

    Parameters
    ----------
    config : SurveyConfig
        Parsed native survey configuration.
    background : object
        Prepared background providing exact growth, geometry and normalization
        metadata.
    template : object
        Prepared smooth/wiggle template with z_ref and h_fid metadata.
    readers : mapping or None, optional
        Density/SNR adapters keyed by field ID; None loads the configured
        resources.
    prepare : bool, optional
        Reserved preparation flag; only the default False is accepted.

    Returns
    -------
    survey : PreparedSurvey
        Field registry, nonempty bin specifications and immutable provenance.

    Raises
    ------
    ValueError
        If configuration or injected inputs disagree, a numerical domain is
        invalid, or prepare=True is requested.

    Notes
    -----
    This function never invokes CAMB. It validates prepared background/template consistency, integrates bin geometry, constructs quadrature, and samples density/SNR inputs. Fixed covariance preparation is performed by Forecast.prepare.
    """

    if not isinstance(config, SurveyConfig):
        raise ValueError("require parsed SurveyConfig")
    if prepare:
        raise ValueError("forecast preparation belongs to the Step-4 Forecast facade")
    if background is None or template is None:
        raise ValueError(
            "prepare_survey requires injected background and template; it never runs CAMB"
        )
    h_fid = getattr(background, "h_fid", None)
    if h_fid is None:
        raise ValueError("background must provide h_fid")
    required_background = (
        "template_growth_redshift",
        "damping_reference_redshift",
        "sigma8_damping_reference",
    )
    missing_background = [
        attribute
        for attribute in required_background
        if not hasattr(background, attribute)
    ]
    if missing_background:
        raise ValueError(
            "background must provide finite template-growth and damping-reference "
            "metadata: " + ", ".join(missing_background)
        )
    template_growth_redshift = _require_finite_scalar(
        background.template_growth_redshift,
        "background template_growth_redshift",
    )
    damping_reference_redshift = _require_finite_scalar(
        background.damping_reference_redshift,
        "background damping_reference_redshift",
    )
    sigma8_damping_reference = _require_finite_scalar(
        background.sigma8_damping_reference,
        "background sigma8_damping_reference",
        positive=True,
    )
    for attribute in ("z_ref", "h_fid"):
        if not hasattr(template, attribute):
            raise ValueError(
                f"template must provide {attribute} for INI consistency validation"
            )
    _require_scalar_match(
        template.z_ref,
        config.cosmology["template_redshift"],
        "template z_ref/configured template_redshift",
    )
    _require_scalar_match(
        template.h_fid,
        h_fid,
        "template h_fid/background h_fid",
    )
    _require_scalar_match(
        template_growth_redshift,
        config.cosmology["template_redshift"],
        "background template_growth_redshift/configured template_redshift",
    )
    _require_scalar_match(
        damping_reference_redshift,
        config.cosmology["damping_reference_redshift"],
        "background damping_reference_redshift/configured damping_reference_redshift",
    )
    # Bind only the configured target/nuisance parameterization.
    fields = config.observed_fields
    mode = config.model.get("mode", "bao")
    full_shape = mode == "full_shape"
    bao_marginalized = mode == "bao_marginalized"
    registry = ParameterRegistry(
        [
            Parameter(f"{name}_{index}", 1.0, "target", step=0.001)
            for index in range(len(config.bins))
            for name in ("ap", "at")
        ]
    )
    if full_shape:
        from .full_shape import make_registry

        registry = make_registry(config, background)
    if bao_marginalized:
        from .bao_marginalized import make_registry

        registry = make_registry(config, background)
    with ExitStack() as stack:
        reader_map = _normalise_readers(config, config.fields, readers, stack)
        densities = {name: item["density"] for name, item in reader_map.items()}
        snrs = {
            name: item["snr"]
            for name, item in reader_map.items()
            if item["snr"] is not None
        }
        specs = []
        bin_provenance = []
        integrated = config.forest_integration is not None
        zq_limits = {}
        for bin_config in config.bins:
            if (full_shape or bao_marginalized) and not bin_config.selected_pairs:
                bin_provenance.append(
                    {
                        "index": bin_config.index,
                        "bounds": [bin_config.z_min, bin_config.z_max],
                        "z_eval": bin_config.z_eval,
                        "selected_pairs": [],
                        "required_pairs": [],
                        "status": "excluded",
                    }
                )
                continue
            z_eval = bin_config.z_eval
            geometry = prepare_geometry(
                bin_config.z_min,
                bin_config.z_max,
                z_eval=z_eval,
                area_deg2=float(config.survey["area_deg2"]),
                h_fid=float(h_fid),
                z_order=int(config.numerical["z_order"]),
                hubble=background.hubble_parameter,
                transverse_distance=background.transverse_comoving_distance,
            )

            # Include every category boundary explicitly in the observed k grid.
            k_min = float(config.numerical["k_min"])
            field_kinds = {field.id: field.kind for field in fields}
            category_cuts = []
            for left, right in bin_config.selected_pairs:
                kinds = (field_kinds[left], field_kinds[right])
                suffix = (
                    "forest_forest"
                    if kinds == ("forest", "forest")
                    else "galaxy_galaxy"
                    if kinds == ("galaxy", "galaxy")
                    else "galaxy_forest"
                )
                category_cuts.append(
                    float(
                        config.numerical.get(
                            "k_max_" + suffix, config.numerical["k_max"]
                        )
                    )
                )
            k_edges = np.linspace(
                k_min, max(category_cuts), int(config.numerical["k_intervals"]) + 1
            )
            k_edges = np.unique(np.concatenate((k_edges, category_cuts)))
            grid = gauss_legendre_grid(
                k_edges,
                k_order=int(config.numerical["k_order"]),
                mu_order=int(config.numerical["mu_order"]),
                h_fid=float(h_fid),
            )

            # Growth and damping normalizations retain distinct reference redshifts.
            selection = PairSelection(fields, bin_config.selected_pairs)
            sigma8, growth_rate = _background_values(background, z_eval)
            sigma8_template, _ = _background_values(background, float(template.z_ref))
            biases, betas, widths = {}, {}, {}
            for field_config in config.fields:
                field = field_config.observed
                if field_config.bias_redshifts is None:
                    from .models.biases import analytic_density_bias

                    bias = analytic_density_bias(z_eval, field_config.tracer)
                else:
                    bias_function = linear_tabulated_bias(
                        field_config.bias_redshifts, field_config.bias_values
                    )
                    bias = bias_function(z_eval)
                biases[field.id] = float(bias)
                if field.kind == "forest":
                    betas[field.id] = float(config.model["forest_beta"])
                    reconstruction = float(config.model["forest_reconstruction_factor"])
                else:
                    reconstruction = float(config.survey["reconstruction_factor"])
                sigma_transverse = (
                    float(config.model["damping_amplitude"])
                    * sigma8
                    / sigma8_damping_reference
                    / np.sqrt(reconstruction)
                )
                widths[field.id] = (
                    (1 + growth_rate) * sigma_transverse,
                    sigma_transverse,
                )
            if full_shape:
                from .full_shape import make_model

                model, binding = make_model(
                    config,
                    bin_config,
                    registry,
                    template,
                    fields,
                    biases,
                    betas,
                    widths,
                    growth_rate,
                    (sigma8 / sigma8_template) ** 2,
                )
            elif bao_marginalized:
                from .bao_marginalized import make_model

                model, binding = make_model(
                    config,
                    bin_config,
                    registry,
                    template,
                    fields,
                    biases,
                    betas,
                    widths,
                    growth_rate,
                    (sigma8 / sigma8_template) ** 2,
                )
            else:
                model = KaiserModel(
                    template,
                    fields,
                    biases=biases,
                    betas=betas,
                    widths=widths,
                    f=growth_rate,
                    local_names=("ap", "at"),
                    wiggle=Scaling("ap_at", ap="ap", at="at"),
                    z=z_eval,
                    growth=(sigma8 / sigma8_template) ** 2,
                )
                binding = BoundParameters(
                    registry,
                    ("ap", "at"),
                    {name: f"{name}_{bin_config.index - 1}" for name in ("ap", "at")},
                )
            p3d = PreparedP3D(
                registry,
                selection,
                [
                    P3DProvider(
                        "native accuracy BAO", model, binding, selection.required_pairs
                    )
                ],
            )

            # Instrument response uses the foreground absorption wavelength.
            wavelength = float(config.survey["lya_rest_angstrom"]) * (1 + z_eval)
            responses = {
                item.observed.id: (
                    InstrumentResponse(
                        pixel_width_angstrom_to_velocity(
                            item.pixel_width_angstrom, lambda_obs_angstrom=wavelength
                        ),
                        resolving_power_fwhm_to_sigma(
                            float(config.survey["resolution"])
                        ),
                    )
                    if item.observed.kind == "forest"
                    else InstrumentResponse(0, 0)
                )
                for item in config.fields
            }

            # Central mode: forest density and SNR use the representative
            # background-source redshift. Integrated mode keeps it only for the
            # (unchanged) magnitude partition and integrates the forest sources
            # inside the bin below.
            z_sources = {
                item.observed.id: (
                    wavelength
                    / np.sqrt(item.min_rest_frame_lya * item.max_rest_frame_lya)
                    - 1
                    if item.observed.kind == "forest"
                    else z_eval
                )
                for item in config.fields
            }
            magnitude_min = float(config.survey["min_band_mag"])
            magnitude_max = float(config.survey["max_band_mag"])
            partition = breakpoints(
                densities, snrs, z_sources, magnitude_min, magnitude_max
            )
            magnitudes, quadrature = composite(
                partition, int(config.numerical["magnitude_order"])
            )

            # Sample active populations once on the shared magnitude quadrature.
            active = set(np.unique(selection.selected_pairs).tolist())
            forests, galaxies = {}, {}
            integration_record = {}
            for field_index, item in enumerate(config.fields):
                if field_index not in active:
                    continue
                if item.observed.kind == "forest" and integrated:
                    field_id = item.observed.id
                    if field_id not in zq_limits:
                        zq_limits[field_id] = _integrated_zq_limits(
                            item,
                            densities[field_id],
                            config.forest_integration["fields"][field_id],
                        )
                    source, density_provenance, snr_provenance = _integrated_source(
                        item,
                        densities[field_id],
                        snrs.get(field_id),
                        responses[field_id],
                        bin_config,
                        integration=config.forest_integration,
                        zq_limits=zq_limits[field_id],
                        magnitudes=magnitudes,
                        quadrature=quadrature,
                        lya_rest=float(config.survey["lya_rest_angstrom"]),
                    )
                    options = dict(
                        method=config.input_policies["weighting_method"],
                        rtol=float(config.numerical["weight_rtol"]),
                        min_updates=int(config.input_policies["weighting_min_updates"]),
                        stable_steps=int(
                            config.input_policies["weighting_stable_steps"]
                        ),
                        max_updates=int(config.input_policies["weighting_max_updates"]),
                    )
                    forests[field_id] = ForestInput(
                        options,
                        default_p1d,
                        BoundParameters(registry, (), {}),
                        registry.fiducials,
                        auxiliary_coordinates=(
                            float(config.input_policies["weighting_reference_k_t_deg"]),
                            float(
                                config.input_policies[
                                    "weighting_reference_k_p_velocity"
                                ]
                            ),
                        ),
                        provenance={
                            "density_policy": density_provenance,
                            "snr_policy": snr_provenance,
                            "input_policies": dict(config.input_policies),
                            "weighting_status": "pending_bin_preparation",
                            "magnitude_bounds": [magnitude_min, magnitude_max],
                            "forest_source_integration": source.info,
                        },
                        integrated=source,
                    )
                    integration_record[field_id] = source.info
                    continue
                sample = _reader_samples(
                    item,
                    densities[item.observed.id],
                    snrs.get(item.observed.id),
                    z_sources[item.observed.id],
                    magnitudes,
                    wavelength,
                )
                if item.observed.kind == "forest":
                    options = dict(
                        z_source=z_sources[item.observed.id],
                        magnitudes=magnitudes,
                        quadrature=quadrature,
                        rho=density_per_velocity(
                            sample["density"], z_source=z_sources[item.observed.id]
                        ),
                        variance=sample["variance"],
                        length_velocity=SPEED_LIGHT_KMS
                        * np.log(item.max_rest_frame_lya / item.min_rest_frame_lya),
                        method=config.input_policies["weighting_method"],
                        rtol=float(config.numerical["weight_rtol"]),
                        min_updates=int(config.input_policies["weighting_min_updates"]),
                        stable_steps=int(
                            config.input_policies["weighting_stable_steps"]
                        ),
                        max_updates=int(config.input_policies["weighting_max_updates"]),
                    )
                    forests[item.observed.id] = ForestInput(
                        options,
                        default_p1d,
                        BoundParameters(registry, (), {}),
                        registry.fiducials,
                        auxiliary_coordinates=(
                            float(config.input_policies["weighting_reference_k_t_deg"]),
                            float(
                                config.input_policies[
                                    "weighting_reference_k_p_velocity"
                                ]
                            ),
                        ),
                        provenance={
                            "density_policy": dict(sample["density_diagnostics"]),
                            "snr_policy": dict(sample["snr_diagnostics"]),
                            "input_policies": dict(config.input_policies),
                            "weighting_status": "pending_bin_preparation",
                            "magnitude_bounds": [magnitude_min, magnitude_max],
                            "z_source": z_sources[item.observed.id],
                        },
                    )
                else:
                    galaxies[item.observed.id] = local_galaxy_density(
                        sample["density"], quadrature, geometry
                    )
            specs.append(
                BinSpec(
                    f"desi2_accuracy-{bin_config.index}",
                    geometry,
                    grid,
                    p3d,
                    responses,
                    forests=forests,
                    galaxies=galaxies,
                    independent_sampling=True,
                )
            )
            bin_provenance.append(
                {
                    "index": bin_config.index,
                    "bounds": [bin_config.z_min, bin_config.z_max],
                    "z_eval": z_eval,
                    "selected_pairs": [list(pair) for pair in selection.selected_pairs],
                    "required_pairs": [list(pair) for pair in selection.required_pairs],
                    "magnitude_partition": partition,
                    "magnitude_order": int(config.numerical["magnitude_order"]),
                    "field_ids": [field.id for field in fields],
                }
            )
            if integrated:
                bin_provenance[-1]["forest_source_integration"] = (
                    _integration_bin_record(specs[-1].id, integration_record)
                )
        provenance_record = {
            "config": config.provenance,
            "registry_ids": list(registry.ids),
            "fields": [field.id for field in fields],
            "bins": bin_provenance,
            "background": type(background).__name__,
            "template": type(template).__name__,
        }
        if integrated:
            provenance_record["forest_source_integration"] = {
                **dict(config.forest_integration),
                "resolved_zq_limits": {
                    field_id: list(limits) for field_id, limits in zq_limits.items()
                },
            }
        provenance = freeze(provenance_record)
    return PreparedSurvey(config, fields, registry, tuple(specs), provenance)


def prepare_ini(source=None, *, background, template, readers=None):
    """Parse an INI and assemble native survey bin specifications.

    Parameters
    ----------
    source : path-like or None, optional
        Native INI path; None selects the bundled accuracy recipe.
    background : object
        Prepared background providing exact growth, geometry and normalization
        metadata.
    template : object
        Prepared smooth/wiggle template with z_ref and h_fid metadata.
    readers : mapping or None, optional
        Density/SNR adapters keyed by field ID; None loads the configured
        resources.

    Returns
    -------
    survey : PreparedSurvey
        Parsed configuration and prepared model/noise bin specifications.

    Raises
    ------
    ValueError
        If the configuration or supplied preparation inputs are invalid.

    Notes
    -----
    Reads configured inputs when readers are omitted; background and template must already be prepared.
    """

    return prepare_survey(
        parse_survey_ini(source),
        background=background,
        template=template,
        readers=readers,
    )


__all__ = [
    "BinConfig",
    "FieldConfig",
    "PreparedSurvey",
    "SCHEMA_NAME",
    "SCHEMA_VERSION",
    "SurveyConfig",
    "UnsupportedSchemaError",
    "parse_ini",
    "parse_survey_ini",
    "prepare_ini",
    "prepare_survey",
]
