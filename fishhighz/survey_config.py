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

from .accuracy import CONTROLS, INI_DEFAULTS, REVISION
from .fields import ObservedField, PairSelection
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
    return MappingProxyType(dict(value))


def _float_list(value, name, minimum=1):
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
    values = tuple(x.strip() for x in value.replace(",", " ").split() if x.strip())
    if not values:
        raise ValueError(f"{name}: require a nonempty list")
    return values


def _section_options(parser, section, allowed, required=()):
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

    @property
    def observed_fields(self):
        return tuple(field.observed for field in self.fields)

    @property
    def provenance(self):
        return freeze(
            {
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
        )


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
        """Alias making the preparation boundary explicit."""

        return self.bins


def _read_parser(source):
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


def parse_survey_ini(source=None):
    """Parse and strictly validate a native FishHighz INI.

    The default source is the bundled DESI-2 accuracy recipe.  A lyaforecast
    INI is rejected as an unsupported schema; no option translation is made.
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

    prescription = None
    if parser.defaults():
        raise ValueError("[DEFAULT] options are unsupported; use explicit sections")
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
        if parser["prescription"].get("revision", REVISION) != REVISION:
            raise ValueError(
                f"[prescription] unsupported revision; expected {REVISION}"
            )
        prescription = _plain_mapping({"name": "accuracy", "revision": REVISION})
        for section, defaults in INI_DEFAULTS.items():
            if section not in parser:
                parser.add_section(section)
            for key, value in defaults.items():
                if key not in parser[section]:
                    parser[section][key] = value
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
        tuple(CONTROLS) + ("weight_rtol",),
        tuple(CONTROLS) + ("weight_rtol",),
    )

    if (
        parser["model"]["parameterization"] != "ap_at"
        or parser["model"]["parameter_names"].replace(" ", "") != "ap,at"
    ):
        raise ValueError(
            "[model] native preparation currently requires ap_at with parameter_names=ap,at"
        )
    fixed_labels = {
        "parameterization": "ap_at",
        "smooth_scaling": "identity",
        "wiggle_scaling": "ap_at",
        "damping": "mixed_squared_width",
        "growth": "camb_sigma8_ratio",
        "reconstruction": "per_field",
    }
    for key, expected in fixed_labels.items():
        if parser["model"][key] != expected:
            raise ValueError(f"[model] {key}={parser['model'][key]!r} is unsupported")
    fixed_policies = {
        "density_semantics": "cell_count_per_deg2",
        "density_width_policy": "legacy_first_spacing",
        "density_redshift_normalization": "target_density",
        "density_interpolation": "RectBivariateSpline_kx2_ky2_s0",
        "density_negative_policy": "floor_negative",
        "snr_smoothing": "legacy",
        "snr_interpolation": "linear_RegularGridInterpolator",
        "snr_bright_policy": "clamp_to_brightest_tabulated_magnitude",
        "magnitude_partition": "density_knots_support_snr_nodes_negative_roots",
        "weighting_method": "early_lyaforecast",
        "weighting_reference": "fiducial_auto_p3d_and_p1d_times_response_squared",
        "sampling_noise": "independent_diagonal",
    }
    for key, expected in fixed_policies.items():
        if parser["input policies"][key] != expected:
            raise ValueError(
                f"[input policies] {key}={parser['input policies'][key]!r} is unsupported"
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

    edges = _float_list(parser["survey"]["z_edges"], "[survey] z_edges", 2)
    try:
        area = float(parser["survey"]["area_deg2"])
        min_mag = float(parser["survey"]["min_band_mag"])
        max_mag = float(parser["survey"]["max_band_mag"])
        nbin = parser["survey"].getint("num_z_bins")
        resolution = float(parser["survey"]["resolution"])
        reconstruction_factor = float(parser["survey"]["reconstruction_factor"])
        lya_rest = float(parser["survey"]["lya_rest_angstrom"])
        numerical = {
            key: float(parser["numerical"][key])
            for key in ("k_min", "k_max", "ap_step", "weight_rtol")
        }
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
    if nbin < 1 or not np.isfinite(resolution) or resolution <= 0:
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
    if len(edges) != nbin + 1 or any(b <= a for a, b in zip(edges[:-1], edges[1:])):
        raise ValueError(
            "[survey] z_edges must contain num_z_bins+1 strictly increasing edges"
        )
    if parser["survey"]["evaluation_redshift"] != "geometric_1plusz":
        raise ValueError("[survey] evaluation_redshift must be geometric_1plusz")
    if parser["survey"]["band"] != "r":
        raise ValueError("[survey] band must be the adopted r-band input")
    if len(set(field_names)) != len(field_names):
        raise ValueError("[fields] ids must be unique")

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
        for token in selected.split(","):
            pair = tuple(x.strip() for x in token.split("x"))
            if len(pair) != 2 or any(x not in ids for x in pair):
                raise ValueError(f"[{section}] invalid selected pair {token!r}")
            pair = tuple(sorted(pair, key=field_order.index))
            if pair in pairs:
                raise ValueError(f"[{section}] duplicate selected pair {pair}")
            pairs.append(pair)
        z_eval = float(np.sqrt((1 + z_min) * (1 + z_max)) - 1)
        bins.append(BinConfig(index, z_min, z_max, z_eval, tuple(pairs)))

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
    )


def parse_ini(source=None):
    """Short alias for :func:`parse_survey_ini`."""

    return parse_survey_ini(source)


def _resolve_resource(value, config):
    value = str(value)
    if value.startswith("package:"):
        return "package", value[len("package:") :]
    if config.path is None:
        return "package", value
    return "external", resolve_input_path(value, config.path)


def _reader_samples(field, density, snr, z_source, magnitudes, wavelength):
    if hasattr(density, "sample"):
        d = density.sample(z_source, magnitudes)
        density_values, density_provenance = d["values"], d["provenance"]
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
            s = snr.sample(
                z_source=z_source,
                magnitudes=magnitudes,
                wavelength=wavelength,
                pixel_width_angstrom=field.pixel_width_angstrom,
                exposure_count=field.num_exposures,
            )
            result["variance"] = np.asarray(s["values"])
            result["snr_diagnostics"] = s["provenance"]
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


def _background_values(background, z):
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
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite numeric") from error
    if not np.isfinite(result) or (positive and result <= 0):
        qualifier = "positive finite" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _build_readers(config, fields, stack):
    """Build adopted density/SNR adapters while package paths are materialized."""

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
    """Assemble native :class:`~fishhighz.survey.BinSpec` objects from a parsed configuration.

    ``background`` and ``template`` are caller-prepared objects.  This function
    never invokes CAMB.  ``readers`` may inject deterministic density/SNR
    adapters; when omitted, the configured Step-1 resources are read through
    the strict readers and explicitly named accuracy adapters.  ``prepare=True``
    additionally runs the existing fixed-bin preparation, returning it as
    ``PreparedSurvey.prepared_bins`` is intentionally deferred to Step 4 and
    therefore currently rejected.
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
    fields = config.observed_fields
    registry = ParameterRegistry(
        [
            Parameter(f"{name}_{index}", 1.0, "target", step=0.001)
            for index in range(len(config.bins))
            for name in ("ap", "at")
        ]
    )
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
        for bin_config in config.bins:
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
            grid = gauss_legendre_grid(
                np.linspace(
                    float(config.numerical["k_min"]),
                    float(config.numerical["k_max"]),
                    int(config.numerical["k_intervals"]) + 1,
                ),
                k_order=int(config.numerical["k_order"]),
                mu_order=int(config.numerical["mu_order"]),
                h_fid=float(h_fid),
            )
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
            lo = float(config.survey["min_band_mag"])
            hi = float(config.survey["max_band_mag"])
            partition = breakpoints(densities, snrs, z_sources, lo, hi)
            magnitudes, quadrature = composite(
                partition, int(config.numerical["magnitude_order"])
            )
            active = set(np.unique(selection.selected_pairs).tolist())
            forests, galaxies = {}, {}
            for field_index, item in enumerate(config.fields):
                if field_index not in active:
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
                            "magnitude_bounds": [lo, hi],
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
        provenance = freeze(
            {
                "config": config.provenance,
                "registry_ids": list(registry.ids),
                "fields": [field.id for field in fields],
                "bins": bin_provenance,
                "background": type(background).__name__,
                "template": type(template).__name__,
            }
        )
    return PreparedSurvey(config, fields, registry, tuple(specs), provenance)


def prepare_ini(source=None, *, background, template, readers=None):
    """Parse ``source`` and assemble its native :class:`~fishhighz.survey.BinSpec` objects."""

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
