"""Lazy INI-driven forecast facade over the FishHighz numerical API."""

from __future__ import annotations

import hashlib
import json
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType

import numpy as np

from .cosmology import prepare_camb
from .fields import PairSelection
from .forecast import PreparedBin, prepare_bin, run_forecast
from .models.external import P3DProvider, PreparedP3D
from .models.templates import load_template
from .resources import bundled_path, resolve_input_path
from .results import FisherResult
from .survey import BinSpec, freeze
from .survey_config import (
    PreparedSurvey,
    SurveyConfig,
    _normalise_readers,
    parse_survey_ini,
    prepare_survey,
)


def _plain(value):
    """Return JSON-compatible values from immutable FishHighz metadata."""

    if isinstance(value, MappingProxyType) or hasattr(value, "items"):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [_plain(item) for item in value.tolist()]
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _constraints(result):
    """Extract two-parameter diagnostics without hiding rank failures."""

    if result.diagnostics.rank != len(result.registry.ids):
        return "unavailable", None, None, None
    errors = result.marginalized_errors()
    correlation = result.correlations()[0, 1]
    return "available", float(errors[0]), float(errors[1]), float(correlation)


@dataclass(frozen=True)
class SpectrumConstraint:
    """One selected or explicitly excluded spectrum constraint."""

    bin_index: int
    bin_id: str
    bounds: tuple[float, float]
    z_eval: float
    kind: str
    pair: tuple[str, str] | None
    parameter_ids: tuple[str, ...]
    status: str
    fisher: FisherResult | None
    sigma_ap: float | None
    sigma_at: float | None
    correlation: float | None
    target_ids: tuple[str, ...] = ()
    target_covariance: object = None
    target_errors: object = None
    target_correlations: object = None
    target_fiducials: tuple[float, ...] = ()
    sigma8_fid: float | None = None

    @property
    def result(self):
        """Alias for the retained FisherResult."""

        return self.fisher

    @property
    def available(self):
        return self.status == "available"


@dataclass(frozen=True)
class PreparedForecast:
    """Prepared background, template, native survey and fixed joint bins."""

    config: SurveyConfig
    survey: PreparedSurvey
    background: object
    template: object
    bins: tuple[PreparedBin, ...]
    input_identity: object = field(default_factory=lambda: MappingProxyType({}))

    @property
    def prepared_bins(self):
        return self.bins

    @property
    def bin_specs(self):
        return self.survey.bins


@dataclass(frozen=True)
class SurveyResult:
    """Public result retaining joint and own-covariance spectra."""

    config: SurveyConfig
    prepared: PreparedForecast
    individual: tuple[SpectrumConstraint, ...]
    joint: tuple[SpectrumConstraint, ...]
    excluded: tuple[SpectrumConstraint, ...]
    convergence: object
    resolved_settings: object
    combined: FisherResult

    @property
    def individual_spectra(self):
        return self.individual

    @property
    def joint_constraints(self):
        return self.joint

    @property
    def constraints(self):
        return self.individual + self.joint

    @property
    def fisher_results(self):
        return tuple(item.fisher for item in self.constraints)

    def save(self, path):
        """Save identities/statuses as JSON and retained matrices as NPZ."""

        output = Path(path).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=False)
        if self.config.model.get("mode", "bao") == "full_shape":
            from .full_shape import save_full_shape

            return save_full_shape(self, output)
        if self.config.model.get("mode", "bao") == "bao_marginalized":
            from .bao_marginalized import save_bao_marginalized

            return save_bao_marginalized(self, output)
        records = []
        for item in (*self.constraints, *self.excluded):
            records.append(
                {
                    "bin_index": item.bin_index,
                    "bin_id": item.bin_id,
                    "bounds": list(item.bounds),
                    "z_eval": item.z_eval,
                    "kind": item.kind,
                    "pair": None if item.pair is None else list(item.pair),
                    "parameter_ids": list(item.parameter_ids),
                    "status": item.status,
                    "sigma_ap": item.sigma_ap,
                    "sigma_at": item.sigma_at,
                    "correlation": item.correlation,
                }
            )
        settings = {
            "schema": {"name": "fishhighz-forecast-result", "version": 1},
            "resolved_settings": _plain(self.resolved_settings),
            "convergence": _plain(self.convergence),
            "records": records,
            "result_counts": {
                "individual": len(self.individual),
                "joint": len(self.joint),
                "excluded": len(self.excluded),
            },
        }
        (output / "settings.json").write_text(
            json.dumps(settings, indent=2, allow_nan=False) + "\n"
        )
        np.savez_compressed(
            output / "results.npz",
            individual_fisher=np.asarray(
                [item.fisher.data_fisher for item in self.individual], dtype=float
            ),
            joint_fisher=np.asarray(
                [item.fisher.data_fisher for item in self.joint], dtype=float
            ),
            combined_fisher=np.asarray(self.combined.data_fisher, dtype=float),
        )
        return output


def _resource_path(config, value, stack):
    value = str(value)
    if value.startswith("package:") or config.path is None:
        return stack.enter_context(bundled_path(value.removeprefix("package:")))
    return resolve_input_path(value, config.path)


def _input_identifier(config, value):
    """Keep package identities independent of temporary extraction paths."""

    if str(value).startswith("package:"):
        return str(value)
    if config.path is None:
        return "package:" + str(value)
    return str(resolve_input_path(value, config.path))


def _object_identity(value, source):
    """Retain known input digests without opening an injected object's files."""

    result = {"source": source, "type": type(value).__name__, "sha256": None}
    provenance = getattr(getattr(value, "reader", value), "provenance", {})
    if "sha256" in provenance:
        result["sha256"] = provenance["sha256"]
        if "paths" in provenance:
            result["identifiers"] = provenance["paths"]
        elif "path" in provenance:
            result["identifier"] = provenance["path"]
    if getattr(value, "source_sha256", None) is not None:
        result.update(identifier=value.source_path, sha256=value.source_sha256)
    return result


def _reader_identity(config, readers, source):
    result = {}
    for item in config.fields:
        values = {}
        for kind in ("density", "snr"):
            reader = readers[item.observed.id][kind]
            if reader is None:
                continue
            identity = _object_identity(reader, source)
            if source == "configured":
                identifier = _input_identifier(config, getattr(item, kind))
                if kind == "snr":
                    identity["identifiers"] = [
                        identifier.rstrip("/") + "/" + Path(path).name
                        for path in identity["identifiers"]
                    ]
                else:
                    identity["identifier"] = identifier
            values[kind] = identity
        result[item.observed.id] = values
    return result


def _default_background(config, stack, identity):
    path = _resource_path(config, config.cosmology["camb_ini"], stack)
    identity.update(
        source="configured",
        identifier=_input_identifier(config, config.cosmology["camb_ini"]),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )
    return prepare_camb(
        ini=path,
        template_growth_redshift=float(config.cosmology["template_redshift"]),
        damping_reference_redshift=float(
            config.cosmology["damping_reference_redshift"]
        ),
        redshifts=[item.z_eval for item in config.bins],
    )


def _default_template(config, background, stack):
    return load_template(
        _resource_path(config, config.cosmology["template"], stack),
        h_fid=float(background.h_fid),
        z_ref=float(config.cosmology["template_redshift"]),
    )


def _pair_spec(spec, pair):
    """Construct one spectrum with a new, one-spectrum covariance closure."""

    fields = spec.p3d.selection.fields
    pair = tuple(pair)
    selection = PairSelection(fields, [pair])
    required = {tuple(value) for value in selection.required_pairs.tolist()}
    providers = []
    for route in spec.p3d.routes:
        declared = [
            tuple(value) for value in route.pairs.tolist() if tuple(value) in required
        ]
        if not declared:
            continue
        original = route.provider
        providers.append(
            P3DProvider(
                original.label,
                original.model,
                original.parameters,
                declared,
                jacobian=original.jacobian,
                analytic_ids=original.analytic_ids,
            )
        )
    p3d = PreparedP3D(spec.p3d.registry, selection, providers)
    active_ids = {fields[index].id for index in np.unique(selection.selected_pairs)}
    forests = {
        name: value
        for name, value in (spec.forests or {}).items()
        if name in active_ids
    }
    galaxies = {
        name: value
        for name, value in (spec.galaxies or {}).items()
        if name in active_ids
    }
    full_noise = None
    if spec.full_noise is not None:
        original_pairs = [
            tuple(value) for value in spec.p3d.selection.required_pairs.tolist()
        ]
        columns = [
            original_pairs.index(tuple(value))
            for value in selection.required_pairs.tolist()
        ]
        full_noise = spec.full_noise[:, columns]
        forests = galaxies = None
        independent_sampling = None
    else:
        independent_sampling = spec.independent_sampling
    return BinSpec(
        spec.id + ":" + "x".join(fields[index].id for index in pair),
        spec.geometry,
        spec.grid,
        p3d,
        spec.responses,
        forests=forests,
        galaxies=galaxies,
        independent_sampling=independent_sampling,
        full_noise=full_noise,
    )


def _convergence(prepared_bins):
    def record(convergence):
        if convergence is None:
            return {"status": "not_applicable"}
        if hasattr(convergence, "items"):
            return _plain(dict(convergence))
        if hasattr(convergence, "__dict__"):
            return _plain(vars(convergence))
        return {"status": str(convergence)}

    return freeze(
        {
            item.id: {
                name: record(weight.convergence)
                for name, weight in item.weights.items()
            }
            for item in prepared_bins
        }
    )


class Forecast:
    """Lazy native-INI forecast facade with deterministic injection hooks."""

    def __init__(
        self,
        source=None,
        *,
        background=None,
        template=None,
        readers=None,
        background_factory=None,
        template_factory=None,
        readers_factory=None,
    ):
        if (
            source is not None
            and str(source) == "desi2_accuracy.ini"
            and not Path(source).exists()
        ):
            source = None
        self.config = parse_survey_ini(source)
        self._background = background
        self._template = template
        self._readers = readers
        self._background_factory = background_factory
        self._template_factory = template_factory
        self._readers_factory = readers_factory
        self._prepared = None

    def prepare(self):
        """Prepare background, native survey bins and fixed covariance state."""

        if self._prepared is not None:
            return self._prepared
        with ExitStack() as stack:
            identities = {"ini": self.config.source_identity}
            background_identity = {}
            background = self._background
            if background is None:
                background = (
                    self._background_factory(self.config)
                    if self._background_factory is not None
                    else _default_background(self.config, stack, background_identity)
                )
            identities["background"] = background_identity or _object_identity(
                background,
                "injected" if self._background is not None else "injected_factory",
            )
            template = self._template
            if template is None:
                template = (
                    self._template_factory(self.config, background)
                    if self._template_factory is not None
                    else _default_template(self.config, background, stack)
                )
            template_source = (
                "injected"
                if self._template is not None
                else "injected_factory"
                if self._template_factory is not None
                else "configured"
            )
            identities["template"] = _object_identity(template, template_source)
            if template_source == "configured":
                identities["template"]["identifier"] = _input_identifier(
                    self.config, self.config.cosmology["template"]
                )
            readers = self._readers
            if readers is None and self._readers_factory is not None:
                readers = self._readers_factory(self.config)
            reader_source = (
                "injected"
                if self._readers is not None
                else "injected_factory"
                if self._readers_factory is not None
                else "configured"
            )
            readers = _normalise_readers(
                self.config, self.config.fields, readers, stack
            )
            identities["fields"] = _reader_identity(self.config, readers, reader_source)
            survey = prepare_survey(
                self.config,
                background=background,
                template=template,
                readers=readers,
            )
            if self.config.model.get("mode", "bao") in (
                "full_shape",
                "bao_marginalized",
            ):
                from .full_shape import validate_template_coverage

                for spec in survey.bins:
                    validate_template_coverage(spec)
            bins = tuple(prepare_bin(spec) for spec in survey.bins)
        self._prepared = PreparedForecast(
            self.config, survey, background, template, bins, freeze(identities)
        )
        return self._prepared

    def run(
        self, *, batch_size=2048, step_scale=None, numerical=False, individuals=True
    ):
        """Run joint and own-covariance individual constraints."""

        if self.config.model.get("mode", "bao") == "full_shape":
            from .full_shape import run_full_shape

            return run_full_shape(
                self,
                batch_size=batch_size,
                step_scale=1.0 if step_scale is None else step_scale,
                numerical=numerical,
                individuals=individuals,
            )
        if self.config.model.get("mode", "bao") == "bao_marginalized":
            from .bao_marginalized import run_bao_marginalized

            return run_bao_marginalized(
                self,
                batch_size=batch_size,
                step_scale=1.0 if step_scale is None else step_scale,
                numerical=numerical,
                individuals=individuals,
            )
        if not individuals:
            raise ValueError("individuals=False requires full_shape mode")
        prepared = self.prepare()
        if step_scale is None:
            step_scale = float(self.config.numerical["ap_step"]) / 0.001
        forecast = run_forecast(
            prepared.bins,
            batch_size=batch_size,
            step_scale=step_scale,
            numerical=numerical,
        )
        individual, joint = [], []
        fields = prepared.survey.fields
        for index, (spec, bin_run) in enumerate(
            zip(prepared.survey.bins, forecast.bins)
        ):
            parameter_ids = tuple(f"{name}_{index}" for name in ("ap", "at"))
            reduced = bin_run.result.fix_except(parameter_ids)
            status, sigma_ap, sigma_at, correlation = _constraints(reduced)
            item = self.config.bins[index]
            joint.append(
                SpectrumConstraint(
                    index,
                    spec.id,
                    (item.z_min, item.z_max),
                    item.z_eval,
                    "joint",
                    None,
                    parameter_ids,
                    status,
                    reduced,
                    sigma_ap,
                    sigma_at,
                    correlation,
                )
            )
            for pair_values in spec.p3d.selection.selected_pairs.tolist():
                pair = tuple(int(value) for value in pair_values)
                own = prepare_bin(_pair_spec(spec, pair))
                own_run = run_forecast(
                    [own],
                    batch_size=batch_size,
                    step_scale=step_scale,
                    numerical=numerical,
                )
                own_result = own_run.bins[0].result.fix_except(parameter_ids)
                status, sigma_ap, sigma_at, correlation = _constraints(own_result)
                individual.append(
                    SpectrumConstraint(
                        index,
                        spec.id,
                        (item.z_min, item.z_max),
                        item.z_eval,
                        "individual",
                        (fields[pair[0]].id, fields[pair[1]].id),
                        parameter_ids,
                        status,
                        own_result,
                        sigma_ap,
                        sigma_at,
                        correlation,
                    )
                )
        selected_by_bin = {
            index: {tuple(pair) for pair in item.selected_pairs}
            for index, item in enumerate(self.config.bins)
        }
        all_pairs = PairSelection(fields).selected_pairs.tolist()
        excluded = []
        for index, spec in enumerate(prepared.survey.bins):
            item = self.config.bins[index]
            for pair_values in all_pairs:
                pair = tuple(int(value) for value in pair_values)
                names = (fields[pair[0]].id, fields[pair[1]].id)
                if names in selected_by_bin[index]:
                    continue
                excluded.append(
                    SpectrumConstraint(
                        index,
                        spec.id,
                        (item.z_min, item.z_max),
                        item.z_eval,
                        "individual",
                        names,
                        tuple(f"{name}_{index}" for name in ("ap", "at")),
                        "excluded",
                        None,
                        None,
                        None,
                        None,
                    )
                )
        resolved = freeze(
            {
                "config": self.config.provenance,
                "inputs": prepared.input_identity,
                "background": type(prepared.background).__name__,
                "template": type(prepared.template).__name__,
                "bin_ids": [item.id for item in prepared.bins],
                "parameter_order": list(prepared.survey.registry.ids),
                "individual_covariance": "one independently prepared covariance per selected spectrum",
                "joint_covariance": "full selected-spectrum covariance per bin",
                "run": {
                    "batch_size": batch_size,
                    "step_scale": float(step_scale),
                    "numerical": bool(numerical),
                },
                "result_counts": {
                    "individual": len(individual),
                    "joint": len(joint),
                    "excluded": len(excluded),
                },
            }
        )
        return SurveyResult(
            self.config,
            prepared,
            tuple(individual),
            tuple(joint),
            tuple(excluded),
            _convergence(prepared.bins),
            resolved,
            forecast.combined,
        )


__all__ = [
    "Forecast",
    "PreparedForecast",
    "SpectrumConstraint",
    "SurveyResult",
]
