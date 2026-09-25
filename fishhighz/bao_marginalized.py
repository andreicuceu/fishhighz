"""Bin-local BAO dilation with fixed galaxy growth and tracer nuisances."""

import numpy as np

from .full_shape import _interval_result, validate_template_coverage
from .models.biases import analytic_density_bias, linear_tabulated_bias
from .models.external import BoundParameters
from .models.kaiser import KaiserModel, Scaling
from .parameters import Parameter, ParameterRegistry

TARGETS = ("alpha_iso", "phi")
NUISANCES = ("b_lya", "b_qso", "b_lbg", "b_lae", "beta_lya")


def active_names(config, index, pairs=None):
    """Two BAO coordinates and nuisances represented in selected spectra."""
    pairs = config.bins[index].selected_pairs if pairs is None else pairs
    ids = {name for pair in pairs for name in pair}
    tracers = {item.tracer for item in config.fields if item.observed.id in ids}
    return TARGETS + tuple(
        name
        for name in NUISANCES
        if name.removeprefix("b_").removeprefix("beta_") in tracers
    )


def make_registry(config, background):
    """Construct independent target and nuisance parameters for nonempty bins."""
    parameters = []
    for index, bin_config in enumerate(config.bins):
        if not bin_config.selected_pairs:
            continue
        values = dict.fromkeys(TARGETS, 1.0)
        values["beta_lya"] = float(config.model["forest_beta"])
        for item in config.fields:
            bias = (
                analytic_density_bias(bin_config.z_eval, item.tracer)
                if item.bias_redshifts is None
                else linear_tabulated_bias(item.bias_redshifts, item.bias_values)(
                    bin_config.z_eval
                )
            )
            key = "b_" + item.tracer
            if key in values and not np.isclose(values[key], bias, rtol=1e-12, atol=0):
                raise ValueError("shared forest bias requires equal fiducials")
            values[key] = float(bias)
        for name in active_names(config, index):
            step = (
                float(config.numerical["scale_step"])
                if name in TARGETS
                else float(config.numerical["relative_step"])
                * (abs(values[name]) or 1.0)
            )
            parameters.append(
                Parameter(
                    f"{name}_{index}",
                    values[name],
                    "target" if name in TARGETS else "nuisance",
                    step=step,
                )
            )
    return ParameterRegistry(parameters)


def make_model(
    config, bin_config, registry, template, fields, biases, betas, widths, f, growth
):
    """Bind selected tracer nuisances while keeping the shared galaxy f fixed."""
    index = bin_config.index - 1
    active = active_names(config, index)
    bindings = {name: f"{name}_{index}" for name in TARGETS}
    biases, betas = dict(biases), dict(betas)
    for field_index, item in enumerate(config.fields):
        name = "b_" + item.tracer
        if name in active:
            local = f"bias_field_{field_index}"
            biases[item.observed.id] = local
            bindings[local] = f"{name}_{index}"
        if item.observed.kind == "forest" and "beta_lya" in active:
            local = f"beta_field_{field_index}"
            betas[item.observed.id] = local
            bindings[local] = f"beta_lya_{index}"
    model = KaiserModel(
        template,
        fields,
        biases=biases,
        betas=betas,
        widths=widths,
        f=f,
        local_names=tuple(bindings),
        smooth=Scaling("ap_at", ap=1, at=1),
        wiggle=Scaling("alpha_iso_phi", alpha_iso="alpha_iso", phi="phi"),
        z=bin_config.z_eval,
        growth=growth,
    )
    return model, BoundParameters(registry, tuple(bindings), bindings)


def reported_constraint(result, target_ids):
    """Marginalize active nuisances only for a full-rank Fisher matrix."""
    if result.diagnostics.rank != len(result.registry.ids):
        return "unavailable", None, None, None
    covariance = result.marginalized_covariance(target_ids)
    errors = np.sqrt(np.diag(covariance))
    return (
        "available",
        covariance,
        errors,
        covariance / errors[:, None] / errors[None, :],
    )


def run_bao_marginalized(
    forecast, *, batch_size, step_scale, numerical, individuals=True
):
    """Compute selected joint and individual constraints with fixed covariance."""
    from .fields import PairSelection
    from .public import SpectrumConstraint, SurveyResult, _convergence
    from .results import combine_results
    from .survey import freeze
    from .survey_config import _background_values

    prepared = forecast.prepare()
    config = forecast.config
    for spec in prepared.survey.bins:
        validate_template_coverage(spec, step_scale=step_scale, numerical=numerical)
    joint, individual, excluded, joint_results = [], [], [], []
    fields = prepared.survey.fields
    active_bins = iter(zip(prepared.survey.bins, prepared.bins))
    selected_by_bin = {}
    growth_by_bin = {}
    for index, bin_config in enumerate(config.bins):
        sigma8, growth = _background_values(prepared.background, bin_config.z_eval)
        growth_by_bin[index] = growth
        if bin_config.selected_pairs:
            spec, prepared_bin = next(active_bins)
            selected = {
                tuple(fields[i].id for i in pair)
                for pair in spec.p3d.selection.selected_pairs.tolist()
            }
        else:
            spec = prepared_bin = None
            selected = set()
        selected_by_bin[index] = sorted(selected)
        targets = tuple(f"{name}_{index}" for name in TARGETS)

        def record(result, pair, status=None):
            ids = tuple(
                f"{name}_{index}"
                for name in (
                    ()
                    if spec is None
                    else active_names(config, index, None if pair is None else [pair])
                )
            )
            covariance = errors = correlations = None
            if result is not None:
                result = result.fix_except(ids)
                status, covariance, errors, correlations = reported_constraint(
                    result, targets
                )
            return SpectrumConstraint(
                index,
                f"desi2_accuracy-{bin_config.index}" if spec is None else spec.id,
                (bin_config.z_min, bin_config.z_max),
                bin_config.z_eval,
                "joint" if pair is None else "individual",
                pair,
                ids,
                status,
                result,
                None,
                None,
                None,
                targets,
                covariance,
                errors,
                correlations,
                (1.0, 1.0),
                sigma8,
            )

        if prepared_bin is None:
            joint.append(record(None, None, "excluded"))
        else:
            result = _interval_result(
                prepared_bin,
                config,
                batch_size=batch_size,
                step_scale=step_scale,
                numerical=numerical,
            )
            joint_results.append(result)
            joint.append(record(result, None))
        if individuals and prepared_bin is not None:
            for pair_values in spec.p3d.selection.selected_pairs.tolist():
                pair = tuple(pair_values)
                pair_names = tuple(fields[i].id for i in pair)
                own_result = _interval_result(
                    prepared_bin,
                    config,
                    pairs=[pair],
                    batch_size=batch_size,
                    step_scale=step_scale,
                    numerical=numerical,
                )
                individual.append(record(own_result, pair_names))
        for pair in PairSelection(fields).selected_pairs.tolist():
            pair_names = tuple(fields[i].id for i in pair)
            if pair_names not in selected:
                excluded.append(record(None, pair_names, "excluded"))
    combined = combine_results(joint_results)
    resolved = freeze(
        {
            "config": config.provenance,
            "inputs": prepared.input_identity,
            "parameter_order": prepared.survey.registry.ids,
            "selected_pairs_by_bin": selected_by_bin,
            "target_names": TARGETS,
            "quadrature": dict(config.numerical),
            "growth_convention": "shared fixed fiducial f(z); no growth parameter or prior",
            "fixed_growth_by_bin": growth_by_bin,
            "covariance": "fiducial; individual spectra use their own covariance",
            "run": {
                "batch_size": batch_size,
                "step_scale": step_scale,
                "numerical": numerical,
                "individuals": individuals,
            },
        }
    )
    return SurveyResult(
        config,
        prepared,
        tuple(individual),
        tuple(joint),
        tuple(excluded),
        _convergence(prepared.bins),
        resolved,
        combined,
    )


def save_bao_marginalized(result, output):
    """Write named Fisher and marginalized target matrices in schema version 1."""
    from .full_shape import save_full_shape

    return save_full_shape(
        result, output, schema_name="fishhighz-bao-marginalized-result"
    )
