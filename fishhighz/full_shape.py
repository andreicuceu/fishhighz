"""Independent component dilation and bin-local full-shape parameters."""

import numpy as np

from .models.biases import analytic_density_bias, linear_tabulated_bias
from .models.external import BoundParameters
from .models.kaiser import KaiserModel, Scaling
from .parameters import Parameter, ParameterRegistry

TARGETS = ("alpha_w", "phi_w", "alpha_s", "phi_s", "f")
ISO_TARGETS = ("alpha_iso_w", "phi_w", "alpha_iso_s", "phi_s", "f")


def basis_targets(config):
    """Ordered full-shape parameters for the configured dilation basis."""
    return (
        ISO_TARGETS
        if config.model.get("parameterization", "alpha_phi") == "alpha_iso_phi"
        else TARGETS
    )


def target_names(config):
    """Physical targets retained by this full-shape calculation."""
    return (
        basis_targets(config)[:4]
        if config.model.get("target_set", "full") == "dilation_only"
        else basis_targets(config)
    )


def active_names(config, index, pairs=None):
    """Retain all targets and only nuisances of selected physical tracers."""
    pairs = config.bins[index].selected_pairs if pairs is None else pairs
    ids = {name for pair in pairs for name in pair}
    tracers = {item.tracer for item in config.fields if item.observed.id in ids}
    return (
        target_names(config)
        + tuple("b_" + name for name in ("lya", "qso", "lbg", "lae") if name in tracers)
        + (("beta_lya",) if "lya" in tracers else ())
    )


def make_registry(config, background):
    from .survey_config import _background_values

    parameters = []
    for index, bin_config in enumerate(config.bins):
        if not bin_config.selected_pairs:
            continue
        _, f = _background_values(background, bin_config.z_eval)
        values = dict.fromkeys(basis_targets(config), 1.0)
        values.update(f=f, beta_lya=float(config.model["forest_beta"]))
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
                if name in basis_targets(config)[:4]
                else float(config.numerical["relative_step"])
                * (abs(values[name]) or 1.0)
            )
            parameters.append(
                Parameter(
                    f"{name}_{index}",
                    values[name],
                    "target" if name in target_names(config) else "nuisance",
                    step=step,
                )
            )
    return ParameterRegistry(parameters)


def make_model(
    config, bin_config, registry, template, fields, biases, betas, widths, f, growth
):
    index = bin_config.index - 1
    active = active_names(config, index)
    bindings = {name: f"{name}_{index}" for name in target_names(config)}
    basis = config.model.get("parameterization", "alpha_phi")
    amplitude = "alpha_iso" if basis == "alpha_iso_phi" else "alpha"
    smooth_name, wiggle_name = basis_targets(config)[2], basis_targets(config)[0]
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
        f="f" if "f" in bindings else f,
        local_names=tuple(bindings),
        smooth=Scaling(basis, **{amplitude: smooth_name, "phi": "phi_s"}),
        wiggle=Scaling(basis, **{amplitude: wiggle_name, "phi": "phi_w"}),
        z=bin_config.z_eval,
        growth=growth,
    )
    return model, BoundParameters(registry, tuple(bindings), bindings)


def reported_constraint(result, target_ids, sigma8):
    """Marginalize active nuisances, then apply the fixed-sigma8 Jacobian."""
    if result.diagnostics.rank != len(result.registry.ids):
        return "unavailable", None, None, None
    covariance = result.marginalized_covariance(target_ids)
    jacobian = np.array([1.0, 1.0, 1.0, 1.0, sigma8][: len(target_ids)])
    covariance = covariance * jacobian[:, None] * jacobian[None, :]
    errors = np.sqrt(np.diag(covariance))
    return (
        "available",
        covariance,
        errors,
        covariance / errors[:, None] / errors[None, :],
    )


def validate_template_coverage(spec, *, step_scale=None, numerical=False):
    """Validate the complete observed domain using the actual derivative schedule.

    Angular extrema are at mu=0 and mu=1, so each independently rescaled
    component requires [k_min/max(ap,at), k_max/min(ap,at)] template support.
    A None step_scale validates the fiducial only during survey preparation.
    """
    from .derivatives import _schedule
    from .kernels.kaiser import _resolve, _scales
    from .models.kaiser import _BASES

    theta = spec.p3d.registry.fiducials
    schedules = (
        (None,) * len(spec.p3d.routes)
        if step_scale is None
        else _schedule(spec.p3d, theta, None, step_scale, numerical)
    )
    for route, schedule in zip(spec.p3d.routes, schedules):
        model = route.provider.model
        binding = route.provider.parameters.binding.local_to_global
        states = [("fiducial", theta)]
        for index, stencil in () if schedule is None else schedule:
            if stencil is None:
                continue
            for point in stencil.points:
                perturbed = theta.copy()
                perturbed[index] = point
                states.append(
                    (
                        f"parameter {spec.p3d.registry.ids[index]}, {stencil.method}, value={point}",
                        perturbed,
                    )
                )
        for context, state in states:
            local = state[binding]
            values = _resolve(model.fixed, model.destinations, model.sources, local)
            n_fields = len(model.fields)
            for component_index, component in enumerate(("smooth", "wiggle")):
                coordinates = values[
                    2 * n_fields + 1 + 2 * component_index : 2 * n_fields
                    + 3
                    + 2 * component_index
                ]
                basis = model.bases[component_index]
                if coordinates[0] <= 0 or coordinates[1] <= (
                    -1 if basis == "alpha_iso_epsilon" else 0
                ):
                    raise ValueError(f"{component}, {context}: invalid {basis} scaling")
                with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                    ap, at, _ = _scales(coordinates, tuple(_BASES).index(basis))
                    mapped = (
                        spec.grid.k_min / max(ap, at),
                        spec.grid.k_max / min(ap, at),
                    )
                if (
                    not np.all(np.isfinite(mapped))
                    or mapped[0] <= 0
                    or mapped[0] < model.template.k[0]
                    or mapped[1] > model.template.k[-1]
                ):
                    raise ValueError(
                        f"{component}, {context}: fixed observed cuts "
                        f"{(spec.grid.k_min, spec.grid.k_max)} map to {mapped} "
                        f"outside template domain {model.template.domain}; pad template for all derivative stencils"
                    )


def _category_cut(config, fields, pair):
    """Observed cutoff of one selected pair in h_fid/Mpc."""
    kinds = (fields[pair[0]].kind, fields[pair[1]].kind)
    suffix = (
        "forest_forest"
        if kinds == ("forest", "forest")
        else "galaxy_galaxy"
        if kinds == ("galaxy", "galaxy")
        else "galaxy_forest"
    )
    return float(config.numerical.get("k_max_" + suffix, config.numerical["k_max"]))


def _selected_nodes(prepared, pairs, lower, upper):
    """Restrict measured spectra and nodes, retaining their covariance closure."""
    from dataclasses import replace

    from .covariance import gaussian_covariance
    from .fields import PairSelection
    from .fisher import factor_covariance
    from .models.external import P3DProvider, PreparedP3D

    original = prepared.p3d.selection
    selected = PairSelection(original.fields, pairs)
    node = np.flatnonzero((prepared.k >= lower) & (prepared.k < upper))
    if not len(node):
        raise ValueError("category interval contains no integration nodes")
    if len(node) == len(prepared.k) and np.array_equal(
        selected.selected_pairs, original.selected_pairs
    ):
        return prepared
    required = [tuple(pair) for pair in selected.required_pairs.tolist()]
    old_required = [tuple(pair) for pair in original.required_pairs.tolist()]
    columns = [old_required.index(pair) for pair in required]
    routes = []
    for route in prepared.p3d.routes:
        declared = [
            tuple(pair) for pair in route.pairs.tolist() if tuple(pair) in required
        ]
        if declared:
            provider = route.provider
            routes.append(
                P3DProvider(
                    provider.label,
                    provider.model,
                    provider.parameters,
                    declared,
                    jacobian=provider.jacobian,
                    analytic_ids=provider.analytic_ids,
                )
            )
    p3d = PreparedP3D(prepared.p3d.registry, selected, routes)
    total = prepared.total[np.ix_(node, columns)]
    factors = factor_covariance(
        gaussian_covariance(total, prepared.modes[node], selected)
    )
    return replace(
        prepared,
        p3d=p3d,
        k=prepared.k[node],
        mu=prepared.mu[node],
        modes=prepared.modes[node],
        response=prepared.response[node],
        products=prepared.products[np.ix_(node, columns)],
        noise=prepared.noise[np.ix_(node, columns)],
        power=prepared.power[np.ix_(node, columns)],
        total=total,
        factors=factors,
    )


def _interval_result(
    prepared, config, *, pairs=None, batch_size, step_scale, numerical
):
    """Sum disjoint measured intervals before any nuisance marginalization."""
    from .forecast import run_bin
    from .results import combine_results

    selection = prepared.p3d.selection
    fields = selection.fields
    available = [tuple(pair) for pair in selection.selected_pairs.tolist()]
    pairs = available if pairs is None else list(pairs)
    cuts = {_category_cut(config, fields, pair) for pair in pairs}
    bounds = (float(config.numerical["k_min"]), *sorted(cuts))
    contributions = []
    for lower, upper in zip(bounds[:-1], bounds[1:]):
        active = [
            pair for pair in pairs if _category_cut(config, fields, pair) >= upper
        ]
        if active:
            subset = _selected_nodes(prepared, active, lower, upper)
            contributions.append(
                run_bin(
                    subset,
                    batch_size=batch_size,
                    step_scale=step_scale,
                    numerical=numerical,
                ).result
            )
    return combine_results(contributions)


def run_full_shape(forecast, *, batch_size, step_scale, numerical, individuals=True):
    from .fields import PairSelection
    from .public import SpectrumConstraint, SurveyResult, _convergence
    from .results import combine_results
    from .survey import freeze
    from .survey_config import _background_values

    prepared = forecast.prepare()
    config = forecast.config
    for spec in prepared.survey.bins:
        validate_template_coverage(spec, step_scale=step_scale, numerical=numerical)
    joint, individual, excluded = [], [], []
    fields = prepared.survey.fields
    active_bins = iter(zip(prepared.survey.bins, prepared.bins))
    joint_results = []
    selected_by_bin = {}
    for index, bin_config in enumerate(config.bins):
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
        bin_config = config.bins[index]
        sigma8, f = _background_values(prepared.background, bin_config.z_eval)
        target_basis = target_names(config)
        targets = tuple(f"{name}_{index}" for name in target_basis)
        reported_ids = (
            targets[:-1] + (f"fsigma8_{index}",) if "f" in target_basis else targets
        )

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
                    result, targets, sigma8
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
                reported_ids,
                covariance,
                errors,
                correlations,
                tuple([1.0] * 4 + ([f * sigma8] if "f" in target_basis else [])),
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
    grid_by_bin = {}
    for spec in prepared.survey.bins:
        grid = spec.grid
        selected_pairs = [
            tuple(pair) for pair in spec.p3d.selection.selected_pairs.tolist()
        ]
        cuts = sorted({_category_cut(config, fields, pair) for pair in selected_pairs})
        intervals = []
        for lower, upper in zip((grid.k_min, *cuts[:-1]), cuts):
            active = [
                tuple(fields[i].id for i in pair)
                for pair in selected_pairs
                if _category_cut(config, fields, pair) >= upper
            ]
            intervals.append(
                {
                    "bounds": (lower, upper),
                    "radial_nodes": int(
                        np.count_nonzero((grid.k >= lower) & (grid.k < upper))
                    ),
                    "selected_pairs": active,
                }
            )
        grid_by_bin[spec.id] = {
            "k_nodes": grid.k.tolist(),
            "k_weights": grid.w_k.tolist(),
            "mu_nodes": grid.mu.tolist(),
            "mu_weights": grid.w_mu.tolist(),
            "category_intervals": intervals,
        }
    resolved = freeze(
        {
            "config": config.provenance,
            "inputs": prepared.input_identity,
            "parameter_order": prepared.survey.registry.ids,
            "requested_k_min": float(config.numerical["k_min"]),
            "requested_k_max": float(config.numerical["k_max"]),
            "effective_k_max": {
                name: float(
                    config.numerical.get("k_max_" + name, config.numerical["k_max"])
                )
                for name in ("galaxy_galaxy", "galaxy_forest", "forest_forest")
            },
            "selected_pairs_by_bin": selected_by_bin,
            "grid_by_bin": grid_by_bin,
            "target_names": target_names(config),
            "quadrature": dict(config.numerical),
            "growth_convention": "f_sigma8_fid; template normalization fixed",
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


def save_full_shape(result, output, *, schema_name="fishhighz-full-shape-result"):
    import json

    from .public import _plain

    arrays = {"combined_fisher": result.combined.data_fisher}
    records = []
    for index, item in enumerate((*result.constraints, *result.excluded)):
        prefix = f"record_{index:04d}"
        record = {
            key: _plain(getattr(item, key))
            for key in (
                "bin_index",
                "bin_id",
                "bounds",
                "z_eval",
                "kind",
                "pair",
                "parameter_ids",
                "status",
                "target_ids",
                "target_fiducials",
                "sigma8_fid",
            )
        }
        record["arrays"] = {}
        if item.fisher is not None:
            diag = item.fisher.diagnostics
            record["parameters"] = [vars(p) for p in item.fisher.registry.parameters]
            record["rank"] = {
                "rank": diag.rank,
                "dimension": len(diag.ids),
                "condition": diag.condition if np.isfinite(diag.condition) else None,
                "tolerance": diag.tolerance,
                "basis": diag.basis,
            }
            for name, value in (
                ("fisher", item.fisher.data_fisher),
                ("scales", diag.scales),
                ("eigenvalues", diag.eigenvalues),
                ("null_directions", diag.null_directions),
                ("target_covariance", item.target_covariance),
                ("target_errors", item.target_errors),
                ("target_correlations", item.target_correlations),
            ):
                if value is not None:
                    key = prefix + "_" + name
                    arrays[key] = value
                    record["arrays"][name] = key
        records.append(record)
    settings = {
        "schema": {"name": schema_name, "version": 1},
        "records": records,
        "resolved_settings": _plain(result.resolved_settings),
        "convergence": _plain(result.convergence),
        "combined_parameters": [vars(p) for p in result.combined.registry.parameters],
        "result_counts": {
            "individual": len(result.individual),
            "joint": len(result.joint),
            "excluded": len(result.excluded),
        },
    }
    (output / "settings.json").write_text(
        json.dumps(settings, indent=2, allow_nan=False) + "\n"
    )
    np.savez_compressed(output / "results.npz", **arrays)
    return output
