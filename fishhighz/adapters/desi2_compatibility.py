"""DESI Run-2 15x2pt compatibility examples from a live lyaforecast run.

This adapter is deliberately case-specific.  It observes the arrays prepared by
the installed reference implementation, then uses FishHighz to assemble the
selected Wick covariance and Fisher matrices.  It does not translate arbitrary
lyaforecast configurations or use a captured validation bundle.
"""

import configparser
import contextlib
import io
from pathlib import Path

import numpy as np

from ..kernels.full_sum_weights import METHODS, solve
from ..validation.cases import bins, recipe, selection
from ..validation.compatibility_weights import WeightInputs, forest_noise
from ..validation.numerics import contract, legacy_jacobian, summaries, wick
from ..validation.profile_definitions import (
    REFERENCE,
    REVISION,
    STOPPING,
    forecast_selection,
)
from ..validation.reference_capture import imported_reference, resolved_resources

CASE = "lya_qso_lbg_lae_15x2pt"
PROFILES = ("full-compatibility", "fixed-compatibility")


def _configuration(reference, work_directory):
    """Validate the one supported INI and make an output-local effective copy."""
    root = Path(reference).resolve()
    source = root / "examples" / "desi2" / f"{CASE}.ini"
    config = configparser.ConfigParser()
    config.optionxform = str
    if config.read(source) != [str(source)]:
        raise ValueError(f"missing DESI Run-2 configuration: {source}")
    if {section: dict(config[section]) for section in config.sections()} != recipe(
        CASE
    ):
        raise ValueError(
            "installed 15x2pt configuration differs from the supported recipe"
        )
    work = Path(work_directory).resolve()
    work.mkdir(parents=True, exist_ok=False)
    config["output"]["filename"] = str(work / "outputs" / "forecast")
    effective = work / "effective.ini"
    with effective.open("w") as stream:
        config.write(stream)
    return root, source, config, effective


def _capture(reference, work_directory):
    """Prepare only selected spectra from one initialized NewForecast model."""
    import lyaforecast.forecast_new as module

    root, source, config, effective = _configuration(reference, work_directory)
    identity = imported_reference(root)
    identity["resources"] = resolved_resources(root, config)
    original = selection(CASE)
    field_ids = [field.id for field in original.fields]
    captured = []
    with contextlib.redirect_stdout(io.StringIO()):
        forecast = module.NewForecast(effective)
        power = forecast.power_spectrum
        for index, z_mean in enumerate(forecast.survey.z_bin_centres):
            selected = forecast_selection(CASE, index)
            required_pairs = [tuple(pair) for pair in selected.required_pairs.tolist()]
            selected_pairs = [tuple(pair) for pair in selected.selected_pairs.tolist()]
            required_labels = [
                f"{field_ids[i]}_{field_ids[j]}" for i, j in required_pairs
            ]
            selected_labels = [
                f"{field_ids[i]}_{field_ids[j]}" for i, j in selected_pairs
            ]
            signals = {}
            totals = {}
            pair_inputs = {}
            lo, hi = forecast.survey.z_bin_edges[:, index]
            for label in required_labels:
                covariance = forecast._covariance[label]
                covariance(lo, hi)
                covariance.compute_eff_density_and_noise(label)
                signals[label] = np.array(
                    [
                        power.compute_p3d_hmpc_smooth(
                            z_mean,
                            power.k,
                            mu,
                            covariance.pix_width_kms,
                            covariance.pix_res_kms,
                            label,
                        )
                        for mu in power.mu
                    ]
                )
                totals[label] = np.array(
                    [
                        covariance.compute_total_power(power.k, mu, label)
                        for mu in power.mu
                    ]
                )
                weights = covariance._weights
                context = {
                    name: getattr(covariance, name, None)
                    for name in (
                        "_z_mean",
                        "_zq",
                        "_pix_kms",
                        "_res_kms",
                        "forest_length",
                        "_distance_to_velocity",
                        "_angle_to_distance",
                        "_aliasing_weights",
                        "_effective_noise_power",
                        "_w_lya",
                        "_volume",
                    )
                }
                context["magnitudes"] = weights.maglist.copy()
                if weights._lya_tracer is not None:
                    context.update(
                        density=weights._get_dn_dkmsdm(
                            weights._zq, weights.maglist, weights._lya_tracer
                        ),
                        variance=weights._get_pix_var_m(),
                        auxiliary_signal=float(weights._p3d_w),
                        auxiliary_p1d=float(weights._p1d_w),
                    )
                pair_inputs[label] = context
            k = np.tile(power.k, len(power.mu))
            mu_flat = np.repeat(power.mu, len(power.k))
            widths = []
            for label in selected_labels:
                reconstruction = (
                    1.0 if "lya" in label else forecast.reconstruction_factor
                )
                transverse = (
                    3.26
                    * forecast.cosmo.growth_factor_ratios[index]
                    / np.sqrt(reconstruction)
                )
                widths.append(
                    [
                        (1 + forecast.cosmo.growth_rate_zbins[index]) * transverse,
                        transverse,
                    ]
                )
            observed_j = np.concatenate(
                [
                    legacy_jacobian(
                        np.stack([signals[name][i] for name in selected_labels]),
                        power.k,
                        mu,
                        widths=widths,
                    )
                    for i, mu in enumerate(power.mu)
                ]
            )
            captured.append(
                dict(
                    bin=index,
                    mean_z=float(z_mean),
                    covariance_z=float(pair_inputs[required_labels[0]]["_z_mean"]),
                    k=k,
                    mu=mu_flat,
                    modes=np.tile(covariance.num_modes, len(power.mu)),
                    total=np.column_stack(
                        [totals[label].reshape(-1) for label in required_labels]
                    ),
                    observed_j=observed_j,
                    required_pairs=required_pairs,
                    selected_pairs=selected_pairs,
                    pair_inputs=pair_inputs,
                )
            )
    if [row["bin"] for row in captured] != list(range(6)):
        raise ValueError("lyaforecast did not prepare all six DESI Run-2 bins")
    return captured, field_ids, original, identity, source


def _fixed_total(row, required_pairs, fields):
    """Replace only forest-auto noise using converged early-lyaforecast weights."""
    if [tuple(pair) for pair in row["required_pairs"]] != list(required_pairs):
        raise ValueError("captured required-pair order differs from selected forecast")
    total = row["total"].copy()
    states = {}
    for column, (i, j) in enumerate(required_pairs):
        field = fields[i]
        if i != j or field.kind != "forest":
            continue
        name = f"{field.id}_{field.id}"
        context = row["pair_inputs"][name]
        inputs = WeightInputs.from_pair(context)
        if len(inputs.magnitudes) != 107:
            raise ValueError("fixed compatibility requires the original 107-node grid")
        state = solve(inputs, METHODS["early_lyaforecast"], **STOPPING)
        states[field.id] = state
        if state["status"] != "converged":
            raise ValueError(f"{name}: {state['status']}; no qualified forecast")
        baseline = forest_noise(
            context,
            row["k"],
            row["mu"],
            context["_aliasing_weights"][-1],
            context["_effective_noise_power"][-1],
        )
        changed = forest_noise(
            context, row["k"], row["mu"], *state["coefficients"][-2:]
        )
        total[:, column] = total[:, column] - baseline + changed
    return total, states


def run_desi2_compatibility(*, reference, profile, work_directory):
    """Return six selected-bin compatibility forecasts and provenance.

    Parameters
    ----------
    reference
        Installed lyaforecast checkout containing ``examples/desi2`` and its
        package resources.  The imported package must resolve to this checkout.
    profile
        ``"full-compatibility"`` or ``"fixed-compatibility"``.
    work_directory
        New directory used for the effective INI and lyaforecast log.
    """
    if profile not in PROFILES:
        raise ValueError(f"profile must be one of {PROFILES}")
    captured, field_ids, original, identity, source = _capture(
        reference, work_directory
    )
    records = []
    convergence = {}
    for row in captured:
        index = row["bin"]
        selected = forecast_selection(CASE, index)
        required = [tuple(pair) for pair in selected.required_pairs.tolist()]
        selected_pairs = [tuple(pair) for pair in selected.selected_pairs.tolist()]
        if required != row["required_pairs"] or selected_pairs != row["selected_pairs"]:
            raise ValueError("prepared and requested selection differ")
        if profile == "fixed-compatibility":
            total, states = _fixed_total(row, required, original.fields)
            convergence[str(index + 1)] = {
                name: {
                    "status": state["status"],
                    "updates": state["updates"],
                    "candidate": state["candidate"],
                    "A_deg2": float(state["coefficients"][-2]),
                    "P_pixel_deg2_km_s": float(state["coefficients"][-1]),
                }
                for name, state in states.items()
            }
        else:
            total = row["total"]
        covariance = wick(
            total, row["modes"], required, selected_pairs, len(original.fields)
        )
        fisher, individual = contract(covariance, row["observed_j"])
        result = summaries(fisher, individual)
        lo, hi = bins(CASE)[index]
        for pair_index, pair in enumerate(selected_pairs):
            available = bool(np.all(result["pair_constrained"][pair_index]))
            records.append(
                dict(
                    bin_index=index,
                    bounds=(lo, hi),
                    redshift=row["mean_z"],
                    kind="individual",
                    pair=tuple(field_ids[i] for i in pair),
                    pair_indices=pair,
                    status="available" if available else "unavailable_rank_deficient",
                    fisher=result["pair_fisher"][pair_index],
                    sigma_ap=(
                        float(result["pair_errors"][pair_index, 0])
                        if available
                        else None
                    ),
                    sigma_at=(
                        float(result["pair_errors"][pair_index, 1])
                        if available
                        else None
                    ),
                    correlation=(
                        float(result["pair_correlation"][pair_index, 0, 1])
                        if available
                        else None
                    ),
                )
            )
        available = bool(np.all(result["constrained"]))
        records.append(
            dict(
                bin_index=index,
                bounds=(lo, hi),
                redshift=row["mean_z"],
                kind="joint",
                pair=None,
                pair_indices=None,
                status="available" if available else "unavailable_rank_deficient",
                fisher=result["fisher"],
                sigma_ap=float(result["errors"][0]) if available else None,
                sigma_at=float(result["errors"][1]) if available else None,
                correlation=float(result["correlation"][0, 1]) if available else None,
            )
        )
    settings = {
        "case": CASE,
        "profile": profile,
        "recipe_revision": REVISION,
        "parameter_order": ["ap", "at"],
        "selection": "bin1-qso-only; bins2-6-all",
        "bin1_excluded_pairs": [
            [field_ids[i], field_ids[j]]
            for i, j in selection(CASE).selected_pairs
            if [i, j] not in forecast_selection(CASE, 0).selected_pairs.tolist()
        ],
        "individual_result_count": sum(r["kind"] == "individual" for r in records),
        "joint_result_count": sum(r["kind"] == "joint" for r in records),
        "bin_redshifts": [
            {
                "bin_index": row["bin"],
                "mean_z": row["mean_z"],
                "covariance_z": row["covariance_z"],
            }
            for row in captured
        ],
        "input": {
            "reference": str(Path(reference).resolve()),
            "configuration": str(source.resolve()),
            "identity": identity,
        },
        "forest_weighting": {
            "method": "legacy"
            if profile == "full-compatibility"
            else "early_lyaforecast",
            "magnitude_measure": "original rectangular 107-node full endpoints",
            "input_policy": "literal signed compatibility",
            "iterations": 3 if profile == "full-compatibility" else None,
            "stopping": STOPPING if profile == "fixed-compatibility" else None,
            "representative_mode": REFERENCE,
            "convergence": convergence if profile == "fixed-compatibility" else None,
        },
        "covariance": "selected FishHighz Wick covariance; joint Fisher is not a sum of individual Fishers",
    }
    return settings, records
