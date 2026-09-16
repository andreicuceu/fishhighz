"""Offline supplied-array peak/operator swaps using saved fixed model settings.

This does not solve CAMB or rerun a survey forecast. It evaluates the recorded
intrinsic model on common legacy nodes to distinguish peak extraction from the
full derivative operator. The separate policy/physical one-setting diagnostics
are in the primary bundle. Array replacement chains are order dependent.
"""

import argparse
import json
from itertools import islice
from pathlib import Path

import numpy as np

from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS
from fishhighz.response import pixel_width_angstrom_to_velocity, velocity_response
from fishhighz.validation.attribution import interpolate
from fishhighz.validation.cases import recipe, selection
from fishhighz.validation.evidence import digest
from fishhighz.validation.numerics import (
    change,
    contract,
    legacy_peak,
    relative,
    summaries,
    wick,
)
from fishhighz.validation.plots import records, values


def legacy_noise(task, inputs, k, mu, adapters):
    """Literal saved-coefficient noise; galaxy interpolation keeps legacy signs.

    No density normalization, floor extension, or overlap prescription is added.
    Arrays here are validation-only upstream inputs, not strict reader defaults.
    """
    from fishhighz.models.p1d import default_p1d

    sel = selection(task["case"])
    config = recipe(task["case"])
    mags = np.linspace(
        float(config["survey"]["min_band_mag"]),
        float(config["survey"]["max_band_mag"]),
        int(config["survey"]["num mag bins"]),
    )
    result = np.zeros((len(k), len(sel.required_pairs)))
    for column, (i, j) in enumerate(sel.required_pairs):
        if i != j:
            continue
        field = sel.fields[i]
        row = inputs[field.id + "_" + field.id]
        av = row["_distance_to_velocity"]
        ddeg = row["_angle_to_distance"]
        z = row["_z_mean"]
        if field.kind == "forest":
            q = k * mu / av
            window = velocity_response(
                q,
                pixel_width_velocity=row["_pix_kms"],
                gaussian_sigma_velocity=row["_res_kms"],
            )
            result[:, column] = (
                (
                    row["_aliasing_weights"][-1] * default_p1d([], z, q) * window**2
                    + row["_effective_noise_power"][-1]
                )
                * ddeg**2
                / av
            )
        else:
            sample = adapters[field.id].sample(z, mags)
            raw = np.where(
                sample["provenance"]["masks"]["density_floor"], 1e-20, sample["raw"]
            )
            density = np.cumsum(raw * (1 + z) / 299800 * (mags[1] - mags[0]))[-1]
            result[:, column] = ddeg**2 / (av * density)
    return result


def negative_floor_noise(task, inputs, k, mu, adapters, original):
    """Isolate the negative-density extension in the literal three-update recipe.

    Formulas adapted from lyaforecast weights.py/covariance.py (GPLv3). This
    validation-only supplied-array diagnostic preserves signed legacy arithmetic;
    it never changes strict readers or the accepted positive-density engine.
    """
    from fishhighz.models.p1d import default_p1d

    sel = selection(task["case"])
    config = recipe(task["case"])
    result = original.copy()
    checks = {}
    for column, (i, j) in enumerate(sel.required_pairs):
        if i != j:
            continue
        field = sel.fields[i]
        row = inputs[field.id + "_" + field.id]
        av = row["_distance_to_velocity"]
        ddeg = row["_angle_to_distance"]
        z = row["_z_mean"]
        if field.kind == "forest":
            density = np.asarray(row["density"])
            var = np.asarray(row["variance"])
            m = np.asarray(row["magnitudes"])
            dm = m[1] - m[0]
            pixel = row["_pix_kms"]
            length = row["forest_length"]
            signal = row["auxiliary_signal"]
            alias = row["auxiliary_p1d"]

            def coefficients(rho):
                # Keep literal operation order, including signed/zero cumulative sums.
                with np.errstate(all="ignore"):
                    v1 = alias / pixel
                    w = v1 / (v1 + var)
                    for _ in range(3):
                        neff = np.cumsum(rho * w * dm) * (length / pixel)
                        w = signal / (signal + var / neff)
                    i1 = np.cumsum(rho * w * dm)[-1]
                    i2 = np.cumsum(rho * w**2 * dm)[-1]
                    i3 = np.cumsum(rho * w**2 * var * dm)[-1]
                    a = i2 / (i1**2 * length)
                    p = i3 * pixel / (i1**2 * length)
                if not np.all(np.isfinite(w)) or not np.isfinite(a + p):
                    raise ValueError("nonfinite literal signed-density diagnostic")
                return w, a, p

            w, a, p = coefficients(density)
            error = max(
                relative(w, row["_w_lya"]),
                relative(
                    [a, p],
                    [row["_aliasing_weights"][-1], row["_effective_noise_power"][-1]],
                ),
            )
            if error > 5e-12:
                raise ValueError(f"literal weight reconstruction {field.id}: {error}")
            negative = density < 0
            checks[field.id] = dict(
                negative_nodes=int(np.count_nonzero(negative)),
                literal_weight_relative=error,
            )
            if np.any(negative):
                _, a, p = coefficients(
                    np.where(negative, 1e-20 * (1 + row["_zq"]) / 299800, density)
                )
                q = k * mu / av
                window = velocity_response(
                    q,
                    pixel_width_velocity=pixel,
                    gaussian_sigma_velocity=row["_res_kms"],
                )
                result[:, column] = (
                    (a * default_p1d([], z, q) * window**2 + p) * ddeg**2 / av
                )
        else:
            mags = np.linspace(
                float(config["survey"]["min_band_mag"]),
                float(config["survey"]["max_band_mag"]),
                int(config["survey"]["num mag bins"]),
            )
            sample = adapters[field.id].sample(z, mags)
            raw = np.where(
                sample["provenance"]["masks"]["density_floor"], 1e-20, sample["raw"]
            )
            checks[field.id] = dict(negative_nodes=int(np.count_nonzero(raw < 0)))
            if np.any(raw < 0):
                density = np.cumsum(
                    np.where(raw < 0, 1e-20, raw)
                    * (1 + z)
                    / 299800
                    * (mags[1] - mags[0])
                )[-1]
                result[:, column] = ddeg**2 / (av * density)
    return result, checks


def density_adapters(task, reference_origin):
    from fishhighz.adapters.legacy_compat import LegacyDensity
    from fishhighz.adapters.legacy_inputs import DensityReader

    config = recipe(task["case"])
    sel = selection(task["case"])
    result = {}
    for field, key in zip(sel.fields, [s for s in config if s.startswith("tracer ")]):
        if field.kind == "forest":
            continue
        t = config[key]
        reader = DensityReader(
            Path(reference_origin) / "resources/data" / t["dn dz"],
            semantics="cell_count_per_deg2",
            target_density=float(t["target density"]),
            z_norm_min=2.15 if field.id == "qso" else None,
            width_policy="legacy_first_spacing",
            label=field.id,
        )
        result[field.id] = LegacyDensity(reader, "floor_negative")
    return result


def reconstruct_weights(accuracy, settings, task, model):
    """Recover final arrays from saved samples; independently check recorded noise."""
    from fishhighz.kernels.weights import _integrals, _iterate
    from fishhighz.models.p1d import default_p1d
    from fishhighz.weights import density_per_velocity

    config = recipe(task["case"])
    sel = selection(task["case"])
    outputs = {}
    checks = {}
    av = settings["geometry"]["a_v"]
    ddeg = settings["geometry"]["d_deg"]
    z = settings["model"]["z_eval"]
    tracers = [config[s] for s in config if s.startswith("tracer ")]
    for i, field in enumerate(sel.fields):
        if field.kind != "forest" or field.id not in settings["samples"]:
            continue
        row = settings["samples"][field.id]
        rho = density_per_velocity(row["density"], z_source=row["z_source"])
        masses = rho * np.asarray(row["quadrature"])
        variance = np.asarray(row["variance"])
        pixel = pixel_width_angstrom_to_velocity(
            float(tracers[i]["pix_width_ang"]),
            lambda_obs_angstrom=LYA_REST_ANGSTROM * (1 + z),
        )
        sigma = SPEED_LIGHT_KMS / (
            float(config["survey"]["resolution"]) * 2 * np.sqrt(2 * np.log(2))
        )
        q = 0.00035
        kp = av * q
        k = np.hypot(kp, 2.4 / ddeg)
        mu = kp / k
        waux = velocity_response(
            [q], pixel_width_velocity=pixel, gaussian_sigma_velocity=sigma
        )[0]
        signal = float(
            model([1.0, 1.0], z, [k], [mu], [(i, i)])[0, 0] * waux**2 * av / ddeg**2
        )
        alias = float(default_p1d([], z, [q])[0] * waux**2)
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            weight, changes = _iterate(
                masses,
                variance,
                row["length_velocity"],
                pixel,
                signal,
                alias,
                settings["controls"]["iterations"],
            )
            i1, i2, i3, a, pixel_power = _integrals(
                masses, weight, variance, row["length_velocity"], pixel
            )
        qgrid = accuracy["k"] * accuracy["mu"] / av
        window = velocity_response(
            qgrid, pixel_width_velocity=pixel, gaussian_sigma_velocity=sigma
        )
        predicted = (
            (a * default_p1d([], z, qgrid) * window**2 + pixel_power) * ddeg**2 / av
        )
        column = next(
            c for c, pair in enumerate(sel.required_pairs) if tuple(pair) == (i, i)
        )
        discrepancy = relative(predicted, accuracy["noise"][:, column])
        if discrepancy > 5e-12:
            raise ValueError(f"{field.id}: weight/noise reconstruction {discrepancy}")
        prefix = f"field_{i}_"
        outputs.update(
            {
                prefix + name: value
                for name, value in dict(
                    weights=weight,
                    changes=changes,
                    rho=rho,
                    masses=masses,
                    I1=i1,
                    I2=i2,
                    I3=i3,
                    coefficients=np.array([a, pixel_power, signal, alias]),
                ).items()
            }
        )
        checks[field.id] = dict(
            noise_relative=discrepancy,
            A=float(a),
            P_pixel=float(pixel_power),
            auxiliary_signal=signal,
            auxiliary_alias=alias,
            interpretation="Reconstructed from saved final samples and fixed model; original volume not recomputed",
        )
    return outputs, checks


def attribution(legacy, accuracy, report, task, template, adapters=None):
    from fishhighz.derivatives import evaluate_derivatives
    from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
    from fishhighz.models.kaiser import KaiserModel, Scaling
    from fishhighz.parameters import Parameter, ParameterRegistry

    settings = report["settings"]
    model_settings = settings["model"]
    sel = selection(task["case"])
    k, mu = legacy["k"], legacy["mu"]
    ka = np.unique(k)
    ma = np.unique(mu)
    model = KaiserModel(
        template,
        sel.fields,
        biases=model_settings["biases"],
        betas=model_settings["betas"],
        widths=model_settings["widths"],
        f=model_settings["f"] if any(f.kind == "galaxy" for f in sel.fields) else None,
        local_names=("ap", "at"),
        wiggle=Scaling("ap_at", ap="ap", at="at"),
        z=model_settings["z_eval"],
        growth=model_settings["G"],
    )
    registry = ParameterRegistry(
        [
            Parameter(n, 1.0, "target", step=report["final_controls"]["step"])
            for n in ("ap", "at")
        ]
    )
    prepared = PreparedP3D(
        registry,
        sel,
        [
            P3DProvider(
                "recorded accuracy",
                model,
                BoundParameters(registry, ("ap", "at"), {"ap": "ap", "at": "at"}),
                sel.required_pairs,
            )
        ],
    )
    config = recipe(task["case"])
    tracers = [config[s] for s in config if s.startswith("tracer ")]
    transfer = np.ones((len(k), len(sel.fields)))
    for i, field in enumerate(sel.fields):
        if field.kind == "forest":
            pixel = pixel_width_angstrom_to_velocity(
                float(tracers[i]["pix_width_ang"]),
                lambda_obs_angstrom=LYA_REST_ANGSTROM * (1 + model_settings["z_eval"]),
            )
            sigma = SPEED_LIGHT_KMS / (
                float(config["survey"]["resolution"]) * 2 * np.sqrt(2 * np.log(2))
            )
            transfer[:, i] = velocity_response(
                k * mu / settings["geometry"]["a_v"],
                pixel_width_velocity=pixel,
                gaussian_sigma_velocity=sigma,
            )
    i, j = sel.required_pairs.T
    products = transfer[:, i] * transfer[:, j]
    derivatives = evaluate_derivatives(
        prepared, registry.fiducials, model_settings["z_eval"], k, mu
    )
    full_j = (products[:, :, None] * derivatives.jacobian)[:, sel.selected_to_required]
    components = template.evaluate(k)
    factors = np.column_stack(
        [
            model_settings["biases"][f.id] * (1 + model_settings["betas"][f.id] * mu**2)
            if f.kind == "forest"
            else model_settings["biases"][f.id] + model_settings["f"] * mu**2
            for f in sel.fields
        ]
    )
    rsd = factors[:, i] * factors[:, j] * model_settings["G"] * products
    widths = np.array([model_settings["widths"][f.id] for f in sel.fields])
    squared = (widths[i] ** 2 + widths[j] ** 2) / 2
    damping = np.exp(
        -0.5
        * k[:, None] ** 2
        * (
            mu[:, None] ** 2 * squared[None, :, 0]
            + (1 - mu[:, None] ** 2) * squared[None, :, 1]
        )
    )
    # template.evaluate columns are smooth, signed wiggle.
    wiggle = rsd * components[:, 1, None] * damping
    undamped = rsd * components.sum(axis=1)[:, None]
    peak = (
        np.concatenate(
            [legacy_peak(row.T, ka).T for row in undamped.reshape(len(ma), len(ka), -1)]
        )
        * damping
    )

    def backward(signal):
        grid = signal.reshape(len(ma), len(ka), -1)
        d = np.zeros_like(grid)
        d[:, 1:] = np.diff(grid, axis=1) / ((ka[1] - ka[0]) / ka[1:])[None, :, None]
        return (
            d.reshape(len(k), -1)[:, sel.selected_to_required, None]
            * np.column_stack((mu**2, 1 - mu**2))[:, None, :]
        )

    jw = backward(wiggle)
    jp = backward(peak)
    # GL uses mu-fastest; convert explicitly to the interpolation helper's mu-major order.
    ak = np.unique(accuracy["k"])
    am = np.unique(accuracy["mu"])

    def aligned(name):
        a = accuracy[name]
        return a.reshape(len(ak), len(am), *a.shape[1:]).swapaxes(0, 1).reshape(a.shape)

    total = interpolate(aligned("total"), ak, am, k, mu)
    accurate_noise = interpolate(aligned("noise"), ak, am, k, mu)
    if adapters is None:
        adapters = density_adapters(
            task, report["provenance"]["reference"]["reference_origin"]
        )
    old_noise = legacy_noise(task, report["legacy_pair_inputs"], k, mu, adapters)
    measured_j = interpolate(aligned("observed_j"), ak, am, k, mu)
    modes = (
        legacy["modes"] * settings["grid"]["volume"] / float(report["legacy_volume"])
    )
    model_checks = {}
    if task["kind"] == "real_bao":
        native = model(
            [1.0, 1.0],
            model_settings["z_eval"],
            accuracy["k"],
            accuracy["mu"],
            sel.required_pairs,
        )
        model_checks["intrinsic_relative"] = relative(native, accuracy["intrinsic"])
        # Reconstruct the fixed response and observed derivative on the actual GL nodes.
        observed = np.ones((len(accuracy["k"]), len(sel.fields)))
        for index, field in enumerate(sel.fields):
            if field.kind == "forest":
                pixel = pixel_width_angstrom_to_velocity(
                    float(tracers[index]["pix_width_ang"]),
                    lambda_obs_angstrom=LYA_REST_ANGSTROM
                    * (1 + model_settings["z_eval"]),
                )
                sigma = SPEED_LIGHT_KMS / (
                    float(config["survey"]["resolution"]) * 2 * np.sqrt(2 * np.log(2))
                )
                observed[:, index] = velocity_response(
                    accuracy["k"] * accuracy["mu"] / settings["geometry"]["a_v"],
                    pixel_width_velocity=pixel,
                    gaussian_sigma_velocity=sigma,
                )
        model_checks["response_relative"] = relative(observed, accuracy["response"])
        prod = observed[:, i] * observed[:, j]
        jac = (
            prod[:, :, None]
            * evaluate_derivatives(
                prepared,
                registry.fiducials,
                model_settings["z_eval"],
                accuracy["k"],
                accuracy["mu"],
            ).jacobian
        )[:, sel.selected_to_required]
        model_checks["observed_j_relative"] = relative(jac, accuracy["observed_j"])
        model_checks["total_relative"] = relative(
            prod * native + accuracy["noise"], accuracy["total"]
        )
        if any(value > 5e-12 for value in model_checks.values()):
            raise ValueError(f"recorded model/response/J mismatch: {model_checks}")
    outputs, weight_checks = reconstruct_weights(accuracy, settings, task, model)
    outputs.update(legacy_noise=old_noise, accuracy_noise=accurate_noise)
    stages = []

    def save(name, t, jac, n, held):
        c = wick(t, n, sel.required_pairs, sel.selected_pairs, len(sel.fields))
        f, p = contract(c, jac, independent=True)
        ix = len(stages)
        outputs.update(
            {
                f"total_{ix}": t,
                f"jacobian_{ix}": jac,
                f"modes_{ix}": n,
                **{f"{key}_{ix}": a for key, a in summaries(f, p).items()},
            }
        )
        stages.append(
            dict(
                name=name,
                held_fixed=held,
                values=values(f),
                fisher_relative_previous=relative(f, outputs[f"fisher_{ix - 1}"])
                if ix
                else 0.0,
            )
        )

    save(
        "legacy",
        legacy["total"],
        legacy["observed_j"],
        legacy["modes"],
        "captured compatibility endpoint",
    )
    save("volume", legacy["total"], legacy["observed_j"], modes, "legacy T and J")
    save(
        "per-field accuracy noise",
        legacy["total"] - old_noise + accurate_noise,
        legacy["observed_j"],
        modes,
        "legacy observed signal and J; integrated modes; replace only diagonal aliasing/pixel/Poisson noise",
    )
    save(
        "accuracy observed signal",
        total,
        legacy["observed_j"],
        modes,
        "legacy J; integrated volume and accuracy noise; replace only observed signal, including interpolation/model/response changes",
    )
    save(
        "template model with polynomial peak",
        total,
        jp,
        modes,
        "accuracy T and integrated modes; legacy peak/operator applied to undamped accuracy mean with accepted widths",
    )
    save(
        "signed PK-PKSB peak",
        total,
        jw,
        modes,
        "same T, modes, RSD, response, widths and backward operator; only peak extraction changes",
    )
    save(
        "full wiggle mapping derivative",
        total,
        full_j,
        modes,
        "same T and modes; replace backward operator including its treatment of response with accepted full mapping at fixed response",
    )
    save(
        "final quadrature",
        accuracy["total"],
        accuracy["observed_j"],
        accuracy["modes"],
        "accuracy endpoint; includes common-node interpolation residual",
    )
    output = dict(
        context=task,
        stages=stages,
        weight_reconstruction=weight_checks,
        model_reconstruction=model_checks,
        endpoint_relative=relative(outputs["fisher_7"], accuracy["fisher"]),
        common_node_j_interpolation_relative=relative(measured_j, full_j),
        interpretation="Sequential supplied-array effects are order dependent and not independent additive physical errors. Model/redshift/noise groups may remain unresolved; see separate one-setting diagnostics.",
    )
    # A separate branch isolates only the legacy forest/galaxy cross-width rule.
    # Autos, template, growth, response, noise, nodes and volume remain fixed.
    old_squared = np.array(
        [
            widths[ii] ** 2
            if sel.fields[ii].kind == "forest"
            else widths[jj] ** 2
            if sel.fields[jj].kind == "forest"
            else widths[ii] ** 2
            for ii, jj in sel.required_pairs
        ]
    )
    delta = old_squared - squared
    smooth = (
        components[:, 0, None] * factors[:, i] * factors[:, j] * model_settings["G"]
    )
    lookup = {tuple(pair): index for index, pair in enumerate(sel.required_pairs)}

    def cross_model(theta, z, knodes, munodes, pairs):
        if not np.array_equal(knodes, k) or not np.array_equal(munodes, mu):
            raise ValueError("cross-width diagnostic uses fixed recorded nodes")
        columns = [lookup[tuple(pair)] for pair in pairs]
        ratio = np.exp(
            -0.5
            * (
                (k * mu / theta[0])[:, None] ** 2 * delta[None, columns, 0]
                + (k * np.sqrt(1 - mu**2) / theta[1])[:, None] ** 2
                * delta[None, columns, 1]
            )
        )
        base = model(theta, z, knodes, munodes, pairs)
        return smooth[:, columns] + (base - smooth[:, columns]) * ratio

    try:
        cp = PreparedP3D(
            registry,
            sel,
            [
                P3DProvider(
                    "legacy cross-width rule",
                    cross_model,
                    BoundParameters(registry, ("ap", "at"), {"ap": "ap", "at": "at"}),
                    sel.required_pairs,
                )
            ],
        )
        cj = (
            products[:, :, None]
            * evaluate_derivatives(
                cp, registry.fiducials, model_settings["z_eval"], k, mu
            ).jacobian
        )[:, sel.selected_to_required]
        ct = total + products * (
            cross_model([1.0, 1.0], model_settings["z_eval"], k, mu, sel.required_pairs)
            - model([1.0, 1.0], model_settings["z_eval"], k, mu, sel.required_pairs)
        )
        cf, cs = contract(
            wick(ct, modes, sel.required_pairs, sel.selected_pairs, len(sel.fields)),
            cj,
            independent=True,
        )
        outputs.update(
            cross_rule_total=ct,
            cross_rule_jacobian=cj,
            cross_rule_modes=modes,
            **{
                "cross_rule_" + name: value for name, value in summaries(cf, cs).items()
            },
        )
        output["cross_rule"] = dict(
            values=values(cf),
            reference_values=values(outputs["fisher_6"]),
            change=change(cf, outputs["fisher_6"], cs, outputs["pair_fisher_6"], 1, 1),
            held_fixed="Only forest-galaxy cross damping changes to the forest auto width; accepted autos, template, growth, response, noise, nodes and volume held fixed",
        )
    except (ValueError, np.linalg.LinAlgError) as error:
        output["cross_rule"] = dict(error=str(error))
    try:
        floor_noise, floor_checks = negative_floor_noise(
            task, report["legacy_pair_inputs"], k, mu, adapters, old_noise
        )
        ft = legacy["total"] - old_noise + floor_noise
        ff, fp = contract(
            wick(
                ft,
                legacy["modes"],
                sel.required_pairs,
                sel.selected_pairs,
                len(sel.fields),
            ),
            legacy["observed_j"],
            independent=True,
        )
        outputs.update(
            negative_floor_total=ft,
            negative_floor_jacobian=legacy["observed_j"],
            negative_floor_modes=legacy["modes"],
            **{
                "negative_floor_" + name: value
                for name, value in summaries(ff, fp).items()
            },
        )
        output["negative_density_floor"] = dict(
            values=values(ff),
            reference_values=values(outputs["fisher_0"]),
            change=change(ff, outputs["fisher_0"], fp, outputs["pair_fisher_0"], 1, 1),
            checks=floor_checks,
            held_fixed="Only negative density values become 1e-20 before legacy velocity conversion; literal three updates, auxiliary signal/P1D, SNR, response, mean/J and modes held fixed. Validation-only signed-array diagnostic, not a physical alternate survey.",
        )
    except (ValueError, np.linalg.LinAlgError) as error:
        output["negative_density_floor"] = dict(error=str(error))
    return outputs, output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle")
    p.add_argument("--template")
    p.add_argument("--output", required=True)
    p.add_argument(
        "--replot", help="saved attribution directory; no template or model evaluation"
    )
    args = p.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    if args.replot:
        path = Path(args.replot) / "manifest.json"
        manifest = json.loads(path.read_text())
        verify_saved(path.parent, manifest)
        plot_saved(manifest, out)
        (out / "data-source.json").write_text(
            json.dumps(dict(path=str(path.resolve()), sha256=digest(path)), indent=2)
            + "\n"
        )
        return
    if not args.bundle or not args.template:
        p.error("calculation requires --bundle and --template")
    template = None
    adapters = {}
    pending = {}
    results = []
    primary_count = len(
        json.loads((Path(args.bundle) / "manifest.json").read_text())["records"]
    )
    # Attribution uses primary endpoints; the full plot reader checks diagnostics.
    for record, a in islice(records(args.bundle), primary_count):
        task = record["task"]
        if task["kind"] == "diagnostic":
            continue
        key = (task["case"], task["bin"])
        if task["profile"] == "compatibility":
            if a is not None:
                pending[key] = (
                    a,
                    record["report"]["settings"]["grid"]["volume"],
                    record["report"]["settings"]["pair_inputs"],
                )
            continue
        if a is None or key not in pending:
            results.append(dict(context=task, error="missing primary endpoint"))
            continue
        legacy, volume, pair_inputs = pending.pop(key)
        try:
            report = dict(record["report"])
            report["legacy_volume"] = volume
            report["legacy_pair_inputs"] = pair_inputs
            if task["case"] not in adapters:
                adapters[task["case"]] = density_adapters(
                    task, report["provenance"]["reference"]["reference_origin"]
                )
            if template is None:
                from fishhighz.models.templates import load_template

                template = load_template(
                    args.template, h_fid=report["settings"]["grid"]["h_fid"]
                )
            arrays, row = attribution(
                legacy, a, report, task, template, adapters[task["case"]]
            )
            name = f"{task['case']}-{task['bin']}.npz"
            np.savez_compressed(out / name, **arrays)
            row.update(array=name, sha256=digest(out / name))
            results.append(row)
        except (ValueError, np.linalg.LinAlgError) as error:
            results.append(dict(context=task, error=str(error)))
    manifest = dict(
        source_bundle=str(Path(args.bundle).resolve()),
        source_sha256=digest(Path(args.bundle) / "manifest.json"),
        script_sha256=digest(__file__),
        template_sha256=digest(args.template),
        records=results,
    )
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    plot_saved(manifest, out)


def verify_saved(root, manifest):
    """Independent matrix check before rendering saved attribution data."""
    from fishhighz.validation.schema import validate_request

    root = Path(root).resolve()
    for row in manifest["records"]:
        task = row["context"]
        validate_request(task)
        if "array" not in row:
            continue
        path = (root / row["array"]).resolve()
        if path.parent != root or digest(path) != row["sha256"]:
            raise ValueError("attribution data hash/path")
        with np.load(path, allow_pickle=False) as f:
            arrays = {n: f[n] for n in f.files}
        if any(
            a.dtype.kind not in "fiu" or not np.all(np.isfinite(a))
            for a in arrays.values()
        ):
            raise ValueError("nonfinite/nonreal attribution arrays")
        if [stage["name"] for stage in row["stages"]] != [
            "legacy",
            "volume",
            "per-field accuracy noise",
            "accuracy observed signal",
            "template model with polynomial peak",
            "signed PK-PKSB peak",
            "full wiggle mapping derivative",
            "final quadrature",
        ]:
            raise ValueError("incorrect attribution stage inventory")
        for n in range(len(row["stages"])):
            c = wick(
                arrays[f"total_{n}"],
                arrays[f"modes_{n}"],
                task["required_pairs"],
                task["selected_pairs"],
                len(task["fields"]),
            )
            fisher, pairs = contract(c, arrays[f"jacobian_{n}"], independent=True)
            if (
                relative(fisher, arrays[f"fisher_{n}"]) > 5e-12
                or relative(pairs, arrays[f"pair_fisher_{n}"]) > 5e-12
            ):
                raise ValueError("attribution matrix mismatch")
            for name, value in summaries(fisher, pairs).items():
                if (
                    arrays[name + f"_{n}"].shape != value.shape
                    or relative(arrays[name + f"_{n}"], value) > 5e-12
                ):
                    raise ValueError("attribution derived quantities mismatch")
            np.testing.assert_allclose(
                values(fisher), row["stages"][n]["values"], rtol=5e-12, atol=0
            )
            previous = relative(fisher, arrays[f"fisher_{n - 1}"]) if n else 0.0
            if not np.isclose(
                row["stages"][n]["fisher_relative_previous"],
                previous,
                rtol=5e-12,
                atol=1e-15,
            ):
                raise ValueError("attribution stage difference mismatch")
        for prefix, name in [
            ("cross_rule", "cross_rule"),
            ("negative_floor", "negative_density_floor"),
        ]:
            if "values" not in row.get(name, {}):
                continue
            c = wick(
                arrays[prefix + "_total"],
                arrays[prefix + "_modes"],
                task["required_pairs"],
                task["selected_pairs"],
                len(task["fields"]),
            )
            fisher, pairs = contract(c, arrays[prefix + "_jacobian"], independent=True)
            if (
                relative(fisher, arrays[prefix + "_fisher"]) > 5e-12
                or relative(pairs, arrays[prefix + "_pair_fisher"]) > 5e-12
            ):
                raise ValueError("attribution branch matrix mismatch")
            for quantity, value in summaries(fisher, pairs).items():
                if (
                    arrays[prefix + "_" + quantity].shape != value.shape
                    or relative(arrays[prefix + "_" + quantity], value) > 5e-12
                ):
                    raise ValueError("attribution branch derived quantities mismatch")
            np.testing.assert_allclose(
                values(fisher), row[name]["values"], rtol=5e-12, atol=0
            )
            reference = 6 if prefix == "cross_rule" else 0
            np.testing.assert_allclose(
                values(arrays[f"fisher_{reference}"]),
                row[name]["reference_values"],
                rtol=5e-12,
                atol=0,
            )
            expected_change = change(
                fisher,
                arrays[f"fisher_{reference}"],
                pairs,
                arrays[f"pair_fisher_{reference}"],
                1,
                1,
            )
            if set(row[name]["change"]) != set(expected_change) or any(
                not np.isclose(row[name]["change"][key], value, rtol=5e-12, atol=1e-15)
                for key, value in expected_change.items()
            ):
                raise ValueError("attribution branch difference mismatch")


def plot_saved(manifest, out):
    """Plot saved numbers without loading a template or evaluating a model."""
    results = manifest["records"]
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = {}
    for case in sorted({r["context"]["case"] for r in results}):
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        axes = axes.ravel()
        for row in results:
            if row["context"]["case"] != case or "stages" not in row:
                continue
            v = np.array([s["values"] for s in row["stages"]])
            for j, ax in enumerate(axes[:2]):
                ax.plot(
                    100 * (v[:, j] / v[0, j] - 1),
                    marker="o",
                    ms=3,
                    label="bin " + str(row["context"]["bin"]),
                )
        for j, ax in enumerate(axes[:2]):
            ax.set_xticks(
                range(8),
                [
                    "legacy",
                    "volume",
                    "noise",
                    "signal",
                    "model + poly",
                    "PK-PKSB",
                    "full J",
                    "final grid",
                ],
                rotation=40,
            )
            ax.set_ylabel(("ap" if j == 0 else "at") + " change from legacy (%)")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=7)
        cross_rows = [
            r
            for r in results
            if r["context"]["case"] == case and "values" in r.get("cross_rule", {})
        ]
        for j, label in enumerate(("ap", "at")):
            axes[2].plot(
                [r["context"]["bin"] for r in cross_rows],
                [
                    100
                    * (
                        r["cross_rule"]["values"][j]
                        / r["cross_rule"]["reference_values"][j]
                        - 1
                    )
                    for r in cross_rows
                ],
                marker="o",
                ms=3,
                label=label,
            )
        axes[2].set_xlabel("bin")
        axes[2].set_ylabel("legacy cross-width rule / accepted (%)")
        axes[2].legend()
        axes[2].grid(alpha=0.2)
        floor_rows = [
            r
            for r in results
            if r["context"]["case"] == case
            and "values" in r.get("negative_density_floor", {})
        ]
        for j, label in enumerate(("ap", "at")):
            axes[3].plot(
                [r["context"]["bin"] for r in floor_rows],
                [
                    100
                    * (
                        r["negative_density_floor"]["values"][j]
                        / r["negative_density_floor"]["reference_values"][j]
                        - 1
                    )
                    for r in floor_rows
                ],
                marker="o",
                ms=3,
                label=label,
            )
        axes[3].set_xlabel("bin")
        axes[3].set_ylabel("negative floor / signed legacy (%)")
        axes[3].legend()
        axes[3].grid(alpha=0.2)
        fig.suptitle(
            case + " | controlled peak/operator swaps; order-dependent array chain"
        )
        unavailable = {
            name: [
                r["context"]["bin"]
                for r in results
                if r["context"]["case"] == case
                and (
                    ("stages" not in r)
                    if name == "chain"
                    else "values" not in r.get(name, {})
                )
            ]
            for name in ("chain", "cross_rule", "negative_density_floor")
        }
        fig.text(
            0.5,
            0.01,
            "Unavailable diagnostic bins: "
            + str({k: v for k, v in unavailable.items() if v}),
            ha="center",
            fontsize=8,
        )
        fig.tight_layout(rect=(0, 0.04, 1, 0.96))
        for suffix in ("png", "pdf"):
            path = out / f"{case}.{suffix}"
            fig.savefig(path, dpi=140)
            files[path.name] = digest(path)
        plt.close(fig)
    (out / "plots-manifest.json").write_text(json.dumps(files, indent=2) + "\n")


if __name__ == "__main__":
    main()
