"""15x2pt saved-sample weighting diagnosis and fixed-weight magnitude refinement.

No CAMB solve, new weighting prescription or changed arithmetic guard. The
fixed-weight diagnostic uses the original initial weight function, held fixed
as a function of magnitude while its integrals are refined.
"""

import argparse
import configparser
import json
from pathlib import Path

import numpy as np

from fishhighz.adapters.legacy_compat import LegacyDensity, LegacySNR, plain
from fishhighz.adapters.legacy_inputs import DensityReader, SNRReader
from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS
from fishhighz.kernels.weights import _integrals
from fishhighz.models.kaiser import KaiserModel, Scaling
from fishhighz.models.p1d import default_p1d
from fishhighz.models.templates import load_template
from fishhighz.response import pixel_width_angstrom_to_velocity, velocity_response
from fishhighz.validation.accuracy import AccuracyRecipe
from fishhighz.validation.cases import recipe, selection
from fishhighz.validation.evidence import digest, record_report
from fishhighz.validation.numerics import change, contract, relative, summaries, wick
from fishhighz.validation.weight_diagnosis import diagnose
from fishhighz.weights import density_per_velocity

CASE = "lya_qso_lbg_lae_15x2pt"


def sample_reader(root):
    """Prepare only the retained raw interpolation policies; no external P3D."""
    r = AccuracyRecipe.__new__(AccuracyRecipe)
    r.root = Path(root).resolve()
    r.case = CASE
    r.selection = selection(CASE)
    r.config = configparser.ConfigParser()
    r.config.read_dict(recipe(CASE))
    r.densities, r.snrs, r.tracers = {}, {}, {}
    for field, key in zip(
        r.selection.fields, [s for s in r.config.sections() if s.startswith("tracer ")]
    ):
        t = r.config[key]
        r.tracers[field.id] = dict(t)
        forest = field.kind == "forest"
        reader = DensityReader(
            r.root / "lyaforecast/resources/data" / t["dn dz"],
            semantics="cell_count_per_deg2",
            target_density=t.getfloat("target density"),
            z_norm_min=2.15 if forest or field.id == "qso" else None,
            magnitude_bounds=(t.getfloat("min_band_mag"), t.getfloat("max_band_mag"))
            if forest
            else None,
            width_policy="legacy_first_spacing",
            label=field.id,
        )
        r.densities[field.id] = LegacyDensity(reader, "floor_negative")
        if forest:
            paths = sorted(
                (r.root / "lyaforecast/resources/data" / t["snr-file-dir"]).glob(
                    "*.dat"
                )
            )
            r.snrs[field.id] = LegacySNR(
                SNRReader(paths, smoothing="legacy", label=field.id)
            )
    r._partitions, r._samples = {}, {}
    return r


def run(bundle, reference, template_path, output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    root = Path(bundle).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    records = [
        r
        for r in manifest["records"]
        if r["task"]["case"] == CASE and r["task"]["profile"] == "accuracy"
    ]
    if [r["task"]["bin"] for r in records] != list(range(6)):
        raise ValueError("diagnosis requires all six 15x2pt bins")
    reader = sample_reader(reference)
    sel = reader.selection
    template = None
    result = dict(
        source_sha256=digest(root / "manifest.json"),
        template_sha256=digest(template_path),
        interpretation="Fixed initial weight function isolates magnitude quadrature; coupled finite iterations are diagnostic and unresolved, not alternate primary surveys",
        bins=[],
    )
    for rec in records:
        index = rec["task"]["bin"]
        report = record_report(root, rec)
        assert digest(root / rec["arrays"]) == rec["sha256"]
        with np.load(root / rec["arrays"], allow_pickle=False) as data:
            a = {n: data[n] for n in data.files}
        s = report["settings"]
        if template is None:
            template = load_template(template_path, h_fid=s["grid"]["h_fid"])
        ms = s["model"]
        model = KaiserModel(
            template,
            sel.fields,
            biases=ms["biases"],
            betas=ms["betas"],
            widths=ms["widths"],
            f=ms["f"],
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
            z=ms["z_eval"],
            growth=ms["G"],
        )
        av, ddeg = s["geometry"]["a_v"], s["geometry"]["d_deg"]
        saved_arrays = {}
        rows = []
        for order in (4, 8, 16, 32, 64):
            sampled, m, q = reader.samples(index, order)
            if order == report["final_controls"]["magnitude_order"]:
                for field in sel.fields:
                    for key in ("density", "quadrature", "magnitudes"):
                        np.testing.assert_allclose(
                            sampled[field.id][key],
                            s["samples"][field.id][key],
                            rtol=5e-12,
                            atol=0,
                        )
            noise_by_count = {count: a["noise"].copy() for count in (0, 3, 6, 12, 24)}
            available = {count: True for count in noise_by_count}
            fields = []
            for fi, field in enumerate(sel.fields):
                if field.kind != "forest":
                    continue
                row = sampled[field.id]
                pixel = pixel_width_angstrom_to_velocity(
                    float(reader.tracers[field.id]["pix_width_ang"]),
                    lambda_obs_angstrom=LYA_REST_ANGSTROM * (1 + ms["z_eval"]),
                )
                sigma = SPEED_LIGHT_KMS / (
                    float(reader.config["survey"]["resolution"])
                    * 2
                    * np.sqrt(2 * np.log(2))
                )
                kp = 0.00035
                k = np.hypot(av * kp, 2.4 / ddeg)
                mu = av * kp / k
                transfer = velocity_response(
                    [kp], pixel_width_velocity=pixel, gaussian_sigma_velocity=sigma
                )[0]
                signal = float(
                    model(
                        [1.0, 1.0],
                        ms["z_eval"],
                        np.array([k]),
                        np.array([mu]),
                        np.array([[fi, fi]]),
                    )[0, 0]
                    * transfer**2
                    * av
                    / ddeg**2
                )
                alias = float(default_p1d([], ms["z_eval"], [kp])[0] * transfer**2)
                masses = (
                    density_per_velocity(row["density"], z_source=row["z_source"]) * q
                )
                variance = row["variance"]
                wa, wr = diagnose(
                    masses, variance, row["length_velocity"], pixel, signal, alias
                )
                prefix = f"order{order}_field{fi}_"
                saved_arrays.update({prefix + n: v for n, v in wa.items()})
                saved_arrays.update(
                    {
                        prefix + n: np.asarray(v)
                        for n, v in dict(
                            magnitudes=m, measure=q, masses=masses, variance=variance
                        ).items()
                    }
                )
                w0 = 1 / (1 + pixel * variance / alias)
                _, _, _, aa, pp = _integrals(
                    masses, w0, variance, row["length_velocity"], pixel
                )
                coefficients = {
                    0: (aa, pp),
                    **{
                        r["iterations"]: (r["A"], r["P_pixel"])
                        for r in wr["rows"]
                        if r["available"]
                    },
                }
                mode = a["k"] * a["mu"] / av
                w = velocity_response(
                    mode, pixel_width_velocity=pixel, gaussian_sigma_velocity=sigma
                )
                p1d = default_p1d([], ms["z_eval"], mode)
                pair = list(map(tuple, sel.required_pairs)).index((fi, fi))
                for count in noise_by_count:
                    if count not in coefficients:
                        available[count] = False
                    else:
                        aa, pp = coefficients[count]
                        noise_by_count[count][:, pair] = (
                            (aa * p1d * w * w + pp) * ddeg * ddeg / av
                        )
                fields.append(
                    dict(
                        field=field.id,
                        signal=signal,
                        alias=alias,
                        fixed_coefficients=list(coefficients[0]),
                        **wr,
                    )
                )
            for count in noise_by_count:
                entry = dict(
                    order=order,
                    iterations=count,
                    fixed_initial_weights=count == 0,
                    available=available[count],
                    fields=fields,
                )
                if available[count]:
                    total = a["observed_signal"] + noise_by_count[count]
                    f, pf = contract(
                        wick(
                            total,
                            a["modes"],
                            sel.required_pairs,
                            sel.selected_pairs,
                            len(sel.fields),
                        ),
                        a["observed_j"],
                        independent=True,
                    )
                    key = f"order{order}_iterations{count}_"
                    saved_arrays.update(
                        {key + n: v for n, v in summaries(f, pf).items()}
                    )
                    entry["fisher_key"] = key + "fisher"
                    entry["pair_fisher_key"] = key + "pair_fisher"
                    if (
                        order == report["final_controls"]["magnitude_order"]
                        and count == report["final_controls"]["iterations"]
                    ):
                        entry["primary_fisher_relative"] = relative(f, a["fisher"])
                        entry["primary_noise_relative"] = relative(
                            noise_by_count[count], a["noise"]
                        )
                        assert (
                            entry["primary_fisher_relative"] < 5e-12
                            and entry["primary_noise_relative"] < 5e-12
                        )
                rows.append(entry)
            reader._samples.clear()
        comparisons = []
        for count in (0, 3, 6, 12, 24):
            valid = [r for r in rows if r["iterations"] == count and r["available"]]
            for low, high in zip(valid[:-1], valid[1:]):
                comparisons.append(
                    dict(
                        iterations=count,
                        orders=[low["order"], high["order"]],
                        **change(
                            saved_arrays[low["fisher_key"]],
                            saved_arrays[high["fisher_key"]],
                            saved_arrays[low["pair_fisher_key"]],
                            saved_arrays[high["pair_fisher_key"]],
                            1,
                            1,
                        ),
                    )
                )
        path = out / f"bin-{index}.npz"
        np.savez_compressed(path, **saved_arrays)
        result["bins"].append(
            dict(
                bin=index,
                partition=s["partition"],
                rows=rows,
                comparisons=comparisons,
                arrays=path.name,
                sha256=digest(path),
            )
        )
        (out / "diagnosis.json").write_text(
            json.dumps(plain(result), indent=2, allow_nan=False)
        )
        print(
            f"bin {index}: fixed and coupled magnitude/weight diagnosis saved",
            flush=True,
        )
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("bundle", "reference", "template", "output"):
        p.add_argument("--" + name, required=True)
    args = p.parse_args()
    run(args.bundle, args.reference, args.template, args.output)
