"""Direct scientific S4 controls; validation only, no production profile changes.

Legacy peak and covariance algebra adapted from lyaforecast (GPLv3).
"""

import argparse
import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np

from fishhighz.adapters.legacy_compat import plain
from fishhighz.geometry import SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.p1d import default_p1d
from fishhighz.response import velocity_response
from fishhighz.validation.accuracy import AccuracyRecipe, background
from fishhighz.validation.adaptive_weights import adaptive_weights
from fishhighz.validation.cases import bins
from fishhighz.validation.compatibility_weights import WeightInputs
from fishhighz.validation.numerics import contract, relative, summaries, wick
from fishhighz.validation.profile_definitions import forecast_selection

CASE = "lya_qso_lbg_lae_15x2pt"
FIXED = dict(
    magnitude="rect107",
    negative=False,
    fourier=False,
    volume=False,
    redshift=False,
    response=False,
    growth=False,
    power=False,
    covariance_damping=False,
    cross_damping=False,
    peak=False,
    operator=False,
    field_response=False,
    constants=False,
    stopping=False,
)
ACCURACY = {k: ("gauss" if k == "magnitude" else True) for k in FIXED}


class Direct:
    def __init__(self, recipe, index, fixed, accurate, full_report):
        self.recipe, self.index = recipe, index
        self.fixed, self.accurate = fixed, accurate
        self.sel = forecast_selection(CASE, index)
        self.fields = self.sel.fields
        self.active = np.unique(self.sel.required_pairs)
        self.pairs = self.sel.required_pairs
        self.names = [f"{self.fields[i].id}_{self.fields[j].id}" for i, j in self.pairs]
        self.rows = full_report["settings"]["pair_inputs"]
        self.zc = recipe.z(index)
        self.zm = sum(bins(CASE)[index]) / 2
        self.ka = np.unique(fixed["k"])
        self.dk = self.ka[1] - self.ka[0]
        self.av = recipe.cosmo.velocity_from_distance(self.zc)
        self.dd = recipe.cosmo.distance_from_degrees(self.zc)
        self.volume_fixed = float(full_report["settings"]["grid"]["volume"])
        lo, hi = bins(CASE)[index]
        self.geometry = prepare_geometry(
            lo,
            hi,
            z_eval=self.zc,
            area_deg2=float(recipe.config["survey"]["survey_area"]),
            h_fid=recipe.h,
            z_order=32,
            hubble=recipe.cosmo.results.hubble_parameter,
            transverse_distance=recipe.cosmo.results.comoving_radial_distance,
        )
        self.pk_k, _, p = recipe.cosmo.results.get_matter_power_spectrum(
            minkh=float(recipe.config["power spectrum"]["k_min_hmpc"]),
            maxkh=float(recipe.config["power spectrum"]["k_max_hmpc"]),
            npoints=1000,
        )
        self.pk = p[0]
        self.cache = {}
        self.peak_cache = {}

    def growth(self, z, c):
        r = self.recipe
        if c["growth"]:
            return (
                r._growth(z)[0]
                / r._growth(r.template.z_ref if c["power"] else r.cosmo.z_ref)[0]
            ) ** 2
        return (
            (1 + (r.template.z_ref if c["power"] else r.cosmo.z_ref)) / (1 + z)
        ) ** 2

    def linear(self, k, z, c):
        raw = (
            self.recipe.template.evaluate(k).sum(axis=1)
            if c["power"]
            else np.interp(k, self.pk_k, self.pk)
        )
        return raw * self.growth(z, c)

    def rsd(self, mu, z):
        b = self.recipe.external.bias
        factors = {}
        for i in self.active:
            f = self.fields[i]
            t = self.recipe.tracers[f.id]["tracer"]
            factors[i] = b._get_density_bias(z, t) * (1 + b._get_beta_rsd(z, t) * mu**2)
        return np.column_stack([factors[i] * factors[j] for i, j in self.pairs])

    def response(self, k, mu, z, c):
        speed = SPEED_LIGHT_KMS if c["constants"] else 299800.0
        div = 2 * np.sqrt(2 * np.log(2)) if c["response"] else 1
        wave = 1215.67 * (1 + self.zc)
        av = self.recipe.cosmo.velocity_from_distance(z)
        own = {}
        for i in self.active:
            f = self.fields[i]
            if f.kind == "forest":
                pixel = float(self.recipe.tracers[f.id]["pix_width_ang"]) * speed / wave
                own[i] = velocity_response(
                    k * mu / av,
                    pixel_width_velocity=pixel,
                    gaussian_sigma_velocity=speed
                    / (float(self.recipe.config["survey"]["resolution"]) * div),
                )
            else:
                own[i] = np.ones(len(k))
        result = []
        for (i, j), name in zip(self.pairs, self.names):
            if c["field_response"]:
                result.append(own[i] * own[j])
                continue
            n = int(self.fields[i].kind == "forest") + int(
                self.fields[j].kind == "forest"
            )
            if not n:
                result.append(np.ones(len(k)))
                continue
            row = self.rows[name]
            pixel = row["_pix_kms"] * speed / 299800.0
            sigma = row["_res_kms"] * speed / 299800.0 / div
            result.append(
                velocity_response(
                    k * mu / av,
                    pixel_width_velocity=pixel,
                    gaussian_sigma_velocity=sigma,
                )
                ** n
            )
        return np.column_stack(result)

    def widths(self, z, c):
        sigma, f = self.recipe._growth(z)
        st = 3.26 * sigma / self.recipe.cosmo.sigma8
        recon = float(self.recipe.config["survey"]["reconstruction factor"])
        widths = []
        for i, j in self.pairs:
            if c["cross_damping"]:
                factor = (
                    (1 if self.fields[i].kind == "forest" else 1 / recon)
                    + (1 if self.fields[j].kind == "forest" else 1 / recon)
                ) / 2
            else:
                factor = (
                    1
                    if any(self.fields[a].kind == "forest" for a in (i, j))
                    else 1 / recon
                )
            widths.append([((1 + f) * st) ** 2 * factor, st**2 * factor])
        return np.array(widths)

    def damping(self, k, mu, z, c):
        w = self.widths(z, c)
        return np.exp(
            -0.5
            * k[:, None] ** 2
            * (mu[:, None] ** 2 * w[:, 0] + (1 - mu[:, None] ** 2) * w[:, 1])
        )

    def peak(self, k, mu, z, c):
        """Intrinsic peak, stripping fixed response after legacy observed fit."""
        if c["peak"]:
            return (
                self.recipe.template.evaluate(k)[:, 1, None]
                * self.growth(z, {**c, "power": True})
                * self.rsd(mu, z)
            )
        result = np.empty((len(k), len(self.pairs)))
        x0 = np.log(self.ka)
        center = x0.mean()
        span = np.ptp(x0)
        x = (x0 - center) / span
        weights = np.ones(len(x))
        weights[:3] = 1e8
        for value in np.unique(mu):
            sel = mu == value
            key = (tuple(c.items()), float(value), z)
            if key not in self.peak_cache:
                native_mu = np.full(len(self.ka), value)
                native = (
                    self.linear(self.ka, z, c)[:, None]
                    * self.rsd(native_mu, z)
                    * self.response(self.ka, native_mu, z, c)
                )
                coef = np.stack(
                    [
                        np.polyfit(x, np.log(abs(row) + 1e-12), 8, w=weights)
                        for row in native.T
                    ]
                )
                self.peak_cache[key] = coef
            coef = self.peak_cache[key]
            response = self.response(k[sel], mu[sel], z, c)
            model = self.linear(k[sel], z, c)[:, None] * self.rsd(mu[sel], z)
            smooth = np.column_stack(
                [
                    np.exp(np.polyval(row, (np.log(k[sel]) - center) / span))
                    for row in coef
                ]
            )
            result[sel] = model - np.sign(model) * smooth / response
        return result

    def intrinsic(self, k, mu, z, c):
        result = self.linear(k, z, c)[:, None] * self.rsd(mu, z)
        if c["covariance_damping"]:
            # Keep power-source and decomposition switches separable: add the
            # selected template wiggle's damping correction to the selected P_lin.
            wiggle = (
                self.recipe.template.evaluate(k)[:, 1, None]
                * self.growth(z, {**c, "power": True})
                * self.rsd(mu, z)
            )
            result = result + wiggle * (self.damping(k, mu, z, c) - 1)
        # The legacy forest-auto coordinate roundtrip uses mu=k_parallel/(k+1e-10).
        if not c["field_response"]:
            shifted = self.rsd(mu * k / (k + 1e-10), z)
            original = self.rsd(mu, z)
            for col, (i, j) in enumerate(self.pairs):
                if i == j and self.fields[i].kind == "forest":
                    result[:, col] *= shifted[:, col] / original[:, col]
        return result

    def inputs(self, c):
        key = tuple(
            (k, v)
            for k, v in c.items()
            if k
            not in (
                "fourier",
                "volume",
                "redshift",
                "peak",
                "operator",
                "cross_damping",
            )
        )
        if key in self.cache:
            return self.cache[key]
        r = self.recipe
        speed = SPEED_LIGHT_KMS if c["constants"] else 299800.0
        sampled, m, q = r.samples(self.index, 16)
        if c["magnitude"] != "gauss":
            m = np.linspace(
                float(r.config["survey"]["min_band_mag"]),
                float(r.config["survey"]["max_band_mag"]),
                int(c["magnitude"][4:]),
            )
            q = np.full(len(m), m[1] - m[0])
        noise = {}
        records = {}
        for i in self.active:
            field = self.fields[i]
            row = sampled[field.id]
            zs = row["z_source"]
            d = r.densities[field.id].sample(zs, m)
            rho = (
                np.array(
                    d["values"]
                    if c["negative"]
                    else np.where(
                        d["provenance"]["masks"]["density_floor"], 1e-20, d["raw"]
                    )
                )
                * (1 + zs)
                / speed
            )
            if field.kind != "forest":
                noise[field.id] = self.dd**2 / self.av / np.sum(rho * q)
                records[field.id] = {
                    "density_integral": float(np.sum(rho * q)),
                    "negative_nodes": int(np.sum(rho < 0)),
                }
                continue
            t = r.tracers[field.id]
            s = r.snrs[field.id].sample(
                z_source=zs,
                magnitudes=m,
                wavelength=1215.67 * (1 + self.zc),
                pixel_width_angstrom=float(t["pix_width_ang"]),
                exposure_count=float(t["num exposures"]),
            )
            var = np.asarray(s["values"])
            # SNR reader uses the physical pixel conversion; legacy SNR per Angstrom
            # conversion depends on wavelength/pixel Angstrom only, hence unchanged.
            pixel = float(t["pix_width_ang"]) * speed / (1215.67 * (1 + self.zc))
            length = speed * np.log(
                float(t["max_rest_frame_lya"]) / float(t["min_rest_frame_lya"])
            )
            kk = np.hypot(2.4 / self.dd, 0.00035 * self.av)
            mm = 0.00035 * self.av / kk
            col = self.names.index(f"{field.id}_{field.id}")
            p = (
                float(
                    (
                        self.intrinsic(np.array([kk]), np.array([mm]), self.zc, c)
                        * self.response(np.array([kk]), np.array([mm]), self.zc, c)
                    )[0, col]
                )
                * self.av
                / self.dd**2
            )
            div = 2 * np.sqrt(2 * np.log(2)) if c["response"] else 1
            sigma = speed / (float(r.config["survey"]["resolution"]) * div)
            b = float(
                default_p1d([], self.zc, [0.00035])[0]
                * velocity_response(
                    [0.00035], pixel_width_velocity=pixel, gaussian_sigma_velocity=sigma
                )[0]
                ** 2
            )
            inp = WeightInputs(m, rho, q, var, length, pixel, p, b)
            solved = adaptive_weights(
                inp,
                "sum_historical",
                context=f"{field.id}_{field.id}",
                rtol=1e-5 if c["stopping"] else 1e-4,
            )
            if solved["status"] != "converged":
                raise ValueError((field.id, solved["status"]))
            a, pp = solved["coefficients"][-2:]
            noise[field.id] = (a, pp, pixel, sigma)
            records[field.id] = dict(
                P=p,
                B=b,
                A=a,
                P_pixel=pp,
                updates=solved["updates"],
                status=solved["status"],
                weights=solved["weights"],
                magnitudes=m,
                quadrature=q,
                rho=rho,
                variance=var,
                length=length,
                pixel=pixel,
                negative_nodes=int(np.sum(rho < 0)),
                snr_counts=s["provenance"]["counts"],
                confirmation=solved["confirmation"],
            )
        self.cache[key] = (noise, records)
        return noise, records

    def evaluate(self, c, *, refine=False):
        if c["fourier"]:
            grid = gauss_legendre_grid(
                np.linspace(0.01, 0.5, 257 if refine else 129),
                k_order=4,
                mu_order=64 if refine else 32,
                h_fid=self.recipe.h,
            )
            k, mu, measure = grid.k_flat, grid.mu_flat, grid.q_mode
        else:
            k, mu = self.fixed["k"], self.fixed["mu"]
            measure = self.fixed["modes"] / self.volume_fixed
        volume = self.geometry.volume if c["volume"] else self.volume_fixed
        modes = measure * volume
        noise, records = self.inputs(c)
        total = self.intrinsic(k, mu, self.zc, c) * self.response(k, mu, self.zc, c)
        for col, (i, j) in enumerate(self.pairs):
            if i != j:
                continue
            name = self.fields[i].id
            if self.fields[i].kind == "forest":
                a, pp, pixel, sigma = noise[name]
                response = velocity_response(
                    k * mu / self.av,
                    pixel_width_velocity=pixel,
                    gaussian_sigma_velocity=sigma,
                )
                total[:, col] += (
                    (a * default_p1d([], self.zc, k * mu / self.av) * response**2 + pp)
                    * self.dd**2
                    / self.av
                )
            else:
                total[:, col] += noise[name]
        z = self.zc if c["redshift"] else self.zm
        response = self.response(k, mu, z, c)
        if c["operator"]:
            step = 1.25e-4 if refine else 2.5e-4
            deriv = []
            for axis in range(2):
                values = []
                for sign in (1, -1):
                    ap = 1 + sign * step if axis == 0 else 1.0
                    at = 1 + sign * step if axis == 1 else 1.0
                    kp = k * mu / ap
                    kt = k * np.sqrt(1 - mu**2) / at
                    km = np.hypot(kp, kt)
                    mm = kp / km
                    values.append(
                        self.peak(km, mm, z, c)
                        * self.damping(km, mm, z, c)
                        / (ap * at**2)
                    )
                deriv.append((values[0] - values[1]) / (2 * step) * response)
            jac = np.stack(deriv, axis=-1)
        else:
            previous = np.maximum(k - self.dk, self.ka[0])
            current = self.peak(k, mu, z, c) * self.damping(k, mu, z, c) * response
            older = (
                self.peak(previous, mu, z, c)
                * self.damping(previous, mu, z, c)
                * self.response(previous, mu, z, c)
            )
            derivative = (current - older) / (self.dk / k)[:, None]
            derivative[k <= self.ka[0] + self.dk * 0.000001] = 0
            jac = (
                derivative[:, :, None] * np.column_stack((mu**2, 1 - mu**2))[:, None, :]
            )
        covariance = wick(
            total, modes, self.pairs, self.sel.selected_pairs, len(self.fields)
        )
        fisher, single = contract(covariance, jac)
        result = summaries(fisher, single)
        result.update(
            k=k,
            mu=mu,
            modes=modes,
            total=total,
            observed_j=jac,
            selected_covariance=covariance,
            required_pairs=self.pairs,
            selected_pairs=self.sel.selected_pairs,
        )
        return result, dict(
            controls=c,
            volume=volume,
            mean_z=z,
            covariance_z=self.zc,
            weights=plain(records),
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bins", type=int, nargs="+", default=list(range(6)))
    parser.add_argument("--endpoints-only", action="store_true")
    parser.add_argument("--witness", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    bundle = Path(".validation/s2/profiles-r1")
    manifest = json.loads((bundle / "manifest.json").read_text())
    cosmo, template = background(
        "../lyaforecast", "../vega/vega/models/Planck18/Planck18_z_2.406.fits"
    )
    recipe = AccuracyRecipe("../lyaforecast", CASE, cosmo, template, {})
    results = []
    for index in args.bins:
        saved = {}
        reports = {}
        for entry in manifest["records"]:
            if entry["task"]["bin"] != index:
                continue
            name = entry["task"]["profile"]
            saved[name] = dict(np.load(bundle / entry["arrays"]))
            reports[name] = json.loads((bundle / entry["report_file"]).read_text())
        direct = Direct(
            recipe,
            index,
            saved["fixed-compatibility"],
            saved["accuracy"],
            reports["full-compatibility"],
        )
        trials = [("fixed", FIXED), ("accuracy", ACCURACY)]
        if not args.endpoints_only:
            for key in FIXED:
                trials.extend(
                    [
                        (f"forward-{key}", {**FIXED, key: ACCURACY[key]}),
                        (f"reverse-{key}", {**ACCURACY, key: FIXED[key]}),
                    ]
                )
            trials += [
                ("forward-rect425", {**FIXED, "magnitude": "rect425"}),
                ("reverse-rect425", {**ACCURACY, "magnitude": "rect425"}),
            ]
        if args.witness:
            trials = [
                (label, c)
                for label, c in trials
                if label
                in (
                    "forward-cross_damping",
                    "reverse-cross_damping",
                    "forward-operator",
                    "reverse-operator",
                )
            ]
        completed = {}

        def run_trial(label, c, refine=False):
            start = time.monotonic()
            arrays, report = direct.evaluate(c, refine=refine)
            name = f"bin{index + 1}-{label}"
            operands = (
                "k",
                "mu",
                "modes",
                "total",
                "observed_j",
                "selected_covariance",
            )
            retained = {k: v for k, v in arrays.items() if k not in operands}
            sample = np.unique(np.linspace(0, len(arrays["k"]) - 1, 7, dtype=int))
            retained.update({f"sample_{k}": arrays[k][sample] for k in operands})
            if args.witness or label in ("fixed", "accuracy"):
                retained.update(
                    {k: arrays[k] for k in operands if k != "selected_covariance"}
                )
            np.savez_compressed(args.output / f"{name}.npz", **retained)
            report["refinement"] = refine
            (args.output / f"{name}.json").write_text(
                json.dumps(plain(report), indent=2, allow_nan=False) + "\n"
            )
            record = dict(
                bin=index,
                label=label,
                arrays=f"{name}.npz",
                report=f"{name}.json",
                errors=arrays["errors"].tolist(),
                pair_errors=arrays["pair_errors"].tolist(),
            )
            if label in ("fixed", "accuracy"):
                old = saved["fixed-compatibility" if label == "fixed" else "accuracy"]
                record["closure"] = {
                    key: relative(arrays[key], old[key])
                    for key in ("fisher", "pair_fisher", "errors", "pair_errors")
                }
                record["closure_max_error"] = float(
                    np.max(
                        abs(
                            np.r_[
                                arrays["errors"] / old["errors"] - 1,
                                (
                                    arrays["pair_errors"] / old["pair_errors"] - 1
                                ).ravel(),
                            ]
                        )
                    )
                )
                if record["closure_max_error"] > 1e-4:
                    raise ValueError(("endpoint closure", record))
            results.append(record)
            completed[label] = record
            (args.output / "summary.json").write_text(
                json.dumps(plain(results), indent=2, allow_nan=False) + "\n"
            )
            print(
                index + 1,
                label,
                record.get("closure", record["errors"]),
                f"{time.monotonic() - start:.2f}s",
                flush=True,
            )
            return record

        for label, c in trials:
            run_trial(label, c)
        if not args.endpoints_only and not args.witness:

            def error_vector(row):
                return np.r_[row["errors"], np.ravel(row["pair_errors"])]

            f0 = error_vector(completed["fixed"])
            a0 = error_vector(completed["accuracy"])
            differences = {}
            for key in FIXED:
                forward = np.log(error_vector(completed["forward-" + key]) / f0)
                reverse = np.log(error_vector(completed["reverse-" + key]) / a0)
                differences[key] = float(np.max(abs(forward + reverse)))
            # Three strongest endpoint-dependent effects implicate at most three
            # combinations, plus the two predeclared estimator interactions.
            strongest = sorted(differences, key=differences.get, reverse=True)[:3]
            pairs = set(combinations(sorted(strongest), 2))
            pairs.update([("fourier", "operator"), ("operator", "peak")])
            for x, y in sorted(pairs):
                for direction, base, target in [
                    ("forward", FIXED, ACCURACY),
                    ("reverse", ACCURACY, FIXED),
                ]:
                    c = {**base, x: target[x], y: target[y]}
                    run_trial(f"{direction}-{x}+{y}", c)
            # Combined numerical refinement of each GL control, including all
            # retained weak spectra; fixed primary remains a reproduction recipe.
            for label, c in trials:
                if c["fourier"]:
                    fine = run_trial(label + "-refined", c, refine=True)
                    coarse = completed[label]
                    fine["max_error_refinement"] = float(
                        np.max(abs(error_vector(fine) / error_vector(coarse) - 1))
                    )
            (args.output / "summary.json").write_text(
                json.dumps(plain(results), indent=2, allow_nan=False) + "\n"
            )


if __name__ == "__main__":
    main()
