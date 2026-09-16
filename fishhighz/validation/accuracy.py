"""Converged existing-model accuracy profile and explicitly labeled input studies."""

import configparser
from collections import OrderedDict
from pathlib import Path

import numpy as np

from ..adapters.legacy_compat import LegacyDensity, LegacySNR, plain
from ..adapters.legacy_inputs import DensityReader, SNRReader
from ..derivatives import evaluate_derivatives
from ..forecast import prepare_bin
from ..geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS, prepare_geometry
from ..grids import gauss_legendre_grid
from ..models.external import BoundParameters, P3DProvider, PreparedP3D
from ..models.kaiser import KaiserModel, Scaling
from ..models.p1d import default_p1d
from ..models.templates import load_template
from ..noise import local_galaxy_density
from ..parameters import Parameter, ParameterRegistry
from ..response import (
    InstrumentResponse,
    pixel_width_angstrom_to_velocity,
    velocity_response,
)
from ..survey import BinSpec, ForestInput
from ..weights import density_per_velocity
from .cases import CASE_IDS, bins, recipe, selection, verify_inventory
from .numerics import contract, relative
from .reference_capture import imported_reference, resolved_resources
from .schema import _assemble

DEFAULT = dict(
    k_intervals=128,
    mu_order=32,
    z_order=32,
    magnitude_order=16,
    iterations=12,
    step=2.5e-4,
)
FIXED_REFERENCE = dict(
    convention="intrinsic_p1d_times_field_response_squared",
    q_star=0.00035,
    q_star_units="s/km",
)


def _forest_input(
    field,
    row,
    magnitudes,
    quadrature,
    response,
    registry,
    z_eval,
    *,
    method,
    iterations=None,
    policy="primary",
):
    """Construct the accuracy-profile forest input and declared weight metadata."""
    common = dict(
        z_source=row["z_source"],
        magnitudes=magnitudes,
        quadrature=quadrature,
        rho=density_per_velocity(row["density"], z_source=row["z_source"]),
        variance=row["variance"],
        length_velocity=row["length_velocity"],
    )
    reference = None
    if method == "inverse_variance":
        if iterations is not None:
            raise ValueError("iterations are inapplicable to inverse_variance accuracy")
        q_star = FIXED_REFERENCE["q_star"]
        intrinsic_p1d = float(default_p1d([], z_eval, [q_star])[0])
        response_factor = float(
            velocity_response(
                [q_star],
                pixel_width_velocity=response.pixel_width_velocity,
                gaussian_sigma_velocity=response.gaussian_sigma_velocity,
            )[0]
        )
        b_star = intrinsic_p1d * response_factor**2
        if not np.isfinite(b_star) or b_star <= 0:
            raise ValueError(
                f"{field.id}: fixed reference B_star is not positive finite"
            )
        options = dict(common, method=method, alias=b_star)
        auxiliary_coordinates = None
        reference = dict(
            **FIXED_REFERENCE,
            intrinsic_p1d=intrinsic_p1d,
            intrinsic_p1d_units="km/s",
            response_factor=response_factor,
            B_star=b_star,
            B_star_units="km/s",
        )
        iteration_status = {"applicable": False, "status": "inapplicable"}
    elif method == "legacy":
        if iterations is None:
            raise ValueError("legacy accuracy requires an explicit iteration count")
        options = dict(common, method=method, iterations=iterations)
        auxiliary_coordinates = (2.4, 0.00035)
        iteration_status = {"applicable": True, "count": iterations}
    else:
        raise ValueError("accuracy weight method must be legacy or inverse_variance")
    weighting = dict(
        method=method,
        iterations=iteration_status,
        reference=reference,
    )
    source = ForestInput(
        options,
        default_p1d,
        BoundParameters(registry, (), {}),
        registry.fiducials,
        auxiliary_coordinates=auxiliary_coordinates,
        provenance=dict(
            policy=policy,
            density=row["density_diagnostics"],
            snr=row["snr_diagnostics"],
            weighting=weighting,
        ),
    )
    return source, weighting


def background(root, template_path):
    """One caller-prepared CAMB background at every actual evaluation redshift."""
    import camb
    from astropy.io import fits
    from lyaforecast.cosmoCAMB import CosmoCamb

    root = Path(root).resolve()
    imported_reference(root)
    with fits.open(template_path) as f:
        zt = float(f[1].header["ZREF"])
    zs = sorted(
        {
            2.3,
            zt,
            1.8,
            4.5,
            *[
                float(np.sqrt((1 + a) * (1 + b)) - 1)
                for c in CASE_IDS
                for a, b in bins(c)
            ],
            *[(a + b) / 2 for c in CASE_IDS for a, b in bins(c)],
        }
    )
    ini = root / "lyaforecast/resources/camb_configs/Planck18.ini"
    parsed = camb.read_ini(str(ini))
    cosmo = CosmoCamb(str(ini), z_ref=2.3, z_centres=zs)
    template = load_template(template_path, h_fid=parsed.H0 / 100)
    return cosmo, template


def composite(partition, order):
    """Ordered Gauss-Legendre magnitude nodes/weights on a fixed partition."""
    x, w = np.polynomial.legendre.leggauss(order)
    lo, hi = np.asarray(partition[:-1]), np.asarray(partition[1:])
    return ((lo[:, None] + hi[:, None]) / 2 + (hi - lo)[:, None] * x / 2).ravel(), (
        (hi - lo)[:, None] * w / 2
    ).ravel()


def breakpoints(densities, snrs, z_queries, lo, hi):
    """Partition spline/support/SNR boundaries and negative-interpolant roots.

    Quadratic pieces at fixed z are fitted on each existing spline interval only
    to locate zeros; this does not replace or modify the density interpolation.
    """
    points = [lo, hi]
    for name, density in densities.items():
        r = density.reader
        knots = np.unique(np.r_[r.magnitudes, density._spline.get_knots()[1], lo, hi])
        knots = knots[(knots >= lo) & (knots <= hi)]
        points.extend(knots)
        for a, b in zip(knots[:-1], knots[1:]):
            if a < r.magnitudes[0] or b > r.magnitudes[-1]:
                continue
            y = density._spline.ev(np.full(3, z_queries[name]), [a, (a + b) / 2, b])
            # Local t in [0,1], coefficients of the exact quadratic segment.
            c = y[0]
            aa = 2 * (y[2] - 2 * y[1] + y[0])
            bb = y[2] - y[0] - aa
            roots = np.roots([aa, bb, c]) if aa != 0 else ([-c / bb] if bb != 0 else [])
            for t in roots:
                if np.isreal(t) and 1e-12 < float(np.real(t)) < 1 - 1e-12:
                    points.append(a + (b - a) * float(np.real(t)))
    for snr in snrs.values():
        points.extend(snr.reader.magnitudes)
    p = np.unique(np.asarray(points))
    p = p[(p >= lo) & (p <= hi)]
    # Coalesce only numerically indistinguishable duplicate boundaries.
    return p[
        np.r_[True, np.diff(p) > 64 * np.finfo(float).eps * np.maximum(1, abs(p[1:]))]
    ]


class AccuracyRecipe:
    """Prepare source snapshots once and reuse fixed state during finite differences."""

    def __init__(
        self,
        root,
        case,
        cosmo,
        template,
        provenance,
        *,
        weight_method="inverse_variance",
    ):
        from lyaforecast.power_spectrum import PowerSpectrum
        from scipy.interpolate import interp1d

        self.root = Path(root).resolve()
        verify_inventory(self.root / "examples/desi2")
        self.case = case
        self.config = configparser.ConfigParser()
        self.config.read_dict(recipe(case))
        self.selection = selection(case)
        self.cosmo = cosmo
        self.template = template
        self.h = template.h_fid
        self.external = PowerSpectrum(self.config, cosmo, {})
        self.provenance = provenance
        if weight_method not in ("inverse_variance", "legacy"):
            raise ValueError(
                "accuracy weight method must be legacy or inverse_variance"
            )
        self.weight_method = weight_method
        resolved_resources(self.root, self.config)
        self.registry = ParameterRegistry(
            [
                Parameter(f"{name}_{i}", 1, "target", step=0.001)
                for i in range(len(bins(case)))
                for name in ("ap", "at")
            ]
        )
        self.densities = {}
        self.snrs = {}
        self.tracers = {}
        for field, key in zip(
            self.selection.fields,
            [s for s in self.config.sections() if s.startswith("tracer ")],
        ):
            t = self.config[key]
            self.tracers[field.id] = dict(t)
            if "bias z" in t:
                self.external.bias.set_density_bias_func(
                    t["tracer"],
                    interp1d(
                        np.fromstring(t["bias z"], sep=" "),
                        np.fromstring(t["bias val"], sep=" "),
                        bounds_error=False,
                        fill_value="extrapolate",
                    ),
                )
            forest = field.kind == "forest"
            reader = DensityReader(
                self.root / "lyaforecast/resources/data" / t["dn dz"],
                semantics="cell_count_per_deg2",
                target_density=t.getfloat("target density"),
                z_norm_min=2.15 if forest or field.id == "qso" else None,
                magnitude_bounds=(
                    t.getfloat("min_band_mag"),
                    t.getfloat("max_band_mag"),
                )
                if forest
                else None,
                width_policy="legacy_first_spacing",
                label=field.id,
            )
            self.densities[field.id] = LegacyDensity(reader, "floor_negative")
            if forest:
                paths = sorted(
                    (self.root / "lyaforecast/resources/data" / t["snr-file-dir"]).glob(
                        "*.dat"
                    )
                )
                self.snrs[field.id] = LegacySNR(
                    SNRReader(paths, smoothing="legacy", label=field.id)
                )
        self._partitions = {}
        self._samples = {}
        self._prepared = OrderedDict()

    def z(self, index):
        lo, hi = bins(self.case)[index]
        return float(np.sqrt((1 + lo) * (1 + hi)) - 1)

    def _growth(self, z):
        indices = np.flatnonzero(self.cosmo.z_bins == z)
        if len(indices) != 1:
            raise ValueError("CAMB must be prepared at exact requested redshift")
        i = indices[0]
        return float(self.cosmo.sigma8_zbins[i]), float(self.cosmo.growth_rate_zbins[i])

    def model(self, index, *, mean_z=None, growth="camb", reconstruction=True):
        z = self.z(index) if mean_z is None else float(mean_z)
        sigma, f = self._growth(z)
        sigma_template, _ = self._growth(self.template.z_ref)
        g = (
            (sigma / sigma_template) ** 2
            if growth == "camb"
            else ((1 + self.template.z_ref) / (1 + z)) ** 2
        )
        biases = {}
        betas = {}
        widths = {}
        for field in self.selection.fields:
            name = self.tracers[field.id]["tracer"]
            biases[field.id] = float(self.external.bias._get_density_bias(z, name))
            if field.kind == "forest":
                betas[field.id] = float(self.external.bias._get_beta_rsd(z, name))
            r = (
                1
                if field.kind == "forest" or not reconstruction
                else self.config["survey"].getfloat("reconstruction factor")
            )
            st = 3.26 * sigma / self.cosmo.sigma8 / np.sqrt(r)
            widths[field.id] = ((1 + f) * st, st)
        model = KaiserModel(
            self.template,
            self.selection.fields,
            biases=biases,
            betas=betas,
            widths=widths,
            f=f if any(f.kind == "galaxy" for f in self.selection.fields) else None,
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
            z=z,
            growth=g,
        )
        binding = BoundParameters(
            self.registry, ("ap", "at"), {n: f"{n}_{index}" for n in ("ap", "at")}
        )
        provider = model
        if mean_z is not None:
            # Supplied-mean diagnostic: keep geometric noise/response coordinates
            # while explicitly evaluating the fixed intrinsic model at mean_z.
            def provider(theta, requested_z, k, mu, pairs):
                if requested_z != self.z(index):
                    raise ValueError("diagnostic geometry redshift mismatch")
                return model(theta, z, k, mu, pairs)

        p3d = PreparedP3D(
            self.registry,
            self.selection,
            [
                P3DProvider(
                    "accuracy BAO", provider, binding, self.selection.required_pairs
                )
            ],
        )
        return p3d, dict(
            z_eval=z,
            biases=biases,
            betas=betas,
            f=f,
            G=g,
            growth=growth,
            sigma8=sigma,
            sigma8_template=sigma_template,
            sigma8_damping_reference=float(self.cosmo.sigma8),
            widths=widths,
        )

    def responses(self, index, *, resolution="fwhm"):
        wave = LYA_REST_ANGSTROM * (1 + self.z(index))
        responses = {}
        divisor = 2 * np.sqrt(2 * np.log(2)) if resolution == "fwhm" else 1
        for f in self.selection.fields:
            t = self.tracers[f.id]
            responses[f.id] = (
                InstrumentResponse(
                    pixel_width_angstrom_to_velocity(
                        float(t["pix_width_ang"]), lambda_obs_angstrom=wave
                    ),
                    SPEED_LIGHT_KMS
                    / (float(self.config["survey"]["resolution"]) * divisor),
                )
                if f.kind == "forest"
                else InstrumentResponse(0, 0)
            )
        return responses

    def samples(self, index, order, *, policy="primary", rectangular=False):
        key = (index, order, policy, rectangular)
        if key in self._samples:
            return self._samples[key]
        z = self.z(index)
        wave = LYA_REST_ANGSTROM * (1 + z)
        z_queries = {
            f.id: wave
            / np.sqrt(
                float(self.tracers[f.id]["min_rest_frame_lya"])
                * float(self.tracers[f.id]["max_rest_frame_lya"])
            )
            - 1
            if f.kind == "forest"
            else z
            for f in self.selection.fields
        }
        lo = float(self.config["survey"]["min_band_mag"])
        hi = float(self.config["survey"]["max_band_mag"])
        if index not in self._partitions:
            self._partitions[index] = breakpoints(
                self.densities, self.snrs, z_queries, lo, hi
            )
        if rectangular:
            m = np.linspace(lo, hi, int(self.config["survey"]["num mag bins"]))
            q = np.full(len(m), m[1] - m[0])
        else:
            m, q = composite(self._partitions[index], order)
        result = {}
        floor = {"floor_low": 1e-22, "floor_high": 1e-18}.get(policy, 1e-20)
        width_factor = {"width_minus": 0.9, "width_plus": 1.1}.get(policy, 1.0)
        if policy not in (
            "primary",
            "floor_low",
            "floor_high",
            "width_minus",
            "width_plus",
            "remove_bright",
            "remove_sentinel",
        ):
            raise ValueError("unknown explicit sensitivity policy")
        active = np.unique(self.selection.selected_pairs)
        for i, f in enumerate(self.selection.fields):
            if i not in active:
                continue
            d = self.densities[f.id].sample(z_queries[f.id], m)
            values = d["values"].copy()
            masks = d["provenance"]["masks"]
            values[
                np.asarray(masks["density_floor"])
                | np.asarray(masks["negative_density"])
            ] = floor
            values /= width_factor
            row = dict(
                magnitudes=m,
                quadrature=q,
                density=values,
                z_source=z_queries[f.id],
                density_diagnostics=plain(d["provenance"]),
            )
            if f.kind == "forest":
                t = self.tracers[f.id]
                s = self.snrs[f.id].sample(
                    z_source=z_queries[f.id],
                    magnitudes=m,
                    wavelength=wave,
                    pixel_width_angstrom=float(t["pix_width_ang"]),
                    exposure_count=float(t["num exposures"]),
                )
                remove = np.zeros(len(m), dtype=bool)
                if policy == "remove_bright":
                    remove = np.asarray(s["provenance"]["masks"]["bright_clamp"])
                if policy == "remove_sentinel":
                    remove = np.asarray(s["provenance"]["masks"]["out_of_range"])
                row["density"][remove] = 0
                row.update(
                    variance=s["values"],
                    snr_diagnostics=plain(s["provenance"]),
                    removed=remove,
                    length_velocity=SPEED_LIGHT_KMS
                    * np.log(
                        float(t["max_rest_frame_lya"]) / float(t["min_rest_frame_lya"])
                    ),
                )
            result[f.id] = row
        self._samples[key] = (result, m, q)
        return result, m, q

    def prepare(
        self,
        index,
        controls,
        *,
        policy="primary",
        resolution="fwhm",
        growth="camb",
        rectangular=False,
        reconstruction=True,
        arithmetic_mean=False,
    ):
        if self.weight_method == "inverse_variance" and "iterations" in controls:
            raise ValueError("iterations are inapplicable to inverse_variance accuracy")
        if self.weight_method == "legacy" and "iterations" not in controls:
            raise ValueError("legacy accuracy requires an explicit iteration count")
        key = (
            index,
            tuple(sorted((k, v) for k, v in controls.items() if k != "step")),
            self.weight_method,
            policy,
            resolution,
            growth,
            rectangular,
            reconstruction,
            arithmetic_mean,
        )
        if key in self._prepared:
            return self._prepared[key]
        lo, hi = bins(self.case)[index]
        z = self.z(index)
        geometry = prepare_geometry(
            lo,
            hi,
            z_eval=z,
            area_deg2=float(self.config["survey"]["survey_area"]),
            h_fid=self.h,
            z_order=controls["z_order"],
            hubble=self.cosmo.results.hubble_parameter,
            transverse_distance=self.cosmo.results.comoving_radial_distance,
        )
        grid = gauss_legendre_grid(
            np.linspace(0.01, 0.5, controls["k_intervals"] + 1),
            k_order=4,
            mu_order=controls["mu_order"],
            h_fid=self.h,
        )
        p3d, model_settings = self.model(
            index,
            growth=growth,
            reconstruction=reconstruction,
            mean_z=(lo + hi) / 2 if arithmetic_mean else None,
        )
        responses = self.responses(index, resolution=resolution)
        sampled, m, q = self.samples(
            index, controls["magnitude_order"], policy=policy, rectangular=rectangular
        )
        forests = {}
        forest_weighting = {}
        galaxies = {}
        for field in self.selection.fields:
            if field.id not in sampled:
                continue
            row = sampled[field.id]
            if field.kind == "forest":
                forests[field.id], forest_weighting[field.id] = _forest_input(
                    field,
                    row,
                    m,
                    q,
                    responses[field.id],
                    self.registry,
                    z,
                    method=self.weight_method,
                    iterations=controls.get("iterations"),
                    policy=policy,
                )
            else:
                galaxies[field.id] = local_galaxy_density(row["density"], q, geometry)
        prepared = prepare_bin(
            BinSpec(
                f"{self.case}-{index}",
                geometry,
                grid,
                p3d,
                responses,
                forests=forests,
                galaxies=galaxies,
                independent_sampling=True,
            )
        )
        settings = dict(
            profile="accuracy",
            parameters=[f"ap_{index}", f"at_{index}"],
            bounds=[lo, hi],
            fields=[f.id for f in self.selection.fields],
            grid=dict(
                kind="gauss_legendre",
                volume=geometry.volume,
                h_fid=self.h,
                k_intervals=controls["k_intervals"],
                k_order=4,
                mu_order=controls["mu_order"],
            ),
            controls=dict(controls),
            model=model_settings,
            resolution=resolution,
            policy=policy,
            rectangular=rectangular,
            reconstruction=reconstruction,
            arithmetic_mean=arithmetic_mean,
            partition=self._partitions[index].tolist(),
            magnitude_nodes=m.tolist(),
            magnitude_weights=q.tolist(),
            samples=plain(sampled),
            forest_weighting=dict(
                method=self.weight_method,
                reference=(
                    FIXED_REFERENCE
                    if self.weight_method == "inverse_variance"
                    else None
                ),
                iterations=(
                    {"applicable": False, "status": "inapplicable"}
                    if self.weight_method == "inverse_variance"
                    else {"applicable": True, "count": controls["iterations"]}
                ),
                forests=forest_weighting,
            ),
            galaxy_nbar=galaxies,
            geometry=dict(a_v=geometry.a_v, d_deg=geometry.d_deg),
            source_width_policy="legacy_first_spacing; unknown physical cells",
            noise_ownership="per-field independent sampling; pixel/Poisson unsmoothed",
            response_ownership="observed J, response already applied exactly once",
        )
        result = (prepared, settings)
        self._prepared[key] = result
        while len(self._prepared) > 3:
            self._prepared.popitem(last=False)
        return result

    def evaluate(self, task, controls, **options):
        index = task["bin"]
        prepared, settings = self.prepare(index, controls, **options)
        settings = {**settings, "controls": dict(controls)}
        active = [self.registry.ids.index(n) for n in task["parameters"]]
        derivatives = evaluate_derivatives(
            prepared.p3d,
            prepared.theta,
            self.z(index),
            prepared.k,
            prepared.mu,
            step_scale=controls["step"] / 0.001,
        )
        selected = self.selection.selected_to_required
        j = (prepared.products[:, :, None] * derivatives.jacobian)[:, selected][
            :, :, active
        ]
        arrays, report = _assemble(
            task,
            prepared.total,
            j,
            settings,
            factors=prepared.factors,
            extra=dict(
                intrinsic=prepared.power,
                response=prepared.response,
                noise=prepared.noise,
                observed_signal=prepared.products * prepared.power,
            ),
        )
        # Independent direct NumPy solve on separately reconstructed Wick C.
        # This never calls production factorization or reuses its contraction.
        f, single = contract(arrays["selected_covariance"], j, independent=True)
        if (
            relative(f, arrays["fisher"]) > 5e-12
            or relative(single, arrays["pair_fisher"]) > 5e-12
        ):
            raise ValueError("direct NumPy and evidence Fisher disagree")
        report["provenance"] = self.provenance
        return arrays, report

    def study(self, task):
        """Run the bounded controller with this prepared scientific recipe."""
        from .study import study

        return study(self, task, payload_factory=OrderedDict)

    def sensitivity(self, task, controls, policy):
        """Recompute physical derived inputs for one explicit artificial policy test."""
        options = {}
        if policy == "growth_eds":
            options = dict(growth="eds")
        elif policy == "legacy_resolution":
            options = dict(resolution="legacy")
        elif policy == "rectangular_magnitude":
            options = dict(rectangular=True)
        elif policy == "no_reconstruction":
            options = dict(reconstruction=False)
        elif policy == "arithmetic_mean":
            options = dict(arithmetic_mean=True)
        elif policy == "three_weights":
            if self.weight_method != "legacy":
                raise ValueError(
                    "three_weights is a legacy cumulative diagnostic and is "
                    "inapplicable to inverse_variance accuracy"
                )
            controls = {**controls, "iterations": 3}
        a, r = self.evaluate(
            task,
            controls,
            policy="primary"
            if policy
            in (
                "growth_eds",
                "legacy_resolution",
                "rectangular_magnitude",
                "no_reconstruction",
                "arithmetic_mean",
                "three_weights",
            )
            else policy,
            **options,
        )
        r.update(
            sensitivity=dict(
                policy=policy,
                weight_method=self.weight_method,
                primary_floor=1e-20,
                diagnostic_floor={"floor_low": 1e-22, "floor_high": 1e-18}.get(
                    policy, 1e-20
                ),
                cell_width_factor={"width_minus": 0.9, "width_plus": 1.1}.get(
                    policy, 1.0
                ),
                artificial=True,
                support_changed=policy in ("remove_bright", "remove_sentinel"),
                interpretation="Not a calibrated systematic error or alternative primary survey",
            )
        )
        return a, r
