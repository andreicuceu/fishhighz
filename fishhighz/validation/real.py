"""Opt-in local DESI-2 recipes; reference imports occur only in caller setup.

No NewForecast invocation, INI translation, raw asset changes or CAMB rebuilds
inside forecast loops. All original settings are verified against explicit recipes.
"""

import configparser
import importlib.metadata
import sys
from pathlib import Path

import numpy as np

from ..adapters.legacy_compat import (
    LegacyDensity,
    LegacySNR,
    plain,
    sample_legacy_forest,
)
from ..adapters.legacy_inputs import DensityReader, SNRReader
from ..adapters.lyaforecast import IntrinsicP3D
from ..derivatives import evaluate_derivatives
from ..forecast import prepare_bin, run_bin
from ..geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS, prepare_geometry
from ..grids import gauss_legendre_grid
from ..models.external import BoundParameters, P3DProvider, PreparedP3D
from ..models.kaiser import KaiserModel, Scaling
from ..models.p1d import default_p1d
from ..models.templates import load_template
from ..noise import local_galaxy_density
from ..parameters import Parameter, ParameterRegistry
from ..response import InstrumentResponse, pixel_width_angstrom_to_velocity
from ..survey import BinSpec, ForestInput
from .cases import bins, recipe, selection, verify_inventory
from .evidence import digest


class RealRecipe:
    """Prepare explicit case assets/background once, then sample fixed bin inputs.

    reference_root and template_path are explicit external read-only paths.
    The current Python must already provide lyaforecast/CAMB; no installation or
    neighboring path injection occurs here. negative_policy is mandatory.
    """

    def __init__(
        self, reference_root, template_path, case, *, negative_policy, bin_indices
    ):
        import camb
        import lyaforecast
        from lyaforecast.cosmoCAMB import CosmoCamb
        from lyaforecast.power_spectrum import PowerSpectrum
        from lyaforecast.tracer import Tracer

        self.root = Path(reference_root).resolve()
        self.inventory = verify_inventory(self.root / "examples/desi2")
        self.case, self.settings = case, recipe(case)
        self.selection = selection(case)
        self.config = configparser.ConfigParser()
        self.config.read_dict(self.settings)
        self.bin_indices = tuple(bin_indices)
        bounds = bins(case)
        self.z = {
            i: np.sqrt((1 + bounds[i][0]) * (1 + bounds[i][1])) - 1
            for i in self.bin_indices
        }
        # Two or more explicit background nodes for reference linear f interpolation.
        zs = sorted(set([1.8, 4.5, *self.z.values()]))
        self.cosmo = CosmoCamb(
            str(self.root / "lyaforecast/resources/camb_configs/Planck18.ini"),
            z_ref=2.3,
            z_centres=zs,
        )
        self.h = self.cosmo._pars.H0 / 100
        self.template = load_template(template_path, h_fid=self.h)
        self.external = PowerSpectrum(self.config, self.cosmo, {})
        self.tracers, self.densities, self.snrs = {}, {}, {}
        self.resources = [
            Path(template_path),
            self.root / "lyaforecast/resources/camb_configs/Planck18.ini",
        ]
        for field, key in zip(
            self.selection.fields, [s for s in self.settings if s.startswith("tracer ")]
        ):
            t = self.config[key]
            tracer = Tracer(t)
            self.tracers[field.id] = tracer
            if tracer.bias_func is not None:
                self.external.bias.set_density_bias_func(
                    tracer.simple_name, tracer.bias_func
                )
            path = self.root / "lyaforecast/resources/data" / t["dn dz"]
            self.resources.append(path)
            forest = field.kind == "forest"
            reader = DensityReader(
                path,
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
            self.densities[field.id] = LegacyDensity(reader, negative_policy)
            if forest:
                paths = sorted(
                    (self.root / "lyaforecast/resources/data" / t["snr-file-dir"]).glob(
                        "*.dat"
                    )
                )
                self.resources.extend(paths)
                self.snrs[field.id] = LegacySNR(
                    SNRReader(paths, smoothing="legacy", label=field.id)
                )
        self.registry = ParameterRegistry(
            [
                Parameter(f"{name}_{i}", 1, "target", step=1e-3)
                for i in self.bin_indices
                for name in ("ap", "at")
            ]
        )
        self.provenance = dict(
            reference_origin=lyaforecast.__file__,
            camb_origin=camb.__file__,
            python=sys.executable,
            versions={
                n: importlib.metadata.version(n)
                for n in ("numpy", "scipy", "astropy", "camb", "lyaforecast")
            },
            resources={str(p.resolve()): digest(p) for p in self.resources},
            sources={
                str(p.resolve()): digest(p)
                for p in (self.root / "lyaforecast").glob("*.py")
            },
            original=self.settings,
            h_fid=self.h,
            template=dict(
                path=self.template.source_path,
                sha256=self.template.source_sha256,
                z_ref=self.template.z_ref,
                h_template=self.template.h_template,
                metadata=plain(self.template.metadata),
            ),
            negative_policy=negative_policy,
            compatibility=True,
            cosmology_note="Vega template header and CAMB Planck18.ini are distinct inputs; no equivalence assumed",
            bin_registry=list(self.registry.ids),
        )
        self._samples = {}

    def prepare(self, index, *, k_intervals=128, mu_order=32, z_order=32):
        """Reprepare changed grids, retaining raw readers and physical input policies."""
        lo, hi = bins(self.case)[index]
        z = self.z[index]
        c = self.cosmo
        g = prepare_geometry(
            lo,
            hi,
            z_eval=z,
            area_deg2=self.config["survey"].getfloat("survey_area"),
            h_fid=self.h,
            z_order=z_order,
            hubble=c.results.hubble_parameter,
            transverse_distance=c.results.comoving_radial_distance,
        )
        grid = gauss_legendre_grid(
            np.linspace(0.01, 0.5, k_intervals + 1),
            k_order=4,
            mu_order=mu_order,
            h_fid=self.h,
        )
        fields = self.selection.fields
        f = float(self.external.bias._growth_rate_func(z))
        ratio = float(np.interp(z, c.z_bins, c.sigma8_zbins) / c.sigma8)
        biases, betas, widths = {}, {}, {}
        for field in fields:
            name = (
                field.physical_tracer
                if hasattr(field, "physical_tracer")
                else self.tracers[field.id].simple_name
            )
            biases[field.id] = float(self.external.bias._get_density_bias(z, name))
            if field.kind == "forest":
                betas[field.id] = float(self.external.bias._get_beta_rsd(z, name))
            reconstruction = (
                1
                if field.kind == "forest"
                else self.config["survey"].getfloat("reconstruction factor")
            )
            transverse = 3.26 * ratio / np.sqrt(reconstruction)
            widths[field.id] = ((1 + f) * transverse, transverse)
        growth = ((1 + self.template.z_ref) / (1 + z)) ** 2
        model = KaiserModel(
            self.template,
            fields,
            biases=biases,
            betas=betas,
            widths=widths,
            f=f if any(t.kind == "galaxy" for t in fields) else None,
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
            z=z,
            growth=growth,
        )
        bound = BoundParameters(
            self.registry,
            ("ap", "at"),
            {name: f"{name}_{index}" for name in ("ap", "at")},
        )
        p3d = PreparedP3D(
            self.registry,
            self.selection,
            [P3DProvider("wiggle BAO", model, bound, self.selection.required_pairs)],
        )
        wave = LYA_REST_ANGSTROM * (1 + z)
        responses = {}
        for field in fields:
            t = self.tracers[field.id]
            responses[field.id] = (
                InstrumentResponse(
                    pixel_width_angstrom_to_velocity(
                        t.pix_ang, lambda_obs_angstrom=wave
                    ),
                    SPEED_LIGHT_KMS / self.config["survey"].getfloat("resolution"),
                )
                if field.kind == "forest"
                else InstrumentResponse(0, 0)
            )
        if index not in self._samples:
            magnitudes = np.linspace(
                self.config["survey"].getfloat("min_band_mag"),
                self.config["survey"].getfloat("max_band_mag"),
                self.config["survey"].getint("num mag bins"),
            )
            quadrature = np.full(len(magnitudes), magnitudes[1] - magnitudes[0])
            forests, galaxies, sampled_meta = {}, {}, {}
            for field in fields:
                if fields.index(field) not in np.unique(self.selection.selected_pairs):
                    continue
                t = self.tracers[field.id]
                if field.kind == "forest":
                    source_z = wave / np.sqrt(t.lrmin * t.lrmax) - 1
                    sampled = sample_legacy_forest(
                        self.densities[field.id],
                        self.snrs[field.id],
                        g,
                        responses[field.id],
                        z_source=source_z,
                        magnitudes=magnitudes,
                        pixel_width_angstrom=t.pix_ang,
                        exposure_count=t.num_exp,
                    )
                    forests[field.id] = ForestInput(
                        dict(
                            z_source=source_z,
                            magnitudes=magnitudes,
                            quadrature=quadrature,
                            rho=sampled["rho"],
                            variance=sampled["variance"],
                            length_velocity=SPEED_LIGHT_KMS * np.log(t.lrmax / t.lrmin),
                            method="legacy",
                            iterations=3,
                        ),
                        default_p1d,
                        BoundParameters(self.registry, (), {}),
                        self.registry.fiducials,
                        auxiliary_coordinates=(2.4, 0.00035),
                        provenance=sampled["provenance"],
                    )
                    sampled_meta[field.id] = sampled["provenance"]
                else:
                    sampled = self.densities[field.id].sample(z, magnitudes)
                    galaxies[field.id] = local_galaxy_density(
                        sampled["values"], quadrature, g
                    )
                    sampled_meta[field.id] = plain(sampled["provenance"])
            self._samples[index] = (
                forests,
                galaxies,
                sampled_meta,
                magnitudes,
                quadrature,
            )
        forests, galaxies, metadata, magnitudes, quadrature = self._samples[index]
        prepared = prepare_bin(
            BinSpec(
                f"{self.case}-{index}",
                g,
                grid,
                p3d,
                responses,
                forests=forests,
                galaxies=galaxies,
                independent_sampling=True,
            )
        )
        settings = dict(
            bin=index,
            bounds=[lo, hi],
            z_eval=z,
            arithmetic_label=(lo + hi) / 2,
            k_intervals=k_intervals,
            k_order=4,
            k_domain=[0.01, 0.5],
            mu_order=mu_order,
            z_order=z_order,
            biases=biases,
            betas=betas,
            f=f,
            growth_G=growth,
            sigma8_ratio=ratio,
            widths=widths,
            magnitude_nodes=magnitudes,
            magnitude_weights=quadrature,
            magnitude_rule="rectangular including both endpoints",
            response_convention="explicit legacy c/R velocity sigma",
            c=SPEED_LIGHT_KMS,
            samples=metadata,
            noise_independence=True,
            weights_iterations=3,
            forest_lengths={
                name: float(source.weight_options["length_velocity"])
                for name, source in forests.items()
            },
            galaxies_nbar=galaxies,
            primary="wiggle ap/at only; full remapping/Q/RSD/damping, fixed identity smooth; no priors",
            parameter_ids=list(self.registry.ids),
            volume=g.volume,
            modes_normalization="V*k^2*w_k*w_mu/(2*pi^2)",
        )
        return prepared, plain(settings)

    def study(self, task):
        """One-bin independent k/mu/z/step refinements with explicit pass metrics."""
        index = task["bin"]
        variants = [
            ("base", 128, 32, 32, 1.0),
            ("k32", 32, 32, 32, 1.0),
            ("k64", 64, 32, 32, 1.0),
            ("mu8", 128, 8, 32, 1.0),
            ("mu16", 128, 16, 32, 1.0),
            ("z8", 128, 32, 8, 1.0),
            ("z16", 128, 32, 16, 1.0),
            ("step_half", 128, 32, 32, 0.5),
            ("step_quarter", 128, 32, 32, 0.25),
        ]
        arrays, reports = {}, {}
        base, settings = self.prepare(index)
        active = [self.registry.ids.index(f"{n}_{index}") for n in ("ap", "at")]
        for name, k, mu, zorder, step in variants:
            prepared = (
                base
                if name in ("base", "step_half", "step_quarter")
                else self.prepare(index, k_intervals=k, mu_order=mu, z_order=zorder)[0]
            )
            result = run_bin(prepared, batch_size=2048, step_scale=step).result
            matrix = result.data_fisher[np.ix_(active, active)]
            errors = np.sqrt(np.diag(np.linalg.inv(matrix)))
            arrays[name + "_fisher"] = matrix
            arrays[name + "_errors"] = errors
            reports[name] = dict(
                volume=prepared.geometry.volume,
                conditioning=float(np.linalg.cond(matrix)),
                grid=[k, 4, mu, zorder],
                step=step * 1e-3,
            )
        metrics = {}
        for axis, first, last in [
            ("k", "k64", "base"),
            ("mu", "mu16", "base"),
            ("z", "z16", "base"),
            ("step", "step_half", "step_quarter"),
        ]:
            a, b = arrays[first + "_fisher"], arrays[last + "_fisher"]
            ferr = float(np.linalg.norm(a - b) / np.linalg.norm(b))
            eerr = float(
                np.max(abs(arrays[first + "_errors"] / arrays[last + "_errors"] - 1))
            )
            volume = abs(reports[first]["volume"] / reports[last]["volume"] - 1)
            metrics[axis] = dict(
                fisher_relative=ferr,
                error_relative=eerr,
                volume_relative=volume,
                passed=ferr <= 1e-3 and eerr <= 5e-3 and volume <= 1e-6,
            )
        derivative = evaluate_derivatives(
            base.p3d, base.theta, base.geometry.z_eval, base.k, base.mu, step_scale=0.25
        )
        # Independent Gaussian covariance construction from the full field matrix.
        n = len(self.selection.fields)
        matrix_total = np.zeros((len(base.k), n, n))
        for pair, column in zip(self.selection.required_pairs, base.total.T):
            i, j = pair
            matrix_total[:, i, j] = matrix_total[:, j, i] = column
        selected = self.selection.selected_pairs
        independent_cov = np.empty((len(base.k), len(selected), len(selected)))
        for a, (i, j) in enumerate(selected):
            for b, (m, n) in enumerate(selected):
                independent_cov[:, a, b] = (
                    matrix_total[:, i, m] * matrix_total[:, j, n]
                    + matrix_total[:, i, n] * matrix_total[:, j, m]
                ) / base.modes
        observed_jac = (base.products[:, :, None] * derivative.jacobian)[
            :, self.selection.selected_to_required
        ][:, :, active]
        solved = np.linalg.solve(independent_cov, observed_jac)
        independent_fisher = np.einsum("nsi,nsj->ij", observed_jac, solved)
        np.testing.assert_allclose(
            independent_cov,
            base.factors @ base.factors.swapaxes(-1, -2),
            rtol=5e-12,
            atol=0,
        )
        np.testing.assert_allclose(
            independent_fisher, arrays["step_quarter_fisher"], rtol=5e-12, atol=0
        )
        independent_metrics = dict(
            covariance_relative=float(
                np.linalg.norm(
                    independent_cov - base.factors @ base.factors.swapaxes(-1, -2)
                )
                / np.linalg.norm(independent_cov)
            ),
            fisher_relative=float(
                np.linalg.norm(independent_fisher - arrays["step_quarter_fisher"])
                / np.linalg.norm(independent_fisher)
            ),
        )
        arrays["independent_covariance"] = independent_cov
        arrays["independent_fisher"] = independent_fisher
        arrays.update(
            fisher=arrays["step_quarter_fisher"],
            errors=arrays["step_quarter_errors"],
            k=base.k,
            mu=base.mu,
            modes=base.modes,
            response=base.response,
            noise=base.noise,
            intrinsic=base.power,
            observed_signal=base.products * base.power,
            total=base.total,
            factors=base.factors,
            jacobian=derivative.jacobian,
            selected_pairs=self.selection.selected_pairs,
            required_pairs=self.selection.required_pairs,
        )
        report = dict(
            independent_metrics=independent_metrics,
            settings=settings,
            refinements=reports,
            metrics=metrics,
            provenance=self.provenance,
            passed=all(m["passed"] for m in metrics.values()),
            coverage="one selected real bin; no magnitude or weight convergence claim",
        )
        return arrays, report

    def amplitude(self, index):
        """Actual intrinsic object plus A wrapper, independent full-field trace oracle."""
        registry = ParameterRegistry([Parameter("A", 1, "target", step=0.001)])
        base, _ = self.prepare(index, k_intervals=2, mu_order=3, z_order=8)
        routes = {
            (i, j): f"{self.selection.fields[i].id}_{self.selection.fields[j].id}"
            for i, j in self.selection.required_pairs
        }
        bridge = IntrinsicP3D(
            self.external,
            routes=routes,
            n_fields=len(self.selection.fields),
            h_source=self.h,
            h_fid=self.h,
            k_domain=(0.01, 0.5),
            z_domain=(1.8, 4.5),
        )

        def model(theta, z, k, mu, pairs):
            return theta[0] * bridge([], z, k, mu, pairs)

        bound = BoundParameters(registry, ("A",), {"A": "A"})
        p3d = PreparedP3D(
            registry,
            self.selection,
            [
                P3DProvider(
                    "external amplitude", model, bound, self.selection.required_pairs
                )
            ],
        )
        intrinsic = bridge(
            [], base.geometry.z_eval, base.k, base.mu, self.selection.required_pairs
        )
        direct = np.column_stack(
            [
                self.external.compute_p3d_hmpc(
                    base.geometry.z_eval, base.k, base.mu, routes[tuple(p)]
                )
                for p in self.selection.required_pairs
            ]
        )
        np.testing.assert_allclose(intrinsic, direct, rtol=5e-12, atol=0)
        np.testing.assert_array_equal(
            intrinsic,
            bridge(
                [], base.geometry.z_eval, base.k, base.mu, self.selection.required_pairs
            ),
        )
        # Reuse identical supplied per-field response and noise; P1D is unchanged.
        from ..covariance import gaussian_covariance
        from ..fisher import factor_covariance, fisher_from_factors

        signal = base.products * intrinsic
        covariance = gaussian_covariance(
            signal + base.noise, base.modes, self.selection
        )
        deriv = evaluate_derivatives(
            p3d, registry.fiducials, base.geometry.z_eval, base.k, base.mu
        )
        fisher = fisher_from_factors(
            (base.products[:, :, None] * deriv.jacobian)[
                :, self.selection.selected_to_required
            ],
            factor_covariance(covariance),
        )[0, 0]
        oracle = 0.0
        n = len(self.selection.fields)
        for row, noise, modes in zip(signal, base.noise, base.modes):
            s = np.zeros((n, n))
            t = np.zeros((n, n))
            for (i, j), power, nn in zip(self.selection.required_pairs, row, noise):
                s[i, j] = s[j, i] = power
                t[i, j] = t[j, i] = power + nn
            x = np.linalg.solve(t, s)
            oracle += modes * np.trace(x @ x) / 2
        np.testing.assert_allclose(fisher, oracle, rtol=5e-12, atol=0)
        return dict(
            fisher=float(fisher),
            oracle=float(oracle),
            relative=float(abs(fisher / oracle - 1)),
            direct_max_relative=float(np.max(abs(intrinsic / direct - 1))),
            response="each field applied once, fixed noise reused",
            p1d="independent default_p1d unchanged",
            setup=self.provenance,
            amplitude_only=True,
        )
