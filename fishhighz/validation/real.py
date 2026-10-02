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
        """Load one explicit DESI-2 recipe and its read-only reference inputs.

        Parameters
        ----------
        reference_root : str or pathlib.Path
            Root of the verified lyaforecast reference checkout.
        template_path : str or pathlib.Path
            Path to the Vega-format K/PK/PKSB FITS template.
        case : str
            Identifier of one of the seven original DESI-2 validation
            configurations.
        negative_policy : str
            Explicit treatment of negative interpolated source densities.
        bin_indices : sequence of int or None
            Zero-based bins to include; None uses the suite-defined bin selection.

        Notes
        -----
        Evaluates the reference cosmology, loads source readers and the supplied template, and records source provenance.
        """
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
        evaluation_redshifts = sorted(set([1.8, 4.5, *self.z.values()]))
        self.cosmo = CosmoCamb(
            str(self.root / "lyaforecast/resources/camb_configs/Planck18.ini"),
            z_ref=2.3,
            z_centres=evaluation_redshifts,
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
            tracer_settings = self.config[key]
            tracer = Tracer(tracer_settings)
            self.tracers[field.id] = tracer
            if tracer.bias_func is not None:
                self.external.bias.set_density_bias_func(
                    tracer.simple_name, tracer.bias_func
                )
            path = self.root / "lyaforecast/resources/data" / tracer_settings["dn dz"]
            self.resources.append(path)
            forest = field.kind == "forest"
            reader = DensityReader(
                path,
                semantics="cell_count_per_deg2",
                interpolation="spline",
                target_density=tracer_settings.getfloat("target density"),
                z_norm_min=2.15 if forest or field.id == "qso" else None,
                magnitude_bounds=(
                    tracer_settings.getfloat("min_band_mag"),
                    tracer_settings.getfloat("max_band_mag"),
                )
                if forest
                else None,
                width_policy="legacy_first_spacing",
                label=field.id,
            )
            self.densities[field.id] = LegacyDensity(reader, negative_policy)
            if forest:
                paths = sorted(
                    (
                        self.root
                        / "lyaforecast/resources/data"
                        / tracer_settings["snr-file-dir"]
                    ).glob("*.dat")
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
        """Reprepare changed grids, retaining raw readers and physical input policies.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.
        k_intervals : int
            Number of equal comoving k intervals between 0.01 and 0.5 h/Mpc, each
            with four Gaussian nodes. Default is ``128``.
        mu_order : int
            Gauss-Legendre order on the dimensionless interval 0 <= mu <= 1. Default
            is ``32``.
        z_order : int
            Gauss-Legendre order for integrating the comoving bin volume. Default is
            ``32``.

        Returns
        -------
        prepared : PreparedBin
            Fiducial powers, responses, fixed source weights, noise and covariance.
        settings : dict
            Redshift, quadrature controls, physical parameters and sampled-input
            metadata.

        Notes
        -----
        Caches source samples by bin while preparing the requested Fourier and volume quadrature.
        """
        lo, hi = bins(self.case)[index]
        forest_redshift = self.z[index]
        cosmology = self.cosmo
        geometry = prepare_geometry(
            lo,
            hi,
            z_eval=forest_redshift,
            area_deg2=self.config["survey"].getfloat("survey_area"),
            h_fid=self.h,
            z_order=z_order,
            hubble=cosmology.results.hubble_parameter,
            transverse_distance=cosmology.results.comoving_radial_distance,
        )

        grid = gauss_legendre_grid(
            np.linspace(0.01, 0.5, k_intervals + 1),
            k_order=4,
            mu_order=mu_order,
            h_fid=self.h,
        )

        fields = self.selection.fields
        growth_rate = float(self.external.bias._growth_rate_func(forest_redshift))
        ratio = float(
            np.interp(forest_redshift, cosmology.z_bins, cosmology.sigma8_zbins)
            / cosmology.sigma8
        )
        biases, betas, widths = {}, {}, {}
        for field in fields:
            name = (
                field.physical_tracer
                if hasattr(field, "physical_tracer")
                else self.tracers[field.id].simple_name
            )
            biases[field.id] = float(
                self.external.bias._get_density_bias(forest_redshift, name)
            )
            if field.kind == "forest":
                betas[field.id] = float(
                    self.external.bias._get_beta_rsd(forest_redshift, name)
                )
            reconstruction = (
                1
                if field.kind == "forest"
                else self.config["survey"].getfloat("reconstruction factor")
            )
            transverse = 3.26 * ratio / np.sqrt(reconstruction)
            widths[field.id] = ((1 + growth_rate) * transverse, transverse)

        growth = ((1 + self.template.z_ref) / (1 + forest_redshift)) ** 2
        model = KaiserModel(
            self.template,
            fields,
            biases=biases,
            betas=betas,
            widths=widths,
            f=growth_rate
            if any(tracer.kind == "galaxy" for tracer in fields)
            else None,
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
            z=forest_redshift,
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

        observed_wavelength = LYA_REST_ANGSTROM * (1 + forest_redshift)
        responses = {}
        for field in fields:
            tracer = self.tracers[field.id]
            responses[field.id] = (
                InstrumentResponse(
                    pixel_width_angstrom_to_velocity(
                        tracer.pix_ang, lambda_obs_angstrom=observed_wavelength
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
                tracer = self.tracers[field.id]
                if field.kind == "forest":
                    source_z = (
                        observed_wavelength / np.sqrt(tracer.lrmin * tracer.lrmax) - 1
                    )
                    sampled = sample_legacy_forest(
                        self.densities[field.id],
                        self.snrs[field.id],
                        geometry,
                        responses[field.id],
                        z_source=source_z,
                        magnitudes=magnitudes,
                        pixel_width_angstrom=tracer.pix_ang,
                        exposure_count=tracer.num_exp,
                    )
                    forests[field.id] = ForestInput(
                        dict(
                            z_source=source_z,
                            magnitudes=magnitudes,
                            quadrature=quadrature,
                            rho=sampled["rho"],
                            variance=sampled["variance"],
                            length_velocity=SPEED_LIGHT_KMS
                            * np.log(tracer.lrmax / tracer.lrmin),
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
                    sampled = self.densities[field.id].sample(
                        forest_redshift, magnitudes
                    )
                    galaxies[field.id] = local_galaxy_density(
                        sampled["values"], quadrature, geometry
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
            bin=index,
            bounds=[lo, hi],
            z_eval=forest_redshift,
            arithmetic_label=(lo + hi) / 2,
            k_intervals=k_intervals,
            k_order=4,
            k_domain=[0.01, 0.5],
            mu_order=mu_order,
            z_order=z_order,
            biases=biases,
            betas=betas,
            f=growth_rate,
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
            volume=geometry.volume,
            modes_normalization="V*k^2*w_k*w_mu/(2*pi^2)",
        )
        return prepared, plain(settings)

    def study(self, task):
        """One-bin independent k/mu/z/step refinements with explicit pass metrics.

        Parameters
        ----------
        task : dict
            Declared case, bin, selected field pairs, parameter order and validation
            thresholds.

        Returns
        -------
        arrays : dict of str to ndarray
            Fiducial and refined Fisher/error arrays plus direct covariance and
            derivative evidence.
        report : dict
            Independent reconstruction checks, finite-refinement metrics and input
            provenance.

        Notes
        -----
        Runs only the fixed one-bin k, mu, volume and derivative-step variants listed in this method.
        """
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
            lower_fisher, upper_fisher = (
                arrays[first + "_fisher"],
                arrays[last + "_fisher"],
            )
            fisher_relative_change = float(
                np.linalg.norm(lower_fisher - upper_fisher)
                / np.linalg.norm(upper_fisher)
            )
            error_relative_change = float(
                np.max(abs(arrays[first + "_errors"] / arrays[last + "_errors"] - 1))
            )
            volume = abs(reports[first]["volume"] / reports[last]["volume"] - 1)
            metrics[axis] = dict(
                fisher_relative=fisher_relative_change,
                error_relative=error_relative_change,
                volume_relative=volume,
                passed=fisher_relative_change <= 1e-3
                and error_relative_change <= 5e-3
                and volume <= 1e-6,
            )

        derivative = evaluate_derivatives(
            base.p3d, base.theta, base.geometry.z_eval, base.k, base.mu, step_scale=0.25
        )
        # Independent Gaussian covariance construction from the full field matrix.
        n_fields = len(self.selection.fields)
        matrix_total = np.zeros((len(base.k), n_fields, n_fields))
        for pair, column in zip(self.selection.required_pairs, base.total.T):
            i, j = pair
            matrix_total[:, i, j] = matrix_total[:, j, i] = column
        selected = self.selection.selected_pairs
        independent_covariance = np.empty((len(base.k), len(selected), len(selected)))
        for a, (i, j) in enumerate(selected):
            for b, (m, n) in enumerate(selected):
                independent_covariance[:, a, b] = (
                    matrix_total[:, i, m] * matrix_total[:, j, n]
                    + matrix_total[:, i, n] * matrix_total[:, j, m]
                ) / base.modes
        observed_jacobian = (base.products[:, :, None] * derivative.jacobian)[
            :, self.selection.selected_to_required
        ][:, :, active]
        solved = np.linalg.solve(independent_covariance, observed_jacobian)
        independent_fisher = np.einsum("nsi,nsj->ij", observed_jacobian, solved)
        np.testing.assert_allclose(
            independent_covariance,
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
                    independent_covariance
                    - base.factors @ base.factors.swapaxes(-1, -2)
                )
                / np.linalg.norm(independent_covariance)
            ),
            fisher_relative=float(
                np.linalg.norm(independent_fisher - arrays["step_quarter_fisher"])
                / np.linalg.norm(independent_fisher)
            ),
        )
        arrays["independent_covariance"] = independent_covariance
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
        """Actual intrinsic object plus A wrapper, independent full-field trace oracle.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.

        Returns
        -------
        comparison : dict
            Scalar amplitude information from both contractions, their relative
            discrepancy and direct external-spectrum checks.

        Notes
        -----
        Uses a tiny grid with the same fixed response and noise. Asserts agreement with the full-field Gaussian trace expression.
        """
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
            """Apply a single amplitude parameter to the external intrinsic spectra.

            Parameters
            ----------
            theta : array_like, shape (n_parameter,)
                Model parameter values in the declared local parameter order.
            z : float
                Dimensionless evaluation redshift.
            k : array_like
                Comoving Fourier wavenumbers in h/Mpc; array shape follows the model or
                paired grid.
            mu : array_like
                Dimensionless line-of-sight direction cosines aligned with the Fourier
                grid.
            pairs : array_like, shape (n_pair, 2)
                Ordered pairs of integer field indices.

            Returns
            -------
            power : ndarray, shape (n_cell, n_pair)
                External intrinsic spectra in (Mpc/h)^3 multiplied by the dimensionless
                amplitude.
            """
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
        derivatives = evaluate_derivatives(
            p3d, registry.fiducials, base.geometry.z_eval, base.k, base.mu
        )
        fisher = fisher_from_factors(
            (base.products[:, :, None] * derivatives.jacobian)[
                :, self.selection.selected_to_required
            ],
            factor_covariance(covariance),
        )[0, 0]
        oracle = 0.0
        n_fields = len(self.selection.fields)
        for row, noise, modes in zip(signal, base.noise, base.modes):
            signal_matrix = np.zeros((n_fields, n_fields))
            total_matrix = np.zeros((n_fields, n_fields))
            for (i, j), power, noise_power in zip(
                self.selection.required_pairs, row, noise
            ):
                signal_matrix[i, j] = signal_matrix[j, i] = power
                total_matrix[i, j] = total_matrix[j, i] = power + noise_power
            inverse_weighted_signal = np.linalg.solve(total_matrix, signal_matrix)
            oracle += (
                modes * np.trace(inverse_weighted_signal @ inverse_weighted_signal) / 2
            )
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
