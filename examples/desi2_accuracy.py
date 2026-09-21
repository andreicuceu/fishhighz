#!/usr/bin/env python3
"""Run the six-bin DESI Run-2 15x2pt forecast with the accuracy recipe.

This tutorial keeps the FishHighz calculation visible.  A small local input
class reads the original lyaforecast survey tables; fields, tracer model,
geometry, Fourier quadrature, noise, covariance, derivatives, and Fisher
contractions are assembled here through the public FishHighz API.

Example
-------
Run from an installed FishHighz environment, supplying an installed
lyaforecast checkout and the Vega-format linear-power template.  The environment
needs ``fishhighz[templates,survey]``, lyaforecast, and CAMB::

    python examples/desi2_accuracy.py \
        --reference ../lyaforecast \
        --template /path/to/Planck18.fits \
        --output accuracy-desi2

The output directory must not already exist.  It receives ``settings.json``
and ``results.npz`` and the same individual and joint constraints are printed.
"""

from __future__ import annotations

import argparse
import configparser
import hashlib
from pathlib import Path

import numpy as np

try:
    from ._desi2_results import create_output, print_results, write_results
except ImportError:  # Direct ``python examples/desi2_accuracy.py`` execution.
    from _desi2_results import create_output, print_results, write_results

from fishhighz.adapters.legacy_compat import LegacyDensity, LegacySNR, plain
from fishhighz.adapters.legacy_inputs import DensityReader, SNRReader
from fishhighz.covariance import gaussian_covariance
from fishhighz.derivatives import evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import fisher_matrix
from fishhighz.forecast import prepare_bin, run_forecast
from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import (
    BoundParameters,
    P3DProvider,
    PreparedP3D,
)
from fishhighz.models.kaiser import KaiserModel, Scaling
from fishhighz.models.p1d import default_p1d
from fishhighz.models.templates import load_template
from fishhighz.noise import local_galaxy_density
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import (
    InstrumentResponse,
    pixel_width_angstrom_to_velocity,
    resolving_power_fwhm_to_sigma,
)
from fishhighz.results import FisherResult
from fishhighz.survey import BinSpec, ForestInput
from fishhighz.validation.accuracy import breakpoints, composite
from fishhighz.validation.cases import bins, verify_inventory
from fishhighz.validation.profile_definitions import REFERENCE, REVISION, STOPPING
from fishhighz.validation.reference_capture import (
    imported_reference,
    resolved_resources,
)
from fishhighz.weights import density_per_velocity

CASE = "lya_qso_lbg_lae_15x2pt"

# These are the final S2/S3 controls.  ``k_intervals`` counts subintervals;
# each has four Gauss--Legendre nodes.  The tighter weight tolerance is an
# intentional override of the package study default (1e-4).
CONTROLS = {
    "k_intervals": 128,
    "k_order": 4,
    "k_min": 0.01,
    "k_max": 0.5,
    "mu_order": 32,
    "z_order": 32,
    "magnitude_order": 16,
    "ap_step": 2.5e-4,
    "weight_rtol": 1e-5,
}


class SurveyInputs:
    """Example-specific owner of raw DESI-2 density and SNR adapters."""

    def __init__(self, reference, config, cosmology, template):
        from lyaforecast.power_spectrum import PowerSpectrum
        from scipy.interpolate import interp1d

        self.reference = Path(reference)
        self.config = config
        self.cosmology = cosmology
        self.template = template

        # ObservedField separates the two forest samples from their shared Ly-a
        # physical model.  Their backgrounds, noise, and responses remain
        # distinct.  The original field order fixes every pair index below.
        self.fields = (
            ObservedField("lya(qso)", "forest", "lya", background="qso"),
            ObservedField("qso", "galaxy", "qso"),
            ObservedField("lbg", "galaxy", "lbg"),
            ObservedField("lae", "galaxy", "lae"),
            ObservedField("lya(lbg)", "forest", "lya", background="lbg"),
        )

        # One global registry contains an independent AP pair for each bin.
        # Its object identity is shared by every P3D and P1D binding.
        self.registry = ParameterRegistry(
            [
                Parameter(f"{name}_{index}", 1.0, "target", step=0.001)
                for index in range(6)
                for name in ("ap", "at")
            ]
        )

        # The lyaforecast PowerSpectrum object is used only to reproduce its
        # tabulated bias and beta evolution.  FishHighz owns the actual P3D,
        # response, covariance, derivatives, and Fisher calculation.
        self.legacy_power = PowerSpectrum(config, cosmology, {})
        self.tracers = {}
        self.densities = {}
        self.snrs = {}
        tracer_sections = [s for s in config.sections() if s.startswith("tracer ")]
        for field, section in zip(self.fields, tracer_sections):
            tracer = config[section]
            self.tracers[field.id] = dict(tracer)
            if "bias z" in tracer:
                self.legacy_power.bias.set_density_bias_func(
                    tracer["tracer"],
                    interp1d(
                        np.fromstring(tracer["bias z"], sep=" "),
                        np.fromstring(tracer["bias val"], sep=" "),
                        bounds_error=False,
                        fill_value="extrapolate",
                    ),
                )

            # DensityReader validates the raw rectangular table and applies the
            # declared target normalization.  The legacy-first-spacing policy
            # is explicit because the historical files lack physical cell edges.
            forest = field.kind == "forest"
            density_reader = DensityReader(
                reference / "lyaforecast/resources/data" / tracer["dn dz"],
                semantics="cell_count_per_deg2",
                target_density=tracer.getfloat("target density"),
                z_norm_min=2.15 if forest or field.id == "qso" else None,
                magnitude_bounds=(
                    tracer.getfloat("min_band_mag"),
                    tracer.getfloat("max_band_mag"),
                )
                if forest
                else None,
                width_policy="legacy_first_spacing",
                label=field.id,
            )
            # LegacyDensity makes spline extension and the accuracy-profile
            # negative-interpolant floor explicit; strict reader defaults are
            # not changed globally.
            self.densities[field.id] = LegacyDensity(density_reader, "floor_negative")
            if forest:
                snr_paths = sorted(
                    (
                        reference
                        / "lyaforecast/resources/data"
                        / tracer["snr-file-dir"]
                    ).glob("*.dat")
                )
                # SNRReader reconstructs the original smoothed DESI tables;
                # LegacySNR labels the bright clamp and out-of-range sentinel.
                self.snrs[field.id] = LegacySNR(
                    SNRReader(snr_paths, smoothing="legacy", label=field.id)
                )
        self.partitions = {}

    def z_eval(self, index):
        z_min, z_max = bins(CASE)[index]
        return float(np.sqrt((1 + z_min) * (1 + z_max)) - 1)

    def selection(self, index):
        # PairSelection builds the full Wick closure even when only three bin-1
        # means are selected.  Bins 2--6 request all 15 upper-triangle pairs.
        selected = None
        if index == 0:
            selected = (
                ("lya(qso)", "lya(qso)"),
                ("lya(qso)", "qso"),
                ("qso", "qso"),
            )
        return PairSelection(self.fields, selected)

    def growth(self, redshift):
        matches = np.flatnonzero(self.cosmology.z_bins == redshift)
        if len(matches) != 1:
            raise ValueError("CAMB must be prepared at the exact requested redshift")
        i = matches[0]
        return (
            float(self.cosmology.sigma8_zbins[i]),
            float(self.cosmology.growth_rate_zbins[i]),
        )

    def p3d(self, index):
        z_eval = self.z_eval(index)
        sigma8, growth_rate = self.growth(z_eval)
        sigma8_template, _ = self.growth(self.template.z_ref)
        power_growth = (sigma8 / sigma8_template) ** 2
        biases, betas, widths = {}, {}, {}
        for field in self.fields:
            tracer = self.tracers[field.id]["tracer"]
            biases[field.id] = float(
                self.legacy_power.bias._get_density_bias(z_eval, tracer)
            )
            if field.kind == "forest":
                betas[field.id] = float(
                    self.legacy_power.bias._get_beta_rsd(z_eval, tracer)
                )
            reconstruction = (
                1.0
                if field.kind == "forest"
                else self.config["survey"].getfloat("reconstruction factor")
            )
            sigma_transverse = (
                3.26 * sigma8 / self.cosmology.sigma8 / np.sqrt(reconstruction)
            )
            widths[field.id] = (
                (1 + growth_rate) * sigma_transverse,
                sigma_transverse,
            )

        # KaiserModel receives every scientific assumption explicitly: forest
        # b(1+beta mu^2), galaxy b+f mu^2, fixed damping widths, template growth,
        # and wiggle-only ap/at scaling.  Smooth Scaling defaults to identity.
        model = KaiserModel(
            self.template,
            self.fields,
            biases=biases,
            betas=betas,
            widths=widths,
            f=growth_rate,
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
            z=z_eval,
            growth=power_growth,
        )
        binding = BoundParameters(
            self.registry,
            ("ap", "at"),
            {name: f"{name}_{index}" for name in ("ap", "at")},
        )
        selection = self.selection(index)
        # P3DProvider declares exact pair ownership and parameter dependencies;
        # PreparedP3D verifies complete covariance closure before any model call.
        p3d = PreparedP3D(
            self.registry,
            selection,
            [P3DProvider("accuracy BAO", model, binding, selection.required_pairs)],
        )
        return p3d, {
            "biases": biases,
            "betas": betas,
            "growth_rate": growth_rate,
            "template_power_growth": power_growth,
            "sigma8": sigma8,
            "sigma8_template": sigma8_template,
            "widths": widths,
            "ap_scaling": "wiggle only",
        }

    def responses(self, index):
        wavelength = LYA_REST_ANGSTROM * (1 + self.z_eval(index))
        result = {}
        for field in self.fields:
            tracer = self.tracers[field.id]
            # InstrumentResponse uses full pixel width and Gaussian one-sigma.
            # The survey resolving power specifies a FWHM, converted once here.
            result[field.id] = (
                InstrumentResponse(
                    pixel_width_angstrom_to_velocity(
                        float(tracer["pix_width_ang"]),
                        lambda_obs_angstrom=wavelength,
                    ),
                    resolving_power_fwhm_to_sigma(
                        float(self.config["survey"]["resolution"])
                    ),
                )
                if field.kind == "forest"
                else InstrumentResponse(0, 0)
            )
        return result

    def samples(self, index):
        z_eval = self.z_eval(index)
        wavelength = LYA_REST_ANGSTROM * (1 + z_eval)
        z_sources = {
            field.id: wavelength
            / np.sqrt(
                float(self.tracers[field.id]["min_rest_frame_lya"])
                * float(self.tracers[field.id]["max_rest_frame_lya"])
            )
            - 1
            if field.kind == "forest"
            else z_eval
            for field in self.fields
        }
        lo = float(self.config["survey"]["min_band_mag"])
        hi = float(self.config["survey"]["max_band_mag"])
        if index not in self.partitions:
            self.partitions[index] = breakpoints(
                self.densities, self.snrs, z_sources, lo, hi
            )
        magnitudes, quadrature = composite(
            self.partitions[index], CONTROLS["magnitude_order"]
        )
        active = np.unique(self.selection(index).selected_pairs)
        sampled = {}
        for field_index, field in enumerate(self.fields):
            if field_index not in active:
                continue
            density = self.densities[field.id].sample(z_sources[field.id], magnitudes)
            row = {
                "magnitudes": magnitudes,
                "quadrature": quadrature,
                "density": np.asarray(density["values"]),
                "z_source": z_sources[field.id],
                "density_diagnostics": plain(density["provenance"]),
            }
            if field.kind == "forest":
                tracer = self.tracers[field.id]
                snr = self.snrs[field.id].sample(
                    z_source=z_sources[field.id],
                    magnitudes=magnitudes,
                    wavelength=wavelength,
                    pixel_width_angstrom=float(tracer["pix_width_ang"]),
                    exposure_count=float(tracer["num exposures"]),
                )
                row.update(
                    variance=np.asarray(snr["values"]),
                    snr_diagnostics=plain(snr["provenance"]),
                    length_velocity=SPEED_LIGHT_KMS
                    * np.log(
                        float(tracer["max_rest_frame_lya"])
                        / float(tracer["min_rest_frame_lya"])
                    ),
                )
            sampled[field.id] = row
        return sampled, magnitudes, quadrature


def _background(reference, template_path):
    """Prepare the supplied lyaforecast CAMB model at every exact required z."""
    import camb
    from astropy.io import fits
    from lyaforecast.cosmoCAMB import CosmoCamb

    reference_identity = imported_reference(reference)
    with fits.open(template_path) as hdus:
        template_redshift = float(hdus[1].header["ZREF"])
    redshifts = sorted(
        {
            template_redshift,
            *[
                float(np.sqrt((1 + z_min) * (1 + z_max)) - 1)
                for z_min, z_max in bins(CASE)
            ],
        }
    )
    ini = reference / "lyaforecast/resources/camb_configs/Planck18.ini"
    parsed = camb.read_ini(str(ini))
    cosmology = CosmoCamb(str(ini), z_ref=2.3, z_centres=redshifts)
    template = load_template(template_path, h_fid=parsed.H0 / 100)
    reference_identity["template"] = {
        "path": str(template_path),
        "sha256": hashlib.sha256(Path(template_path).read_bytes()).hexdigest(),
        "z_ref": template_redshift,
        "h_fid": template.h_fid,
    }
    return cosmology, template, reference_identity


def _forest_input(sample, registry):
    """Translate one sampled forest into the public ``ForestInput`` record."""
    # density_per_velocity converts the sampled angular density into the
    # velocity-density convention required by the forest-weight recurrence.
    weight_options = {
        "z_source": sample["z_source"],
        "magnitudes": sample["magnitudes"],
        "quadrature": sample["quadrature"],
        "rho": density_per_velocity(sample["density"], z_source=sample["z_source"]),
        "variance": sample["variance"],
        "length_velocity": sample["length_velocity"],
        "method": "early_lyaforecast",
        # Adaptive stopping checks amplitude, shape, the full vector, A, and
        # P_pixel.  No capped result or fallback can enter a prepared bin.
        **{**STOPPING, "rtol": CONTROLS["weight_rtol"]},
    }
    # BoundParameters states that default_p1d has no free parameters.  The
    # registry identity must still match the P3D registry used in every bin.
    p1d_parameters = BoundParameters(registry, (), {})
    return ForestInput(
        weight_options,
        default_p1d,
        p1d_parameters,
        registry.fiducials,
        # The representative angular/velocity mode is converted independently
        # in each bin, so its comoving k and mu correctly vary with redshift.
        auxiliary_coordinates=(REFERENCE["k_t_deg"], REFERENCE["k_p_velocity"]),
        provenance={
            "density_policy": "legacy interpolation; negative values floored at 1e-20",
            "snr_policy": "legacy smoothing, bright clamp, and out-of-range sentinel",
            "density": sample["density_diagnostics"],
            "snr": sample["snr_diagnostics"],
        },
    )


def _prepare_bins(reference, template_path):
    """Prepare all fixed fiducial survey state and explanatory settings."""
    # Verify that imports and raw resources come from the requested installed
    # checkout.  The example reads the original 15x2pt INI rather than a captured
    # baseline or historical forecast output.
    inventories = verify_inventory(reference / "examples/desi2")
    config = configparser.ConfigParser()
    config.read(reference / "examples/desi2" / f"{CASE}.ini")
    resources = resolved_resources(reference, config)

    # _background reads Planck18.ini with lyaforecast/CAMB and asks CAMB for
    # every exact geometric-bin redshift plus the template reference redshift.
    # load_template retains the supplied signed Vega K/PK/PKSB decomposition.
    cosmology, template, provenance = _background(reference, template_path)
    provenance["resources"] = resources
    provenance["ini"] = inventories[CASE]
    inputs = SurveyInputs(reference, config, cosmology, template)
    prepared_bins = []
    bin_settings = []
    for index, (z_min, z_max) in enumerate(bins(CASE)):
        # The evaluation redshift is geometric in 1+z.  Forest density/SNR
        # samples use the corresponding foreground/source redshift convention:
        # z_source follows the geometric centre of the forest wavelength range,
        # while discrete tracers are sampled at z_eval.
        z_eval = inputs.z_eval(index)

        # prepare_geometry integrates c*D_M^2/H through the finite bin with 32
        # redshift nodes and records a_v and d_deg for velocity/angular units.
        geometry = prepare_geometry(
            z_min,
            z_max,
            z_eval=z_eval,
            area_deg2=float(config["survey"]["survey_area"]),
            h_fid=template.h_fid,
            z_order=CONTROLS["z_order"],
            hubble=cosmology.results.hubble_parameter,
            transverse_distance=cosmology.results.comoving_radial_distance,
        )

        # gauss_legendre_grid integrates 128 fixed observed-k intervals with
        # order four and mu in [0,1] at order 32.  Its q_mode already includes
        # the Fourier measure k^2 dk dmu/(2*pi^2).
        grid = gauss_legendre_grid(
            np.linspace(
                CONTROLS["k_min"],
                CONTROLS["k_max"],
                CONTROLS["k_intervals"] + 1,
            ),
            k_order=CONTROLS["k_order"],
            mu_order=CONTROLS["mu_order"],
            h_fid=template.h_fid,
        )

        # model returns a PreparedP3D with one KaiserModel provider.  Smooth
        # power remains fixed in observed coordinates; ap_i/at_i rescale the
        # wiggle component only.  Damping widths are fixed per field, and mixed
        # pairs use the mean of the two squared auto widths in KaiserModel.
        p3d, model_settings = inputs.p3d(index)

        # responses converts the instrumental resolving-power FWHM to Gaussian
        # sigma exactly once and the full pixel width from Angstrom to km/s.
        responses = inputs.responses(index)

        # samples evaluates the explicit legacy density/SNR policies on an
        # order-16 composite Gauss--Legendre rule.  Its partition is the union
        # of density spline pieces, support boundaries, SNR magnitude nodes,
        # and negative-interpolant roots; the latter are floored explicitly.
        sampled, magnitude_nodes, magnitude_weights = inputs.samples(index)

        forests = {}
        galaxies = {}
        for field in inputs.fields:
            if field.id not in sampled:
                continue
            sample = sampled[field.id]
            if field.kind == "forest":
                forests[field.id] = _forest_input(sample, inputs.registry)
            else:
                # local_galaxy_density integrates dn/(dz dm) over magnitude and
                # converts the angular surface density to comoving nbar.
                galaxies[field.id] = local_galaxy_density(
                    sample["density"], magnitude_weights, geometry
                )

        # BinSpec is the complete immutable public input contract.  Independent
        # sampling means the per-field forest pixel noise and galaxy shot noise
        # form diagonal noise; all signal cross-powers remain in the covariance.
        spec = BinSpec(
            f"{CASE}-{index + 1}",
            geometry,
            grid,
            p3d,
            responses,
            forests=forests,
            galaxies=galaxies,
            independent_sampling=True,
        )

        # prepare_bin evaluates weights, P1D, the fiducial P3D, response, total
        # powers, mode counts, full selected Gaussian covariance, and its fixed
        # Cholesky factors exactly once.  AP derivatives reuse this fixed state.
        prepared = prepare_bin(spec)
        prepared_bins.append(prepared)
        bin_settings.append(
            {
                "index": index + 1,
                "bounds": [z_min, z_max],
                "z_eval": z_eval,
                "volume": geometry.volume,
                "model": model_settings,
                "selected_pairs": prepared.diagnostics["selected_pairs"],
                "required_pairs": prepared.diagnostics["required_pairs"],
                "magnitude_partition": inputs.partitions[index],
                "magnitude_nodes": magnitude_nodes,
                "magnitude_weights": magnitude_weights,
                "forest_weights": {
                    name: {
                        "convergence": plain(weight.convergence),
                        "A": weight.A,
                        "P_pixel": weight.P_pixel,
                        "representative_k": weight.auxiliary.k,
                        "representative_mu": weight.auxiliary.mu,
                        "P": weight.signal,
                        "B": weight.alias,
                    }
                    for name, weight in prepared.weights.items()
                },
            }
        )
    return inputs, prepared_bins, bin_settings, provenance


def _constraints(result):
    """Use the public result diagnostics, returning an explicit rank status."""
    if result.diagnostics.rank != len(result.registry.ids):
        return "unavailable", None, None, None
    sigma = result.marginalized_errors()
    correlation = result.correlations()[0, 1]
    return "available", float(sigma[0]), float(sigma[1]), float(correlation)


def _forecast_records(inputs, prepared_bins):
    """Run joint bins once and contract every individual selected spectrum."""
    # Every prepared bin shares the same 12-parameter registry: ap_i/at_i are
    # independent between redshift bins.  run_forecast checks that identity,
    # evaluates mean derivatives at step 2.5e-4, and sums independent-bin
    # information.  No prior or covariance derivative is included.
    forecast = run_forecast(
        prepared_bins,
        step_scale=CONTROLS["ap_step"] / 0.001,
        batch_size=2048,
    )
    records = []
    for index, (prepared, bin_run) in enumerate(zip(prepared_bins, forecast.bins)):
        active = [inputs.registry.ids.index(f"{name}_{index}") for name in ("ap", "at")]
        parameter_ids = [f"ap_{index}", f"at_{index}"]
        # Each bin run has zero rows for the other five independent bins.
        # fix_except explicitly removes only those absent-bin AP parameters;
        # marginalized_errors then inverts the retained two-parameter result.
        joint_result = bin_run.result.fix_except(parameter_ids)
        status, sigma_ap, sigma_at, correlation = _constraints(joint_result)
        records.append(
            {
                "bin_index": index,
                "bounds": [prepared.geometry.z_min, prepared.geometry.z_max],
                "redshift": prepared.geometry.z_eval,
                "kind": "joint",
                "pair": None,
                "pair_indices": None,
                "status": status,
                "parameter_ids": parameter_ids,
                "fisher": joint_result.data_fisher,
                "sigma_ap": sigma_ap,
                "sigma_at": sigma_at,
                "correlation": correlation,
            }
        )

        # evaluate_derivatives exposes the required-pair intrinsic Jacobian.
        # Multiplying by the frozen pair response gives derivatives of the
        # observed mean; selected_to_required then restores selected order.
        derivative = evaluate_derivatives(
            prepared.p3d,
            prepared.theta,
            prepared.geometry.z_eval,
            prepared.k,
            prepared.mu,
            step_scale=CONTROLS["ap_step"] / 0.001,
        )
        selected = prepared.p3d.selection.selected_to_required
        observed_jacobian = (prepared.products[:, :, None] * derivative.jacobian)[
            :, selected
        ][:, :, active]

        # gaussian_covariance reconstructs the public Wick covariance from the
        # fixed total powers and mode counts.  Slicing one observable at a time
        # produces its marginal variance, not a block from a joint inverse.
        covariance = gaussian_covariance(
            prepared.total,
            prepared.modes,
            prepared.p3d.selection,
        )
        fields = prepared.p3d.selection.fields
        for pair_index, pair in enumerate(prepared.p3d.selection.selected_pairs):
            individual_fisher = fisher_matrix(
                observed_jacobian[:, pair_index : pair_index + 1],
                covariance[:, pair_index : pair_index + 1, pair_index : pair_index + 1],
            )
            individual_registry = ParameterRegistry(
                [inputs.registry.parameters[i] for i in active]
            )
            individual_result = FisherResult(individual_registry, individual_fisher)
            status, sigma_ap, sigma_at, correlation = _constraints(individual_result)
            records.append(
                {
                    "bin_index": index,
                    "bounds": [prepared.geometry.z_min, prepared.geometry.z_max],
                    "redshift": prepared.geometry.z_eval,
                    "kind": "individual",
                    "pair": [fields[pair[0]].id, fields[pair[1]].id],
                    "pair_indices": [int(pair[0]), int(pair[1])],
                    "status": status,
                    "parameter_ids": parameter_ids,
                    "fisher": individual_result.data_fisher,
                    "sigma_ap": sigma_ap,
                    "sigma_at": sigma_at,
                    "correlation": correlation,
                }
            )
    return records, forecast.combined.data_fisher


def run(*, reference, template, output):
    """Run and save the DESI Run-2 accuracy forecast."""
    reference = Path(reference).resolve(strict=True)
    template = Path(template).resolve(strict=True)
    # Fail immediately if the requested output already exists, before CAMB or
    # any forecast work starts.  An interrupted new run may leave this directory
    # as an explicit record rather than overwriting earlier scientific output.
    output = create_output(output)
    inputs, prepared_bins, bin_settings, provenance = _prepare_bins(reference, template)
    records, combined_fisher = _forecast_records(inputs, prepared_bins)
    full_selection = PairSelection(inputs.fields)
    selected_bin1 = {tuple(pair) for pair in inputs.selection(0).selected_pairs}
    bin1_excluded_pairs = [
        [full_selection.fields[i].id, full_selection.fields[j].id]
        for i, j in full_selection.selected_pairs
        if (int(i), int(j)) not in selected_bin1
    ]
    settings = {
        "case": CASE,
        "profile": "accuracy",
        "recipe_revision": REVISION,
        "selection": "bin1-qso-only; bins2-6-all",
        "reference": str(reference),
        "template": str(template),
        "input_provenance": provenance,
        "controls": CONTROLS,
        "weighting": {
            "method": "early_lyaforecast",
            "reference": REFERENCE,
            "stopping": {**STOPPING, "rtol": CONTROLS["weight_rtol"]},
        },
        "parameter_order": list(inputs.registry.ids),
        "bins": bin_settings,
        "combined_fisher": combined_fisher,
        "bin1_excluded_pairs": bin1_excluded_pairs,
        "result_counts": {"individual": 78, "joint": 6},
    }
    write_results(output, settings=plain(settings), records=records)
    print_results(records)
    return settings, records


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        required=True,
        help="installed lyaforecast checkout containing examples/desi2 and resources",
    )
    parser.add_argument(
        "--template", required=True, help="Vega-format K/PK/PKSB FITS template"
    )
    parser.add_argument(
        "--output", required=True, help="new directory for settings.json/results.npz"
    )
    return parser


if __name__ == "__main__":
    run(**vars(_parser().parse_args()))
