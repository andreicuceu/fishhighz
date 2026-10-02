"""Converged existing-model accuracy profile and explicitly labeled input studies."""

import configparser
from collections import OrderedDict
from pathlib import Path

import numpy as np

from ..accuracy import ADAPTIVE, REFERENCE, REVISION, STOPPING
from ..adapters.legacy_compat import LegacyDensity, LegacySNR, plain
from ..adapters.legacy_inputs import DensityReader, SNRReader
from ..derivatives import evaluate_derivatives
from ..forecast import prepare_bin
from ..geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS, prepare_geometry
from ..grids import gauss_legendre_grid
from ..magnitude import breakpoints, composite
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
from .profile_definitions import forecast_selection, identity
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
    weight_rtol=1e-4,
):
    """Construct the accuracy-profile forest input and declared weight metadata.

    Parameters
    ----------
    field : ObservedField
        Observed forest or galaxy field and its physical tracer identity.
    row : dict
        Sampled density, pixel variance, source redshift, forest length and
        input-policy diagnostics.
    magnitudes : array_like, shape (n_magnitude,)
        Apparent-magnitude quadrature nodes in mag.
    quadrature : array_like, shape (n_magnitude,)
        Magnitude integration weights in mag.
    response : InstrumentResponse
        Per-field pixel width and Gaussian resolution in km/s.
    registry : ParameterRegistry
        Global parameter definitions, fiducials and local bindings.
    z_eval : float
        Dimensionless forest evaluation redshift.
    method : str
        Forest-weight prescription: legacy, inverse_variance, early_lyaforecast
        or mcdonald, as applicable.
    iterations : int
        Number of completed nonlinear weight updates. Default is ``None``.
    policy : str
        Explicit source-input policy or sensitivity variant. Default is
        ``'primary'``.
    weight_rtol : float
        Dimensionless relative tolerance for adaptive forest-weight convergence.
        Default is ``0.0001``.

    Returns
    -------
    source : ForestInput
        Frozen per-field forest-weight and independent P1D inputs.
    weighting : dict
        Declared weight method, applicable iteration controls and auxiliary-
        reference metadata.

    Raises
    ------
    ValueError :
        If inputs, declared identities or numerical validation conditions are
        inconsistent.
    """
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
    elif method in ADAPTIVE:
        options = dict(
            common,
            method=method,
            iterations=iterations,
            **{**STOPPING, "rtol": weight_rtol},
        )
        auxiliary_coordinates = (REFERENCE["k_t_deg"], REFERENCE["k_p_velocity"])
        reference = REFERENCE
        iteration_status = dict(
            applicable=True, stopping={**STOPPING, "rtol": weight_rtol}
        )
    elif method == "legacy":
        if iterations is None:
            raise ValueError("legacy accuracy requires an explicit iteration count")
        options = dict(common, method=method, iterations=iterations)
        auxiliary_coordinates = (2.4, 0.00035)
        iteration_status = {"applicable": True, "count": iterations}
    else:
        raise ValueError(
            "accuracy weight method must be legacy, inverse_variance, "
            "early_lyaforecast or mcdonald"
        )

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
    """One caller-prepared CAMB background at every actual evaluation redshift.

    Parameters
    ----------
    root : str or pathlib.Path
        Root directory of the input checkout or saved evidence bundle.
    template_path : str or pathlib.Path
        Path to the Vega-format K/PK/PKSB FITS template.

    Returns
    -------
    cosmology : CosmoCamb
        CAMB background prepared at all requested geometric/arithmetic bin
        redshifts and the template redshift.
    template : PowerTemplate
        Prepared Vega wiggle/no-wiggle template in comoving units.

    Notes
    -----
    Reads the caller-selected template and verified reference cosmology configuration, then evaluates CAMB once.
    """
    import camb
    from astropy.io import fits
    from lyaforecast.cosmoCAMB import CosmoCamb

    root = Path(root).resolve()
    imported_reference(root)
    with fits.open(template_path) as template_file:
        template_redshift = float(template_file[1].header["ZREF"])
    evaluation_redshifts = sorted(
        {
            2.3,
            template_redshift,
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
    cosmo = CosmoCamb(str(ini), z_ref=2.3, z_centres=evaluation_redshifts)
    template = load_template(template_path, h_fid=parsed.H0 / 100)
    return cosmo, template


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
        weight_method="early_lyaforecast",
        recipe_revision=REVISION,
    ):
        """Prepare a verified case, source readers and fixed parameter bindings.

        Parameters
        ----------
        root : str or pathlib.Path
            Root directory of the input checkout or saved evidence bundle.
        case : str
            Identifier of one of the seven original DESI-2 validation
            configurations.
        cosmo : CosmoCamb
            Background cosmology already evaluated at all required redshifts.
        template : PowerTemplate
            Prepared wiggle and smooth power template in comoving units.
        provenance : dict
            Input paths, hashes and software identities retained with the
            calculation.
        weight_method : str
            Forest-weight prescription: legacy, inverse_variance, early_lyaforecast
            or mcdonald, as applicable. Default is ``'early_lyaforecast'``.
        recipe_revision : str or None
            Declared revision of the physical recipe; None retains historical
            request semantics. Default is ``REVISION``.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.

        Notes
        -----
        Loads density/SNR readers, configures tracer bias interpolation and initializes bounded preparation caches.
        """
        from lyaforecast.power_spectrum import PowerSpectrum
        from scipy.interpolate import interp1d

        self.root = Path(root).resolve()
        verify_inventory(self.root / "examples/desi2")
        self.case = case
        self.config = configparser.ConfigParser()
        self.config.read_dict(recipe(case))
        self.selection = selection(case)
        self.recipe_revision = recipe_revision
        self.cosmo = cosmo
        self.template = template
        self.h = template.h_fid
        self.external = PowerSpectrum(self.config, cosmo, {})
        self.provenance = provenance
        if weight_method not in ("inverse_variance", "legacy", *ADAPTIVE):
            raise ValueError(
                "accuracy weight method must be legacy, inverse_variance, "
                "early_lyaforecast or mcdonald"
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
            tracer_settings = self.config[key]
            self.tracers[field.id] = dict(tracer_settings)
            if "bias z" in tracer_settings:
                self.external.bias.set_density_bias_func(
                    tracer_settings["tracer"],
                    interp1d(
                        np.fromstring(tracer_settings["bias z"], sep=" "),
                        np.fromstring(tracer_settings["bias val"], sep=" "),
                        bounds_error=False,
                        fill_value="extrapolate",
                    ),
                )
            forest = field.kind == "forest"
            reader = DensityReader(
                self.root / "lyaforecast/resources/data" / tracer_settings["dn dz"],
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
            self.densities[field.id] = LegacyDensity(reader, "floor_negative")
            if forest:
                paths = sorted(
                    (
                        self.root
                        / "lyaforecast/resources/data"
                        / tracer_settings["snr-file-dir"]
                    ).glob("*.dat")
                )
                self.snrs[field.id] = LegacySNR(
                    SNRReader(paths, smoothing="legacy", label=field.id)
                )
        self._partitions = {}
        self._samples = {}
        self._prepared = OrderedDict()

    def z(self, index):
        """Return the forest evaluation redshift for one bin.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.

        Returns
        -------
        redshift : float
            Dimensionless geometric mean in 1 + z of the bin edges.
        """
        lo, hi = bins(self.case)[index]
        return float(np.sqrt((1 + lo) * (1 + hi)) - 1)

    def _growth(self, z):
        """Read sigma8 and growth rate at an exactly prepared CAMB redshift.

        Parameters
        ----------
        z : float
            Dimensionless evaluation redshift.

        Returns
        -------
        sigma8 : float
            Linear matter fluctuation amplitude at the exact redshift.
        growth_rate : float
            Dimensionless logarithmic growth rate at the same redshift.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.
        """
        indices = np.flatnonzero(self.cosmo.z_bins == z)
        if len(indices) != 1:
            raise ValueError("CAMB must be prepared at exact requested redshift")
        i = indices[0]
        return float(self.cosmo.sigma8_zbins[i]), float(self.cosmo.growth_rate_zbins[i])

    def model(self, index, *, mean_z=None, growth="camb", reconstruction=True):
        """Construct the intrinsic wiggle-dilation model with fixed fiducial physics.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.
        mean_z : float or None
            Alternative model redshift; None uses the geometric forest redshift.
            Default is ``None``.
        growth : str
            Growth prescription; camb uses the prepared sigma8 ratio, while other
            values use the specified redshift scaling. Default is ``'camb'``.
        reconstruction : bool
            Whether to retain the configured galaxy BAO reconstruction factor.
            Default is ``True``.

        Returns
        -------
        power : PreparedP3D
            Bound BAO model on the bin selection, with fixed per-field damping
            widths.
        settings : dict
            Evaluation redshift, biases, RSD factors, growth normalization and
            damping lengths.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.
        """
        selected = self.bin_selection(index)
        model_redshift = self.z(index) if mean_z is None else float(mean_z)
        sigma8, growth_rate = self._growth(model_redshift)
        sigma_template, _ = self._growth(self.template.z_ref)
        growth_factor = (
            (sigma8 / sigma_template) ** 2
            if growth == "camb"
            else ((1 + self.template.z_ref) / (1 + model_redshift)) ** 2
        )

        # Bias, RSD and damping remain fixed while the BAO dilation varies.
        biases = {}
        betas = {}
        widths = {}
        for field in self.selection.fields:
            name = self.tracers[field.id]["tracer"]
            biases[field.id] = float(
                self.external.bias._get_density_bias(model_redshift, name)
            )
            if field.kind == "forest":
                betas[field.id] = float(
                    self.external.bias._get_beta_rsd(model_redshift, name)
                )
            reconstruction_factor = (
                1
                if field.kind == "forest" or not reconstruction
                else self.config["survey"].getfloat("reconstruction factor")
            )
            sigma_transverse = (
                3.26 * sigma8 / self.cosmo.sigma8 / np.sqrt(reconstruction_factor)
            )
            widths[field.id] = ((1 + growth_rate) * sigma_transverse, sigma_transverse)

        model = KaiserModel(
            self.template,
            self.selection.fields,
            biases=biases,
            betas=betas,
            widths=widths,
            f=growth_rate
            if any(field.kind == "galaxy" for field in self.selection.fields)
            else None,
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
            z=model_redshift,
            growth=growth_factor,
        )

        binding = BoundParameters(
            self.registry, ("ap", "at"), {n: f"{n}_{index}" for n in ("ap", "at")}
        )
        provider = model
        if mean_z is not None:
            # Supplied-mean diagnostic: keep geometric noise/response coordinates
            # while explicitly evaluating the fixed intrinsic model at mean_z.
            def provider(theta, requested_z, k, mu, pairs):
                """Evaluate the mean-redshift diagnostic with fixed geometry coordinates.

                Parameters
                ----------
                theta : array_like, shape (n_parameter,)
                    Model parameter values in the declared local parameter order.
                requested_z : float
                    Dimensionless geometry redshift expected by the prepared provider.
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
                    Intrinsic pair spectra in (Mpc/h)^3 at the separately declared model
                    redshift.

                Raises
                ------
                ValueError :
                    If inputs, declared identities or numerical validation conditions are
                    inconsistent.
                """
                if requested_z != self.z(index):
                    raise ValueError("diagnostic geometry redshift mismatch")
                return model(theta, model_redshift, k, mu, pairs)

        p3d = PreparedP3D(
            self.registry,
            selected,
            [P3DProvider("accuracy BAO", provider, binding, selected.required_pairs)],
        )
        return p3d, dict(
            z_eval=model_redshift,
            biases=biases,
            betas=betas,
            f=growth_rate,
            G=growth_factor,
            growth=growth,
            sigma8=sigma8,
            sigma8_template=sigma_template,
            sigma8_damping_reference=float(self.cosmo.sigma8),
            widths=widths,
        )

    def responses(self, index, *, resolution="fwhm"):
        """Construct each field response at the bin observed wavelength.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.
        resolution : str
            Instrument-resolution convention; fwhm converts the quoted width to a
            Gaussian sigma. Default is ``'fwhm'``.

        Returns
        -------
        responses : dict of str to InstrumentResponse
            Velocity-space pixel and Gaussian widths for each observed field; galaxy
            widths are zero.
        """
        observed_wavelength = LYA_REST_ANGSTROM * (1 + self.z(index))
        responses = {}
        divisor = 2 * np.sqrt(2 * np.log(2)) if resolution == "fwhm" else 1
        for f in self.selection.fields:
            tracer_settings = self.tracers[f.id]
            responses[f.id] = (
                InstrumentResponse(
                    pixel_width_angstrom_to_velocity(
                        float(tracer_settings["pix_width_ang"]),
                        lambda_obs_angstrom=observed_wavelength,
                    ),
                    SPEED_LIGHT_KMS
                    / (float(self.config["survey"]["resolution"]) * divisor),
                )
                if f.kind == "forest"
                else InstrumentResponse(0, 0)
            )
        return responses

    def bin_selection(self, index):
        """Apply the selection belonging to this recipe revision.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.

        Returns
        -------
        selection : PairSelection
            Revised bin-dependent selection, or the explicitly retained historical
            selection.
        """
        return (
            forecast_selection(self.case, index)
            if getattr(self, "recipe_revision", None)
            else self.selection
        )

    def samples(self, index, order, *, policy="primary", rectangular=False):
        """Sample source populations on the declared magnitude measure.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.
        order : int
            Gauss-Legendre quadrature order within each magnitude interval.
        policy : str
            Explicit source-input policy or sensitivity variant. Default is
            ``'primary'``.
        rectangular : bool
            Whether to use the legacy uniform magnitude measure instead of composite
            Gaussian quadrature. Default is ``False``.

        Returns
        -------
        samples : dict
            Per-field density, pixel variance and provenance at each magnitude node.
        magnitudes : ndarray, shape (n_magnitude,)
            Source magnitude nodes in mag.
        quadrature : ndarray, shape (n_magnitude,)
            Magnitude integration weights in mag.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.

        Notes
        -----
        Caches partitions and sampled arrays by bin, quadrature order and explicit policy. Returned arrays belong to that cache.
        """
        key = (index, order, policy, rectangular)
        if key in self._samples:
            return self._samples[key]
        forest_redshift = self.z(index)
        observed_wavelength = LYA_REST_ANGSTROM * (1 + forest_redshift)
        z_queries = {
            f.id: observed_wavelength
            / np.sqrt(
                float(self.tracers[f.id]["min_rest_frame_lya"])
                * float(self.tracers[f.id]["max_rest_frame_lya"])
            )
            - 1
            if f.kind == "forest"
            else forest_redshift
            for f in self.selection.fields
        }
        lo = float(self.config["survey"]["min_band_mag"])
        hi = float(self.config["survey"]["max_band_mag"])
        if index not in self._partitions:
            self._partitions[index] = breakpoints(
                self.densities, self.snrs, z_queries, lo, hi
            )

        if rectangular:
            magnitude_grid = np.linspace(
                lo, hi, int(self.config["survey"]["num mag bins"])
            )
            magnitude_weights = np.full(
                len(magnitude_grid), magnitude_grid[1] - magnitude_grid[0]
            )
        else:
            magnitude_grid, magnitude_weights = composite(
                self._partitions[index], order
            )
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

        # Restrict source preparation before constructing noise and covariance.
        active = np.unique(self.bin_selection(index).selected_pairs)
        for i, f in enumerate(self.selection.fields):
            if i not in active:
                continue
            density_sample = self.densities[f.id].sample(
                z_queries[f.id], magnitude_grid
            )
            values = density_sample["values"].copy()
            masks = density_sample["provenance"]["masks"]
            values[
                np.asarray(masks["density_floor"])
                | np.asarray(masks["negative_density"])
            ] = floor
            values /= width_factor
            row = dict(
                magnitudes=magnitude_grid,
                quadrature=magnitude_weights,
                density=values,
                z_source=z_queries[f.id],
                density_diagnostics=plain(density_sample["provenance"]),
            )
            if f.kind == "forest":
                tracer_settings = self.tracers[f.id]
                noise_sample = self.snrs[f.id].sample(
                    z_source=z_queries[f.id],
                    magnitudes=magnitude_grid,
                    wavelength=observed_wavelength,
                    pixel_width_angstrom=float(tracer_settings["pix_width_ang"]),
                    exposure_count=float(tracer_settings["num exposures"]),
                )
                removal_mask = np.zeros(len(magnitude_grid), dtype=bool)
                if policy == "remove_bright":
                    removal_mask = np.asarray(
                        noise_sample["provenance"]["masks"]["bright_clamp"]
                    )
                if policy == "remove_sentinel":
                    removal_mask = np.asarray(
                        noise_sample["provenance"]["masks"]["out_of_range"]
                    )
                row["density"][removal_mask] = 0
                row.update(
                    variance=noise_sample["values"],
                    snr_diagnostics=plain(noise_sample["provenance"]),
                    removed=removal_mask,
                    length_velocity=SPEED_LIGHT_KMS
                    * np.log(
                        float(tracer_settings["max_rest_frame_lya"])
                        / float(tracer_settings["min_rest_frame_lya"])
                    ),
                )
            result[f.id] = row
        self._samples[key] = (result, magnitude_grid, magnitude_weights)
        return result, magnitude_grid, magnitude_weights

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
        """Prepare one bin with the selected model, source measure and noise prescription.

        Parameters
        ----------
        index : int
            Zero-based redshift-bin index.
        controls : dict
            Quadrature orders, grid subdivisions, derivative step and applicable
            forest-weight convergence controls.
        policy : str
            Explicit source-input policy or sensitivity variant. Default is
            ``'primary'``.
        resolution : str
            Instrument-resolution convention; fwhm converts the quoted width to a
            Gaussian sigma. Default is ``'fwhm'``.
        growth : str
            Growth prescription; camb uses the prepared sigma8 ratio, while other
            values use the specified redshift scaling. Default is ``'camb'``.
        rectangular : bool
            Whether to use the legacy uniform magnitude measure instead of composite
            Gaussian quadrature. Default is ``False``.
        reconstruction : bool
            Whether to retain the configured galaxy BAO reconstruction factor.
            Default is ``True``.
        arithmetic_mean : bool
            Evaluate the intrinsic model at the arithmetic bin centre while
            retaining the geometric noise coordinates. Default is ``False``.

        Returns
        -------
        prepared : PreparedBin
            Fixed fiducial powers, responses, weights, noise and covariance factors.
        settings : dict
            Physical and numerical controls, sampled inputs and weight diagnostics.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.

        Notes
        -----
        Caches up to three prepared states. Derivative step is excluded from the preparation key because fiducial weights and covariance remain fixed.
        """
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
        forest_redshift = self.z(index)
        geometry = prepare_geometry(
            lo,
            hi,
            z_eval=forest_redshift,
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
        sampled, magnitude_grid, magnitude_weights = self.samples(
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
                    magnitude_grid,
                    magnitude_weights,
                    responses[field.id],
                    self.registry,
                    forest_redshift,
                    method=self.weight_method,
                    iterations=controls.get("iterations"),
                    policy=policy,
                    weight_rtol=controls.get("weight_rtol", 1e-4),
                )
            else:
                galaxies[field.id] = local_galaxy_density(
                    row["density"], magnitude_weights, geometry
                )

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
        if self.weight_method in ADAPTIVE:
            for name, weight in prepared.weights.items():
                forest_weighting[name].update(
                    result=plain(weight.convergence),
                    weights=weight.weights.tolist(),
                    A=weight.A,
                    P_pixel=weight.P_pixel,
                    auxiliary=dict(
                        k=weight.auxiliary.k,
                        mu=weight.auxiliary.mu,
                        P=weight.signal,
                        B=weight.alias,
                    ),
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
            magnitude_nodes=magnitude_grid.tolist(),
            magnitude_weights=magnitude_weights.tolist(),
            samples=plain(sampled),
            forest_weighting=dict(
                method=self.weight_method,
                reference=(
                    FIXED_REFERENCE
                    if self.weight_method == "inverse_variance"
                    else REFERENCE
                    if self.weight_method in ADAPTIVE
                    else None
                ),
                iterations=(
                    {"applicable": False, "status": "inapplicable"}
                    if self.weight_method == "inverse_variance"
                    else {
                        "applicable": True,
                        "stopping": {
                            **STOPPING,
                            "rtol": controls.get("weight_rtol", 1e-4),
                        },
                    }
                    if self.weight_method in ADAPTIVE
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
        if getattr(self, "recipe_revision", None):
            settings["recipe_identity"] = identity("accuracy", self.weight_method)
        result = (prepared, settings)
        self._prepared[key] = result
        while len(self._prepared) > 3:
            self._prepared.popitem(last=False)
        return result

    def evaluate(self, task, controls, **options):
        """Differentiate the observed BAO mean and verify its Fisher contraction.

        Parameters
        ----------
        task : dict
            Declared case, bin, selected field pairs, parameter order and validation
            thresholds.
        controls : dict
            Quadrature orders, grid subdivisions, derivative step and applicable
            forest-weight convergence controls.
        options : dict
            Keyword options forwarded to prepare for explicit sensitivity variants.

        Returns
        -------
        arrays : dict of str to ndarray
            Bound numerical evidence and independent Fisher comparison operands.
        report : dict
            Exact settings, input provenance and numerical metadata.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.
        """
        if task.get("recipe_revision") != getattr(self, "recipe_revision", None):
            raise ValueError("accuracy task and recipe revision differ")
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

        # Apply each field response once, then select the requested observables.
        selected = prepared.p3d.selection.selected_to_required
        observed_jacobian = (prepared.products[:, :, None] * derivatives.jacobian)[
            :, selected
        ][:, :, active]
        arrays, report = _assemble(
            task,
            prepared.total,
            observed_jacobian,
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
        independent_fisher, single = contract(
            arrays["selected_covariance"], observed_jacobian, independent=True
        )
        if (
            relative(independent_fisher, arrays["fisher"]) > 5e-12
            or relative(single, arrays["pair_fisher"]) > 5e-12
        ):
            raise ValueError("direct NumPy and evidence Fisher disagree")
        report["provenance"] = self.provenance
        return arrays, report

    def study(self, task):
        """Run the bounded controller with this prepared scientific recipe.

        Parameters
        ----------
        task : dict
            Declared case, bin, selected field pairs, parameter order and validation
            thresholds.

        Returns
        -------
        arrays : dict of str to ndarray
            Final payload and recorded numerical-refinement operands.
        report : dict
            Bounded convergence results, actual trials and unresolved controls.
        """
        from .study import study

        return study(self, task, payload_factory=OrderedDict)

    def sensitivity(self, task, controls, policy):
        """Recompute physical derived inputs for one explicit artificial policy test.

        Parameters
        ----------
        task : dict
            Declared case, bin, selected field pairs, parameter order and validation
            thresholds.
        controls : dict
            Quadrature orders, grid subdivisions, derivative step and applicable
            forest-weight convergence controls.
        policy : str
            Explicit source-input policy or sensitivity variant.

        Returns
        -------
        arrays : dict of str to ndarray
            Evidence recomputed under the specified artificial policy.
        report : dict
            Changed physical settings and explicit sensitivity interpretation.

        Raises
        ------
        ValueError :
            If inputs, declared identities or numerical validation conditions are
            inconsistent.
        """
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
        arrays, report = self.evaluate(
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
        report.update(
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
        return arrays, report
