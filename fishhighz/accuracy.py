"""Named numerical settings for the qualified DESI-2 accuracy prescription.

The values in this module are preparation settings, not a new estimator.  They
are kept outside ``fishhighz.validation`` so a native survey preparation can
record the exact S2--S4 recipe without importing historical evidence code.
"""

REVISION = "early-lyaforecast-2026-09-18"
ADAPTIVE = ("early_lyaforecast", "mcdonald")
STOPPING = {"rtol": 1e-4, "min_updates": 3, "stable_steps": 3, "max_updates": 96}
REFERENCE = {
    "k_t_deg": 2.4,
    "k_p_velocity": 0.00035,
    "convention": "fiducial_auto_p3d_and_p1d_times_response_squared",
}
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
FIXED_REFERENCE = {
    "convention": "intrinsic_p1d_times_field_response_squared",
    "q_star": 0.00035,
    "q_star_units": "s/km",
}


# Native INI defaults expand the qualified prescription without changing its
# scientific inputs. Strings preserve the explicit native-schema representation.
INI_DEFAULTS = {
    "model": {
        "parameterization": "ap_at",
        "parameter_names": "ap, at",
        "smooth_scaling": "identity",
        "wiggle_scaling": "ap_at",
        "damping": "mixed_squared_width",
        "growth": "camb_sigma8_ratio",
        "reconstruction": "per_field",
        "forest_beta": "1.45",
        "damping_amplitude": "3.26",
        "forest_reconstruction_factor": "1.0",
    },
    "input policies": {
        "density_semantics": "cell_count_per_deg2",
        "density_width_policy": "legacy_first_spacing",
        "density_redshift_normalization": "target_density",
        "density_negative_policy": "floor_negative",
        "density_interpolation": "RectBivariateSpline_kx2_ky2_s0",
        "snr_smoothing": "legacy",
        "snr_interpolation": "linear_RegularGridInterpolator",
        "snr_bright_policy": "clamp_to_brightest_tabulated_magnitude",
        "snr_clamp": "1e-10",
        "snr_sentinel": "1e20",
        "magnitude_partition": "density_knots_support_snr_nodes_negative_roots",
        "weighting_method": "early_lyaforecast",
        "weighting_reference": "fiducial_auto_p3d_and_p1d_times_response_squared",
        "weighting_reference_k_t_deg": "2.4",
        "weighting_reference_k_p_velocity": "0.00035",
        "weighting_rtol": "1e-5",
        "weighting_min_updates": "3",
        "weighting_stable_steps": "3",
        "weighting_max_updates": "96",
        "sampling_noise": "independent_diagonal",
    },
    "numerical": {
        "k_intervals": "128",
        "k_order": "4",
        "k_min": "0.01",
        "k_max": "0.5",
        "mu_order": "32",
        "z_order": "32",
        "magnitude_order": "16",
        "ap_step": "2.5e-4",
        "weight_rtol": "1e-5",
    },
    "survey": {
        "band": "r",
        "lya_rest_angstrom": "1215.67",
        "evaluation_redshift": "geometric_1plusz",
    },
}


def accuracy_settings():
    """Return plain copies of the adopted settings for provenance."""

    return {
        "revision": REVISION,
        "adaptive_methods": list(ADAPTIVE),
        "stopping": dict(STOPPING),
        "reference": dict(REFERENCE),
        "controls": dict(CONTROLS),
        "fixed_reference": dict(FIXED_REFERENCE),
    }


__all__ = [
    "ADAPTIVE",
    "CONTROLS",
    "FIXED_REFERENCE",
    "REFERENCE",
    "REVISION",
    "STOPPING",
    "accuracy_settings",
]
