"""Named numerical settings for the qualified DESI-2 accuracy prescription.

The values in this module are preparation settings, not a new estimator.  They
are kept outside ``fishhighz.validation`` so a native survey preparation can
record the exact S2--S4 recipe without importing historical evidence code.
"""

# Historical recipe identity. It keys the saved S2--S4 validation evidence
# (profiles, trials, schemas) and must not change; the spline density treatment
# of that evidence is selected by this revision in the native INI path as well.
REVISION = "early-lyaforecast-2026-09-18"

# Native-INI prescription identity from 2026-10-01: same weighting, response and
# quadrature recipe, but piecewise-constant source-density cells by default. The
# density tables and redshift edges are survey inputs set by the INI, not by the
# revision (the bundled INI switched to the SRD v2 inputs on the same date).
NATIVE_REVISION = "early-lyaforecast-2026-10-01"

# Integrated forest-source prescription identity from 2026-10-02: the native
# recipe (same weighting, response, magnitude quadrature and piecewise-constant
# density cells), but every forest field integrates over the source redshifts
# z_q whose forest overlaps the bin and over the forest pixels inside the bin
# slice, instead of using one representative source redshift. Selected by
# [input policies] forest_source_integration = integrated or by this revision.
INTEGRATED_REVISION = "early-lyaforecast-integrated-2026-10-02"
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


# Paired density interpolation and magnitude-partition policies. Which member is
# the default depends on the prescription revision (REVISION_DENSITY_INTERPOLATION);
# the other remains selectable through [input policies] density_interpolation.
DENSITY_INTERPOLATION_POLICIES = {
    "piecewise_constant_cells": (
        "piecewise_constant",
        "density_cell_edges_support_snr_nodes",
    ),
    "RectBivariateSpline_kx2_ky2_s0": (
        "spline",
        "density_knots_support_snr_nodes_negative_roots",
    ),
}

# Default density interpolation of each supported [prescription] revision. The
# historical revision reproduces the pre-2026-10-01 quadratic-spline treatment of
# the source densities; the native revision uses piecewise-constant cells.
REVISION_DENSITY_INTERPOLATION = {
    NATIVE_REVISION: "piecewise_constant_cells",
    INTEGRATED_REVISION: "piecewise_constant_cells",
    REVISION: "RectBivariateSpline_kx2_ky2_s0",
}

# Native INI defaults expand the qualified prescription without changing its
# scientific inputs. Strings preserve the explicit native-schema representation.
# density_interpolation is absent here: it follows the prescription revision
# (REVISION_DENSITY_INTERPOLATION). magnitude_partition is absent too: it follows
# density_interpolation unless stated explicitly.
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
        "snr_smoothing": "legacy",
        "snr_interpolation": "linear_RegularGridInterpolator",
        "snr_bright_policy": "clamp_to_brightest_tabulated_magnitude",
        "snr_clamp": "1e-10",
        "snr_sentinel": "1e20",
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

# Gauss-Legendre orders of the integrated forest-source quadrature, used when the
# [numerical] forest_* keys are absent in integrated mode. They are deliberately
# not part of INI_DEFAULTS: the central mode must expand to exactly the same
# values as before, so compact and expanded INIs and all existing provenance stay
# unchanged. The orders are provisional until the refinement study fixes them; on
# the DESI-2 QSO and LBG S/N tables (16, 16, 4) with breakpoints at the S/N source
# redshifts is within 3.5e-5 of a (64, 64, 8) reference for N1 and N3/N1^2, whereas
# (16, 16, 1) leaves 0.15-0.22 per cent.
#   forest_zq_order      : order per panel in ln(1+z_q)
#   forest_lambda_order  : order per wavelength panel inside the bin slice
#   forest_lambda_panels : number of equal wavelength panels per source redshift
INTEGRATED_DEFAULTS = {
    "forest_zq_order": 16,
    "forest_lambda_order": 16,
    "forest_lambda_panels": 4,
}


def accuracy_settings():
    """Return independent copies of the adopted accuracy settings.

    Returns
    -------
    settings : dict
        Historical validation recipe revision (``REVISION``), weighting
        controls, reference modes and quadrature settings for provenance. The
        native-INI prescription revision is ``NATIVE_REVISION`` and is recorded
        by ``parse_survey_ini``, not here.
    """

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
    "INTEGRATED_DEFAULTS",
    "INTEGRATED_REVISION",
    "NATIVE_REVISION",
    "REFERENCE",
    "REVISION",
    "REVISION_DENSITY_INTERPOLATION",
    "STOPPING",
    "accuracy_settings",
]
