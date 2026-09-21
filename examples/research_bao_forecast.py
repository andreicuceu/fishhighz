"""Standalone synthetic joint-BAO forecast with the research-baseline APIs.

This example is intentionally small and is not a DESI-2 forecast.  It generates
its template, geometry and survey arrays in memory.  Install
``fishhighz[templates]``: SciPy is used only to prepare the synthetic spline;
the forecast itself uses NumPy.
"""

import json

import numpy as np

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.forecast import prepare_bin, run_forecast
from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.models.kaiser import KaiserModel, Scaling
from fishhighz.models.p1d import default_p1d
from fishhighz.models.templates import prepare_template
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import (
    InstrumentResponse,
    pixel_width_angstrom_to_velocity,
    resolving_power_fwhm_to_sigma,
)
from fishhighz.survey import BinSpec, ForestInput


def _geometry():
    """Use a simple, internally consistent Einstein-de Sitter background."""
    hubble_0 = 70.0
    z_min, z_max = 2.3, 2.5
    z_eval = np.sqrt((1 + z_min) * (1 + z_max)) - 1

    def hubble(redshift):
        return hubble_0 * (1 + redshift) ** 1.5

    def transverse_distance(redshift):
        return 2 * SPEED_LIGHT_KMS / hubble_0 * (1 - 1 / np.sqrt(1 + redshift))

    return prepare_geometry(
        z_min,
        z_max,
        z_eval=z_eval,
        area_deg2=1000.0,
        h_fid=0.7,
        z_order=8,
        hubble=hubble,
        transverse_distance=transverse_distance,
    )


def _template():
    """Prepare a smooth positive toy spectrum with a damped BAO-like ripple."""
    k = np.geomspace(0.005, 1.0, 500)
    smooth = 6000 * (k / 0.1) ** -2.1 / (1 + (k / 0.25) ** 2)
    wiggle = 0.06 * smooth * np.sin(105 * k) * np.exp(-((k / 0.45) ** 2))
    return prepare_template(
        k,
        smooth + wiggle,
        smooth,
        z_ref=np.sqrt(3.3 * 3.5) - 1,
        h_template=0.7,
        h_fid=0.7,
        metadata={"interpretation": "synthetic BAO-like template"},
    )


def _prepare_selection(registry, geometry, template, selected_pairs):
    fields = (
        ObservedField("forest", "forest", "lya", background="qso"),
        ObservedField("galaxy", "galaxy", "qso"),
    )
    selection = PairSelection(fields, selected_pairs)
    model = KaiserModel(
        template,
        fields,
        biases={"forest": -0.14, "galaxy": 3.2},
        betas={"forest": 1.4},
        widths={"forest": (6.0, 3.0), "galaxy": (4.0, 2.0)},
        f=1.0,
        local_names=("ap", "at"),
        # Smooth power stays at identity.  AP derivatives act on the full
        # wiggle component, including its Kaiser factors, volume and damping.
        smooth=Scaling("ap_at", ap=1.0, at=1.0),
        wiggle=Scaling("ap_at", ap="ap", at="at"),
        z=geometry.z_eval,
    )
    binding = BoundParameters(
        registry,
        model.local_names,
        {"ap": "ap", "at": "at"},
    )
    p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("synthetic_kaiser", model, binding, selection.required_pairs)],
    )

    wavelength = LYA_REST_ANGSTROM * (1 + geometry.z_eval)
    forest_response = InstrumentResponse(
        pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength),
        # R is lambda/FWHM_lambda; the helper converts that physical FWHM to sigma.
        resolving_power_fwhm_to_sigma(3200.0),
    )
    responses = {
        "forest": forest_response,
        "galaxy": InstrumentResponse(0.0, 0.0),
    }
    forests = {}
    if any(0 in pair for pair in selection.selected_pairs):
        forests["forest"] = ForestInput(
            {
                "method": "early_lyaforecast",
                "z_source": 2.8,
                "magnitudes": [20.0, 21.0, 22.0, 23.0],
                "quadrature": [0.5, 1.0, 1.0, 0.5],
                "rho": [0.30, 0.24, 0.15, 0.07],
                "variance": [1.0, 2.0, 5.0, 12.0],
                "length_velocity": SPEED_LIGHT_KMS * np.log(1200.0 / 1040.0),
            },
            default_p1d,
            BoundParameters(registry, (), {}),
            registry.fiducials,
            auxiliary_coordinates=(2.4, 0.00035),
            provenance={"interpretation": "synthetic normalized forest samples"},
        )
    galaxies = (
        {"galaxy": 2.0e-4}
        if any(1 in pair for pair in selection.selected_pairs)
        else {}
    )
    grid = gauss_legendre_grid(
        [0.03, 0.08, 0.16, 0.25], k_order=4, mu_order=6, h_fid=0.7
    )
    return prepare_bin(
        BinSpec(
            "synthetic-z2.4",
            geometry,
            grid,
            p3d,
            responses,
            forests=forests,
            galaxies=galaxies,
            independent_sampling=True,
        )
    )


def run():
    """Return joint and individual synthetic AP errors plus preparation metadata."""
    registry = ParameterRegistry(
        [
            Parameter("ap", 1.0, "target", step=0.002),
            Parameter("at", 1.0, "target", step=0.002),
        ]
    )
    geometry = _geometry()
    template = _template()
    selected = (("forest", "forest"), ("forest", "galaxy"), ("galaxy", "galaxy"))
    prepared = _prepare_selection(registry, geometry, template, selected)
    joint_run = run_forecast([prepared], batch_size=24)
    joint = joint_run.combined
    individual_errors = {}
    individual_fishers = []
    for pair in selected:
        one = _prepare_selection(registry, geometry, template, [pair])
        result = run_forecast([one], batch_size=24).combined
        individual_fishers.append(result.data_fisher)
        individual_errors[" x ".join(pair)] = result.marginalized_errors(
            ["ap", "at"]
        ).tolist()
    summed_individual_fisher = np.sum(individual_fishers, axis=0)
    fisher_scale = np.max(np.abs(joint.data_fisher))
    joint_sum_relative_difference = (
        np.max(np.abs(joint.data_fisher - summed_individual_fisher)) / fisher_scale
    )

    forest_weights = prepared.weights["forest"]
    return {
        "interpretation": "synthetic API example; not a DESI-2 forecast",
        "selected_spectra": [" x ".join(pair) for pair in selected],
        "joint_ap_at_errors": joint.marginalized_errors(["ap", "at"]).tolist(),
        "individual_ap_at_errors": individual_errors,
        "joint_ap_at_correlation": float(joint.correlations(["ap", "at"])[0, 1]),
        "joint_data_fisher": joint.data_fisher.tolist(),
        "summed_individual_data_fisher": summed_individual_fisher.tolist(),
        "joint_vs_sum_max_relative_difference": float(joint_sum_relative_difference),
        "joint_rank": joint.diagnostics.rank,
        "weighting": {
            "method": forest_weights.method,
            "auxiliary_coordinates_deg_inv_s_per_km": [2.4, 0.00035],
            "status": forest_weights.convergence["status"],
            "updates": forest_weights.convergence["updates"],
            "candidate": forest_weights.convergence["candidate"],
            "stopping_controls": {
                "rtol": 1e-4,
                "min_updates": 3,
                "stable_steps": 3,
                "max_updates": 96,
            },
        },
        "fixed_preparation": {
            "noise_convention": prepared.diagnostics["noise_convention"],
            "selected_pair_count": prepared.diagnostics["selected_pair_count"],
            "required_pair_count": prepared.diagnostics["required_pair_count"],
            "node_count": prepared.diagnostics["node_count"],
            "derivative_batches": len(joint_run.bins[0].node_slices),
        },
        "damping_widths_mpc_over_h": {
            "forest_parallel_transverse": [6.0, 3.0],
            "galaxy_parallel_transverse": [4.0, 2.0],
            "cross_squared_width_rule": "arithmetic mean of the two auto squared widths",
        },
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
