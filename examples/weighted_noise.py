"""Fixed weighted forest and Poisson noise with an external amplitude model."""

import json

import numpy as np

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.derivatives import evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.geometry import (
    LYA_REST_ANGSTROM,
    SPEED_LIGHT_KMS,
    mode_counts,
    p1d_velocity_to_comoving,
    prepare_geometry,
    wavenumber_comoving_to_velocity,
)
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import (
    BoundParameters,
    P3DProvider,
    PreparedP3D,
    evaluate_p1d,
)
from fishhighz.models.p1d import default_p1d
from fishhighz.noise import forest_noise, galaxy_noise, prepare_noise
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import (
    InstrumentResponse,
    pair_response,
    pixel_width_angstrom_to_velocity,
    prepare_response,
    resolving_power_fwhm_to_sigma,
)
from fishhighz.weights import prepare_forest_weights, sample_auxiliary


def run(area=1000):
    """Compose a three-field amplitude forecast; verify an independent matrix oracle."""
    geometry = prepare_geometry(
        2,
        3,
        z_eval=2.4,
        area_deg2=area,
        h_fid=0.7,
        z_order=8,
        hubble=lambda z: 70 * (1 + z) ** 1.5,
        transverse_distance=lambda z: (
            2 * SPEED_LIGHT_KMS / 70 * (1 - 1 / np.sqrt(1 + z))
        ),
    )
    grid = gauss_legendre_grid([0.02, 0.1, 0.25], k_order=3, mu_order=4, h_fid=0.7)
    fields = [
        ObservedField("fq", "forest", "lya", background="qso"),
        ObservedField("fl", "forest", "lya", background="lbg"),
        ObservedField("g", "galaxy", "galaxy"),
    ]
    selection = PairSelection(fields, [("fl", "g"), ("fq", "fq"), ("fq", "g")])
    wavelength = LYA_REST_ANGSTROM * (1 + geometry.z_eval)
    settings = {
        "fq": InstrumentResponse(
            pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength),
            resolving_power_fwhm_to_sigma(2500),
        ),
        "fl": InstrumentResponse(100, 60),
        "g": InstrumentResponse(0, 0),
    }
    response = prepare_response(
        fields, grid.k_flat, grid.mu_flat, a_v=geometry.a_v, settings=settings
    )
    products = pair_response(response, selection)
    modes = mode_counts(geometry, grid)
    registry = ParameterRegistry([Parameter("amplitude", 1.3, "target", step=0.001)])
    bias = np.array([-0.7, -0.4, 1.2])
    base = np.outer(bias, bias)
    calls = []

    def power(theta, z, k, mu, pairs):
        calls.append(float(theta[0]))
        return theta[0] * (1 + k[:, None]) * base[pairs[:, 0], pairs[:, 1]]

    prepared = PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "amplitude",
                power,
                BoundParameters(registry, ("A",), {"A": "amplitude"}),
                selection.required_pairs,
            )
        ],
    )
    result = evaluate_derivatives(
        prepared, registry.fiducials, geometry.z_eval, grid.k_flat, grid.mu_flat
    )
    i, j = selection.required_pairs.T
    # Once per field/bin, using this same intrinsic callable and independent P1D.
    weight_results = []
    diagonal = {"g": np.full(len(modes), galaxy_noise(1.0))}
    p1d_binding = BoundParameters(registry, (), {})
    q = grid.k_flat * grid.mu_flat / geometry.a_v
    p1d = evaluate_p1d(default_p1d, p1d_binding, registry.fiducials, geometry.z_eval, q)
    for field, rho, variance in [
        (fields[0], [1, 2, 4], [1, 4, 9]),
        (fields[1], [2, 1, 3], [2, 7, 12]),
    ]:
        auxiliary = sample_auxiliary(
            field,
            geometry,
            settings[field.id],
            prepared,
            registry.fiducials,
            p1d_model=default_p1d,
            p1d_parameters=p1d_binding,
            theta_p1d=registry.fiducials,
            k_t_deg=2.4,
            k_p_velocity=0.00035,
        )
        weights = prepare_forest_weights(
            field,
            geometry,
            settings[field.id],
            z_source=3.0,
            magnitudes=[20, 21, 23],
            quadrature=[0.2, 0.5, 1.3],
            rho=np.array(rho) * 0.01,
            variance=variance,
            length_velocity=10000,
            method="legacy",
            iterations=3,
            auxiliary=auxiliary,
        )
        weight_results.append(weights)
        diagonal[field.id] = forest_noise(
            weights, field, geometry, settings[field.id], grid.k_flat, grid.mu_flat, p1d
        ).total
    noise = prepare_noise(
        selection, len(modes), diagonal=diagonal, independent_sampling=True
    )
    total = combine_observed_power(products * result.power, noise)
    covariance = gaussian_covariance(total, modes, selection)
    factors = factor_covariance(covariance)
    frozen = [
        geometry.z_nodes,
        geometry.w_z,
        response,
        modes,
        grid.k_flat,
        grid.mu_flat,
        grid.q_mode,
        factors,
        noise,
        *[p.weights for p in weight_results],
    ]
    snapshots = [a.copy() for a in frozen]
    selected = selection.selected_to_required
    expected_fisher = 0.0
    for node, (k, row, n) in enumerate(zip(grid.k_flat, response, modes)):
        observed = np.outer(row, row) * base * (1 + k)
        noise_matrix = np.diag([diagonal[f.id][node] for f in fields])
        matrix = registry.fiducials[0] * observed + noise_matrix
        pairs = selection.selected_pairs
        c = np.array(
            [
                [
                    (matrix[a, c] * matrix[b, d] + matrix[a, d] * matrix[b, c]) / n
                    for c, d in pairs
                ]
                for a, b in pairs
            ]
        )
        jac = np.array([observed[a, b] for a, b in pairs])
        np.testing.assert_allclose(covariance[node], c, rtol=5e-14, atol=1e-18)
        expected_fisher += jac @ np.linalg.solve(c, jac)
    fisher_values = []
    for scale in (1.0, 0.5, 0.25):
        result = evaluate_derivatives(
            prepared,
            registry.fiducials,
            geometry.z_eval,
            grid.k_flat,
            grid.mu_flat,
            step_scale=scale,
        )
        jacobian = (products[:, :, None] * result.jacobian)[:, selected, :]
        fisher = fisher_from_factors(jacobian, factors)[0, 0]
        np.testing.assert_allclose(fisher, expected_fisher, rtol=3e-12)
        fisher_values.append(float(fisher))
    q = wavenumber_comoving_to_velocity(grid.k_flat * grid.mu_flat, a_v=geometry.a_v)
    before_calls = len(calls)
    p1d_binding = BoundParameters(registry, (), {})
    intrinsic = evaluate_p1d(
        default_p1d, p1d_binding, registry.fiducials, geometry.z_eval, q
    )
    changed = evaluate_p1d(default_p1d, p1d_binding, [2.0], geometry.z_eval, q)
    np.testing.assert_array_equal(intrinsic, changed)
    assert len(calls) == before_calls
    smoothed = p1d_velocity_to_comoving(
        intrinsic * response[:, 0] ** 2, a_v=geometry.a_v
    )
    # An independent external P1D receives its own explicit binding and treatment.
    external_binding = BoundParameters(registry, ("norm",), {"norm": "amplitude"})
    external = evaluate_p1d(
        lambda t, z, k: t[0] / (1 + 100 * k),
        external_binding,
        registry.fiducials,
        geometry.z_eval,
        q,
    )
    external_smoothed = p1d_velocity_to_comoving(
        external * response[:, 0] ** 2, a_v=geometry.a_v
    )
    for a, b in zip(frozen, snapshots):
        np.testing.assert_array_equal(a, b)
    return dict(
        volume=geometry.volume,
        a_v=geometry.a_v,
        modes_sum=float(modes.sum()),
        amplitude_fisher=fisher_values,
        analytic_fisher=float(expected_fisher),
        amplitude_error=float(1 / np.sqrt(expected_fisher)),
        default_smoothed_p1d_comoving=smoothed.tolist(),
        external_smoothed_p1d_comoving=external_smoothed.tolist(),
        noise=noise.tolist(),
        weight_changes=[p.weight_changes.tolist() for p in weight_results],
        weights=[p.weights.tolist() for p in weight_results],
        null_directions=[],  # One amplitude, positive scalar information; no prior.
        interpretation="synthetic populations, not a physical DESI-2 forecast",
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
