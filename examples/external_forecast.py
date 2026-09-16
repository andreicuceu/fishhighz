"""Standalone synthetic forecast; no external cosmology or survey accuracy claim."""

import json

import numpy as np

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.derivatives import check_convergence, evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import (
    BoundParameters,
    P3DProvider,
    PreparedP3D,
    evaluate_p1d,
)
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult, diagonal_prior


# These plain callables can live in any package; no FishHighz types are used.
def power(theta, z, k, mu, pairs):
    amplitude, nuisance = theta
    matrix = np.array([[3.0, -0.6], [-0.6, 2.0]])
    shape = np.exp(amplitude * k) * (1 + nuisance * mu**2)
    return shape[:, None] * matrix[pairs[:, 0], pairs[:, 1]]


def jacobian(theta, z, k, mu, pairs):
    value = power(theta, z, k, mu, pairs)
    return np.stack(
        [value * k[:, None], value * (mu**2 / (1 + theta[1] * mu**2))[:, None]], axis=-1
    )


def run():
    """Return small reviewable analytic/FD and independent P1D evidence."""
    registry = ParameterRegistry(
        [
            Parameter("shared", 0.7, "target", step=0.02),
            Parameter("nuisance", 0.2, "nuisance", step=0.01),
        ]
    )
    selection = PairSelection(
        [
            ObservedField("g", "galaxy", "external"),
            ObservedField("f", "forest", "external", background="qso"),
        ],
        selected=[("f", "f"), ("f", "g")],
    )
    grid = gauss_legendre_grid([0.05, 0.15, 0.3], k_order=3, mu_order=4, h_fid=0.7)
    bound = BoundParameters(registry, ["a", "n"], {"a": "shared", "n": "nuisance"})
    prepared = PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "external",
                power,
                bound,
                [("g", "g"), ("f", "g"), ("f", "f")],
                jacobian=jacobian,
            )
        ],
    )
    args = (prepared, registry.fiducials, 2.4, grid.k_flat, grid.mu_flat)
    analytic = evaluate_derivatives(*args)
    # Supplied synthetic W and N remain fixed. Providers return intrinsic P.
    response = np.column_stack([0.8 + 0.1 * grid.k_flat, 0.6 + 0.1 * grid.mu_flat])
    i, j = selection.required_pairs.T
    products = response[:, i] * response[:, j]
    noise = np.tile([0.7, 0.05, 1.1], (len(grid.k_flat), 1))
    total = combine_observed_power(products * analytic.power, noise)
    covariance = gaussian_covariance(total, 1e6 * grid.q_mode, selection)
    factors = factor_covariance(covariance)  # exactly once
    selected = selection.selected_to_required

    def information(result):
        observed_jacobian = products[:, :, None] * result.jacobian
        return fisher_from_factors(observed_jacobian[:, selected, :], factors)

    exact = information(analytic)
    prior = diagonal_prior(registry, {"nuisance": 0.5})
    exact_errors = FisherResult(
        registry, exact, prior_fisher=prior
    ).marginalized_errors()
    comparisons = []
    for scale in (1.0, 0.5, 0.25):
        result = evaluate_derivatives(*args, numerical=True, step_scale=scale)
        matrix = information(result)
        errors = FisherResult(
            registry, matrix, prior_fisher=prior
        ).marginalized_errors()
        np.testing.assert_allclose(matrix, exact, rtol=2e-5, atol=1e-10)
        np.testing.assert_allclose(errors, exact_errors, rtol=1e-5, atol=1e-10)
        comparisons.append(
            {
                "step_scale": scale,
                "max_fisher_relative_error": float(np.max(np.abs(matrix / exact - 1))),
                "marginalized_errors": errors.tolist(),
                "calls": result.calls[0].model,
            }
        )
    convergence = check_convergence(
        *args, numerical=True, atol=1e-8, rtol=1e-5, refinements=2
    )
    assert convergence.passed
    # Independent velocity-space P1D binding, with a different registry/order.
    registry_1d = ParameterRegistry([Parameter("velocity_amplitude", 4.0, "nuisance")])
    binding_1d = BoundParameters(registry_1d, ["A"], {"A": "velocity_amplitude"})
    velocity_k = np.array([0.0, 0.001, 0.01])  # s/km
    p1d = evaluate_p1d(
        lambda t, z, q: t[0] / (1 + q * 100),
        binding_1d,
        registry_1d.fiducials,
        2.4,
        velocity_k,
    )
    return {
        "analytic_marginalized_errors": exact_errors.tolist(),
        "comparisons": comparisons,
        "convergence_passed": convergence.passed,
        "p1d_km_per_s": p1d.tolist(),
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
