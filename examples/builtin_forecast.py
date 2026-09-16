"""Synthetic built-in BAO and common-AP plus f forecasts; no survey physics."""

import json

import numpy as np

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.derivatives import check_convergence, evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.models.kaiser import KaiserModel, Scaling
from fishhighz.models.templates import prepare_template
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult, diagonal_prior


def run():
    """Return step convergence and identifiable target errors in two toy modes."""
    knots = np.geomspace(0.01, 0.7, 600)
    smooth = 100 / (1 + 2 * knots)
    wiggle = 8 * np.sin(110 * knots) * np.exp(-((knots / 0.4) ** 2))
    template = prepare_template(
        knots, smooth + wiggle, smooth, z_ref=2.4, h_template=0.7, h_fid=0.7
    )
    grid = gauss_legendre_grid([0.04, 0.12, 0.25], k_order=6, mu_order=6, h_fid=0.7)
    if grid.h_fid != template.h_fid:
        raise ValueError("grid and template h_fid must match")
    fields = [
        ObservedField("F", "forest", "toy", background="qso"),
        ObservedField("g", "galaxy", "toy"),
    ]
    selection = PairSelection(fields, [("g", "g"), ("F", "g")])
    response = np.column_stack([0.8 + 0.1 * grid.mu_flat, 0.9 + 0.1 * grid.k_flat])
    i, j = selection.required_pairs.T
    products = response[:, i] * response[:, j]
    noise = np.tile(
        [0.5, 0, 0.7], (len(grid.k_flat), 1)
    )  # supplied positive diagonal noise
    saved = [
        (value, value.copy())
        for value in (
            grid.k_flat,
            grid.mu_flat,
            grid.weights,
            grid.q_mode,
            response,
            noise,
        )
    ]
    reports = {}
    for mode in ("bao", "ap_f"):
        common = mode == "ap_f"
        parameters = [
            Parameter("ap", 1, "target", step=0.001),
            Parameter("at", 1, "target", step=0.001),
            Parameter("bg", 1.7, "nuisance", step=0.001),
            Parameter("bf", -0.3, "nuisance", step=0.001),
            Parameter("beta", 1.2, "nuisance", step=0.001),
        ]
        if common:
            parameters.append(Parameter("f", 0.8, "target", step=0.001))
        registry = ParameterRegistry(parameters)
        bindings = {"w_ap": "ap", "w_at": "at", "bg": "bg", "bf": "bf", "beta": "beta"}
        if common:
            bindings.update({"nw_ap": "ap", "nw_at": "at", "rate": "f"})
        model = KaiserModel(
            template,
            fields,
            biases={"F": "bf", "g": "bg"},
            betas={"F": "beta"},
            widths={"F": (3, 2), "g": (3, 2)},
            f="rate" if common else 0.8,
            local_names=tuple(bindings),
            smooth=Scaling("ap_at", ap="nw_ap", at="nw_at") if common else None,
            wiggle=Scaling("ap_at", ap="w_ap", at="w_at"),
        )
        bound = BoundParameters(registry, model.local_names, bindings)
        collection = PreparedP3D(
            registry,
            selection,
            [P3DProvider(mode, model, bound, selection.required_pairs)],
        )
        args = (collection, registry.fiducials, model.z, grid.k_flat, grid.mu_flat)
        derivatives = evaluate_derivatives(*args)
        total = combine_observed_power(products * derivatives.power, noise)
        covariance = gaussian_covariance(total, 1e6 * grid.q_mode, selection)
        factors = factor_covariance(covariance)  # one fixed factorization per mode
        factor_snapshot = factors.copy()
        # Explicit toy external constraints on nuisances, not inferred from bounds.
        prior = diagonal_prior(registry, {"bg": 0.5, "bf": 0.2, "beta": 0.5})
        matrices, jacobians, errors = [], [], []
        for scale in (1.0, 0.5, 0.25):
            result = (
                derivatives
                if scale == 1
                else evaluate_derivatives(*args, step_scale=scale)
            )
            observed_jacobian = products[:, :, None] * result.jacobian
            matrix = fisher_from_factors(
                observed_jacobian[:, selection.selected_to_required], factors
            )
            forecast = FisherResult(registry, matrix, prior_fisher=prior)
            targets = ["ap", "at", "f"] if common else ["ap", "at"]
            errors.append(forecast.marginalized_errors(targets).tolist())
            matrices.append(matrix)
            jacobians.append(result.jacobian)
        study = check_convergence(*args, step_scale=0.5, atol=1e-5, rtol=2e-3)
        assert study.passed
        np.testing.assert_allclose(matrices[0], matrices[-1], rtol=2e-3, atol=1e-7)
        np.testing.assert_array_equal(factors, factor_snapshot)
        reports[mode] = {
            "targets": targets,
            "target_errors_by_step_scale": errors,
            "step_scales": [1, 0.5, 0.25],
            "jacobian_max_changes": [
                float(np.max(np.abs(jacobians[i] - jacobians[i + 1]))) for i in range(2)
            ],
            "fisher_max_changes": [
                float(np.max(np.abs(matrices[i] - matrices[i + 1]))) for i in range(2)
            ],
            "convergence_passed": study.passed,
            "rank_with_explicit_nuisance_priors": forecast.diagnostics.rank,
        }
    for value, snapshot in saved:
        np.testing.assert_array_equal(value, snapshot)
    return reports


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
