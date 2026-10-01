"""Observed matrix, selection, amplitude and area normalization acceptance."""

import runpy
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from fishhighz.covariance import combine_observed_power, gaussian_covariance
from fishhighz.derivatives import evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.geometry import mode_counts, prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import (
    BoundParameters,
    P3DProvider,
    PreparedP3D,
    evaluate_p1d,
)
from fishhighz.models.p1d import default_p1d
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import (
    InstrumentResponse,
    pair_response,
    prepare_response,
    velocity_response,
)


def test_example():
    """Check example."""
    report = runpy.run_path(
        str(Path(__file__).parents[1] / "examples/survey_primitives.py")
    )["run"]()
    assert_allclose(report["amplitude_fisher"], report["analytic_fisher"], rtol=3e-12)


def test_area_scaling_one_field():
    """Check area scaling one field."""
    grid = gauss_legendre_grid([0.05, 0.2], k_order=3, mu_order=3, h_fid=0.7)
    selection = PairSelection([ObservedField("g", "galaxy", "g")])
    results = []
    for area in (100, 200):
        bin_geometry = prepare_geometry(
            2,
            3,
            z_eval=2.5,
            area_deg2=area,
            h_fid=0.7,
            z_order=4,
            hubble=lambda z: z * 100,
            transverse_distance=lambda z: z * 1000,
        )
        modes = mode_counts(bin_geometry, grid)
        total = np.full((len(modes), 1), 5.0)
        covariance = gaussian_covariance(total, modes, selection)
        assert_allclose(covariance[:, 0, 0], 2 * 25 / modes, rtol=2e-15)
        fisher = fisher_from_factors(
            np.full((len(modes), 1, 1), 2.0), factor_covariance(covariance)
        )[0, 0]
        assert_allclose(fisher, np.sum(4 * modes / 50), rtol=3e-15)
        results.append((modes, covariance, fisher))
    assert_allclose(results[1][0], 2 * results[0][0], rtol=0, atol=0)
    assert_allclose(results[1][1], results[0][1] / 2, rtol=0, atol=0)
    assert_allclose(results[1][2], 2 * results[0][2], rtol=3e-15)
    assert_allclose(
        1 / np.sqrt(results[1][2]), 1 / np.sqrt(results[0][2]) / np.sqrt(2), rtol=2e-15
    )


def test_permutation_closure_slices_and_fixed_preparation():
    """Check permutation closure slices and fixed preparation."""
    fields = [
        ObservedField("a", "forest", "lya", background="q"),
        ObservedField("b", "forest", "lya", background="l"),
        ObservedField("g", "galaxy", "g"),
    ]
    k_grid = np.array([0.1, 0.3, 0.2, 0.5])
    mu_grid = np.array([0, 1, 0.3, 0.8])
    settings = dict(
        a=InstrumentResponse(100, 30),
        b=InstrumentResponse(200, 50),
        g=InstrumentResponse(0, 0),
    )
    response = prepare_response(fields, k_grid, mu_grid, a_v=100, settings=settings)
    registry = ParameterRegistry([Parameter("amplitude", 1.0, "target", step=0.001)])
    biases = np.array([-0.3, -0.5, 2.0])
    signal = np.outer(biases, biases)
    noise = np.diag([1.0, 2.0, 3.0])

    def model(t, z, k, mu, pairs):
        """Evaluate the synthetic spectrum used by the enclosing regression test.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines.
        pairs : ndarray of int, shape (n_pairs, 2)
            Observed-field indices defining the requested spectra.

        Returns
        -------
        power : ndarray of shape (n_nodes, n_pairs)
            Synthetic intrinsic power in (Mpc/h)^3 before response and noise.
        """
        return np.broadcast_to(
            t[0] * signal[pairs[:, 0], pairs[:, 1]], (len(k), len(pairs))
        )

    selections = [
        PairSelection(fields, selected)
        for selected in [
            [("b", "g"), ("a", "a"), ("g", "a")],
            [("a", "g"), ("a", "a"), ("g", "b")],
        ]
    ]
    outputs = []
    for selection in selections:
        product = pair_response(response, selection)
        i, j = selection.required_pairs.T
        prepared = PreparedP3D(
            registry,
            selection,
            [
                P3DProvider(
                    "p",
                    model,
                    BoundParameters(registry, ["A"], {"A": "amplitude"}),
                    selection.required_pairs,
                )
            ],
        )
        result = evaluate_derivatives(prepared, [1], 2.4, k_grid, mu_grid)
        modes = np.array([10, 20, 30, 40.0])
        cov = gaussian_covariance(
            combine_observed_power(
                product * result.power, np.broadcast_to(noise[i, j], product.shape)
            ),
            modes,
            selection,
        )
        factors = factor_covariance(cov)
        snapshots = [x.copy() for x in (response, modes, factors, k_grid, mu_grid)]
        analytic = product * signal[i, j]
        assert_allclose(
            product[:, :, None] * result.jacobian,
            analytic[:, :, None],
            rtol=3e-13,
            atol=1e-15,
        )
        observed = analytic[:, selection.selected_to_required]
        for node in range(len(k_grid)):
            matrix = np.outer(response[node], response[node]) * signal + noise
            oracle = [
                [
                    (
                        matrix[a, c] * matrix[biases, d]
                        + matrix[a, d] * matrix[biases, c]
                    )
                    / modes[node]
                    for c, d in selection.selected_pairs
                ]
                for a, biases in selection.selected_pairs
            ]
            assert_allclose(cov[node], oracle, rtol=3e-15, atol=1e-18)
        fisher = fisher_from_factors(observed[:, :, None], factors)
        for theta in ([0.8], [1.2], [1.0]):
            evaluate_derivatives(prepared, theta, 2.4, k_grid, mu_grid)
        for x, y in zip((response, modes, factors, k_grid, mu_grid), snapshots):
            assert_array_equal(x, y)
        slice_result = evaluate_derivatives(
            prepared, [1], 2.4, k_grid[::2], mu_grid[::2]
        )
        assert_allclose(slice_result.jacobian, result.jacobian[::2], rtol=0, atol=0)
        outputs.append((cov, observed, fisher))
    assert_allclose(outputs[1][0], outputs[0][0][:, ::-1, ::-1], rtol=0, atol=0)
    assert_allclose(outputs[1][1], outputs[0][1][:, ::-1], rtol=0, atol=0)
    assert_allclose(outputs[1][2], outputs[0][2], rtol=3e-15)
    velocity_k = k_grid * mu_grid / 100
    power_1d = evaluate_p1d(
        default_p1d, BoundParameters(registry, (), {}), [1], 2.4, velocity_k
    )
    velocity_w = velocity_response(
        velocity_k, pixel_width_velocity=100, gaussian_sigma_velocity=30
    )
    comoving_w = np.sinc(k_grid * mu_grid / (2 * np.pi)) * np.exp(
        -0.5 * (k_grid * mu_grid * 0.3) ** 2
    )
    assert_allclose(
        power_1d * velocity_w**2 / 100, (power_1d / 100) * comoving_w**2, rtol=5e-14
    )
