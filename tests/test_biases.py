"""Analytic and tabulated bias controls against the lyaforecast formulae."""

import numpy as np

from fishhighz.models.biases import (
    AnalyticBias,
    analytic_beta_rsd,
    analytic_density_bias,
    linear_tabulated_bias,
)


def test_analytic_density_bias_and_beta_evolution():
    z = 2.9
    expected_lya = -0.1352 * ((1 + z) / (1 + 2.33)) ** 2.9
    expected_qso = 3.54 * ((1 + z) / (1 + 2.33)) ** 1.44
    assert analytic_density_bias(z, "lya") == expected_lya
    assert analytic_density_bias(z, "qso") == expected_qso
    np.testing.assert_allclose(analytic_beta_rsd(z, "lya"), 1.45)
    np.testing.assert_allclose(
        analytic_beta_rsd(z, "qso", growth_rate=0.9), 0.9 / expected_qso
    )


def test_tabulated_linear_interpolation_and_endpoint_extrapolation():
    bias = linear_tabulated_bias([2.0, 3.0, 4.0], [10.0, 20.0, 35.0])
    np.testing.assert_allclose(bias([1.0, 2.5, 5.0]), [0.0, 15.0, 50.0])
    assert not bias.redshifts.flags.writeable
    assert not bias.values.flags.writeable


def test_tabulated_bias_registration_preserves_beta_definition():
    model = AnalyticBias(growth_rate=lambda z: np.asarray(z) * 0 + 1.0)
    model.set_density_bias_func("qso", linear_tabulated_bias([2.0, 3.0], [2.0, 4.0]))
    assert model.density_bias(2.5, "qso") == 3.0
    assert model.beta_rsd(2.5, "qso") == 1 / 3
