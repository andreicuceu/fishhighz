"""Analytic AP, covariance and rank checks for the optional full EFT route."""

import numpy as np
from numpy.testing import assert_allclose

from fishhighz.covariance import gaussian_covariance
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.models.eft import component_coordinates


def test_alpha_phi_endpoints_and_jacobian():
    k = np.array([0.1, 0.2])
    mu = np.array([0.0, 1.0])
    alpha, phi = 1.03, 0.96
    km, mm, jac = component_coordinates(k, mu, alpha, phi)
    assert_allclose(
        km, k / np.array([alpha * np.sqrt(phi), alpha / np.sqrt(phi)]), rtol=2e-15
    )
    assert_allclose(mm, mu, atol=2e-15)
    assert_allclose(jac, alpha**-3 * phi**-0.5, rtol=2e-15)


def test_three_spectrum_contraction_equals_field_trace():
    selection = PairSelection(
        [ObservedField("F", "forest", "f", "qso"), ObservedField("G", "galaxy", "g")]
    )
    power = np.array([[3.0, -0.7, 5.0]])
    jac = np.array([[[0.2, 0.4], [-0.1, 0.3], [0.8, -0.2]]])
    counts = np.array([123.0])
    cov = gaussian_covariance(power, counts, selection)
    actual = fisher_from_factors(jac, factor_covariance(cov))
    c = np.array([[3.0, -0.7], [-0.7, 5.0]])
    inverse = np.linalg.inv(c)
    ds = [
        np.array([[jac[0, 0, i], jac[0, 1, i]], [jac[0, 1, i], jac[0, 2, i]]])
        for i in range(2)
    ]
    expected = np.array(
        [[counts[0] / 2 * np.trace(inverse @ a @ inverse @ b) for b in ds] for a in ds]
    )
    assert_allclose(actual, expected, rtol=2e-14, atol=1e-14)
