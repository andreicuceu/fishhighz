import numpy as np
import pytest

from fishhighz.models.physical_ap import physical_ap_coordinates


def test_identity_and_endpoints():
    k = np.array([0.01, 0.1, 0.2])
    mu = np.array([0.0, 0.4, 1.0])
    km, mm, j = physical_ap_coordinates(k, mu, 1.0, 1.0)
    np.testing.assert_allclose(km, k, rtol=2e-16)
    np.testing.assert_allclose(mm, mu, atol=2e-16)
    assert j == 1.0
    km, mm, j = physical_ap_coordinates(k, mu, 1.1, 0.9)
    np.testing.assert_allclose(km[[0, 2]], k[[0, 2]] / [1.1, 0.9])
    np.testing.assert_array_equal(mm[[0, 2]], [0.0, 1.0])
    assert j == 1 / (1.1**2 * 0.9)


def test_uniform_dilation_and_reference_units():
    k = np.geomspace(0.001, 1.0, 21)
    mu = np.linspace(0, 1, 21)
    km, mm, j = physical_ap_coordinates(k, mu, 1.03, 1.03)
    np.testing.assert_allclose(km, k / 1.03, rtol=3e-16)
    np.testing.assert_allclose(mm, mu, atol=2e-16)
    assert j == 1 / (1.03**3)
    physical = physical_ap_coordinates(k * 0.676, mu, 1.03, 1.03)[0]
    # Separate multiply/divide/hypot paths differ by a few float64 ulps.
    np.testing.assert_allclose(physical, 0.676 * km, rtol=8 * np.finfo(float).eps)


@pytest.mark.parametrize("qp,ql", [(0.0, 1.0), (1.0, -1.0), (np.nan, 1.0)])
def test_bad_dilations(qp, ql):
    with pytest.raises(ValueError):
        physical_ap_coordinates(0.1, 0.5, qp, ql)
