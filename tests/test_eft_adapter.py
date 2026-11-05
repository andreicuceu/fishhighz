"""Analytic AP, covariance and rank checks for the optional full EFT route."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fishhighz.covariance import gaussian_covariance
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.models.eft import EFTModelAdapter, component_coordinates


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


@pytest.fixture(scope="module")
def eft_adapter(request):
    """Prepare a real EFT adapter on a small deterministic output grid.

    Parameters
    ----------
    request : pytest.FixtureRequest
        Indirect parameter specifying ``tree`` or ``loop`` order.

    Returns
    -------
    EFTModelAdapter
        Adapter with synthetic power in (Mpc/h)^3 at z=2.3, FFTLog resolution
        128, eight mu samples and twelve output k samples in h/Mpc.
    """
    pytest.importorskip("eftkernels")
    k_grid = np.geomspace(5e-5, 100.0, 512)
    smooth_power = 1e4 * k_grid**0.96 / (1 + (k_grid / 0.15) ** 2) ** 2
    full_power = smooth_power * (
        1 + 0.03 * np.sin(100 * k_grid) * np.exp(-((k_grid / 0.3) ** 2))
    )
    fiducial = dict(
        k=k_grid,
        pk_full=full_power,
        pk_nw=smooth_power,
        growth_rate=0.97,
        rd_mpc_h=100.0,
        z_fiducial=2.3,
    )
    defaults = dict(
        bias_F=-0.14,
        bias_eta_F=-0.21,
        bias_G=2.4,
        b2_F=0.05,
        bG2_F=-0.02,
        bKKPar_F=-0.02,
        bPi2Par_F=0.02,
        bPi3Par_F=0.01,
        bKPi2Par_F=0.01,
        b2_G=0.5,
        bG2_G=-0.3,
        bGamma3_G=0.2,
    )
    return EFTModelAdapter(
        fiducial,
        ("f",),
        defaults,
        order=request.param,
        k_output=np.geomspace(0.02, 0.2, 12),
        num_mu=8,
        nmax=128,
    )


@pytest.mark.parametrize("eft_adapter", ["tree", "loop"], indirect=True)
def test_eft_adapter_native_nodes_and_pairs(eft_adapter):
    """Compare FF, FG and GG with direct EFT components at identity AP.

    Parameters
    ----------
    eft_adapter : EFTModelAdapter
        Real tree or loop adapter using synthetic input power.
    """
    adapter = eft_adapter
    # Interior native nodes avoid endpoint roundoff in the AP coordinate map.
    mu_grid, k_grid = np.meshgrid(adapter.mu, adapter.k[1:-1], indexing="ij")
    growth = adapter.fiducial["growth_rate"]
    pairs = ((1, 1), (1, 0), (0, 0), (0, 1))
    actual = adapter([growth], 2.3, k_grid.ravel(), mu_grid.ravel(), pairs)
    expected = np.column_stack(
        [
            (
                adapter.models[index].compute(adapter.defaults, component="smooth")
                + adapter.models[index].compute(adapter.defaults, component="wiggle")
            )[:, 1:-1].ravel()
            for index in (2, 1, 0, 1)
        ]
    )
    assert actual.shape == (k_grid.size, len(pairs))
    assert actual.dtype == np.float64
    assert np.isfinite(actual).all()
    # Use the high-level EFT spectrum regression tolerances, including the
    # adapter's interpolation at native nodes and the tree K1 factorization.
    assert_allclose(actual, expected, rtol=1e-8, atol=1e-7)
    assert_allclose(actual[:, 1], actual[:, 3], rtol=0, atol=0)


@pytest.mark.parametrize("eft_adapter", ["loop"], indirect=True)
def test_eft_adapter_growth_restoration(eft_adapter):
    """Check that a growth update changes all pairs and restores their means.

    Parameters
    ----------
    eft_adapter : EFTModelAdapter
        Real loop adapter using synthetic input power.
    """
    adapter = eft_adapter
    k_grid = np.array([0.04, 0.08, 0.15])
    mu_grid = np.array([0.2, 0.6, 0.9])
    growth = adapter.fiducial["growth_rate"]
    pairs = ((0, 0), (0, 1), (1, 1))
    fiducial = adapter([growth], 2.3, k_grid, mu_grid, pairs)
    perturbed = adapter([growth + 0.02], 2.3, k_grid, mu_grid, pairs)
    assert np.isfinite(perturbed).all()
    assert np.all(np.max(np.abs(perturbed - fiducial), axis=0) > 1e-7)
    restored = adapter([growth], 2.3, k_grid, mu_grid, pairs)
    assert_allclose(restored, fiducial, rtol=1e-8, atol=1e-7)
