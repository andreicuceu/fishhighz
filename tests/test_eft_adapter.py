"""Analytic AP, covariance and rank checks for the optional full EFT route."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.interpolate import CubicSpline

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
    # Include both boundaries: identity AP must preserve them exactly.
    mu_grid, k_grid = np.meshgrid(adapter.mu, adapter.k, indexing="ij")
    growth = adapter.fiducial["growth_rate"]
    pairs = ((1, 1), (1, 0), (0, 0), (0, 1))
    actual = adapter([growth], 2.3, k_grid.ravel(), mu_grid.ravel(), pairs)
    expected = np.column_stack(
        [
            (
                adapter.models[index].compute(adapter.defaults, component="smooth")
                + adapter.models[index].compute(adapter.defaults, component="wiggle")
            ).ravel()
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
def test_eft_adapter_growth_restoration(eft_adapter, monkeypatch):
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
    fiducial_models = adapter.models
    perturbed = adapter([growth + 0.02], 2.3, k_grid, mu_grid, pairs)
    assert np.isfinite(perturbed).all()
    assert np.all(np.max(np.abs(perturbed - fiducial), axis=0) > 1e-7)
    from eftkernels import EFTModel

    def unexpected_update(*args, **kwargs):
        raise AssertionError(
            "Returning to fiducial growth must reuse its prepared loops"
        )

    monkeypatch.setattr(EFTModel, "update_fiducial", unexpected_update)
    restored = adapter([growth], 2.3, k_grid, mu_grid, pairs)
    assert adapter.models is fiducial_models
    assert_allclose(restored, fiducial, rtol=1e-8, atol=1e-7)


def test_component_isotropic_boundaries():
    """Preserve exact k and mu at identity and isotropic dilations."""
    k_grid = np.array([0.005, 1.12, 1.12])
    mu_grid = np.array([0.0, 0.00029, 1.0])
    for alpha in (1.0, 1.03):
        mapped_k, mapped_mu, jacobian = component_coordinates(
            k_grid, mu_grid, alpha, 1.0
        )
        np.testing.assert_array_equal(mapped_k, k_grid / alpha)
        np.testing.assert_array_equal(mapped_mu, mu_grid)
        assert_allclose(jacobian, alpha**-3, rtol=2e-15)


@pytest.mark.parametrize("eft_adapter", ["tree", "loop"], indirect=True)
@pytest.mark.parametrize("shared_ap", [False, True])
def test_eft_nonunit_components(eft_adapter, shared_ap):
    """Compare dilated components with direct EFT evaluations at mapped nodes."""
    adapter = eft_adapter
    parameters = dict(
        adapter.defaults, alpha_s=1.03, phi_s=0.96, alpha_w=0.98, phi_w=1.04
    )
    if shared_ap:
        parameters.update(alpha_w=parameters["alpha_s"], phi_w=parameters["phi_s"])
    # Use a denser output grid than the small identity-regression fixture. The
    # 2e-6 bound below tests interpolation, not FFTLog or physical accuracy.
    dense = EFTModelAdapter(
        adapter.fiducial,
        tuple(parameters),
        parameters,
        order=adapter.order,
        k_output=np.geomspace(0.02, 0.2, 120),
        num_mu=32,
        nmax=128,
    )
    k_grid = np.array([0.04, 0.075, 0.15])
    mu_grid = np.array([0.15, 0.55, 0.9])
    pairs = ((0, 0), (0, 1), (1, 1))
    actual = dense(list(parameters.values()), 2.3, k_grid, mu_grid, pairs)
    expected = np.zeros_like(actual)
    for component, suffix in (("smooth", "s"), ("wiggle", "w")):
        alpha, phi = parameters["alpha_" + suffix], parameters["phi_" + suffix]
        # Independent analytic alpha/phi map and volume, without using the
        # adapter helper. Dense angular sampling supplies the direct reference.
        parallel_scale, transverse_scale = alpha / np.sqrt(phi), alpha * np.sqrt(phi)
        dilation = np.sqrt(
            mu_grid**2 / parallel_scale**2 + (1 - mu_grid**2) / transverse_scale**2
        )
        mapped_k, mapped_mu = k_grid * dilation, mu_grid / parallel_scale / dilation
        jacobian = 1 / (parallel_scale * transverse_scale**2)
        for column, model in enumerate(dense.models):
            model.set_output_grid(mapped_k, 256)
            power = model.compute(parameters, component=component)
            direct = CubicSpline(model.muk_grid[:, 0], power, axis=0)(mapped_mu)
            expected[:, column] += jacobian * np.diag(direct)
    assert_allclose(actual, expected, rtol=2e-6, atol=1e-7)


@pytest.mark.parametrize("eft_adapter", ["tree", "loop"], indirect=True)
def test_eft_stochastic_polynomial_and_scales(eft_adapter):
    """Map each pair's scaled stochastic polynomial with the smooth AP map."""
    adapter = eft_adapter
    defaults = dict(
        adapter.defaults, alpha_s=1.03, phi_s=0.96, alpha_w=0.98, phi_w=1.04
    )
    scales = {}
    for index, label in enumerate(("FF", "FG", "GG")):
        for suffix, value in (("0", 1.0), ("k2", -0.4), ("mu2", 0.3)):
            name = f"stoch_{label}_{suffix}"
            defaults[name] = value * (index + 1)
            scales[name] = 2.5 + index
    adapter_with_noise = EFTModelAdapter(
        adapter.fiducial,
        tuple(defaults),
        defaults,
        scales,
        order=adapter.order,
        k_output=adapter.k,
        num_mu=len(adapter.mu),
        nmax=128,
        k_natural=0.3,
    )
    no_stochastic = dict(defaults)
    for name in scales:
        no_stochastic[name] = 0.0
    k_grid, mu_grid = np.array([0.04, 0.08, 0.15]), np.array([0.0, 0.6, 1.0])
    pairs = ((0, 0), (1, 0), (1, 1), (0, 1))
    baseline = adapter_with_noise(
        list(no_stochastic.values()), 2.3, k_grid, mu_grid, pairs
    )
    actual = adapter_with_noise(list(defaults.values()), 2.3, k_grid, mu_grid, pairs)
    parallel = k_grid * mu_grid / (defaults["alpha_s"] / np.sqrt(defaults["phi_s"]))
    transverse = (
        k_grid
        * np.sqrt(1 - mu_grid**2)
        / (defaults["alpha_s"] * np.sqrt(defaults["phi_s"]))
    )
    mapped_k = np.hypot(parallel, transverse)
    mapped_mu = parallel / mapped_k
    jacobian = defaults["alpha_s"] ** -3 * defaults["phi_s"] ** -0.5
    expected = np.column_stack(
        [
            jacobian
            * (index + 1)
            * (2.5 + index)
            * (1 + (mapped_k / 0.3) ** 2 * (-0.4 + 0.3 * mapped_mu**2))
            for index in (0, 1, 2, 1)
        ]
    )
    assert_allclose(actual - baseline, expected, rtol=1e-12, atol=1e-10)


@pytest.mark.parametrize("eft_adapter", ["tree", "loop"], indirect=True)
def test_eft_invalid_nodes_and_ap_bounds(eft_adapter):
    """Reject unsupported angular nodes and AP excursions rather than clamp."""
    adapter = eft_adapter
    growth = adapter.fiducial["growth_rate"]
    for mu in (-0.1, 1.01, np.nan):
        with pytest.raises(ValueError, match=r"mu must be on \[0, 1\]"):
            adapter([growth], 2.3, np.array([0.1]), np.array([mu]), ((0, 0),))
    defaults = dict(adapter.defaults, alpha_s=0.9, alpha_w=1.0)
    bounded = EFTModelAdapter(
        adapter.fiducial,
        tuple(defaults),
        defaults,
        order=adapter.order,
        k_output=adapter.k,
        num_mu=len(adapter.mu),
        nmax=128,
    )
    with pytest.raises(ValueError, match="outside padded EFT output grid"):
        bounded(
            list(defaults.values()),
            2.3,
            np.array([adapter.k[-1]]),
            np.array([0.6]),
            ((0, 0),),
        )
    defaults.update(alpha_s=1.0, alpha_w=1.1)
    with pytest.raises(ValueError, match="outside padded EFT output grid"):
        bounded(
            list(defaults.values()),
            2.3,
            np.array([adapter.k[0]]),
            np.array([0.6]),
            ((0, 0),),
        )


def test_missing_template_keys():
    """Fail before any expensive EFT work when the template lacks PKNW."""
    with pytest.raises(ValueError, match="missing required keys:.*pk_nw.*PKNW"):
        EFTModelAdapter(
            dict(k=[], pk_full=[], growth_rate=1, rd_mpc_h=100, z_fiducial=2.3), (), {}
        )


@pytest.mark.parametrize("eft_adapter", ["loop"], indirect=True)
def test_reversed_pair_compute_count(eft_adapter, monkeypatch):
    """Compute and spline-fit each unique pair only once per component."""
    adapter = eft_adapter
    from eftkernels import EFTModel

    from fishhighz.models import eft

    compute_calls, spline_calls = [], []
    original_compute, original_spline = EFTModel.compute, eft.RectBivariateSpline

    def compute(model, *args, **kwargs):
        compute_calls.append(
            (model.tracer1_name, model.tracer2_name, kwargs["component"])
        )
        return original_compute(model, *args, **kwargs)

    def spline(*args, **kwargs):
        spline_calls.append(1)
        return original_spline(*args, **kwargs)

    monkeypatch.setattr(EFTModel, "compute", compute)
    monkeypatch.setattr(eft, "RectBivariateSpline", spline)
    adapter(
        [adapter.fiducial["growth_rate"]],
        2.3,
        np.array([0.08]),
        np.array([0.6]),
        ((0, 1), (1, 0), (0, 1)),
    )
    assert compute_calls == [("F", "G", "smooth"), ("F", "G", "wiggle")]
    assert len(spline_calls) == 2
