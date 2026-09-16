"""Scalar trigonometric and explicit instrument-unit response checks."""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.geometry import SPEED_LIGHT_KMS, width_velocity_to_comoving
from fishhighz.response import (
    InstrumentResponse,
    gaussian_fwhm_velocity_to_sigma,
    gaussian_sigma_angstrom_to_velocity,
    legacy_resolving_power_to_sigma,
    pair_response,
    pixel_width_angstrom_to_velocity,
    prepare_response,
    resolving_power_fwhm_to_sigma,
    velocity_response,
)


def test_scalar_response_limits():
    q = np.array([0.0, 1e-8, 0.01, 2 * np.pi / 100, 3 * np.pi / 100, 20.0])
    w = velocity_response(q, pixel_width_velocity=100, gaussian_sigma_velocity=3)
    oracle = [
        (math.sin(v * 50) / (v * 50) if v else 1) * math.exp(-0.5 * (v * 3) ** 2)
        for v in q
    ]
    assert_allclose(w, oracle, rtol=5e-14, atol=3e-16)
    assert w[0] == 1 and w[-1] == 0 and w[-2] < 0
    assert abs(w[3]) < 1e-15
    small = velocity_response(
        [1e-8], pixel_width_velocity=100, gaussian_sigma_velocity=0
    )
    assert_allclose(small, [1 - (5e-7) ** 2 / 6], rtol=0, atol=2e-16)
    assert_allclose(
        velocity_response(q, pixel_width_velocity=0, gaussian_sigma_velocity=0),
        1,
        rtol=0,
        atol=0,
    )
    assert (
        velocity_response([1e300], pixel_width_velocity=0, gaussian_sigma_velocity=0)[0]
        == 1
    )
    assert (
        velocity_response([1e200], pixel_width_velocity=0, gaussian_sigma_velocity=1)[0]
        == 0
    )


def test_field_order_pairs_slices_and_units():
    fields = [
        ObservedField("f1", "forest", "same", background="qso"),
        ObservedField("g", "galaxy", "g"),
        ObservedField("f2", "forest", "same", background="lbg"),
    ]
    settings = dict(
        f1=InstrumentResponse(100, 20),
        f2=InstrumentResponse(200, 30),
        g=InstrumentResponse(0, 0),
    )
    k, mu = np.array([0.1, 0.3, 0.2, 0.7]), np.array([0.0, 1.0, 0.5, 0.8])
    k.flags.writeable = False
    w = prepare_response(fields, k, mu, a_v=100, settings=settings)
    assert np.all(w[:, 1] == 1) and np.all(w[0] == 1)
    assert not np.array_equal(w[:, 0], w[:, 2])
    with pytest.raises(ValueError):
        w.flags.writeable = True
    selection = PairSelection(fields, [("f2", "g"), ("f1", "f1")])
    products = pair_response(w, selection)
    oracle = np.array(
        [[row[i] * row[j] for i, j in selection.required_pairs] for row in w]
    )
    assert_allclose(products, oracle, rtol=0, atol=0)
    assert_allclose(
        prepare_response(fields, k[::2], mu[::2], a_v=100, settings=settings),
        w[::2],
        rtol=0,
        atol=0,
    )
    assert_allclose(pair_response(w[::2], selection), products[::2], rtol=0, atol=0)
    for n, f in enumerate(fields):
        item = settings[f.id]
        delta, sigma = width_velocity_to_comoving(
            [item.pixel_width_velocity, item.gaussian_sigma_velocity], a_v=100
        )
        oracle = np.sinc(k * mu * delta / (2 * np.pi)) * np.exp(
            -0.5 * (k * mu * sigma) ** 2
        )
        assert_allclose(w[:, n], oracle, rtol=5e-14)


@pytest.mark.parametrize("z", [2.0, 3.0, 4.0])
def test_width_conventions(z):
    wavelength = 1215.67 * (1 + z)
    delta = pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength)
    sigma = gaussian_sigma_angstrom_to_velocity(0.6, lambda_obs_angstrom=wavelength)
    assert_allclose(
        [delta, sigma],
        [299792.458 * 0.8 / wavelength, 299792.458 * 0.6 / wavelength],
        rtol=5e-15,
    )
    r = 2500
    fwhm = 299792.458 / r
    assert_allclose(
        resolving_power_fwhm_to_sigma(r),
        gaussian_fwhm_velocity_to_sigma(fwhm),
        rtol=1e-15,
    )
    assert_allclose(
        gaussian_sigma_angstrom_to_velocity(
            wavelength / r / (2 * np.sqrt(2 * np.log(2))),
            lambda_obs_angstrom=wavelength,
        ),
        resolving_power_fwhm_to_sigma(r),
        rtol=2e-15,
    )
    legacy = legacy_resolving_power_to_sigma(r)
    assert_allclose(legacy / (299800 / r), SPEED_LIGHT_KMS / 299800, rtol=1e-15)
    assert_allclose(
        legacy / resolving_power_fwhm_to_sigma(r),
        2 * np.sqrt(2 * np.log(2)),
        rtol=1e-15,
    )


@pytest.mark.parametrize(
    "q,p,s",
    [
        ([], 1, 1),
        ([1j], 1, 1),
        ([np.nan], 1, 1),
        ([-1], 1, 1),
        ([[1]], 1, 1),
        ([1], -1, 1),
        ([1], 1, np.inf),
        ([1e308], 1e308, 1),
    ],
)
def test_response_invalid(q, p, s):
    with pytest.raises(ValueError):
        velocity_response(q, pixel_width_velocity=p, gaussian_sigma_velocity=s)


@pytest.mark.parametrize(
    "helper", [resolving_power_fwhm_to_sigma, legacy_resolving_power_to_sigma]
)
def test_bad_resolving_power(helper):
    for value in (0, -1, np.nan, np.inf, True, [10], 1e-320):
        with pytest.raises(ValueError):
            helper(value)


def test_bad_settings_shapes_and_widths():
    f = [ObservedField("g", "galaxy", "g")]
    for settings in (
        {},
        {"g": (0, 0)},
        {"g": InstrumentResponse(0, 0), "extra": InstrumentResponse(0, 0)},
    ):
        with pytest.raises(ValueError):
            prepare_response(f, [1], [0], a_v=100, settings=settings)
    for k, mu in [([], []), ([1], [2]), ([1, 2], [0]), ([-1], [0]), ([1e308], [1])]:
        with pytest.raises(ValueError):
            prepare_response(
                f, k, mu, a_v=1e-308, settings={"g": InstrumentResponse(0, 0)}
            )
    for bad in (-1, np.nan, np.inf, True):
        with pytest.raises(ValueError):
            InstrumentResponse(bad, 0)
        with pytest.raises(ValueError):
            gaussian_fwhm_velocity_to_sigma(bad)
    for helper in (
        pixel_width_angstrom_to_velocity,
        gaussian_sigma_angstrom_to_velocity,
    ):
        for width, wave in [(-1, 2), (1, 0), (np.inf, 1), (1e308, 1e-308)]:
            with pytest.raises(ValueError):
                helper(width, lambda_obs_angstrom=wave)
    with pytest.raises(ValueError):
        pair_response([[1, 2]], PairSelection(f))
