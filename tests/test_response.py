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
    """Check scalar response limits."""
    velocity_k = np.array([0.0, 1e-8, 0.01, 2 * np.pi / 100, 3 * np.pi / 100, 20.0])
    response = velocity_response(
        velocity_k, pixel_width_velocity=100, gaussian_sigma_velocity=3
    )
    oracle = [
        (math.sin(v * 50) / (v * 50) if v else 1) * math.exp(-0.5 * (v * 3) ** 2)
        for v in velocity_k
    ]
    assert_allclose(response, oracle, rtol=5e-14, atol=3e-16)
    assert response[0] == 1 and response[-1] == 0 and response[-2] < 0
    assert abs(response[3]) < 1e-15
    small = velocity_response(
        [1e-8], pixel_width_velocity=100, gaussian_sigma_velocity=0
    )
    assert_allclose(small, [1 - (5e-7) ** 2 / 6], rtol=0, atol=2e-16)
    assert_allclose(
        velocity_response(
            velocity_k, pixel_width_velocity=0, gaussian_sigma_velocity=0
        ),
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
    """Check field order pairs slices and units."""
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
    k_grid, mu_grid = np.array([0.1, 0.3, 0.2, 0.7]), np.array([0.0, 1.0, 0.5, 0.8])
    k_grid.flags.writeable = False
    response = prepare_response(fields, k_grid, mu_grid, a_v=100, settings=settings)
    assert np.all(response[:, 1] == 1) and np.all(response[0] == 1)
    assert not np.array_equal(response[:, 0], response[:, 2])
    with pytest.raises(ValueError):
        response.flags.writeable = True
    selection = PairSelection(fields, [("f2", "g"), ("f1", "f1")])
    products = pair_response(response, selection)
    oracle = np.array(
        [[row[i] * row[j] for i, j in selection.required_pairs] for row in response]
    )
    assert_allclose(products, oracle, rtol=0, atol=0)
    assert_allclose(
        prepare_response(fields, k_grid[::2], mu_grid[::2], a_v=100, settings=settings),
        response[::2],
        rtol=0,
        atol=0,
    )
    assert_allclose(
        pair_response(response[::2], selection), products[::2], rtol=0, atol=0
    )
    for n, f in enumerate(fields):
        item = settings[f.id]
        delta, sigma = width_velocity_to_comoving(
            [item.pixel_width_velocity, item.gaussian_sigma_velocity], a_v=100
        )
        oracle = np.sinc(k_grid * mu_grid * delta / (2 * np.pi)) * np.exp(
            -0.5 * (k_grid * mu_grid * sigma) ** 2
        )
        assert_allclose(response[:, n], oracle, rtol=5e-14)


@pytest.mark.parametrize("z", [2.0, 3.0, 4.0])
def test_width_conventions(z):
    """Check width conventions.

    Parameters
    ----------
    z : float
        Dimensionless redshift test input, supplied by pytest parametrization.
    """
    wavelength = 1215.67 * (1 + z)
    delta = pixel_width_angstrom_to_velocity(0.8, lambda_obs_angstrom=wavelength)
    sigma = gaussian_sigma_angstrom_to_velocity(0.6, lambda_obs_angstrom=wavelength)
    assert_allclose(
        [delta, sigma],
        [299792.458 * 0.8 / wavelength, 299792.458 * 0.6 / wavelength],
        rtol=5e-15,
    )
    resolving_power = 2500
    fwhm = 299792.458 / resolving_power
    assert_allclose(
        resolving_power_fwhm_to_sigma(resolving_power),
        gaussian_fwhm_velocity_to_sigma(fwhm),
        rtol=1e-15,
    )
    assert_allclose(
        gaussian_sigma_angstrom_to_velocity(
            wavelength / resolving_power / (2 * np.sqrt(2 * np.log(2))),
            lambda_obs_angstrom=wavelength,
        ),
        resolving_power_fwhm_to_sigma(resolving_power),
        rtol=2e-15,
    )
    legacy = legacy_resolving_power_to_sigma(resolving_power)
    assert_allclose(
        legacy / (299800 / resolving_power), SPEED_LIGHT_KMS / 299800, rtol=1e-15
    )
    assert_allclose(
        legacy / resolving_power_fwhm_to_sigma(resolving_power),
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
    """Check response invalid.

    Parameters
    ----------
    q : list
        Parametrized quadrature or wavenumber input, supplied by pytest
        parametrization.
    p : int or float
        Parametrized power or parameter input, supplied by pytest
        parametrization.
    s : int or float
        Parametrized scale input, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError):
        velocity_response(q, pixel_width_velocity=p, gaussian_sigma_velocity=s)


@pytest.mark.parametrize(
    "helper", [resolving_power_fwhm_to_sigma, legacy_resolving_power_to_sigma]
)
def test_bad_resolving_power(helper):
    """Check bad resolving power.

    Parameters
    ----------
    helper : callable
        Helper callable under examination, supplied by pytest parametrization.
    """
    for value in (0, -1, np.nan, np.inf, True, [10], 1e-320):
        with pytest.raises(ValueError):
            helper(value)


def test_bad_settings_shapes_and_widths():
    """Check bad settings shapes and widths."""
    observed_fields = [ObservedField("g", "galaxy", "g")]
    for settings in (
        {},
        {"g": (0, 0)},
        {"g": InstrumentResponse(0, 0), "extra": InstrumentResponse(0, 0)},
    ):
        with pytest.raises(ValueError):
            prepare_response(observed_fields, [1], [0], a_v=100, settings=settings)
    for k, mu in [([], []), ([1], [2]), ([1, 2], [0]), ([-1], [0]), ([1e308], [1])]:
        with pytest.raises(ValueError):
            prepare_response(
                observed_fields,
                k,
                mu,
                a_v=1e-308,
                settings={"g": InstrumentResponse(0, 0)},
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
        pair_response([[1, 2]], PairSelection(observed_fields))
