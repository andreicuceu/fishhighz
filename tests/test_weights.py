"""Scalar checks of adaptive forest weighting, noise moments and input contracts."""

import numpy as np
import pytest

from fishhighz.fields import ObservedField
from fishhighz.geometry import _immutable, prepare_geometry
from fishhighz.kernels.full_sum_weights import _integrals
from fishhighz.response import InstrumentResponse
from fishhighz.weights import (
    ForestWeights,
    _context,
    density_per_velocity,
    prepare_forest_weights,
)


def geometry(h=0.7, area=100):
    """Prepare a synthetic common-volume redshift-bin geometry.

    Parameters
    ----------
    h : float, optional
        Dimensionless fiducial Hubble parameter H0/(100 km/s/Mpc). Default is
        0.7.
    area : float, optional
        Survey area in deg^2. Default is 100.

    Returns
    -------
    geometry : BinGeometry
        Distances in Mpc/h, volume in (Mpc/h)^3, and velocity conversion in km/s
        per Mpc/h.
    """
    return prepare_geometry(
        2,
        3,
        z_eval=2.4,
        area_deg2=area,
        h_fid=h,
        hubble=lambda z: np.full_like(z, 250),
        transverse_distance=lambda z: np.full_like(z, 5500),
        z_order=4,
    )


FIELD = ObservedField("forest", "forest", "lya", "qso")
RESPONSE = InstrumentResponse(0.5, 0)


def prepare(**kwargs):
    """Prepare the standard three-magnitude adaptive early-lyaforecast fixture.

    Parameters
    ----------
    **kwargs : dict
        Overrides to prepare_forest_weights, including source samples,
        quadrature, variance, length, method, stopping controls, and the fixed
        reference spectra S (signal) and B (alias).

    Returns
    -------
    weights : ForestWeights
        Converged source weights, full-sample integrals, noise coefficients and
        convergence record.
    """
    args = dict(
        z_source=3,
        magnitudes=[20, 21, 23],
        quadrature=[0.2, 0.5, 1.3],
        rho=np.array([1, 2, 4]) / [0.2, 0.5, 1.3],
        variance=[1, 4, 9],
        length_velocity=10,
        method="early_lyaforecast",
        signal=3,
        alias=2,
    )
    args.update(kwargs)
    return prepare_forest_weights(FIELD, geometry(), RESPONSE, **args)


def inject_forest_weights(
    weights,
    *,
    field=FIELD,
    bin_geometry=None,
    response=RESPONSE,
    z_source=3,
    magnitudes=(20, 21, 23),
    quadrature=(0.2, 0.5, 1.3),
    rho=None,
    variance=(1, 4, 9),
    length_velocity=10,
):
    """Construct ForestWeights from prescribed weights, bypassing the solver.

    Downstream noise and forecast tests use this to impose exactly known
    weights. The moments and noise coefficients are evaluated with the same
    kernel as prepare_forest_weights; no convergence record is attached.

    Parameters
    ----------
    weights : array_like of shape (n_magnitude,)
        Dimensionless prescribed forest weights.
    field : ObservedField, optional
        Forest field identity. Default is the module FIELD.
    bin_geometry : BinGeometry, optional
        Forest evaluation geometry. Default None uses ``geometry()``.
    response : InstrumentResponse, optional
        Instrumental response with positive pixel width. Default is RESPONSE.
    z_source : float, optional
        Dimensionless source redshift. Default is 3.
    magnitudes, quadrature : array_like of shape (n_magnitude,), optional
        Source magnitudes and positive magnitude quadrature weights.
    rho : array_like of shape (n_magnitude,), optional
        Source density in deg^-2 (km/s)^-1 mag^-1. Default None gives the
        standard fixture density [1, 2, 4] / quadrature.
    variance : array_like of shape (n_magnitude,), optional
        Dimensionless pixel-noise variance.
    length_velocity : float, optional
        Forest length in km/s. Default is 10.

    Returns
    -------
    prepared : ForestWeights
        Immutable weights with full-sample I1, I2, I3 in deg^-2 (km/s)^-1, A in
        deg^2 and P_pixel in deg^2 km/s.
    """
    if bin_geometry is None:
        bin_geometry = geometry()
    quadrature = np.asarray(quadrature, dtype=float)
    if rho is None:
        rho = np.array([1, 2, 4]) / quadrature
    rho = np.asarray(rho, dtype=float)
    variance = np.asarray(variance, dtype=float)
    weights = np.asarray(weights, dtype=float)

    masses = rho * quadrature
    first_moment, second_moment, noise_moment, aliasing_coefficient, pixel_power = (
        _integrals(
            masses, weights, variance, length_velocity, response.pixel_width_velocity
        )
    )

    result = object.__new__(ForestWeights)
    values = dict(
        context=_context(field, bin_geometry, response),
        z_source=float(z_source),
        magnitudes=np.asarray(magnitudes, dtype=float),
        quadrature=quadrature,
        rho=rho,
        variance=variance,
        length_velocity=float(length_velocity),
        method="prescribed",
        auxiliary=None,
        signal=None,
        alias=None,
        weights=weights,
        I1=float(first_moment),
        I2=float(second_moment),
        I3=float(noise_moment),
        A=float(aliasing_coefficient),
        P_pixel=float(pixel_power),
        convergence=None,
    )
    for name, value in values.items():
        object.__setattr__(
            result, name, _immutable(value) if isinstance(value, np.ndarray) else value
        )
    return result


@pytest.mark.parametrize("method", ["early_lyaforecast", "mcdonald"])
@pytest.mark.parametrize("variance", [0, 7])
def test_single_population(variance, method):
    """Check that one source population gives weight-independent A and P_pixel.

    For a single magnitude sample A = 1/(rho dm L) and
    P_pixel = variance pixel/(rho dm L), whatever the converged weight.

    Parameters
    ----------
    variance : int
        Pixel-noise variance test input, supplied by pytest parametrization.
    method : str
        Adaptive weighting prescription, supplied by pytest parametrization.
    """
    forest_weights = prepare(
        magnitudes=[21],
        quadrature=[0.3],
        rho=[10],
        variance=[variance],
        method=method,
    )
    assert forest_weights.convergence["status"] == "converged"
    assert forest_weights.A == pytest.approx(1 / 30, rel=5e-15)
    assert forest_weights.P_pixel == pytest.approx(variance * 0.5 / 30, rel=5e-15)


def test_moment_scaling():
    """Check weight-amplitude invariance and density/length scaling of A and P_pixel."""
    masses = np.array([1.0, 2.0, 4.0])
    variance = np.array([1.0, 4.0, 9.0])
    source_weights = np.array([0.3, 0.7, 2.0])
    length, pixel = 10.0, 0.5
    *_, aliasing, pixel_power = _integrals(
        masses, source_weights, variance, length, pixel
    )

    # A and P_pixel are invariant under a global rescaling of the weights.
    for scale in [0.001, 100]:
        *_, other_aliasing, other_pixel_power = _integrals(
            masses, source_weights * scale, variance, length, pixel
        )
        np.testing.assert_allclose(
            [aliasing, pixel_power], [other_aliasing, other_pixel_power], rtol=5e-15
        )

    # Both scale inversely with the source density and with the forest length.
    for doubled_masses, doubled_length in [(2 * masses, length), (masses, 2 * length)]:
        *_, other_aliasing, other_pixel_power = _integrals(
            doubled_masses, source_weights, variance, doubled_length, pixel
        )
        np.testing.assert_allclose(
            [aliasing, pixel_power],
            2 * np.array([other_aliasing, other_pixel_power]),
            rtol=5e-15,
        )

    # Uniform weights and variance: A = 1/(sum(rho dm) L), P_pixel/A = v pixel.
    *_, uniform_aliasing, uniform_pixel_power = _integrals(
        masses, np.ones(3), np.full(3, 2.0), length, pixel
    )
    assert uniform_aliasing == pytest.approx(1 / 70, rel=5e-15)
    assert uniform_pixel_power / uniform_aliasing == pytest.approx(1, rel=5e-15)


def test_immutable():
    """Check that the prepared arrays are read-only."""
    forest_weights = prepare()
    for name in (
        "weights",
        "rho",
        "variance",
        "magnitudes",
        "quadrature",
    ):
        with pytest.raises(ValueError):
            getattr(forest_weights, name).flags.writeable = True


@pytest.mark.parametrize("method", ["early_lyaforecast", "mcdonald"])
@pytest.mark.parametrize("r,v", [([0, 2, 4], [9, 0, 1]), ([1, 0, 4], [0, 9, 1])])
def test_zero_density_samples_do_not_contribute(r, v, method):
    """Check that zero-density magnitudes leave weights and coefficients unchanged.

    Parameters
    ----------
    r : list
        Density times magnitude quadrature in deg^-2 (km/s)^-1.
    v : list
        Dimensionless pixel-noise variance at the same magnitude samples.
    method : str
        Adaptive weighting prescription, supplied by pytest parametrization.
    """
    quadrature = np.array([0.2, 0.5, 1.3])
    forest_weights = prepare(rho=np.array(r) / quadrature, variance=v, method=method)
    support_mask = np.array(r) > 0
    supported = prepare(
        magnitudes=np.array([20, 21, 23])[support_mask],
        quadrature=quadrature[support_mask],
        rho=(np.array(r) / quadrature)[support_mask],
        variance=np.array(v)[support_mask],
        method=method,
    )
    np.testing.assert_allclose(
        forest_weights.weights[support_mask], supported.weights, rtol=5e-15
    )
    np.testing.assert_allclose(
        [forest_weights.I1, forest_weights.I2, forest_weights.I3],
        [supported.I1, supported.I2, supported.I3],
        rtol=5e-15,
    )
    np.testing.assert_allclose(
        [forest_weights.A, forest_weights.P_pixel],
        [supported.A, supported.P_pixel],
        rtol=5e-15,
    )


@pytest.mark.parametrize(
    "change",
    [
        dict(z_source=2.4),
        dict(magnitudes=[21, 20, 23]),
        dict(magnitudes=[]),
        dict(quadrature=[1, 0, 1]),
        dict(rho=[0, 0, 0]),
        dict(rho=[-1, 1, 1]),
        dict(variance=[1, np.nan, 1]),
        dict(signal=0),
        dict(alias=-1),
        dict(signal=None),
        dict(method="guess"),
        dict(method="legacy"),
        dict(method="inverse_variance"),
        dict(method="supplied"),
        dict(max_updates=3),
        dict(rho=[1e308] * 3),
    ],
)
def test_invalid(change):
    """Check invalid.

    Parameters
    ----------
    change : dict
        Input override exercising the specified validation boundary, supplied by
        pytest parametrization.
    """
    with pytest.raises(ValueError):
        prepare(**change)


def test_density_conversion_and_context():
    """Check density conversion and context."""
    row = np.array([0, 10, 3.2])
    np.testing.assert_allclose(
        density_per_velocity(row, z_source=3.1), row * 4.1 / 299792.458, rtol=5e-15
    )
    forest_weights = prepare()
    forest_weights.validate_context(FIELD, geometry(area=200), RESPONSE)
    for f, g, response in [
        (FIELD, geometry(h=0.8), RESPONSE),
        (FIELD, geometry(), InstrumentResponse(1, 0)),
        (ObservedField("other", "forest", "lya", "qso"), geometry(), RESPONSE),
    ]:
        with pytest.raises(ValueError, match="mismatch"):
            forest_weights.validate_context(f, g, response)


@pytest.mark.parametrize(
    "row,z",
    [
        ([-1], 3),
        ([np.nan], 3),
        ([1], -1),
        ([1], np.inf),
        ([], 3),
        ([[1]], 3),
        ([1e-323], 0),
    ],
)
def test_density_invalid(row, z):
    """Check density invalid.

    Parameters
    ----------
    row : list
        Density-table row, supplied by pytest parametrization.
    z : int or float
        Dimensionless redshift test input, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError):
        density_per_velocity(row, z_source=z)


@pytest.mark.parametrize(
    "change",
    [
        dict(variance=[-1, 1, 1]),
        dict(quadrature=[1, 1]),
        dict(quadrature=[1e-300] * 3, rho=[1e-300] * 3),
        dict(length_velocity=0),
        dict(magnitudes=[1, np.inf, 3]),
        dict(alias=np.inf),
        dict(signal=np.nan),
    ],
)
def test_additional_invalid_boundaries(change):
    """Check additional invalid boundaries.

    Parameters
    ----------
    change : dict
        Input override exercising the specified validation boundary, supplied by
        pytest parametrization.
    """
    with pytest.raises(ValueError):
        prepare(**change)
