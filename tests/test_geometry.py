"""Independent dimensional, volume, optional curvature and boundary oracles."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fishhighz.geometry import (
    SPEED_LIGHT_KMS as C,
)
from fishhighz.geometry import (
    mode_counts,
    p1d_comoving_to_velocity,
    p1d_velocity_to_comoving,
    prepare_astropy_geometry,
    prepare_geometry,
    wavenumber_comoving_to_velocity,
    wavenumber_velocity_to_comoving,
    width_comoving_to_velocity,
    width_velocity_to_comoving,
)
from fishhighz.grids import IntegrationGrid, gauss_legendre_grid


def geometry(**kwargs):
    args = dict(
        z_eval=2.4,
        area_deg2=1000,
        h_fid=0.7,
        z_order=4,
        hubble=lambda z: np.full_like(z, 200),
        transverse_distance=lambda z: C * z / 200,
    )
    args.update(kwargs)
    return prepare_geometry(2, 3, **args)


@pytest.mark.parametrize("z", [2.0, 2.4, 3.0])
@pytest.mark.parametrize("h", [0.5, 0.7, 1.0])
def test_geometry_units_and_shell(z, h):
    g = geometry(z_eval=z, h_fid=h)
    omega = 1000 * (np.pi / 180) ** 2
    assert_allclose(g.volume, omega * h**3 * (C / 200) ** 3 * (27 - 8) / 3, rtol=5e-15)
    assert_allclose(g.a_v, 200 / ((1 + z) * h), rtol=2e-15)
    assert_allclose(g.d_deg, h * C * z / 200 * np.pi / 180, rtol=2e-15)
    assert g.speed_light_kms == C
    k, p = np.array([0.0, 0.1, 0.3]), np.array([2.0, 5.0, 8.0])
    q = wavenumber_comoving_to_velocity(k, a_v=g.a_v)
    pc = p1d_velocity_to_comoving(p, a_v=g.a_v)
    assert_allclose(k * pc, q * p, rtol=2e-15)
    assert_allclose(wavenumber_velocity_to_comoving(q, a_v=g.a_v), k)
    assert_allclose(p1d_comoving_to_velocity(pc, a_v=g.a_v), p)
    width = width_velocity_to_comoving(p, a_v=g.a_v)
    assert_allclose(width_comoving_to_velocity(width, a_v=g.a_v), p)


def test_h_scaling_and_modes():
    a, b = geometry(h_fid=0.4), geometry(h_fid=0.8)
    assert_allclose(
        [b.volume / a.volume, b.d_deg / a.d_deg, b.a_v / a.a_v], [8, 2, 0.5]
    )
    grid = gauss_legendre_grid([0.01, 0.1, 0.3], k_order=3, mu_order=3, h_fid=0.4)
    original = grid.q_mode.copy()
    assert_allclose(
        mode_counts(a, grid).sum(),
        a.volume * (0.3**3 - 0.01**3) / (6 * np.pi**2),
        rtol=3e-15,
    )
    custom = IntegrationGrid(
        grid.k, grid.w_k, grid.mu, grid.w_mu, k_min=0.01, k_max=0.3, h_fid=0.4
    )
    assert_allclose(mode_counts(a, custom), mode_counts(a, grid), rtol=0, atol=0)
    assert_allclose(grid.q_mode, original, rtol=0, atol=0)
    with pytest.raises(ValueError, match="h_fid"):
        mode_counts(b, grid)


def test_eds_convergence():
    h0 = 70

    def hubble(z):
        return h0 * (1 + z) ** 1.5

    def distance(z):
        return 2 * C / h0 * (1 - 1 / np.sqrt(1 + z))

    exact = (
        1000 * (np.pi / 180) ** 2 * 0.7**3 * (distance(3) ** 3 - distance(2) ** 3) / 3
    )
    errors = [
        abs(
            geometry(hubble=hubble, transverse_distance=distance, z_order=n).volume
            / exact
            - 1
        )
        for n in (4, 8, 16)
    ]
    assert errors[0] > errors[1]
    assert errors[2] <= max(errors[1], 3e-15)
    assert errors[-1] < 1e-10


def test_preparation_ownership():
    calls, outputs = [], []
    state = [200.0]

    def hubble(z):
        assert z.ndim == 1 and not z.flags.writeable
        with pytest.raises(ValueError):
            z.flags.writeable = True
        calls.append(z.copy())
        result = np.full_like(z, state[0])
        outputs.append(result)
        return result

    g = geometry(hubble=hubble)
    assert [len(z) for z in calls] == [4, 1]
    old = g.volume
    for output in outputs:
        output[:] = 9
    state[0] = 400
    grid = gauss_legendre_grid([0.1, 0.2], k_order=2, mu_order=2, h_fid=0.7)
    mode_counts(g, grid)
    assert g.volume == old and len(calls) == 2
    assert np.all(g.hubble_nodes == 200)
    for array in (g.z_nodes, g.w_z, g.hubble_nodes, g.transverse_distance_nodes):
        with pytest.raises(ValueError):
            array.flags.writeable = True
    with pytest.raises(FrozenInstanceError):
        g.volume = 0


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(z_eval=0),
        dict(z_eval=4),
        dict(z_eval=np.nan),
        dict(area_deg2=0),
        dict(area_deg2=42000),
        dict(area_deg2=1e-323),
        dict(h_fid=0),
        dict(h_fid=1e300),
        dict(h_fid=1e-300),
        dict(z_order=0),
        dict(z_order=True),
        dict(z_order=2.5),
        dict(hubble=lambda z: 200),
        dict(hubble=lambda z: np.zeros_like(z)),
        dict(hubble=lambda z: z + 1j),
        dict(hubble=lambda z: np.full_like(z, np.inf)),
        dict(transverse_distance=lambda z: np.ones((len(z), 1))),
        dict(transverse_distance=lambda z: -z),
        dict(transverse_distance=lambda z: np.full_like(z, 1e300)),
    ],
)
def test_geometry_invalid(kwargs):
    with pytest.raises(ValueError):
        geometry(**kwargs)


@pytest.mark.parametrize(
    "lo,hi", [(-1, 2), (2, 2), (3, 2), (2, np.nextafter(2.0, 3.0))]
)
def test_invalid_bin(lo, hi):
    with pytest.raises(ValueError):
        prepare_geometry(
            lo,
            hi,
            z_eval=2,
            area_deg2=1,
            h_fid=0.7,
            z_order=2,
            hubble=lambda z: z,
            transverse_distance=lambda z: z,
        )


@pytest.mark.parametrize(
    "convert",
    [
        wavenumber_comoving_to_velocity,
        wavenumber_velocity_to_comoving,
        p1d_velocity_to_comoving,
        p1d_comoving_to_velocity,
        width_velocity_to_comoving,
        width_comoving_to_velocity,
    ],
)
def test_conversion_errors(convert):
    inverse = convert in (
        wavenumber_velocity_to_comoving,
        p1d_comoving_to_velocity,
        width_comoving_to_velocity,
    )
    extremes = (
        [([1e308], 1e308), ([1e-308], 1e-308)]
        if inverse
        else [([1e308], 1e-308), ([1e-308], 1e308)]
    )
    for value, av in [([np.nan], 2), ([1.0], 0)] + extremes:
        with pytest.raises(ValueError):
            convert(value, a_v=av)


@pytest.mark.parametrize("curvature", [0.0, 0.15, -0.15])
def test_astropy(curvature):
    astropy = pytest.importorskip("astropy.cosmology")
    u = pytest.importorskip("astropy.units")
    model = (
        astropy.FlatLambdaCDM(H0=68, Om0=0.3)
        if curvature == 0
        else astropy.LambdaCDM(H0=68, Om0=0.3, Ode0=0.7 - curvature)
    )
    values = [
        prepare_astropy_geometry(
            model, 2, 3, z_eval=2.4, area_deg2=1000, h_fid=0.53, z_order=n
        )
        for n in (4, 8, 16)
    ]
    g = values[-1]
    assert g.h_fid != model.h
    assert_allclose(
        g.hubble_eval, model.H(2.4).to_value(u.km / u.s / u.Mpc), rtol=1e-15
    )
    assert_allclose(
        g.transverse_distance_eval,
        model.comoving_transverse_distance(2.4).to_value(u.Mpc),
        rtol=1e-15,
    )
    assert_allclose(g.hubble_nodes, model.H(g.z_nodes).to_value(u.km / u.s / u.Mpc))
    assert_allclose(
        g.transverse_distance_nodes,
        model.comoving_transverse_distance(g.z_nodes).to_value(u.Mpc),
    )
    exact = (
        g.solid_angle
        / (4 * np.pi)
        * g.h_fid**3
        * (model.comoving_volume(3) - model.comoving_volume(2)).to_value(u.Mpc**3)
    )
    errors = [abs(item.volume / exact - 1) for item in values]
    assert errors[1] < errors[0]
    assert errors[2] <= max(errors[1], 1e-13)
    assert errors[-1] < 1e-9
    with pytest.raises(ValueError, match="FLRW"):
        prepare_astropy_geometry(
            None, 2, 3, z_eval=2.4, area_deg2=1, h_fid=0.7, z_order=4
        )


def test_zero_lower_edge_and_mode_representability():
    g = prepare_geometry(
        0,
        1,
        z_eval=1,
        area_deg2=4 * np.pi / (np.pi / 180) ** 2,
        h_fid=0.7,
        z_order=3,
        hubble=lambda z: np.full_like(z, 100),
        transverse_distance=lambda z: 299792.458 * z / 100,
    )
    assert g.z_eval == g.z_max and np.all(g.z_nodes > 0)
    assert_allclose(
        g.volume, 4 * np.pi * 0.7**3 * (299792.458 / 100) ** 3 / 3, rtol=3e-15
    )
    enormous = geometry(transverse_distance=lambda z: np.full_like(z, 1e149))
    grid = gauss_legendre_grid([1e5, 2e5], k_order=2, mu_order=2, h_fid=0.7)
    with pytest.raises(ValueError, match="mode counts"):
        mode_counts(enormous, grid)
    tiny = geometry(transverse_distance=lambda z: np.full_like(z, 1e-150))
    grid = gauss_legendre_grid([1e-20, 2e-20], k_order=2, mu_order=2, h_fid=0.7)
    with pytest.raises(ValueError, match="mode counts"):
        mode_counts(tiny, grid)
