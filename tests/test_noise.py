"""Physical noise, normalization and independent matrix/Fisher oracles."""

import runpy
from pathlib import Path

import numpy as np
import pytest
from test_weights import FIELD, RESPONSE, geometry, prepare

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.noise import (
    forest_noise,
    galaxy_noise,
    local_galaxy_density,
    prepare_noise,
)
from fishhighz.response import InstrumentResponse
from fishhighz.weights import prepare_forest_weights


@pytest.mark.parametrize("v", [0, 7])
@pytest.mark.parametrize("h", [0.5, 0.7, 1.0])
def test_single_noise_and_h_units(v, h):
    """Check single noise and h units.

    Parameters
    ----------
    v : int
        Parametrized numeric input, including invalid values where specified,
        supplied by pytest parametrization.
    h : float
        Hubble parameter or finite-difference scale specified by the case,
        supplied by pytest parametrization.
    """
    bin_geometry = geometry(h=h)
    response = InstrumentResponse(100, 0)
    forest_weights = prepare_forest_weights(
        FIELD,
        bin_geometry,
        response,
        z_source=3,
        magnitudes=[20],
        quadrature=[0.2],
        rho=[15],
        variance=[v],
        length_velocity=10,
        method="supplied",
        weights=[0.07],
    )
    velocity_k = np.array([0, 2 * np.pi / 100, 3 * np.pi / 100, 0.002])
    noise = forest_noise(
        forest_weights,
        FIELD,
        bin_geometry,
        response,
        velocity_k * bin_geometry.a_v,
        np.ones(4),
        [2] * 4,
    )
    pixel_response_squared = np.array(
        [1, 0, (1 / (1.5 * np.pi)) ** 2, (np.sin(0.1) / 0.1) ** 2]
    )
    expected = (
        (2 * pixel_response_squared + 100 * v)
        * bin_geometry.d_deg**2
        / bin_geometry.a_v
        / 30
    )
    np.testing.assert_allclose(noise.total, expected, rtol=5e-14, atol=1e-28)
    np.testing.assert_allclose(
        noise.pixel, 100 * v * bin_geometry.d_deg**2 / bin_geometry.a_v / 30, rtol=5e-15
    )
    reference_geometry = geometry(h=0.7)
    np.testing.assert_allclose(
        noise.total,
        (2 * pixel_response_squared + 100 * v)
        * reference_geometry.d_deg**2
        / reference_geometry.a_v
        / 30
        * (h / 0.7) ** 3,
        rtol=5e-14,
        atol=1e-28,
    )
    for value in (noise.total, noise.aliasing, noise.pixel):
        with pytest.raises(ValueError):
            value.flags.writeable = True


def test_galaxy_local_density():
    """Check galaxy local density."""
    bin_geometry = geometry()
    galaxy_density = local_galaxy_density([3, 7, 10], [0.2, 0.5, 1.3], bin_geometry)
    expected = (
        (0.2 * 3 + 0.5 * 7 + 1.3 * 10)
        * (1 + bin_geometry.z_eval)
        / 299792.458
        * bin_geometry.a_v
        / bin_geometry.d_deg**2
    )
    assert galaxy_density == pytest.approx(expected, rel=5e-15)
    assert galaxy_noise(galaxy_density) == pytest.approx(1 / expected, rel=5e-15)
    assert galaxy_noise(2 * galaxy_density) == galaxy_noise(galaxy_density) / 2
    assert (
        local_galaxy_density([3, 7, 10], [0.2, 0.5, 1.3], geometry(area=200))
        == galaxy_density
    )


FIELDS = [
    FIELD,
    ObservedField("f2", "forest", "lya", "lbg"),
    ObservedField("g", "galaxy", "g"),
]


@pytest.mark.parametrize("permutation", [(0, 1, 2), (2, 0, 1), (1, 2, 0)])
def test_packing_psd_permutations_slices(permutation):
    """Check packing psd permutations slices.

    Parameters
    ----------
    permutation : tuple
        Index permutation applied to the fixture, supplied by pytest
        parametrization.
    """
    fields = [FIELDS[i] for i in permutation]
    selection = PairSelection(fields, [("f2", "g"), ("forest", "forest")])
    matrix = np.array([[4, -2, 1], [-2, 5, 0.5], [1, 0.5, 3]])
    matrix = matrix[np.ix_(permutation, permutation)]
    i, j = selection.required_pairs.T
    full = np.array([matrix[i, j] * scale for scale in [1, 2, 3, 4]])
    np.testing.assert_array_equal(prepare_noise(selection, 4, full=full), full)
    np.testing.assert_array_equal(
        prepare_noise(selection, 2, full=full[::2]), full[::2]
    )
    diagonal = {f.id: np.full(4, matrix[a, a]) for a, f in enumerate(fields)}
    generated = prepare_noise(
        selection, 4, diagonal=diagonal, independent_sampling=True
    )
    np.testing.assert_array_equal(
        generated, np.broadcast_to(np.diag(matrix.diagonal())[i, j], full.shape)
    )
    for invalid in [
        np.array([[1, 2, 0], [2, 1, 0], [0, 0, 1]]),
        np.array([[0, 1, 0], [1, 2, 0], [0, 0, 1]]),
    ]:
        with pytest.raises(ValueError, match="known noise"):
            prepare_noise(selection, 1, full=invalid[i, j][None, :])
    for valid in [np.zeros((3, 3)), np.outer([1, -2, 3], [1, -2, 3])]:
        prepare_noise(selection, 1, full=valid[i, j][None, :])


def test_subset_and_replacement():
    """Check subset and replacement."""
    selection = PairSelection(FIELDS, [("forest", "forest")])
    supplied = np.array([[2.0], [3.0]])
    np.testing.assert_array_equal(prepare_noise(selection, 2, full=supplied), supplied)
    np.testing.assert_array_equal(
        prepare_noise(
            selection, 2, diagonal={"forest": [2, 3]}, independent_sampling=True
        ),
        supplied,
    )
    invalid = [
        dict(diagonal={"forest": [2, 3]}),
        dict(full=supplied, diagonal={}),
        dict(full=supplied, independent_sampling=True),
        dict(full=np.ones((2, 2))),
        dict(diagonal={}, independent_sampling=True),
        dict(diagonal={"forest": [2, 3], "g": [1, 1]}, independent_sampling=True),
        dict(diagonal={"forest": 2}, independent_sampling=True),
    ]
    for kw in invalid:
        with pytest.raises(ValueError):
            prepare_noise(selection, 2, **kw)


@pytest.mark.parametrize("p1d", [[-1], [np.inf], [1, 2]])
def test_bad_p1d(p1d):
    """Check bad p1d.

    Parameters
    ----------
    p1d : list
        Intrinsic one-dimensional forest power, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError):
        forest_noise(prepare(), FIELD, geometry(), RESPONSE, [0.1], [0.5], p1d)


def test_zero_p1d_and_supplied_zero_weights():
    """Check zero p1d and supplied zero weights."""
    forest_weights = prepare(
        method="supplied", weights=[0, 1, 0], iterations=None, signal=None, alias=None
    )
    noise = forest_noise(forest_weights, FIELD, geometry(), RESPONSE, [0.1], [0], [0])
    assert noise.aliasing[0] == 0
    assert noise.total[0] == noise.pixel[0]


def test_example_fisher_freezing_and_area():
    """Check example fisher freezing and area."""
    run = runpy.run_path(str(Path(__file__).parents[1] / "examples/weighted_noise.py"))[
        "run"
    ]
    first, second = run(), run(area=2000)
    np.testing.assert_array_equal(first["noise"], second["noise"])
    np.testing.assert_allclose(
        second["amplitude_fisher"], 2 * np.array(first["amplitude_fisher"]), rtol=5e-14
    )
    assert first["null_directions"] == []
