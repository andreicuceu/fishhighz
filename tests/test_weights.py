"""Independent scalar checks of finite cumulative weighting and input contracts."""

import numpy as np
import pytest

from fishhighz.fields import ObservedField
from fishhighz.geometry import prepare_geometry
from fishhighz.response import InstrumentResponse
from fishhighz.weights import density_per_velocity, prepare_forest_weights


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
    """Prepare the standard three-magnitude legacy-weight fixture.

    Parameters
    ----------
    **kwargs : dict
        Overrides to prepare_forest_weights, including source samples, quadrature, variance, length, method, iteration count, and fixed auxiliary spectra.

    Returns
    -------
    weights : ForestWeights
        Source weights, cumulative integrals, and iteration diagnostics.
    """
    args = dict(
        z_source=3,
        magnitudes=[20, 21, 23],
        quadrature=[0.2, 0.5, 1.3],
        rho=np.array([1, 2, 4]) / [0.2, 0.5, 1.3],
        variance=[1, 4, 9],
        length_velocity=10,
        method="legacy",
        iterations=3,
        signal=3,
        alias=2,
    )
    args.update(kwargs)
    return prepare_forest_weights(FIELD, geometry(), RESPONSE, **args)


def oracle(r, v, count):
    """Evaluate the cumulative legacy recurrence with explicit scalar sums.

    Parameters
    ----------
    r : array_like of shape (n_magnitudes,)
        Source density times magnitude quadrature, in deg^-2 (km/s)^-1.
    v : array_like of shape (n_magnitudes,)
        Dimensionless pixel-noise variance for each magnitude sample.
    count : int
        Number of cumulative weight updates.

    Returns
    -------
    reference : tuple
        Final weights, three cumulative integrals, and maximum changes after
        each update.
    """
    source_weights = [2 / (2 + 0.5 * x) if mass else 0 for mass, x in zip(r, v)]
    changes = []
    for _ in range(count):
        next_w = []
        for m in range(len(r)):
            prefix = sum(r[j] * source_weights[j] for j in range(m + 1))
            next_w.append(3 / (3 + v[m] / (prefix * 10 / 0.5)) if r[m] else 0)
        changes.append(max(abs(a - b) for a, b in zip(source_weights, next_w)))
        source_weights = next_w
    integrals = [
        [
            sum(
                r[j] * source_weights[j] ** power * (v[j] if power == 3 else 1)
                for j in range(m + 1)
            )
            for m in range(len(r))
        ]
        for power in (1, 2)
    ]
    integrals.append(
        [
            sum(r[j] * source_weights[j] ** 2 * v[j] for j in range(m + 1))
            for m in range(len(r))
        ]
    )
    return source_weights, integrals, changes


@pytest.mark.parametrize("count", [0, 1, 2, 3, 6])
def test_scalar_prefix_oracle(count):
    """Check scalar prefix oracle.

    Parameters
    ----------
    count : int
        Number of iterations, samples, or records selected by this case,
        supplied by pytest parametrization.
    """
    result = prepare(iterations=count)
    source_weights, integrals, changes = oracle([1, 2, 4], [1, 4, 9], count)
    np.testing.assert_allclose(result.weights, source_weights, rtol=5e-15)
    np.testing.assert_allclose([result.I1, result.I2, result.I3], integrals, rtol=5e-15)
    np.testing.assert_allclose(result.weight_changes, changes, rtol=5e-14, atol=2e-16)
    assert result.A == pytest.approx(
        integrals[1][-1] / integrals[0][-1] ** 2 / 10, rel=5e-15
    )
    assert result.P_pixel == pytest.approx(
        integrals[2][-1] * 0.5 / integrals[0][-1] ** 2 / 10, rel=5e-15
    )


def test_distinguishes_wrong_iterations():
    """Check distinguishes wrong iterations."""
    source_weights = np.array(oracle([1, 2, 4], [1, 4, 9], 0)[0])
    total = 3 / (3 + np.array([1, 4, 9]) / (sum(source_weights * [1, 2, 4]) * 20))
    sweep = source_weights.copy()
    for m in range(3):
        sweep[m] = 3 / (
            3 + [1, 4, 9][m] / (sum(sweep[: m + 1] * np.array([1, 2, 4])[: m + 1]) * 20)
        )
    assert not np.allclose(prepare(iterations=1).weights, total)
    assert not np.allclose(prepare(iterations=1).weights, sweep)


@pytest.mark.parametrize("variance", [0, 7])
@pytest.mark.parametrize("weight", [0.01, 1, 19])
def test_single_population(variance, weight):
    """Check single population.

    Parameters
    ----------
    variance : int
        Pixel-noise variance test input, supplied by pytest parametrization.
    weight : int or float
        Source weight, supplied by pytest parametrization.
    """
    forest_weights = prepare(
        magnitudes=[21],
        quadrature=[0.3],
        rho=[10],
        variance=[variance],
        method="supplied",
        weights=[weight],
        iterations=None,
        signal=None,
        alias=None,
    )
    assert forest_weights.A == pytest.approx(1 / 30, rel=5e-15)
    assert forest_weights.P_pixel == pytest.approx(variance * 0.5 / 30, rel=5e-15)


def test_scaling_and_immutable():
    """Check scaling and immutable."""
    weight_options = dict(
        method="supplied",
        weights=[0.3, 0.7, 2],
        iterations=None,
        signal=None,
        alias=None,
    )
    forest_weights = prepare(**weight_options)
    for scale in [0.001, 100]:
        other = prepare(
            **dict(weight_options, weights=np.array(weight_options["weights"]) * scale)
        )
        np.testing.assert_allclose(
            [forest_weights.A, forest_weights.P_pixel],
            [other.A, other.P_pixel],
            rtol=5e-15,
        )
    other = prepare(**weight_options, rho=np.array([1, 2, 4]) / [0.2, 0.5, 1.3] * 2)
    np.testing.assert_allclose(
        [forest_weights.A, forest_weights.P_pixel],
        2 * np.array([other.A, other.P_pixel]),
        rtol=5e-15,
    )
    other = prepare(**weight_options, length_velocity=20)
    np.testing.assert_allclose(
        [forest_weights.A, forest_weights.P_pixel],
        2 * np.array([other.A, other.P_pixel]),
        rtol=5e-15,
    )
    for name in (
        "weights",
        "I1",
        "I2",
        "I3",
        "rho",
        "variance",
        "magnitudes",
        "quadrature",
        "weight_changes",
    ):
        with pytest.raises(ValueError):
            getattr(forest_weights, name).flags.writeable = True
    uniform = prepare(**dict(weight_options, weights=[1, 1, 1]), variance=[2, 2, 2])
    assert uniform.A == pytest.approx(1 / 70, rel=5e-15)
    assert uniform.P_pixel / uniform.A == pytest.approx(1, rel=5e-15)


@pytest.mark.parametrize("r,v", [([0, 2, 4], [9, 0, 1]), ([1, 0, 4], [0, 9, 1])])
def test_zero_support_prefixes(r, v):
    """Check zero support prefixes.

    Parameters
    ----------
    r : list
        Density times magnitude quadrature in deg^-2 (km/s)^-1.
    v : list
        Dimensionless pixel-noise variance at the same magnitude samples.
    """
    forest_weights = prepare(rho=np.array(r) / [0.2, 0.5, 1.3], variance=v)
    source_weights, ints, _ = oracle(r, v, 3)
    np.testing.assert_allclose(forest_weights.weights, source_weights, rtol=5e-15)
    np.testing.assert_allclose(
        [forest_weights.I1, forest_weights.I2, forest_weights.I3], ints, rtol=5e-15
    )
    assert np.all(forest_weights.weights[np.array(r) == 0] == 0)


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
        dict(iterations=True),
        dict(iterations=-1),
        dict(iterations=1.2),
        dict(signal=0),
        dict(alias=-1),
        dict(weights=[1, 1, 1]),
        dict(method="guess"),
        dict(rho=[1e308] * 3),
        dict(
            method="supplied",
            iterations=None,
            signal=None,
            alias=None,
            weights=[0, 0, 0],
        ),
        dict(
            method="supplied",
            iterations=None,
            signal=None,
            alias=None,
            weights=[1e200] * 3,
        ),
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
        dict(method="supplied", weights=[1, 1, 1]),
        dict(
            method="supplied",
            weights=[1e-300] * 3,
            iterations=None,
            signal=None,
            alias=None,
        ),
        dict(
            method="supplied",
            weights=[-1, 1, 1],
            iterations=None,
            signal=None,
            alias=None,
        ),
        dict(
            method="supplied", weights=[1, 1], iterations=None, signal=None, alias=None
        ),
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
