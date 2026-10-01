"""Portable generated-FITS and independent analytic template checks.

Install .[dev,templates] for development: missing extras must not silently skip
these acceptance tests. No sibling checkout or cosmological asset is required.
"""

import builtins
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from scipy.interpolate import CubicSpline

from fishhighz.models.templates import load_template, prepare_template

K_GRID = np.exp(np.array([-3.0, -2.6, -1.8, -1.1, -0.7, 0.0]))


def functions(k, derivative=0):
    """Evaluate signed polynomial template components in log wavenumber.

    Parameters
    ----------
    k : ndarray of shape (n_nodes,)
        Comoving wavenumbers in h/Mpc.
    derivative : int, optional
        Derivative selector: zero returns power; nonzero returns its first
        derivative with respect to k. Default is 0.

    Returns
    -------
    components : ndarray of shape (n_nodes, 2)
        Smooth and wiggle power in (Mpc/h)^3, or first k derivatives.
    """
    log_k = np.log(k)
    if derivative:
        return np.column_stack(
            (
                (0.3 - 0.4 * log_k + 0.6 * log_k * log_k) / k,
                (1 - 0.2 * log_k - 0.3 * log_k * log_k) / k,
            )
        )
    return np.column_stack(
        (
            -1 + 0.3 * log_k - 0.2 * log_k * log_k + 0.2 * log_k**3,
            0.5 + log_k - 0.1 * log_k * log_k - 0.1 * log_k**3,
        )
    )


def prepared(**kwargs):
    """Prepare the standard signed template fixture.

    Parameters
    ----------
    **kwargs : dict
        prepare_template options, including interpolation policy and template provenance metadata.

    Returns
    -------
    template : PowerTemplate
        Synthetic smooth and wiggle splines.
    """
    parts = functions(K_GRID)
    return prepare_template(
        K_GRID,
        parts.sum(axis=1),
        parts[:, 0],
        z_ref=2.4,
        h_template=0.7,
        h_fid=0.7,
        **kwargs,
    )


def table(*, arrays=None, formats=None, units=None, header=True, omit=()):
    """Build a synthetic Vega-format template FITS table.

    Parameters
    ----------
    arrays : mapping or None, optional
        Numerical column overrides keyed by K, PK, or PKSB. Default is None.
    formats : mapping or None, optional
        FITS column format overrides keyed by column name. Default is None.
    units : mapping or None, optional
        FITS column unit overrides keyed by column name. Default is None.
    header : bool, optional
        Whether to include synthetic cosmological template metadata. Default is
        True.
    omit : sequence of str, optional
        Names of template columns omitted to test validation. Default is ().

    Returns
    -------
    hdu : astropy.io.fits.BinTableHDU
        K, PK, and PKSB columns with requested metadata and test overrides.
    """
    components = functions(K_GRID)
    arrays = arrays or {}
    formats = formats or {}
    units = units or {}
    columns = [
        fits.Column(
            name=name,
            format=formats.get(name, "D"),
            array=arrays.get(name, value),
            unit=units.get(name),
        )
        for name, value in (
            ("K", K_GRID),
            ("PK", components.sum(axis=1)),
            ("PKSB", components[:, 0]),
        )
        if name not in omit
    ]
    result = fits.BinTableHDU.from_columns(columns, name="PK")
    if header:
        result.header["ZREF"] = 2.4
        result.header["H0"] = 70.0
        result.header["F_ZREF"] = 0.9
        result.header["OM"] = 0.3
    return result


def write(tmp_path, hdu=None, *, extras=()):
    """Write a temporary FITS file with one synthetic template extension.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest for generated inputs and results.
    hdu : astropy.io.fits.BinTableHDU or None, optional
        Template table; None builds the standard synthetic table. Default is
        None.
    extras : sequence of astropy.io.fits.hdu.base._BaseHDU, optional
        Additional HDUs appended after the synthetic template. Default is ().

    Returns
    -------
    path : pathlib.Path
        Path to the generated toy.fits file.
    """
    path = tmp_path / "toy.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(np.ones(3), name="UNRELATED"),
            table() if hdu is None else hdu,
            *extras,
        ]
    ).writeto(path)
    return path


def test_fits_snapshot_named_selection_metadata_and_signed_arrays(
    tmp_path, monkeypatch
):
    """Check fits snapshot named selection metadata and signed arrays.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    path = write(tmp_path)
    content = path.read_bytes()
    original_read = Path.read_bytes
    reads = []

    def read_snapshot(self):
        """Read the FITS snapshot, then replace its file to test single-read behavior.

        Returns
        -------
        snapshot : bytes
            Original FITS bytes before the temporary file is replaced.
        """
        reads.append(self)
        value = original_read(self)
        # Change backing file after snapshot; hash and parse must still agree.
        self.write_bytes(b"replacement is not FITS")
        return value

    monkeypatch.setattr(Path, "read_bytes", read_snapshot)
    result = load_template(path, h_fid=0.7)
    assert reads == [path.resolve()]
    assert result.source_sha256 == hashlib.sha256(content).hexdigest()
    assert result.source_path == str(path.resolve())
    assert result.z_ref == 2.4 and result.h_template == result.h_fid == 0.7
    assert result.metadata["H0"] == 70 and result.metadata["OM"] == 0.3
    assert result.component_names == ("smooth", "wiggle")
    assert np.any(result.components[:, 0] < 0) and np.any(result.components[:, 1] < 0)
    np.testing.assert_array_equal(result.k_file, K_GRID)
    np.testing.assert_allclose(
        result.evaluate(K_GRID).sum(axis=1), result.pk_file, atol=2e-14
    )
    for name in (
        "k_file",
        "pk_file",
        "pksb_file",
        "k",
        "full",
        "components",
        "log_k",
        "coefficients",
    ):
        array = getattr(result, name)
        assert (
            array.flags.owndata
            and array.flags.c_contiguous
            and not array.flags.writeable
        )
        assert array.dtype == np.float64
    with pytest.raises(TypeError):
        result.metadata["OM"] = 0.9
    # Source file is now invalid: evaluations must use the prepared arrays only.
    monkeypatch.setattr(fits, "open", lambda *a, **kw: pytest.fail("FITS reopened"))
    monkeypatch.setattr(
        "scipy.interpolate.CubicSpline",
        lambda *a, **kw: pytest.fail("spline prepared again"),
    )
    monkeypatch.setattr(
        CubicSpline, "__call__", lambda *a, **kw: pytest.fail("SciPy evaluation")
    )
    np.testing.assert_allclose(result.evaluate(K_GRID), functions(K_GRID), atol=2e-14)
    np.testing.assert_allclose(
        result.evaluate(K_GRID, derivative=1), functions(K_GRID, 1), atol=3e-13
    )


@pytest.mark.parametrize("kind", ["missing", "ambiguous", "image", "ascii"])
def test_invalid_named_table(tmp_path, kind):
    """Check invalid named table.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    kind : str
        Tracer, input, or calculation classification for this case, supplied by
        pytest parametrization.
    """
    if kind == "missing":
        hdu = fits.ImageHDU(name="OTHER")
    elif kind == "image":
        hdu = fits.ImageHDU(name="PK")
    elif kind == "ascii":
        hdu = fits.TableHDU.from_columns(
            [fits.Column(name="K", format="D", array=K_GRID)], name="PK"
        )
    else:
        hdu = table()
    path = write(tmp_path, hdu, extras=[table()] if kind == "ambiguous" else [])
    with pytest.raises(ValueError, match="unique named PK binary table"):
        load_template(path, h_fid=0.7)


@pytest.mark.parametrize("name", ["K", "PK", "PKSB"])
def test_missing_column(tmp_path, name):
    """Check missing column.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    name : str
        Named quantity or policy under examination, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError, match=f"one {name} column"):
        load_template(write(tmp_path, table(omit=[name])), h_fid=0.7)


@pytest.mark.parametrize(
    "format_,data",
    [
        ("2D", np.ones((6, 2))),
        ("C", np.ones(6, dtype=complex)),
        ("L", np.ones(6, dtype=bool)),
        ("4A", np.array(["bad"] * 6)),
        ("PD()", np.array([np.ones(i + 1) for i in range(6)], dtype=object)),
        ("D", np.array([1, 2, np.nan, 4, 5, 6])),
    ],
)
def test_invalid_fits_column_type_and_shape(tmp_path, format_, data):
    """Check invalid fits column type and shape.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    format_ : str
        FITS column format, supplied by pytest parametrization.
    data : ndarray
        Numerical input fixture, supplied by pytest parametrization.
    """
    hdu = table(arrays={"PKSB": data}, formats={"PKSB": format_})
    with pytest.raises(ValueError, match="PKSB|scalar 1D"):
        load_template(write(tmp_path, hdu), h_fid=0.7)


@pytest.mark.parametrize(
    "k",
    [
        [0.1, 0.2, 0.3],
        [0, 0.2, 0.3, 0.4],
        [0.1, 0.3, 0.2, 0.4],
        [0.1, 0.2, 0.2, 0.4],
        [0.1, 0.2, np.inf, 0.4],
    ],
)
def test_invalid_fits_knots(tmp_path, k):
    """Check invalid fits knots.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    k : list
        Wavenumber test input in the convention stated by the tested function,
        supplied by pytest parametrization.
    """
    hdu = table(arrays={"K": k, "PK": np.ones(len(k)), "PKSB": np.ones(len(k))})
    with pytest.raises(ValueError, match="K|four knots"):
        load_template(write(tmp_path, hdu), h_fid=0.7)


@pytest.mark.parametrize(
    "k,pk,smooth",
    [
        (K_GRID, np.ones(5), np.ones(6)),
        (K_GRID[:, None], np.ones(6), np.ones(6)),
        (K_GRID, np.ones(6, dtype=bool), np.ones(6)),
        (K_GRID, np.ones(6), np.ones(6, dtype=object)),
    ],
)
def test_array_contract(k, pk, smooth):
    """Check array contract.

    Parameters
    ----------
    k : ndarray
        Wavenumber test input in the convention stated by the tested function,
        supplied by pytest parametrization.
    pk : ndarray
        Full three-dimensional template power, supplied by pytest
        parametrization.
    smooth : ndarray
        Smooth template component, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError):
        prepare_template(k, pk, smooth, z_ref=2.4, h_template=0.7, h_fid=0.7)


@pytest.mark.parametrize(
    "units",
    [
        {"K": "h/Mpc", "PK": "(Mpc/h)^3", "PKSB": "Mpc^3/h^3"},
        {"K": "h Mpc**-1", "PK": "Mpc**3 h**-3", "PKSB": "(Mpc/h)**3"},
        {"K": "h Mpc-1", "PK": "Mpc3/h3", "PKSB": "(Mpc/h)3"},
    ],
)
def test_explicit_equivalent_units(tmp_path, units):
    """Check explicit equivalent units.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    units : dict
        Input unit metadata, supplied by pytest parametrization.
    """
    result = load_template(write(tmp_path, table(units=units)), h_fid=0.7)
    np.testing.assert_array_equal(result.components, prepared().components)


@pytest.mark.parametrize(
    "name,unit",
    [("K", "1/Mpc"), ("PK", "Mpc^3"), ("PKSB", "km/s"), ("K", "h/Mpc_wrong")],
)
def test_inconsistent_units(tmp_path, name, unit):
    """Check inconsistent units.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    name : str
        Named quantity or policy under examination, supplied by pytest
        parametrization.
    unit : str
        Input unit string, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError, match=f"{name}: unsupported/inconsistent units"):
        load_template(write(tmp_path, table(units={name: unit})), h_fid=0.7)


def test_required_metadata_and_precedence(tmp_path):
    """Check required metadata and precedence.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    hdu = table(header=False)
    path = write(tmp_path, hdu)
    for kwargs, missing in (({}, "ZREF"), ({"z_ref": 2.4}, "H0")):
        with pytest.raises(ValueError, match=f"missing PK header {missing}"):
            load_template(path, h_fid=0.7, **kwargs)
    result = load_template(path, z_ref=2.4, h_template=0.7, h_fid=0.7)
    assert result.z_ref == 2.4 and "ZREF" not in result.metadata
    # Only the table header counts; a primary value does not silently fill it.
    with fits.open(path, mode="update") as hdus:
        hdus[0].header["H0"] = 70
        hdus[0].header["ZREF"] = 2.4
    with pytest.raises(ValueError, match="missing PK header ZREF"):
        load_template(path, h_fid=0.7)


@pytest.mark.parametrize("kwargs", [{"z_ref": 2.3}, {"h_template": 0.71}])
def test_conflicting_metadata(tmp_path, kwargs):
    """Check conflicting metadata.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    kwargs : dict
        Keyword arguments selecting the parametrized case, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError, match="contradictory"):
        load_template(write(tmp_path), h_fid=0.7, **kwargs)


def test_consistent_metadata_roundoff_and_optional_provenance(tmp_path):
    """Check consistent metadata roundoff and optional provenance.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    result = load_template(
        write(tmp_path),
        h_fid=0.7,
        z_ref=np.nextafter(2.4, 3),
        h_template=np.nextafter(0.7, 1),
    )
    assert result.z_ref == 2.4 and result.h_template == 0.7
    first_template = prepared(metadata={"F_ZREF": 0.9, "SIGMA8_ZREF": 0.3})
    second_template = prepared(metadata={"F_ZREF": 0.1, "SIGMA8_ZREF": 10})
    np.testing.assert_array_equal(
        first_template.evaluate(K_GRID), second_template.evaluate(K_GRID)
    )


@pytest.mark.parametrize(
    "key,value",
    [("ZREF", -1), ("H0", 0), ("H0", -1), ("ZREF", "unknown"), ("H0", True)],
)
def test_bad_header_metadata(tmp_path, key, value):
    """Check bad header metadata.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    key : str
        Dictionary or configuration key under examination, supplied by pytest
        parametrization.
    value : bool or int or str
        Value at the tested validation boundary, supplied by pytest
        parametrization.
    """
    hdu = table()
    hdu.header[key] = value
    with pytest.raises(ValueError, match=key):
        load_template(write(tmp_path, hdu), h_fid=0.7)


def test_conversion_and_derivative_physical_units():
    """Check conversion and derivative physical units."""
    source = prepared()
    converted = prepare_template(
        K_GRID, source.pk_file, source.pksb_file, z_ref=2.4, h_template=0.7, h_fid=0.5
    )
    q_source = np.exp(np.linspace(np.log(K_GRID[0]), np.log(K_GRID[-1]), 51))
    q_fid = q_source * (0.7 / 0.5)
    np.testing.assert_allclose(converted.k * 0.5, K_GRID * 0.7, rtol=2e-16)
    np.testing.assert_allclose(
        converted.full / 0.5**3, source.pk_file / 0.7**3, rtol=4e-16
    )
    for derivative in (0, 1):
        physical_new = converted.evaluate(q_fid, derivative=derivative) / 0.5 ** (
            3 + derivative
        )
        physical_old = source.evaluate(q_source, derivative=derivative) / 0.7 ** (
            3 + derivative
        )
        np.testing.assert_allclose(physical_new, physical_old, rtol=5e-14, atol=3e-12)
    np.testing.assert_allclose(
        converted.evaluate(converted.k).sum(axis=1), converted.full, atol=1e-14
    )


@pytest.mark.parametrize(
    "h_template,h_fid",
    [(1e-300, 1e300), (1e300, 1e-300), (0, 0.7), (0.7, np.inf), (0.7, True)],
)
def test_invalid_or_unrepresentable_conversion(h_template, h_fid):
    """Check invalid or unrepresentable conversion.

    Parameters
    ----------
    h_template : int or float
        Template reference Hubble parameter, supplied by pytest parametrization.
    h_fid : bool or float
        Forecast reference Hubble parameter, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError):
        prepare_template(
            K_GRID, np.ones(6), np.zeros(6), z_ref=0, h_template=h_template, h_fid=h_fid
        )


def test_decomposition_and_coefficients_cannot_overflow():
    """Check decomposition and coefficients cannot overflow."""
    with pytest.raises(ValueError, match="decomposition"):
        prepare_template(
            K_GRID,
            np.full(6, 1e308),
            np.full(6, -1e308),
            z_ref=0,
            h_template=1,
            h_fid=1,
        )
    with pytest.raises(ValueError, match="coefficient"):
        prepare_template(
            K_GRID,
            np.array([1, -1, 1, -1, 1, -1]) * 1e308,
            np.zeros(6),
            z_ref=0,
            h_template=1,
            h_fid=1,
        )
    tiny_grid = 1 + np.arange(4) * np.spacing(1.0)
    huge_grid = tiny_grid * 1e100
    with pytest.raises(ValueError, match="ln"):
        prepare_template(
            huge_grid, np.ones(4), np.ones(4), z_ref=0, h_template=1, h_fid=1
        )


def test_exact_cubic_values_derivatives_endpoints_and_knots():
    """Check exact cubic values derivatives endpoints and knots."""
    result = prepared()
    near = np.concatenate(
        [np.nextafter(K_GRID[1:-1], 0), K_GRID, np.nextafter(K_GRID[1:-1], np.inf)]
    )
    queries = np.concatenate([near, np.exp(np.linspace(-2.99, -0.01, 71))])[::-1]
    for order, atol in ((0, 2e-14), (1, 3e-13)):
        np.testing.assert_allclose(
            result.evaluate(queries, derivative=order),
            functions(queries, order),
            rtol=2e-13,
            atol=atol,
        )
        left = result.evaluate(np.nextafter(K_GRID[1:-1], 0), derivative=order)
        right = result.evaluate(np.nextafter(K_GRID[1:-1], np.inf), derivative=order)
        np.testing.assert_allclose(left, right, rtol=2e-13, atol=atol)


def oscillatory(k, derivative=0):
    """Evaluate a damped oscillatory template and its analytic derivative.

    Parameters
    ----------
    k : ndarray of shape (n_nodes,)
        Comoving wavenumbers in h/Mpc.
    derivative : int, optional
        Derivative selector: zero returns power; nonzero returns its first
        derivative with respect to k. Default is 0.

    Returns
    -------
    components : ndarray of shape (n_nodes, 2)
        Smooth and wiggle power in (Mpc/h)^3, or first k derivatives.
    """
    damping = np.exp(-((k / 0.4) ** 2))
    if derivative:
        return np.column_stack(
            (
                -5000 / (1 + 5 * k) ** 2,
                50
                * damping
                * (110 * np.cos(110 * k) - 2 * k / 0.4**2 * np.sin(110 * k)),
            )
        )
    return np.column_stack((1000 / (1 + 5 * k), 50 * damping * np.sin(110 * k)))


def test_oscillatory_sampling_convergence():
    """Check oscillatory sampling convergence."""
    query = np.linspace(0.020013, 0.499987, 1501)
    errors = []
    for count in (48, 192, 768):
        k_grid = np.geomspace(0.02, 0.5, count)
        components = oscillatory(k_grid)
        template = prepare_template(
            k_grid,
            components.sum(axis=1),
            components[:, 0],
            z_ref=0,
            h_template=0.7,
            h_fid=0.7,
        )
        errors.append(
            [
                np.max(
                    np.abs(
                        template.evaluate(query, derivative=order)
                        - oscillatory(query, order)
                    ),
                    axis=0,
                )
                for order in (0, 1)
            ]
        )
    errors = np.array(errors)  # grid, derivative order, component
    # Absolute errors are meaningful through zero crossings; no division by P_w.
    assert np.all(errors[1:] < errors[:-1] / 2)
    assert errors[0, 0, 1] > 1  # coarse off-knot sampling demonstrably inadequate
    # Fourfold refinement resolves the endpoint-dominated derivative error.
    assert errors[-1, 0, 1] < 1e-3
    assert errors[-1, 1, 1] < 3
    assert errors[-1, 0, 0] < 2e-9
    assert errors[-1, 1, 0] < 4e-6


def test_independent_scipy_coefficients_and_full_sum():
    """Check independent scipy coefficients and full sum."""
    rng = np.random.default_rng(701)
    k_grid = np.exp(np.array([-4.0, -3.7, -2.4, -1.1, -0.4, 0.0]))
    components = rng.normal(size=(6, 2))
    template = prepare_template(
        k_grid, components.sum(axis=1), components[:, 0], z_ref=0, h_template=1, h_fid=1
    )
    reference = CubicSpline(np.log(k_grid), components, axis=0, bc_type="not-a-knot")
    full = CubicSpline(np.log(k_grid), components.sum(axis=1), bc_type="not-a-knot")
    query = np.exp(np.linspace(-3.99, -0.01, 79))
    for order in (0, 1):
        expected = reference(np.log(query), order) / query[:, None] ** order
        np.testing.assert_allclose(
            template.evaluate(query, derivative=order), expected, rtol=3e-13, atol=4e-13
        )
        np.testing.assert_allclose(
            template.evaluate(query, derivative=order).sum(axis=1),
            full(np.log(query), order) / query**order,
            rtol=4e-13,
            atol=4e-13,
        )


def test_redshift_power_amplitude_once_and_immutable():
    """Check redshift power amplitude once and immutable."""
    template = prepared()
    coefficients = template.coefficients.copy()
    for order in (0, 1):
        expected = template.evaluate(K_GRID, derivative=order)
        np.testing.assert_array_equal(
            template.evaluate(K_GRID, z=template.z_ref, growth=1, derivative=order),
            expected,
        )
        np.testing.assert_array_equal(
            template.evaluate(K_GRID, z=3, growth=0.36, derivative=order),
            0.36 * expected,
        )
    np.testing.assert_array_equal(template.coefficients, coefficients)
    with pytest.raises(ValueError, match="explicit growth"):
        template.evaluate(K_GRID, z=np.nextafter(template.z_ref, 3))
    with pytest.raises(ValueError, match="equal 1"):
        template.evaluate(K_GRID, growth=0.5)
    np.testing.assert_allclose(
        template.evaluate(K_GRID, growth=np.nextafter(1.0, 2)),
        template.evaluate(K_GRID),
        rtol=5e-16,
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"z": -1},
        {"z": np.nan},
        {"z": True},
        {"growth": 0},
        {"growth": -1},
        {"growth": np.inf},
        {"growth": 1j},
        {"derivative": 2},
        {"derivative": True},
        {"derivative": 1.0},
    ],
)
def test_bad_evaluation_controls(kwargs):
    """Check bad evaluation controls.

    Parameters
    ----------
    kwargs : dict
        Keyword arguments selecting the parametrized case, supplied by pytest
        parametrization.
    """
    with pytest.raises(ValueError):
        prepared().evaluate(K_GRID, **kwargs)


@pytest.mark.parametrize(
    "query",
    [
        [],
        [[0.1]],
        [0],
        [-1],
        [np.nan],
        [np.inf],
        [True],
        [0.1j],
        np.array([0.1], dtype=object),
        [np.nextafter(K_GRID[0], 0)],
        [np.nextafter(K_GRID[-1], np.inf)],
    ],
)
def test_invalid_queries_and_strict_domain(query):
    """Check invalid queries and strict domain.

    Parameters
    ----------
    query : ndarray or list
        Coordinate query, supplied by pytest parametrization.
    """
    with pytest.raises(ValueError):
        prepared().evaluate(query)


def test_order_slices_ownership_and_finite_arithmetic():
    """Check order slices ownership and finite arithmetic."""
    template = prepared()
    k_grid = np.array(
        [K_GRID[-1], 0.3, K_GRID[0], 0.3, K_GRID[2], 0.3, K_GRID[2], 0.3]
    )[::2]
    k_grid.flags.writeable = False
    for order in (0, 1):
        result = template.evaluate(k_grid, derivative=order)
        joined = np.concatenate(
            [
                template.evaluate(k_grid[:2], derivative=order),
                template.evaluate(k_grid[2:], derivative=order),
            ]
        )
        np.testing.assert_array_equal(result, joined)
        assert (
            result.flags.owndata
            and result.flags.c_contiguous
            and result.dtype == np.float64
        )
        result[:] = 0
        np.testing.assert_allclose(
            template.evaluate(k_grid, derivative=order),
            functions(k_grid, order),
            atol=3e-13,
        )
    with pytest.raises(ValueError, match="requested k range.*template domain"):
        template.evaluate([np.nextafter(K_GRID[0], 0)])
    with pytest.raises(ValueError, match="nonfinite template evaluation"):
        template.evaluate(K_GRID, z=0, growth=1e308)


def test_preparation_copies_caller_arrays():
    """Check preparation copies caller arrays."""
    k_grid = K_GRID.copy()
    parts = functions(k_grid)
    pk, smooth = parts.sum(axis=1), parts[:, 0].copy()
    metadata = {"OM": 0.3}
    template = prepare_template(
        k_grid, pk, smooth, z_ref=0, h_template=1, h_fid=1, metadata=metadata
    )
    k_grid[:], pk[:], smooth[:] = 9, 0, 0
    metadata["OM"] = 0.9
    assert template.metadata["OM"] == 0.3
    np.testing.assert_allclose(template.evaluate(K_GRID), parts, atol=2e-14)


def test_quiet_lazy_module_imports(tmp_path):
    """Check quiet lazy module imports.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    script = """
import sys
import fishhighz
import fishhighz.models.templates
import fishhighz.models.external
import fishhighz.derivatives
import fishhighz.fisher
assert not {'astropy', 'scipy', 'vega', 'camb', 'lyaforecast'} & set(sys.modules)
"""
    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == result.stderr == ""


@pytest.mark.parametrize("missing", ["astropy", "scipy"])
def test_missing_extras_actionable(tmp_path, monkeypatch, missing):
    """Check missing extras actionable.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    missing : str
        Omitted input or metadata item, supplied by pytest parametrization.
    """
    path = write(tmp_path)
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        """Reject the selected optional import and forward all other imports.

        Parameters
        ----------
        name : str
            Name of the artifact, module, or result under examination.
        *args : tuple
            Positional arguments forwarded to the original callable or accepted by
            the test callback.
        **kwargs : dict
            Keyword options forwarded to the original callable or inspected by the
            test callback.

        Returns
        -------
        module : module
            Module returned by the original import operation.

        Raises
        ------
        ModuleNotFoundError
            Deliberately raised to exercise the rejection path in the enclosing
            test.
        """
        if name == missing or name.startswith(missing + "."):
            raise ModuleNotFoundError(missing)
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match=r"install 'fishhighz\[templates\]'"):
        if missing == "astropy":
            load_template(path, h_fid=0.7)
        else:
            prepared()
    # In-memory preparation never needs Astropy.
    if missing == "astropy":
        np.testing.assert_allclose(
            prepared().evaluate(K_GRID), functions(K_GRID), atol=2e-14
        )
