"""On-disk sigma8/fsigma8 cache of the native CAMB preparation."""

from types import SimpleNamespace

import numpy as np
import pytest

from fishhighz.cosmology import prepare_camb


class _Results:
    def __init__(self, redshifts):
        """Initialize the synthetic Results fixture.

        Parameters
        ----------
        redshifts : array_like of shape (n_redshifts,)
            Dimensionless redshifts at which the synthetic background is tabulated.

        Notes
        -----
        Sets the instance state used by the enclosing test; no scientific calculation is run.
        """
        self.redshifts = np.asarray(redshifts, dtype=float)
        self.Params = SimpleNamespace(
            Transfer=SimpleNamespace(PK_redshifts=self.redshifts.tolist())
        )

    def get_sigma8(self):
        """Return the synthetic sigma8.

        Returns
        -------
        values : float or ndarray
            Synthetic sigma8 in dimensionless units, matching the query or stored-
            redshift shape.
        """
        return 0.2 + 0.01 * self.redshifts

    def get_fsigma8(self):
        """Return the synthetic growth rate times sigma8.

        Returns
        -------
        values : float or ndarray
            Synthetic growth rate times sigma8 in dimensionless units, matching the
            query or stored-redshift shape.
        """
        return 0.1 + 0.005 * self.redshifts

    def hubble_parameter(self, redshift):
        """Return the synthetic Hubble parameter.

        Parameters
        ----------
        redshift : float or ndarray
            Dimensionless evaluation redshift; arrays retain their input shape.

        Returns
        -------
        values : float or ndarray
            Synthetic Hubble parameter in km/s/Mpc, matching the query or stored-
            redshift shape.
        """
        return 67.36 + 2.0 * np.asarray(redshift)

    def comoving_radial_distance(self, redshift):
        """Return the synthetic radial comoving distance.

        Parameters
        ----------
        redshift : float or ndarray
            Dimensionless evaluation redshift; arrays retain their input shape.

        Returns
        -------
        values : float or ndarray
            Synthetic radial comoving distance in Mpc, matching the query or stored-
            redshift shape.
        """
        return 1000.0 + 100.0 * np.asarray(redshift)


class _Background(_Results):
    def get_sigma8(self):
        """Return the synthetic sigma8.

        Raises
        ------
        AssertionError
            Deliberately raised to exercise the rejection path in the enclosing
            test.
        """
        raise AssertionError("a background-only surface has no growth arrays")

    get_fsigma8 = get_sigma8


class _FakeCamb:
    __version__ = "fake-1"

    def __init__(self):
        """Initialize the synthetic FakeCamb fixture.

        Notes
        -----
        Sets the instance state used by the enclosing test; no scientific calculation is run.
        """
        self.transfer_solves = 0
        self.background_solves = 0

    def read_ini(self, path):
        """Record a synthetic CAMB INI read and return minimal parameters.

        Parameters
        ----------
        path : pathlib.Path
            Path of the temporary test artifact to read or write.

        Returns
        -------
        parameters : types.SimpleNamespace
            Synthetic H0 and transfer-redshift settings.
        """
        return SimpleNamespace(
            H0=67.36,
            Transfer=SimpleNamespace(PK_redshifts=[2.3], PK_num_redshifts=1),
        )

    def get_results(self, parameters):
        """Record a transfer solve and construct synthetic CAMB results.

        Parameters
        ----------
        parameters : types.SimpleNamespace
            Synthetic CAMB parameters, including the requested transfer redshifts.

        Returns
        -------
        results : _Results
            Synthetic CAMB expansion, distance, and growth arrays.
        """
        self.transfer_solves += 1
        return _Results(parameters.Transfer.PK_redshifts)

    def get_background(self, parameters):
        """Record a background solve without transfer-growth information.

        Parameters
        ----------
        parameters : types.SimpleNamespace
            Synthetic CAMB parameters, including the requested transfer redshifts.

        Returns
        -------
        background : _Background
            Synthetic background-only result.
        """
        self.background_solves += 1
        return _Background(parameters.Transfer.PK_redshifts)


def _prepare(fake, ini, cache_dir, redshifts=(3.0, 2.0)):
    """Prepare the synthetic CAMB background with persistent caching.

    Parameters
    ----------
    fake : _FakeCamb
        Synthetic CAMB implementation recording background and transfer solves.
    ini : pathlib.Path
        Temporary CAMB INI input fixture.
    cache_dir : pathlib.Path
        Directory receiving the prepared-background cache.
    redshifts : array_like of shape (n_redshifts,), optional
        Dimensionless redshifts at which the synthetic background is tabulated.
        Default is (3.0, 2.0).

    Returns
    -------
    background : CAMBBackground
        Prepared expansion, distance, and growth evaluated through the fake CAMB
        module.
    """
    return prepare_camb(
        ini,
        redshifts=redshifts,
        template_growth_redshift=2.4,
        camb_module=fake,
        cache_dir=cache_dir,
    )


def _values(background):
    """Collect tabulated and interpolated synthetic background values.

    Parameters
    ----------
    background : CAMBBackground
        Prepared cosmological background with sampled expansion, distance, and
        growth.

    Returns
    -------
    values : tuple of ndarray
        Redshifts, H in km/s/Mpc, distances in Mpc, dimensionless sigma8 and
        growth rate, and interpolated H/distance.
    """
    nodes = np.array([2.1, 2.7, 3.3])
    return (
        background.redshifts,
        background.hubble_values,
        background.transverse_distance_values,
        background.sigma8_values,
        background.growth_rate_values,
        background.hubble_parameter(nodes),
        background.comoving_radial_distance(nodes),
    )


@pytest.fixture
def ini(tmp_path):
    """Write a minimal temporary CAMB parameter file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest for generated inputs and results.

    Returns
    -------
    path : pathlib.Path
        Path to params.ini in the pytest temporary directory.
    """
    path = tmp_path / "params.ini"
    path.write_text("ombh2 = 0.0224\n")
    return path


def test_hit_reproduces_miss_without_transfer_solve(ini, tmp_path, monkeypatch):
    """Check hit reproduces miss without transfer solve.

    Parameters
    ----------
    ini : pathlib.Path
        Temporary CAMB parameter file supplied by the ini fixture.
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    monkeypatch.delenv("FISHHIGHZ_CAMB_CACHE", raising=False)
    fake = _FakeCamb()
    first = _prepare(fake, ini, tmp_path / "cache")
    second = _prepare(fake, ini, tmp_path / "cache")
    assert fake.transfer_solves == 1 and fake.background_solves == 1
    assert first.cache["status"] == "miss, stored"
    assert second.cache["status"] == "hit"
    assert first.cache["key"] == second.cache["key"]
    for a, b in zip(_values(first), _values(second)):
        np.testing.assert_array_equal(a, b)


def test_key_changes_with_ini_redshifts_and_version(ini, tmp_path, monkeypatch):
    """Check key changes with ini redshifts and version.

    Parameters
    ----------
    ini : pathlib.Path
        Temporary CAMB parameter file supplied by the ini fixture.
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    monkeypatch.delenv("FISHHIGHZ_CAMB_CACHE", raising=False)
    fake = _FakeCamb()
    cache = tmp_path / "cache"
    keys = {_prepare(fake, ini, cache).cache["key"]}
    keys.add(_prepare(fake, ini, cache, redshifts=(3.0, 2.5)).cache["key"])
    ini.write_text("ombh2 = 0.0225\n")
    keys.add(_prepare(fake, ini, cache).cache["key"])
    fake.__version__ = "fake-2"
    keys.add(_prepare(fake, ini, cache).cache["key"])
    assert len(keys) == 4
    assert fake.transfer_solves == 4


def test_corrupt_entry_is_recomputed(ini, tmp_path, monkeypatch):
    """Check corrupt entry is recomputed.

    Parameters
    ----------
    ini : pathlib.Path
        Temporary CAMB parameter file supplied by the ini fixture.
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    monkeypatch.delenv("FISHHIGHZ_CAMB_CACHE", raising=False)
    fake = _FakeCamb()
    cache = tmp_path / "cache"
    reference = _prepare(fake, ini, cache)
    (cache / f"{reference.cache['key']}.npz").write_bytes(b"not an npz")
    again = _prepare(fake, ini, cache)
    assert again.cache["status"] == "miss, stored"
    assert fake.transfer_solves == 2
    np.testing.assert_array_equal(again.sigma8_values, reference.sigma8_values)


def test_disabled_and_implicit_test_surfaces_are_not_cached(ini, tmp_path, monkeypatch):
    """Check disabled and implicit test surfaces are not cached.

    Parameters
    ----------
    ini : pathlib.Path
        Temporary CAMB parameter file supplied by the ini fixture.
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    fake = _FakeCamb()
    monkeypatch.setenv("FISHHIGHZ_CAMB_CACHE", "0")
    assert _prepare(fake, ini, tmp_path / "cache").cache is None
    assert not (tmp_path / "cache").exists()
    monkeypatch.delenv("FISHHIGHZ_CAMB_CACHE")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    assert _prepare(fake, ini, None).cache is None
    assert not (tmp_path / "xdg").exists()
    assert fake.transfer_solves == 2


def test_unwritable_directory_still_returns_results(ini, tmp_path, monkeypatch):
    """Check unwritable directory still returns results.

    Parameters
    ----------
    ini : pathlib.Path
        Temporary CAMB parameter file supplied by the ini fixture.
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    monkeypatch.delenv("FISHHIGHZ_CAMB_CACHE", raising=False)
    blocker = tmp_path / "file"
    blocker.write_text("")
    background = _prepare(_FakeCamb(), ini, blocker / "cache")
    assert background.cache["status"] == "miss, unwritable"
    assert np.all(background.sigma8_values > 0)
