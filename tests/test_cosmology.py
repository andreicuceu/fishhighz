"""Deterministic controls for the optional native CAMB boundary."""

import sys
from types import SimpleNamespace

import numpy as np
import pytest

from fishhighz.cosmology import prepare_camb


class _Results:
    def __init__(self, redshifts, *, disagree=False):
        """Initialize the synthetic Results fixture.

        Parameters
        ----------
        redshifts : array_like of shape (n_redshifts,)
            Dimensionless redshifts at which the synthetic background is tabulated.
        disagree : bool, optional
            Whether to perturb one returned redshift to test adapter consistency
            checks. Default is False.

        Notes
        -----
        Sets the instance state used by the enclosing test; no scientific calculation is run.
        """
        self.redshifts = np.asarray(redshifts, dtype=float)
        returned = self.redshifts.copy()
        if disagree:
            returned[-1] += 0.01
        self.Params = SimpleNamespace(
            Transfer=SimpleNamespace(PK_redshifts=returned.tolist())
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


class _FakeCamb:
    def __init__(self, *, disagree=False):
        """Initialize the synthetic FakeCamb fixture.

        Parameters
        ----------
        disagree : bool, optional
            Whether to perturb one returned redshift to test adapter consistency
            checks. Default is False.

        Notes
        -----
        Sets the instance state used by the enclosing test; no scientific calculation is run.
        """
        self.requests = []
        self.disagree = disagree

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
        self.ini_path = path
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
        assert parameters.Transfer.PK_num_redshifts == len(
            parameters.Transfer.PK_redshifts
        )
        self.requests.append(tuple(parameters.Transfer.PK_redshifts))
        return _Results(parameters.Transfer.PK_redshifts, disagree=self.disagree)


def test_exact_redshift_order_and_named_normalizations():
    """Check exact redshift order and named normalizations."""
    fake = _FakeCamb()
    background = prepare_camb(
        redshifts=[3.0, 2.0],
        template_growth_redshift=2.4,
        damping_reference_redshift=2.3,
        camb_module=fake,
    )
    np.testing.assert_array_equal(background.redshifts, [3.0, 2.0, 2.4, 2.3])
    assert fake.requests == [(3.0, 2.4, 2.3, 2.0)]
    np.testing.assert_array_equal(background.redshifts, [3.0, 2.0, 2.4, 2.3])
    np.testing.assert_allclose(
        background.f,
        background.sigma8_values * 0.5 / background.sigma8_values,
    )
    assert background.sigma8_template == background.sigma8_at(2.4)
    assert background.sigma8_damping_reference == background.sigma8_at(2.3)
    assert background.template_growth_redshift != background.damping_reference_redshift
    assert not background.redshifts.flags.writeable
    with pytest.raises(ValueError, match="exact redshift"):
        background.sigma8_at(2.01)


def test_bulk_mapping_rejects_returned_redshift_metadata_mismatch():
    """Check bulk mapping rejects returned redshift metadata mismatch."""
    with pytest.raises(ValueError, match="differ from the requested"):
        prepare_camb(
            redshifts=[2.0, 3.0],
            template_growth_redshift=2.4,
            damping_reference_redshift=2.3,
            camb_module=_FakeCamb(disagree=True),
        )


def test_template_redshift_alias_cannot_override_named_input():
    """Check template redshift alias cannot override named input."""
    with pytest.raises(ValueError, match="either template_growth_redshift"):
        prepare_camb(
            redshifts=[2.3],
            template_growth_redshift=2.4,
            template_redshift=2.4,
            damping_reference_redshift=2.3,
            camb_module=_FakeCamb(),
        )


def test_background_surface_supports_geometry_arrays():
    """Check background surface supports geometry arrays."""
    background = prepare_camb(
        redshifts=[2.3],
        template_growth_redshift=2.3,
        damping_reference_redshift=2.3,
        camb_module=_FakeCamb(),
    )
    np.testing.assert_allclose(background.hubble_parameter([1.0, 2.0]), [69.36, 71.36])
    np.testing.assert_allclose(
        background.transverse_comoving_distance([1.0, 2.0]), [1100.0, 1200.0]
    )


def test_camb_import_is_lazy_and_missing_extra_is_actionable(monkeypatch):
    """Check camb import is lazy and missing extra is actionable.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    monkeypatch.setitem(sys.modules, "camb", None)
    with pytest.raises(ImportError, match=r"fishhighz\[camb\]"):
        prepare_camb(
            redshifts=[2.3],
            template_growth_redshift=2.4,
            damping_reference_redshift=2.3,
        )


@pytest.mark.parametrize("template_z,damping_z", [(0.0, 2.3), (2.4, 0.0), (0.0, 0.0)])
@pytest.mark.parametrize("angular", [False, True])
def test_zero_redshift_normalizations(monkeypatch, template_z, damping_z, angular):
    # A small analytic background has the physical observer limit D_M(0)=0.
    """Check zero redshift normalizations.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    template_z : float
        Template reference redshift, supplied by pytest parametrization.
    damping_z : float
        Damping reference redshift, supplied by pytest parametrization.
    angular : bool
        Whether the tested spectrum has angular dependence, supplied by pytest
        parametrization.
    """
    method = "angular_diameter_distance" if angular else "comoving_radial_distance"
    monkeypatch.setattr(
        _Results,
        method,
        lambda self, z: 1000.0 * z / (1 + z) if angular else 1000.0 * z,
        raising=False,
    )
    fake = _FakeCamb()
    background = prepare_camb(
        redshifts=[3.0, 2.0],
        template_growth_redshift=template_z,
        damping_reference_redshift=damping_z,
        camb_module=fake,
    )
    expected_order = list(dict.fromkeys([3.0, 2.0, template_z, damping_z]))
    np.testing.assert_array_equal(background.redshifts, expected_order)
    assert fake.requests == [tuple(sorted(expected_order, reverse=True))]
    for z in expected_order:
        assert background.sigma8_at(z) == 0.2 + 0.01 * z
        assert background.growth_rate_at(z) == (0.1 + 0.005 * z) / (0.2 + 0.01 * z)
    assert background.transverse_distance_at(0.0) == 0.0
    assert background.sigma8_template == 0.2 + 0.01 * template_z
    assert background.sigma8_damping_reference == 0.2 + 0.01 * damping_z


@pytest.mark.parametrize(
    "bad_z,distance",
    [(0.0, -1.0), (2.0, -1.0), (2.0, 0.0), (0.0, np.nan), (0.0, np.inf)],
)
def test_invalid_distances_remain_rejected(monkeypatch, bad_z, distance):
    """Check invalid distances remain rejected.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    bad_z : float
        Invalid redshift query, supplied by pytest parametrization.
    distance : float
        Distance prescription, supplied by pytest parametrization.
    """
    monkeypatch.setattr(
        _Results,
        "comoving_radial_distance",
        lambda self, z: distance if z == bad_z else 1000.0 * z,
    )
    with pytest.raises(ValueError, match="distance"):
        prepare_camb(
            redshifts=[2.0],
            template_growth_redshift=0.0,
            damping_reference_redshift=2.3,
            camb_module=_FakeCamb(),
        )
