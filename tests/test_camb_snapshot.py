"""Provider contracts with a deterministic fake solver; real CAMB checked separately."""

import copy
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from fishhighz.models.camb_snapshot import CambSnapshotProvider, evaluate_linear_power


def smooth_copy(k, p):
    return p.copy()


class FakeParams:
    def __init__(self):
        self.omk = 0.0
        self.ombh2 = 0.02237
        self.omch2 = 0.12
        self.H0 = 67.36
        self.omnuh2 = 0.0006
        self.DarkEnergy = SimpleNamespace(w=-1.0, wa=0.0)
        self.InitPower = SimpleNamespace(ns=0.9649, As=2.1e-9)
        self.Transfer = SimpleNamespace(PK_redshifts=[2.3])

    def copy(self):
        return copy.deepcopy(self)

    def __str__(self):
        return "fixed fake neutrinos .0006, w=-1, helium fixed"

    def set_accuracy(self, **kw):
        self.accuracy = kw

    def set_matter_power(self, redshifts, kmax):
        self.Transfer = SimpleNamespace(PK_redshifts=redshifts)
        self.kmax = kmax


class FakeResults:
    def __init__(self, p):
        self.Params = p

    def get_matter_power_interpolator(self, **kw):
        assert not kw["nonlinear"] and not kw["hubble_units"] and not kw["k_hunit"]
        return SimpleNamespace(
            kmin=1e-8,
            kmax=1000.0,
            P=lambda z, k: (
                100 * k**-0.5 / (1 + z) ** 2 * self.Params.InitPower.As / 2.1e-9
            ),
        )

    def get_sigma8(self):
        return 0.8 / (1 + np.array(self.Params.Transfer.PK_redshifts))

    def get_fsigma8(self):
        return 0.95 * self.get_sigma8()

    def hubble_parameter(self, z):
        return 100 * np.sqrt(1 + z)

    def angular_diameter_distance(self, z):
        return 1000 * z / (1 + z)

    def get_derived_params(self):
        return {"rdrag": 147.0}


@pytest.fixture
def fake_camb(monkeypatch):
    calls = []

    def solve(p):
        calls.append(p)
        return FakeResults(p)

    fake = SimpleNamespace(
        __version__="test", get_results=solve, model=SimpleNamespace(NonLinear_none=0)
    )
    monkeypatch.setitem(sys.modules, "camb", fake)
    return calls


def test_units_order_cache_and_immutable(tmp_path, fake_camb):
    k = np.geomspace(4e-5, 110.0, 60)
    z = [3.0, 2.0, 2.7]
    h = 0.6736
    provider = CambSnapshotProvider(
        FakeParams(), k, z, h, smooth_copy, cache_dir=tmp_path
    )
    state = provider.evaluate()
    assert provider.solve_count == 1 and len(fake_camb) == 1
    np.testing.assert_allclose(
        state.pk_linear,
        [100 * (k * h) ** -0.5 / (1 + v) ** 2 * h**3 for v in z],
        rtol=1e-15,
    )
    np.testing.assert_array_equal(state.redshifts, z)
    np.testing.assert_allclose(state.sigma8, 0.8 / (1 + np.array(z)))
    np.testing.assert_allclose(state.growth_rate, 0.95)
    assert fake_camb[0].omnuh2 == 0.0006 and fake_camb[0].DarkEnergy.w == -1.0
    assert provider.evaluate() is state and provider.memory_hits == 1
    with pytest.raises(ValueError):
        state.pk_linear[0, 0] = 0.0
    again = CambSnapshotProvider(FakeParams(), k, z, h, smooth_copy, cache_dir=tmp_path)
    loaded = again.evaluate()
    assert again.disk_hits == 1 and again.solve_count == 0 and len(fake_camb) == 1
    np.testing.assert_array_equal(loaded.pk_linear, state.pk_linear)
    changed = provider.evaluate({"h": 0.72})
    assert changed.key != state.key and changed.h_fid == h and provider.solve_count == 2
    np.testing.assert_array_equal(changed.k_reference, k)


def test_no_unvalidated_neutrino_coordinates(fake_camb):
    p = CambSnapshotProvider(
        FakeParams(), np.geomspace(4e-5, 110.0, 40), [2.0], 0.6736, smooth_copy
    )
    with pytest.raises(ValueError, match="unsupported"):
        p.evaluate({"mnu": 0.1})
    assert not fake_camb


def test_low_k_clamping_is_never_silent():
    interp = SimpleNamespace(
        kmin=1e-4, kmax=10.0, P=lambda z, k: 3 * np.maximum(k, 1e-4) ** 0.965
    )
    k = np.array([2e-5, 1e-4, 1e-3])
    with pytest.raises(ValueError, match="below CAMB"):
        evaluate_linear_power(interp, 2.0, k)
    np.testing.assert_allclose(
        evaluate_linear_power(interp, 2.0, k, low_k_policy="matched_power_law"),
        3 * k**0.965,
        rtol=1e-12,
    )
    with pytest.raises(ValueError, match="upper"):
        evaluate_linear_power(interp, 2.0, [11.0])
