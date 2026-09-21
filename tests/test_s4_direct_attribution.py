"""Independent numerical checks for the bounded S4 diagnostic."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from fishhighz.validation.numerics import legacy_peak

spec = importlib.util.spec_from_file_location(
    "attribute_s4", Path(__file__).parents[1] / "scripts" / "attribute_s4.py"
)
s4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s4)


def test_direct_legacy_peak_recovers_native_signed_estimator():
    direct = s4.Direct.__new__(s4.Direct)
    direct.ka = np.linspace(0.01, 0.5, 500)
    direct.pairs = [(0, 0), (0, 1)]
    direct.peak_cache = {}
    direct.linear = lambda k, z, c: 100 * k**-0.8 * (1 + 0.04 * np.sin(105 * k))
    direct.rsd = lambda mu, z: np.column_stack(
        (1 + 0.3 * mu**2, -0.4 * (1 + 0.8 * mu**2))
    )
    direct.response = lambda k, mu, z, c: np.column_stack(
        (np.exp(-2 * (k * mu) ** 2), np.exp(-((k * mu) ** 2)))
    )
    mu = np.full(500, 0.45)
    observed = (
        direct.linear(direct.ka, 2.3, {})[:, None]
        * direct.rsd(mu, 2.3)
        * direct.response(direct.ka, mu, 2.3, {})
    )
    result = direct.peak(direct.ka, mu, 2.3, {"peak": False}) * direct.response(
        direct.ka, mu, 2.3, {}
    )
    np.testing.assert_allclose(
        result, legacy_peak(observed.T, direct.ka).T, rtol=3e-10, atol=1e-11
    )


def test_cross_width_mean_retains_auto_limits():
    direct = s4.Direct.__new__(s4.Direct)
    direct.fields = [SimpleNamespace(kind="forest"), SimpleNamespace(kind="galaxy")]
    direct.pairs = [(0, 0), (0, 1), (1, 1)]
    direct.recipe = SimpleNamespace(
        _growth=lambda z: (0.3, 0.96),
        cosmo=SimpleNamespace(sigma8=0.31),
        config={"survey": {"reconstruction factor": ".5"}},
    )
    legacy = direct.widths(2.3, {"cross_damping": False})
    revised = direct.widths(2.3, {"cross_damping": True})
    np.testing.assert_array_equal(revised[[0, 2]], legacy[[0, 2]])
    np.testing.assert_allclose(revised[1], (revised[0] + revised[2]) / 2)
    np.testing.assert_array_equal(legacy[1], legacy[0])


def test_response_width_conventions_against_explicit_gaussian_sinc():
    direct = s4.Direct.__new__(s4.Direct)
    direct.active = [0]
    direct.fields = [SimpleNamespace(kind="forest", id="forest")]
    direct.pairs = [(0, 0)]
    direct.names = ["forest_forest"]
    direct.zc = 2.4
    direct.recipe = SimpleNamespace(
        cosmo=SimpleNamespace(velocity_from_distance=lambda z: 100.0),
        tracers={"forest": {"pix_width_ang": "1.0"}},
        config={"survey": {"resolution": "3000"}},
    )
    direct.rows = {
        "forest_forest": {
            "_pix_kms": 299800.0 / (1215.67 * 3.4),
            "_res_kms": 299800.0 / 3000,
        }
    }
    k = np.array([0.03, 0.1, 0.3])
    mu = np.array([0.3, 0.6, 0.9])
    q = k * mu / 100.0
    for physical in (False, True):
        c = {
            **s4.FIXED,
            "response": physical,
            "constants": physical,
            "field_response": physical,
        }
        speed = 299792.458 if physical else 299800.0
        sigma = speed / 3000 / (np.sqrt(8 * np.log(2)) if physical else 1)
        argument = q * speed / (1215.67 * 3.4) / 2
        expected = (np.sin(argument) / argument) ** 2 * np.exp(-((q * sigma) ** 2))
        np.testing.assert_allclose(
            direct.response(k, mu, 2.4, c)[:, 0], expected, rtol=5e-15
        )
