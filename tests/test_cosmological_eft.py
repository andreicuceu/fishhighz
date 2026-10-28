"""Direct AP, null features and updated/fresh EFT state contracts."""

import json

import numpy as np
import pytest
from eftkernels import EFTTracerParameters
from eftkernels.primordial_features import FeatureTemplate

from fishhighz.models.camb_snapshot import CosmologySnapshot
from fishhighz.models.cosmological_eft import CosmologicalEFTAdapter
from fishhighz.models.physical_ap import physical_ap_coordinates


class Provider:
    defaults = {
        "omega_b": 0.02237,
        "omega_c": 0.12,
        "h": 0.6736,
        "n_s": 0.9649,
        "ln_1e10_A_s": np.log(21.0),
    }

    def __init__(self):
        self.cache = {}
        self.solve_count = 0

    def evaluate(self, values=None):
        p = {**self.defaults, **(values or {})}
        key = json.dumps(p, sort_keys=True)
        if key not in self.cache:
            k = np.geomspace(4e-5, 100.1, 300)
            power = (
                2e4
                * k
                / (1 + (k / 0.02) ** 1.5) ** 2
                * np.exp(p["ln_1e10_A_s"] - self.defaults["ln_1e10_A_s"])
            )
            self.cache[key] = CosmologySnapshot(
                key,
                tuple(p.items()),
                np.array([2.5]),
                k,
                power[None],
                power[None],
                np.array([0.97 + 0.1 * (p["omega_c"] - 0.12)]),
                np.array([0.3]),
                np.array([300 * p["h"] / 0.6736]),
                np.array([5000 * (1 + p["h"] - 0.6736)]),
                147.0,
                0.6736,
                "{}",
            )
            self.solve_count += 1
        return self.cache[key]


def objects(p, f):
    return {
        "F": EFTTracerParameters(b1=p["bF"], bEta=p["beF"]),
        "G": EFTTracerParameters(b1=p["bG"], bEta=-1.0, bEta2=1.0, bdeltaeta=-p["bG"]),
    }


def make(provider, **kw):
    defaults = {**provider.defaults, "bF": -0.14, "beF": 0.2, "bG": 3.0}
    if kw.get("feature") is not None:
        defaults.update(a_sin=0.0, a_cos=0.0)
    return CosmologicalEFTAdapter(
        provider,
        2.5,
        ["F", "G", "F"],
        list(defaults),
        defaults,
        objects,
        nmax=128,
        num_mu=8,
        k_output=np.geomspace(0.0001, 0.35, 180),
        **kw,
    )


def test_tree_null_feature_ap_and_shared_forest():
    provider = Provider()
    template = FeatureTemplate("log", 15.0)
    a = make(provider, order="tree", no_ir_resummation=True, feature=template)
    theta = np.array([a.defaults[n] for n in a.names])
    k = np.array([0.02, 0.05, 0.15])
    mu = np.array([0.0, 0.5, 1.0])
    pairs = np.array([[0, 0], [0, 1], [1, 1], [0, 2]])
    null = a(theta, 2.5, k, mu, pairs)
    baseline = make(provider, order="tree", no_ir_resummation=True)
    np.testing.assert_array_equal(
        null, baseline(theta[: len(baseline.names)], 2.5, k, mu, pairs)
    )
    np.testing.assert_array_equal(null[:, 0], null[:, 3])
    shifted = theta.copy()
    shifted[a.names.index("h")] *= 1.01
    p = {**a.defaults, **dict(zip(a.names, shifted))}
    state = provider.evaluate({n: p[n] for n in provider.defaults})
    km, _mm, _j = physical_ap_coordinates(
        k, mu, state.DM[0] / a.fiducial.DM[0], a.fiducial.H[0] / state.H[0]
    )
    output = a(shifted, 2.5, k, mu, pairs)
    for col in (0, 1):
        plus = shifted.copy()
        minus = shifted.copy()
        name = ("a_sin", "a_cos")[col]
        plus[a.names.index(name)] = 0.003
        minus[a.names.index(name)] = -0.003
        numeric = (a(plus, 2.5, k, mu, pairs) - a(minus, 2.5, k, mu, pairs)) / 0.006
        np.testing.assert_allclose(
            numeric,
            output * template.quadratures(km, 0.6736)[col, :, None],
            rtol=2e-12,
            atol=1e-12,
        )


def test_updated_and_fresh_loop_match_and_nuisance_reuses_snapshot():
    p = Provider()
    a = make(p, order="loop", no_ir_resummation=True)
    theta = np.array([a.defaults[n] for n in a.names])
    k = np.array([0.03, 0.08, 0.15])
    mu = np.array([0.2, 0.5, 0.9])
    pairs = np.array([[0, 0], [0, 1], [1, 1]])
    a(theta, 2.5, k, mu, pairs)
    nuisance = theta.copy()
    nuisance[a.names.index("bG")] += 0.01
    a(nuisance, 2.5, k, mu, pairs)
    assert p.solve_count == 1 and a.update_count == 1
    shifted = theta.copy()
    shifted[a.names.index("h")] *= 1.01
    updated = a(shifted, 2.5, k, mu, pairs)
    fresh = make(p, order="loop", no_ir_resummation=True)
    reference = fresh(shifted, 2.5, k, mu, pairs)
    np.testing.assert_allclose(updated, reference, rtol=2e-13, atol=1e-12)
    assert a.update_count == 2 and p.solve_count == 2


def test_scientific_production_guards():
    p = Provider()
    with pytest.raises(ValueError, match="not validated"):
        make(p, order="loop", feature=FeatureTemplate("linear", 200.0))
    with pytest.raises(ValueError, match="not validated"):
        make(p, order="tree", feature=FeatureTemplate("log", 15.0))
    with pytest.raises(ValueError, match="independent AP"):
        CosmologicalEFTAdapter(p, 2.5, ["F"], ["f"], {"f": 0.97}, objects)
