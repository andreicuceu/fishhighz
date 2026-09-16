"""Step 12 independent branch, routing, inventory and BAO oracles."""

import json

import numpy as np
import pytest

from fishhighz.adapters.legacy_compat import LegacyDensity, LegacySNR, plain
from fishhighz.adapters.legacy_inputs import DensityReader, SNRReader
from fishhighz.adapters.lyaforecast import IntrinsicP3D
from fishhighz.validation.cases import CASE_IDS, selection
from fishhighz.validation.evidence import check, execute, requests
from fishhighz.validation.synthetic import evidence_payload, run_case


@pytest.fixture
def density(tmp_path):
    p = tmp_path / "density.txt"
    np.savetxt(
        p,
        [[z, m, (z - 3) ** 2 + (m - 21) ** 2] for z in [2, 3, 4] for m in [20, 21, 22]],
    )
    return DensityReader(
        p, semantics="cell_count_per_deg2", target_density=None, z_norm_min=None
    )


@pytest.fixture
def snr(tmp_path):
    paths = []
    for m in [20, 21, 22]:
        p = tmp_path / f"snr{m}.dat"
        np.savetxt(
            p,
            [[w, 2, 2, 2] for w in [4000, 4500, 5000]],
            header=f"BAND=r MAG={m} EXPTIME=1000 NEXP=4\nWave SN(z=2) SN(z=3) SN(z=4)",
        )
        paths.append(p)
    return SNRReader(paths, smoothing="none")


@pytest.mark.parametrize("z", [1.9, 2, 2.1, 3, 4, 4.1])
def test_density_extension(density, z):
    m = np.array([19.9, 20, 20.5, 21, 22, 22.1])
    sample = LegacyDensity(density, "floor_negative").sample(z, m)
    effective = np.clip(z, 2, 4)
    expected = (effective - 3) ** 2 + (m - 21) ** 2
    expected[[0, -1]] = 1e-20
    np.testing.assert_allclose(sample["values"], expected, rtol=5e-12, atol=1e-15)
    assert sample["provenance"]["counts"]["density_floor"] == 2
    assert sample["provenance"]["counts"]["redshift_extension"] == (
        6 if z < 2 or z > 4 else 0
    )
    assert not sample["values"].flags.writeable


def test_negative_and_small_values(tmp_path):
    # Nonnegative nodes of an exact quadratic with negative inter-node values.
    p = tmp_path / "overshoot"
    np.savetxt(
        p, [[z, m, (m - 20.5) ** 2 - 0.2] for z in [2, 3, 4] for m in [20, 21, 22]]
    )
    d = DensityReader(
        p, semantics="cell_count_per_deg2", target_density=None, z_norm_min=None
    )
    with pytest.raises(ValueError, match="negative"):
        LegacyDensity(d, "reject").sample(3, [20.5])
    s = LegacyDensity(d, "floor_negative").sample(3, [20.5])
    assert s["values"][0] == 1e-20
    np.testing.assert_allclose(s["raw"], [-0.2], rtol=5e-13, atol=0)
    for count in (0.0, 1e-25):
        np.savetxt(p, [[z, m, count] for z in [2, 3, 4] for m in [20, 21, 22]])
        d = DensityReader(
            p, semantics="cell_count_per_deg2", target_density=None, z_norm_min=None
        )
        np.testing.assert_allclose(
            LegacyDensity(d, "floor_negative").sample(3, [21])["values"],
            [count],
            rtol=5e-13,
            atol=0,
        )


@pytest.mark.parametrize("m", [19.9, 20, 21, 22, 22.00001])
@pytest.mark.parametrize("z", [1.99999, 2, 3, 4, 4.00001])
@pytest.mark.parametrize("wave", [3999.999, 4000, 4500, 5000, 5000.001])
def test_snr_branches(snr, m, z, wave):
    s = LegacySNR(snr).sample(
        z_source=z,
        magnitudes=[m],
        wavelength=wave,
        pixel_width_angstrom=2,
        exposure_count=8,
    )
    outside = m > 22 or z < 2 or z > 4 or wave < 4000 or wave > 5000
    assert s["values"][0] == pytest.approx(1e20 if outside else 1 / 16, rel=5e-13)
    assert s["provenance"]["counts"]["out_of_range"] == int(outside)
    assert s["provenance"]["counts"]["bright_clamp"] == int(m < 20 and not outside)
    if outside:
        assert s["values"][0] == 1e20


@pytest.mark.parametrize("pixel,count", [(1e-24, 4), (1, 4), (2, 8)])
def test_post_scale_floor(snr, pixel, count):
    s = LegacySNR(snr).sample(
        z_source=3,
        magnitudes=[21],
        wavelength=4500,
        pixel_width_angstrom=pixel,
        exposure_count=count,
    )
    expected = 1 / max(2 * np.sqrt(pixel) * np.sqrt(count / 4), 1e-10) ** 2
    np.testing.assert_allclose(s["values"], [expected], rtol=5e-13, atol=0)


@pytest.mark.parametrize("bad", [True, 1j, "3", np.nan, np.inf])
def test_invalid_queries(density, snr, bad):
    with pytest.raises(ValueError):
        LegacyDensity(density, "reject").sample(bad, [21])
    with pytest.raises(ValueError):
        LegacySNR(snr).sample(
            z_source=3,
            magnitudes=[bad],
            wavelength=4500,
            pixel_width_angstrom=1,
            exposure_count=4,
        )


@pytest.mark.parametrize(
    "key,value",
    [
        ("pixel_width_angstrom", 0),
        ("pixel_width_angstrom", 1j),
        ("exposure_count", -1),
        ("exposure_time", 2),
        ("wavelength", 0),
        ("z_source", -1),
    ],
)
def test_invalid_exposure(snr, key, value):
    kwargs = dict(
        z_source=3,
        magnitudes=[30],
        wavelength=4500,
        pixel_width_angstrom=1,
        exposure_count=4,
    )
    kwargs[key] = value
    with pytest.raises(ValueError):
        LegacySNR(snr).sample(**kwargs)


class External:
    def __init__(self):
        self.calls = 0

    def compute_p3d_hmpc(self, z, k, mu, corr):
        self.calls += 1
        return (-1 if corr == "cross" else 2) * (z + k + mu)


def bridge(obj=None, **kwargs):
    options = dict(
        routes={(0, 0): "auto", (0, 1): "cross", (1, 1): "auto"},
        n_fields=2,
        h_source=0.5,
        h_fid=1.0,
        k_domain=(0.01, 1),
        z_domain=(2, 4),
    )
    options.update(kwargs)
    return IntrinsicP3D(obj or External(), **options)


def test_external_order_units():
    b = bridge()
    k, mu = np.array([0.2, 0.1]), np.array([0.1, 0.7])
    pairs = [(0, 1), (1, 1), (0, 0), (0, 1)]
    expected = (3 + 2 * k + mu)[:, None] * np.array([-1, 2, 2, -1]) * 8
    for _ in range(2):
        np.testing.assert_allclose(b([], 3, k, mu, pairs), expected, rtol=5e-13, atol=0)


@pytest.mark.parametrize(
    "key,value",
    [
        ("k", [0]),
        ("k", [0.6]),
        ("k", [np.nan]),
        ("mu", [1.1]),
        ("mu", [True]),
        ("z", 1.99),
        ("pairs", [(1, 0)]),
        ("pairs", [(2, 2)]),
        ("pairs", [(True, 0)]),
        ("theta_local", [1]),
    ],
)
def test_external_reject_before_calls(key, value):
    obj = External()
    b = bridge(obj)
    args = dict(theta_local=[], z=3, k=[0.1], mu=[0.2], pairs=[(0, 0)])
    args[key] = value
    with pytest.raises(ValueError):
        b(**args)
    assert obj.calls == 0


@pytest.mark.parametrize("output", [1.0, [np.inf], [1j], [True], [[1.0]]])
def test_external_bad_outputs(output):
    class Bad:
        def compute_p3d_hmpc(self, *args):
            return output

    with pytest.raises(ValueError):
        bridge(Bad())([], 3, [0.1], [0.2], [(0, 0)])


@pytest.mark.parametrize("case,count", list(zip(CASE_IDS, [3, 3, 6, 2, 15, 4, 8])))
def test_selections(case, count):
    assert len(selection(case).selected_pairs) == count
    run_case(case)


def test_full_inventory_and_failures(tmp_path):
    seen = []

    def worker(task):
        seen.append((task["case"], task["bin"]))
        return evidence_payload(task)

    out = tmp_path / "full"
    m = execute(out, suite="full", worker=worker, kind="synthetic_bao")
    assert len(seen) == 39 and len(set(seen)) == 39
    assert check(out)["complete"]
    with pytest.raises(FileExistsError):
        execute(out, suite="full", worker=worker, kind="synthetic_bao")

    def failure(task):
        if task["bin"] == 2:
            raise ValueError("deliberate")
        return worker(task)

    m = execute(tmp_path / "failed", suite="full", worker=failure, kind="synthetic_bao")
    assert len(m["records"]) == 39 and not m["complete"]
    with pytest.raises(ValueError, match="partial"):
        check(tmp_path / "failed")


@pytest.mark.parametrize(
    "change",
    ["case", "bin", "pair", "missing", "nonfinite", "hash", "source", "partial"],
)
def test_checker_corruption(tmp_path, change):
    source = tmp_path / "input"
    source.write_text("original")
    out = tmp_path / "evidence"
    execute(
        out,
        suite="quick",
        inputs=[source],
        worker=evidence_payload,
        kind="synthetic_bao",
    )
    m = json.loads((out / "manifest.json").read_text())
    if change == "case":
        m["records"][0]["task"]["case"] = CASE_IDS[0]
    if change == "bin":
        m["records"][0]["task"]["bin"] = 0
    if change == "pair":
        m["records"][0]["task"]["selected_pairs"].pop()
    if change == "missing":
        m["records"] = []
    if change == "nonfinite":
        np.savez(out / "records-000.npz", fisher=[[np.nan]], errors=[1])
    if change == "hash":
        m["records"][0]["sha256"] = "stale"
    if change == "source":
        source.write_text("changed")
    if change == "partial":
        m["records"][0]["status"] = "failed"
    (out / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError):
        check(out)


def test_plain_diagnostics(density):
    json.dumps(plain(LegacyDensity(density, "reject").sample(2, [20])), allow_nan=False)
    with pytest.raises(ValueError):
        LegacyDensity(density, "implicit")
    assert requests("quick")[0]["case"] == "lya_qso_lbg_lae_15x2pt"


@pytest.mark.parametrize("anisotropic", [False, True])
def test_bao_oracle_and_independent_bins(anisotropic):
    from fishhighz.covariance import gaussian_covariance
    from fishhighz.derivatives import evaluate_derivatives
    from fishhighz.fields import ObservedField, PairSelection
    from fishhighz.fisher import fisher_matrix
    from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
    from fishhighz.models.kaiser import KaiserModel, Scaling
    from fishhighz.models.templates import prepare_template
    from fishhighz.parameters import Parameter, ParameterRegistry
    from fishhighz.results import FisherResult, combine_results

    knots = np.geomspace(0.001, 1, 16000)
    template = prepare_template(
        knots,
        20 + np.sin(30 * knots),
        np.full_like(knots, 20),
        z_ref=2.4,
        h_template=0.7,
        h_fid=0.7,
    )
    field = ObservedField("g", "galaxy", "toy")
    select = PairSelection([field], [("g", "g")])
    registry = ParameterRegistry(
        [
            Parameter(f"{p}{i}", 1, "target", step=1e-5)
            for i in (0, 1)
            for p in ("ap", "at")
        ]
    )
    k = np.array([0.04, 0.07, 0.11, 0.17, 0.23, 0.31])
    mu = np.array([0.1, 0.8, 0.4, 0.6, 0.95, 0.2])
    results = []
    for i in (0, 1):
        f = 0.8 if anisotropic else 0
        sp, st = (5.0, 2.0) if anisotropic else (0.0, 0.0)
        model = KaiserModel(
            template,
            [field],
            biases={"g": 1.0},
            betas={},
            widths={"g": (sp, st)},
            f=f,
            local_names=("ap", "at"),
            wiggle=Scaling("ap_at", ap="ap", at="at"),
        )
        collection = PreparedP3D(
            registry,
            select,
            [
                P3DProvider(
                    "BAO",
                    model,
                    BoundParameters(
                        registry, ("ap", "at"), {"ap": f"ap{i}", "at": f"at{i}"}
                    ),
                    select.required_pairs,
                )
            ],
        )
        d = evaluate_derivatives(collection, registry.fiducials, 2.4, k, mu)
        w = np.sin(30 * k)
        wp = 30 * np.cos(30 * k)
        oracle = np.column_stack((-w - k * mu**2 * wp, -2 * w - k * (1 - mu**2) * wp))
        if anisotropic:
            # Independent closed expression, fourth-order differences; no model calls.
            def power(ap, at):
                kp = k * mu / ap
                kt = k * np.sqrt(1 - mu**2) / at
                kk = np.hypot(kp, kt)
                mm = kp / kk
                return 20 * (1 + f * mu**2) ** 2 + np.sin(30 * kk) * (
                    1 + f * mm**2
                ) ** 2 * np.exp(-0.5 * ((kp * sp) ** 2 + (kt * st) ** 2)) / (ap * at**2)

            h = 1e-4
            oracle = np.column_stack(
                [
                    (
                        -power(1 + 2 * h, 1)
                        + 8 * power(1 + h, 1)
                        - 8 * power(1 - h, 1)
                        + power(1 - 2 * h, 1)
                    )
                    / (12 * h),
                    (
                        -power(1, 1 + 2 * h)
                        + 8 * power(1, 1 + h)
                        - 8 * power(1, 1 - h)
                        + power(1, 1 - 2 * h)
                    )
                    / (12 * h),
                ]
            )
        np.testing.assert_allclose(
            d.jacobian[:, 0, 2 * i : 2 * i + 2], oracle, rtol=3e-7, atol=1e-8
        )
        np.testing.assert_array_equal(
            d.jacobian[:, 0, 2 * (1 - i) : 2 * (1 - i) + 2], 0
        )
        modes = np.arange(1, 7) * 100.0
        covariance = gaussian_covariance(d.power + 2, modes, select)
        expected_variance = 2 * (d.power[:, 0] + 2) ** 2 / modes
        np.testing.assert_allclose(
            covariance[:, 0, 0], expected_variance, rtol=5e-13, atol=0
        )
        fisher = fisher_matrix(d.jacobian, covariance)
        expected = oracle.T @ (oracle / expected_variance[:, None])
        np.testing.assert_allclose(
            fisher[2 * i : 2 * i + 2, 2 * i : 2 * i + 2], expected, rtol=3e-7, atol=0
        )
        result = FisherResult(registry, fisher)
        assert result.diagnostics.rank == 2
        results.append(result)
    combined = combine_results(results)
    assert combined.diagnostics.rank == 4
    np.testing.assert_array_equal(combined.data_fisher[:2, 2:], 0)


def test_exact_pair_order():
    expected = {
        "lbg_lae_3x2pt": [(0, 0), (0, 1), (1, 1)],
        "lya_lbg_lae_3x2pt": [(0, 0), (0, 1), (0, 2)],
        "lya_lbg_lae_6x2pt": [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)],
        "lya_qso_2x2pt": [(0, 0), (0, 1)],
        "lya_qso_lbg_lae_4x2pt": [(0, 0), (0, 1), (0, 2), (0, 3)],
        "lya_qso_lbg_lae_8x2pt": [
            (0, 0),
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 4),
            (2, 4),
            (3, 4),
            (4, 4),
        ],
        "lya_qso_lbg_lae_15x2pt": [
            (0, 0),
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 1),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 2),
            (2, 3),
            (2, 4),
            (3, 3),
            (3, 4),
            (4, 4),
        ],
    }
    assert set(expected) == set(CASE_IDS)
    for case, pairs in expected.items():
        np.testing.assert_array_equal(selection(case).selected_pairs, pairs)


def test_optional_imports_in_subprocess():
    import subprocess
    import sys

    code = """
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, *args):
        if name.split('.')[0] in ('scipy','astropy','camb','lyaforecast','vega'):
            raise ImportError('blocked optional package')
sys.meta_path.insert(0,Block())
from fishhighz.validation.synthetic import run
assert len(run())==7
from fishhighz.adapters.legacy_inputs import DensityReader
try:
    DensityReader('unused',semantics='cell_count_per_deg2',target_density=None,z_norm_min=None)
except ImportError as e:
    assert 'fishhighz[survey]' in str(e)
else:
    raise AssertionError('missing extra accepted')
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_compatibility_owns_spline_snapshots(density, snr):
    d = LegacyDensity(density, "reject")
    s = LegacySNR(snr)
    before = d.sample(2.5, [20.5])["values"].copy()
    # Mutate caller-owned interpolators; the compatibility snapshots stay fixed.
    object.__setattr__(density, "_spline", None)
    snr._interpolator.values[:] = 99
    np.testing.assert_array_equal(d.sample(2.5, [20.5])["values"], before)
    actual = s.sample(
        z_source=3,
        magnitudes=[21],
        wavelength=4500,
        pixel_width_angstrom=1,
        exposure_count=4,
    )
    np.testing.assert_allclose(actual["values"], [0.25], rtol=5e-13, atol=0)
