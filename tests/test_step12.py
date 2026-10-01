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
    """Write and read a quadratic synthetic source-density table.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest for generated inputs and results.

    Returns
    -------
    reader : DensityReader
        Temporary cell-count table with unnormalized source densities.
    """
    density_path = tmp_path / "density.txt"
    np.savetxt(
        density_path,
        [[z, m, (z - 3) ** 2 + (m - 21) ** 2] for z in [2, 3, 4] for m in [20, 21, 22]],
    )
    return DensityReader(
        density_path,
        semantics="cell_count_per_deg2",
        target_density=None,
        z_norm_min=None,
    )


@pytest.fixture
def snr(tmp_path):
    """Write and read three constant synthetic pixel-SNR tables.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory supplied by pytest for generated inputs and results.

    Returns
    -------
    reader : SNRReader
        Reader returning dimensionless SNR without smoothing.
    """
    paths = []
    for m in [20, 21, 22]:
        snr_path = tmp_path / f"snr{m}.dat"
        np.savetxt(
            snr_path,
            [[w, 2, 2, 2] for w in [4000, 4500, 5000]],
            header=f"BAND=r MAG={m} EXPTIME=1000 NEXP=4\nWave SN(z=2) SN(z=3) SN(z=4)",
        )
        paths.append(snr_path)
    return SNRReader(paths, smoothing="none")


@pytest.mark.parametrize("z", [1.9, 2, 2.1, 3, 4, 4.1])
def test_density_extension(density, z):
    """Check density extension.

    Parameters
    ----------
    density : DensityReader
        Synthetic source-density reader supplied by the density fixture.
    z : int or float
        Dimensionless redshift test input, supplied by pytest parametrization.
    """
    magnitude_grid = np.array([19.9, 20, 20.5, 21, 22, 22.1])
    sample = LegacyDensity(density, "floor_negative").sample(z, magnitude_grid)
    effective = np.clip(z, 2, 4)
    expected = (effective - 3) ** 2 + (magnitude_grid - 21) ** 2
    expected[[0, -1]] = 1e-20
    np.testing.assert_allclose(sample["values"], expected, rtol=5e-12, atol=1e-15)
    assert sample["provenance"]["counts"]["density_floor"] == 2
    assert sample["provenance"]["counts"]["redshift_extension"] == (
        6 if z < 2 or z > 4 else 0
    )
    assert not sample["values"].flags.writeable


def test_negative_and_small_values(tmp_path):
    # Nonnegative nodes of an exact quadratic with negative inter-node values.
    """Check negative and small values.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    density_path = tmp_path / "overshoot"
    np.savetxt(
        density_path,
        [[z, m, (m - 20.5) ** 2 - 0.2] for z in [2, 3, 4] for m in [20, 21, 22]],
    )
    density_reader = DensityReader(
        density_path,
        semantics="cell_count_per_deg2",
        target_density=None,
        z_norm_min=None,
    )
    with pytest.raises(ValueError, match="negative"):
        LegacyDensity(density_reader, "reject").sample(3, [20.5])
    density_sample = LegacyDensity(density_reader, "floor_negative").sample(3, [20.5])
    assert density_sample["values"][0] == 1e-20
    np.testing.assert_allclose(density_sample["raw"], [-0.2], rtol=5e-13, atol=0)
    for count in (0.0, 1e-25):
        np.savetxt(
            density_path, [[z, m, count] for z in [2, 3, 4] for m in [20, 21, 22]]
        )
        density_reader = DensityReader(
            density_path,
            semantics="cell_count_per_deg2",
            target_density=None,
            z_norm_min=None,
        )
        np.testing.assert_allclose(
            LegacyDensity(density_reader, "floor_negative").sample(3, [21])["values"],
            [count],
            rtol=5e-13,
            atol=0,
        )


@pytest.mark.parametrize("m", [19.9, 20, 21, 22, 22.00001])
@pytest.mark.parametrize("z", [1.99999, 2, 3, 4, 4.00001])
@pytest.mark.parametrize("wave", [3999.999, 4000, 4500, 5000, 5000.001])
def test_snr_branches(snr, m, z, wave):
    """Check snr branches.

    Parameters
    ----------
    snr : SNRReader
        Synthetic SNR reader supplied by the snr fixture.
    m : int or float
        Parametrized magnitude or size input, supplied by pytest
        parametrization.
    z : int or float
        Dimensionless redshift test input, supplied by pytest parametrization.
    wave : int or float
        Observed wavelength input, supplied by pytest parametrization.
    """
    snr_sample = LegacySNR(snr).sample(
        z_source=z,
        magnitudes=[m],
        wavelength=wave,
        pixel_width_angstrom=2,
        exposure_count=8,
    )
    outside = m > 22 or z < 2 or z > 4 or wave < 4000 or wave > 5000
    assert snr_sample["values"][0] == pytest.approx(
        1e20 if outside else 1 / 16, rel=5e-13
    )
    assert snr_sample["provenance"]["counts"]["out_of_range"] == int(outside)
    assert snr_sample["provenance"]["counts"]["bright_clamp"] == int(
        m < 20 and not outside
    )
    if outside:
        assert snr_sample["values"][0] == 1e20


@pytest.mark.parametrize("pixel,count", [(1e-24, 4), (1, 4), (2, 8)])
def test_post_scale_floor(snr, pixel, count):
    """Check post scale floor.

    Parameters
    ----------
    snr : SNRReader
        Synthetic SNR reader supplied by the snr fixture.
    pixel : int or float
        Pixel width test input, supplied by pytest parametrization.
    count : int
        Number of iterations, samples, or records selected by this case,
        supplied by pytest parametrization.
    """
    snr_sample = LegacySNR(snr).sample(
        z_source=3,
        magnitudes=[21],
        wavelength=4500,
        pixel_width_angstrom=pixel,
        exposure_count=count,
    )
    expected = 1 / max(2 * np.sqrt(pixel) * np.sqrt(count / 4), 1e-10) ** 2
    np.testing.assert_allclose(snr_sample["values"], [expected], rtol=5e-13, atol=0)


@pytest.mark.parametrize("bad", [True, 1j, "3", np.nan, np.inf])
def test_invalid_queries(density, snr, bad):
    """Check invalid queries.

    Parameters
    ----------
    density : DensityReader
        Synthetic source-density reader supplied by the density fixture.
    snr : SNRReader
        Synthetic SNR reader supplied by the snr fixture.
    bad : bool or float or str or complex
        Invalid input exercising the specified rejection path, supplied by
        pytest parametrization.
    """
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
    """Check invalid exposure.

    Parameters
    ----------
    snr : SNRReader
        Synthetic SNR reader supplied by the snr fixture.
    key : str
        Dictionary or configuration key under examination, supplied by pytest
        parametrization.
    value : int or complex
        Value at the tested validation boundary, supplied by pytest
        parametrization.
    """
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
        """Initialize the synthetic External fixture.

        Notes
        -----
        Sets the instance state used by the enclosing test; no scientific calculation is run.
        """
        self.calls = 0

    def compute_p3d_hmpc(self, z, k, mu, corr):
        """Evaluate signed synthetic power and increment the call counter.

        Parameters
        ----------
        z : float
            Dimensionless evaluation redshift.
        k : ndarray of shape (n_nodes,)
            Wavenumbers in the external provider convention, h_source/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines.
        corr : str
            Correlation route: cross gives negative power and auto positive power.

        Returns
        -------
        power : ndarray of shape (n_nodes,)
            Synthetic intrinsic power in (Mpc/h_source)^3.
        """
        self.calls += 1
        return (-1 if corr == "cross" else 2) * (z + k + mu)


def bridge(obj=None, **kwargs):
    """Wrap a synthetic external P3D provider with explicit unit conversion.

    Parameters
    ----------
    obj : object or None, optional
        External object exposing compute_p3d_hmpc; None constructs the synthetic
        provider. Default is None.
    **kwargs : dict
        IntrinsicP3D overrides for pair routes, field count, Hubble conventions, or redshift/wavenumber domains.

    Returns
    -------
    provider : IntrinsicP3D
        Adapter with pair routes and validated redshift/wavenumber domains.
    """
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
    """Check external order units."""
    external_provider = bridge()
    k_grid, mu_grid = np.array([0.2, 0.1]), np.array([0.1, 0.7])
    pairs = [(0, 1), (1, 1), (0, 0), (0, 1)]
    expected = (3 + 2 * k_grid + mu_grid)[:, None] * np.array([-1, 2, 2, -1]) * 8
    for _ in range(2):
        np.testing.assert_allclose(
            external_provider([], 3, k_grid, mu_grid, pairs),
            expected,
            rtol=5e-13,
            atol=0,
        )


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
    """Check external reject before calls.

    Parameters
    ----------
    key : str
        Dictionary or configuration key under examination, supplied by pytest
        parametrization.
    value : float or list
        Value at the tested validation boundary, supplied by pytest
        parametrization.
    """
    obj = External()
    external_provider = bridge(obj)
    args = dict(theta_local=[], z=3, k=[0.1], mu=[0.2], pairs=[(0, 0)])
    args[key] = value
    with pytest.raises(ValueError):
        external_provider(**args)
    assert obj.calls == 0


@pytest.mark.parametrize("output", [1.0, [np.inf], [1j], [True], [[1.0]]])
def test_external_bad_outputs(output):
    """Check external bad outputs.

    Parameters
    ----------
    output : float or list
        Output quantity under examination, supplied by pytest parametrization.
    """

    class Bad:
        def compute_p3d_hmpc(self, *args):
            """Return the deliberately invalid external-provider output.

            Parameters
            ----------
            *args : tuple
                Evaluation redshift, wavenumbers, direction cosines, and correlation
                route; ignored by this callback.

            Returns
            -------
            output : object
                Parametrized invalid power output retained without conversion.
            """
            return output

    with pytest.raises(ValueError):
        bridge(Bad())([], 3, [0.1], [0.2], [(0, 0)])


@pytest.mark.parametrize("case,count", list(zip(CASE_IDS, [3, 3, 6, 2, 15, 4, 8])))
def test_selections(case, count):
    """Check selections.

    Parameters
    ----------
    case : str
        Named forecast or validation case, supplied by pytest parametrization.
    count : int
        Expected number of selected auto- and cross-spectra for the named case.
    """
    assert len(selection(case).selected_pairs) == count
    run_case(case)


def test_full_inventory_and_failures(tmp_path):
    """Check full inventory and failures.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    """
    seen = []

    def worker(task):
        """Construct synthetic evidence for the requested forecast task.

        Parameters
        ----------
        task : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        seen.append((task["case"], task["bin"]))
        return evidence_payload(task)

    out = tmp_path / "full"
    manifest = execute(out, suite="full", worker=worker, kind="synthetic_bao")
    assert len(seen) == 39 and len(set(seen)) == 39
    assert check(out)["complete"]
    with pytest.raises(FileExistsError):
        execute(out, suite="full", worker=worker, kind="synthetic_bao")

    def failure(task):
        """Inject one redshift-bin failure into the synthetic worker.

        Parameters
        ----------
        task : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.

        Returns
        -------
        payload : tuple
            Numerical arrays and evidence report for non-failing bins.

        Raises
        ------
        ValueError
            Deliberately raised to exercise the rejection path in the enclosing
            test.
        """
        if task["bin"] == 2:
            raise ValueError("deliberate")
        return worker(task)

    manifest = execute(
        tmp_path / "failed", suite="full", worker=failure, kind="synthetic_bao"
    )
    assert len(manifest["records"]) == 39 and not manifest["complete"]
    with pytest.raises(ValueError, match="partial"):
        check(tmp_path / "failed")


@pytest.mark.parametrize(
    "change",
    ["case", "bin", "pair", "missing", "nonfinite", "hash", "source", "partial"],
)
def test_checker_corruption(tmp_path, change):
    """Check checker corruption.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    change : str
        Input override exercising the specified validation boundary, supplied by
        pytest parametrization.
    """
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
    manifest = json.loads((out / "manifest.json").read_text())
    if change == "case":
        manifest["records"][0]["task"]["case"] = CASE_IDS[0]
    if change == "bin":
        manifest["records"][0]["task"]["bin"] = 0
    if change == "pair":
        manifest["records"][0]["task"]["selected_pairs"].pop()
    if change == "missing":
        manifest["records"] = []
    if change == "nonfinite":
        np.savez(out / "records-000.npz", fisher=[[np.nan]], errors=[1])
    if change == "hash":
        manifest["records"][0]["sha256"] = "stale"
    if change == "source":
        source.write_text("changed")
    if change == "partial":
        manifest["records"][0]["status"] = "failed"
    (out / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        check(out)


def test_plain_diagnostics(density):
    """Check plain diagnostics.

    Parameters
    ----------
    density : DensityReader
        Synthetic source-density reader supplied by the density fixture.
    """
    json.dumps(plain(LegacyDensity(density, "reject").sample(2, [20])), allow_nan=False)
    with pytest.raises(ValueError):
        LegacyDensity(density, "implicit")
    assert requests("quick")[0]["case"] == "lya_qso_lbg_lae_15x2pt"


@pytest.mark.parametrize("anisotropic", [False, True])
def test_bao_oracle_and_independent_bins(anisotropic):
    """Check bao oracle and independent bins.

    Parameters
    ----------
    anisotropic : bool
        Whether parallel and transverse dilations differ, supplied by pytest
        parametrization.
    """
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
    k_grid = np.array([0.04, 0.07, 0.11, 0.17, 0.23, 0.31])
    mu_grid = np.array([0.1, 0.8, 0.4, 0.6, 0.95, 0.2])
    results = []
    for i in (0, 1):
        growth_rate = 0.8 if anisotropic else 0
        sp, st = (5.0, 2.0) if anisotropic else (0.0, 0.0)
        model = KaiserModel(
            template,
            [field],
            biases={"g": 1.0},
            betas={},
            widths={"g": (sp, st)},
            f=growth_rate,
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
        derivatives = evaluate_derivatives(
            collection, registry.fiducials, 2.4, k_grid, mu_grid
        )
        wiggle_power = np.sin(30 * k_grid)
        wiggle_derivative = 30 * np.cos(30 * k_grid)
        oracle = np.column_stack(
            (
                -wiggle_power - k_grid * mu_grid**2 * wiggle_derivative,
                -2 * wiggle_power - k_grid * (1 - mu_grid**2) * wiggle_derivative,
            )
        )
        if anisotropic:
            # Independent closed expression, fourth-order differences; no model calls.
            def power(ap, at):
                """Evaluate a scalar dilated wiggle spectrum for the BAO derivative check.

                Parameters
                ----------
                ap : float
                    Dimensionless parallel BAO dilation.
                at : float
                    Dimensionless transverse BAO dilation.

                Returns
                -------
                power : ndarray
                    Synthetic power in (Mpc/h)^3 on the enclosing test Fourier nodes.
                """
                k_parallel = k_grid * mu_grid / ap
                k_transverse = k_grid * np.sqrt(1 - mu_grid**2) / at
                mapped_k = np.hypot(k_parallel, k_transverse)
                mapped_mu = k_parallel / mapped_k
                return 20 * (1 + growth_rate * mu_grid**2) ** 2 + np.sin(
                    30 * mapped_k
                ) * (1 + growth_rate * mapped_mu**2) ** 2 * np.exp(
                    -0.5 * ((k_parallel * sp) ** 2 + (k_transverse * st) ** 2)
                ) / (ap * at**2)

            step = 1e-4
            oracle = np.column_stack(
                [
                    (
                        -power(1 + 2 * step, 1)
                        + 8 * power(1 + step, 1)
                        - 8 * power(1 - step, 1)
                        + power(1 - 2 * step, 1)
                    )
                    / (12 * step),
                    (
                        -power(1, 1 + 2 * step)
                        + 8 * power(1, 1 + step)
                        - 8 * power(1, 1 - step)
                        + power(1, 1 - 2 * step)
                    )
                    / (12 * step),
                ]
            )
        np.testing.assert_allclose(
            derivatives.jacobian[:, 0, 2 * i : 2 * i + 2], oracle, rtol=3e-7, atol=1e-8
        )
        np.testing.assert_array_equal(
            derivatives.jacobian[:, 0, 2 * (1 - i) : 2 * (1 - i) + 2], 0
        )
        modes = np.arange(1, 7) * 100.0
        covariance = gaussian_covariance(derivatives.power + 2, modes, select)
        expected_variance = 2 * (derivatives.power[:, 0] + 2) ** 2 / modes
        np.testing.assert_allclose(
            covariance[:, 0, 0], expected_variance, rtol=5e-13, atol=0
        )
        fisher = fisher_matrix(derivatives.jacobian, covariance)
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
    """Check exact pair order."""
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
    """Check optional imports in subprocess."""
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
    """Check compatibility owns spline snapshots.

    Parameters
    ----------
    density : DensityReader
        Synthetic source-density reader supplied by the density fixture.
    snr : SNRReader
        Synthetic SNR reader supplied by the snr fixture.
    """
    density_adapter = LegacyDensity(density, "reject")
    snr_adapter = LegacySNR(snr)
    before = density_adapter.sample(2.5, [20.5])["values"].copy()
    # Mutate caller-owned interpolators; the compatibility snapshots stay fixed.
    object.__setattr__(density, "_spline", None)
    snr._interpolator.values[:] = 99
    np.testing.assert_array_equal(density_adapter.sample(2.5, [20.5])["values"], before)
    actual = snr_adapter.sample(
        z_source=3,
        magnitudes=[21],
        wavelength=4500,
        pixel_width_angstrom=1,
        exposure_count=4,
    )
    np.testing.assert_allclose(actual["values"], [0.25], rtol=5e-13, atol=0)
