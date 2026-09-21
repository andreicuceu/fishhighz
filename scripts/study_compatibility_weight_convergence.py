"""Stage 2: bounded signed trajectories and original-reader grid refinement.

Density/SNR arithmetic adapted from lyaforecast tracer.py/spectrograph.py (GPLv3).
No CAMB, power-spectrum model, or forecast is evaluated.
"""

import argparse
import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from scipy.ndimage import gaussian_filter1d

from fishhighz.adapters.legacy_inputs import SNRReader
from fishhighz.validation.compatibility_weights import VARIANTS, WeightInputs
from fishhighz.validation.evidence import digest
from fishhighz.validation.weight_convergence import (
    classify,
    relative_change,
    trajectory,
)


class RawSampling:
    """Literal compatibility density/SNR preparation for the two backgrounds."""

    def __init__(self, reference):
        self.hashes = {}
        self.density, self.snr = {}, {}
        for pop, filename, target in (
            ("qso", "dn_dzdr_qso_desi_2.dat", 90.0),
            ("lbg", "lbg_matched_dndzdr.txt", 380.0),
        ):
            path = reference / "lyaforecast/resources/data" / filename
            self.hashes[str(path.resolve())] = digest(path)
            z, m, values = np.loadtxt(path, unpack=True)
            values *= (m >= 16.1) & (m <= 26.75)
            values *= target / np.sum(values * (z > 2.15))
            z, m = np.unique(z), np.unique(m)
            values /= (z[1] - z[0]) * (m[1] - m[0])
            self.density[pop] = (
                RectBivariateSpline(z, m, values.reshape(len(z), len(m)), kx=2, ky=2),
                m,
            )
            paths = sorted(
                (
                    reference / "lyaforecast/resources/data" / ("DESI-2-" + pop.upper())
                ).glob("*.dat")
            )
            for path in paths:
                self.hashes[str(path.resolve())] = digest(path)
            self.snr[pop] = SNRReader(paths, smoothing="none")
        a, b = self.snr["qso"], self.snr["lbg"]
        assert np.array_equal(a.magnitudes, b.magnitudes) and np.array_equal(a.z, b.z)
        n = min(len(a.wavelength), len(b.wavelength))
        # Legacy combines equal row indices, retaining primary QSO wavelengths.
        # The two raw wavelength axes are not identical.
        tensors = {
            "qso": (a, a.raw_snr, a.wavelength),
            "lbg": (b, b.raw_snr, b.wavelength),
            "mixed": (
                a,
                np.sqrt(a.raw_snr[:, :, :n] * b.raw_snr[:, :, :n]),
                a.wavelength[:n],
            ),
        }
        self.interpolators = {}
        for key, (reader, raw, wave) in tensors.items():
            self.interpolators[key] = (
                reader.magnitudes,
                reader.z,
                wave,
                RegularGridInterpolator(
                    (reader.magnitudes, reader.z, wave),
                    gaussian_filter1d(raw, 10, axis=2),
                    bounds_error=False,
                    fill_value=None,
                ),
            )

    def sample(self, name, row, m):
        pop = "lbg" if "lya(lbg)" in name else "qso"
        spline, axis = self.density[pop]
        density = spline(row["_zq"], m, grid=False)
        density[(m < axis[0]) | (m > axis[-1])] = 1e-20
        density /= 2.998e5 / (1 + row["_zq"])
        key = "mixed" if name == "lya(qso)_lya(lbg)" else pop
        mags, z, wave, interp = self.interpolators[key]
        wavelength = 1215.67 * (1 + row["_z_mean"])
        pixel = row["_pix_kms"] / (2.998e5 / 1215.67 / (1 + row["_z_mean"]))
        outside = (
            (m > mags[-1])
            | (row["_zq"] < z[0])
            | (row["_zq"] > z[-1])
            | (wavelength < wave[0])
            | (wavelength > wave[-1])
        )
        variance = np.full(m.shape, 1e20)
        inside = ~outside
        snr = interp(
            np.column_stack(
                (
                    np.maximum(m[inside], mags[0]),
                    np.full(sum(inside), row["_zq"]),
                    np.full(sum(inside), wavelength),
                )
            )
        )
        variance[inside] = (
            1 / np.fmax(snr * np.sqrt(pixel) * np.sqrt(4 / 4), 1e-10)
        ) ** 2
        return density, variance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--bin", type=int, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    raw = RawSampling(args.reference)
    path = args.input / f"records-{2 * args.bin:03d}.report.json"
    report = json.loads(path.read_text())
    records = []
    saved = {}
    reader_checks = {}
    for name, row in report["settings"]["pair_inputs"].items():
        if "density" not in row:
            continue
        original = WeightInputs.from_pair(row)
        density, variance = raw.sample(name, row, original.magnitudes)
        checks = dict(
            density=relative_change(density, original.density),
            variance=relative_change(variance, original.variance),
        )
        if max(checks.values()) > 2e-12:
            raise ValueError((name, checks))
        reader_checks[name] = checks
        for nodes in (107, 213, 425):
            m = np.linspace(original.magnitudes[0], original.magnitudes[-1], nodes)
            density, variance = raw.sample(name, row, m)
            inputs = (
                original
                if nodes == 107
                else replace(
                    original,
                    magnitudes=m,
                    density=density,
                    variance=variance,
                    quadrature=np.full(m.shape, m[1] - m[0]),
                )
            )
            prefix = f"{name}/{nodes}"
            saved[prefix + "/magnitudes"] = inputs.magnitudes
            saved[prefix + "/density"] = inputs.density
            saved[prefix + "/variance"] = inputs.variance
            for variant in VARIANTS:
                data = trajectory(inputs, variant)
                if data["failure"] is not None:
                    last = data["weights"][-1]
                    j1 = np.cumsum(inputs.density * last * inputs.quadrature)
                    data["failure"].update(
                        last_finite_prefix_J1_min_abs=float(np.min(np.abs(j1))),
                        last_finite_full_J1=float(j1[-1]),
                        last_finite_weight_min_abs=float(np.min(np.abs(last))),
                    )
                key = prefix + "/" + variant
                for quantity in (
                    "weights",
                    "coefficients",
                    "changes",
                    "residual",
                    "amplitude",
                ):
                    saved[key + "/" + quantity] = data[quantity]
                results = {str(tol): classify(data, tol) for tol in (1e-3, 1e-4)}
                records.append(
                    dict(
                        bin=args.bin,
                        context=name,
                        nodes=nodes,
                        variant=variant,
                        auto=name in ("lya(qso)_lya(qso)", "lya(lbg)_lya(lbg)"),
                        results=results,
                        last_iteration=len(data["weights"]) - 1,
                        failure=data["failure"],
                        negative_density_nodes=int(sum(inputs.density < 0)),
                        signal=inputs.signal,
                        negative_moment_states=np.sum(
                            data["coefficients"][:, :3] < 0, axis=0
                        ).tolist(),
                        zero_moment_states=np.sum(
                            data["coefficients"][:, :3] == 0, axis=0
                        ).tolist(),
                        final_coefficients=data["coefficients"][-1].tolist(),
                        amplitude_ratio=float(
                            data["amplitude"][-1] / data["amplitude"][0]
                        ),
                    )
                )
    np.savez_compressed(args.output / f"bin-{args.bin}.npz", **saved)
    result = dict(
        bin=args.bin,
        report_sha256=digest(path),
        source=str(path.resolve()),
        reader_checks=reader_checks,
        raw_hashes=raw.hashes,
        records=records,
    )
    (args.output / f"bin-{args.bin}.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    print(
        args.bin,
        Counter(r["results"]["0.001"]["status"] for r in records if r["nodes"] == 107),
    )


if __name__ == "__main__":
    main()
