"""Bounded actual Tracer/Spectrograph branch checks, without CAMB/forecast runs."""

import argparse
import configparser
import json
from pathlib import Path

import numpy as np

from fishhighz.adapters.legacy_compat import LegacyDensity, LegacySNR, plain
from fishhighz.adapters.legacy_inputs import DensityReader, SNRReader
from fishhighz.validation.evidence import digest


def run(reference, saved):
    from lyaforecast.spectrograph import Spectrograph
    from lyaforecast.survey import Survey
    from lyaforecast.tracer import Tracer

    root = Path(reference).resolve()
    history = json.loads(Path(saved).read_text())
    hashes = {}
    for old, sha in history["source_resource_hashes"].items():
        p = root / old.split("/lib/lyaforecast/", 1)[1]
        assert digest(p) == sha
        hashes[str(p)] = sha
    c = configparser.ConfigParser()
    c.read(root / "examples/desi2/lya_qso_lbg_lae_15x2pt.ini")
    survey = Survey(c)
    comparisons = []
    extensions = []

    def compare(name, new, old):
        a, b = np.asarray(new), np.asarray(old)
        np.testing.assert_allclose(a, b, rtol=5e-12, atol=0)
        comparisons.append(
            dict(
                name=name,
                new=a.tolist(),
                reference=b.tolist(),
                max_relative=float(np.max(abs(a - b) / np.maximum(abs(b), 1e-300))),
            )
        )

    for key in [s for s in c.sections() if s.startswith("tracer ")]:
        t = Tracer(c[key])
        forest = t.type == "continuous"
        reader = DensityReader(
            root / "lyaforecast/resources/data" / c[key]["dn dz"],
            semantics="cell_count_per_deg2",
            target_density=t.tracer_density,
            z_norm_min=2.15 if forest or t.simple_name == "qso" else None,
            magnitude_bounds=(t.mag_min, t.mag_max) if forest else None,
            width_policy="legacy_first_spacing",
        )
        d = LegacyDensity(reader, "floor_negative")
        mags = np.array(
            [
                reader.magnitudes[0] - 0.01,
                reader.magnitudes[0],
                reader.magnitudes[0] + 0.01,
                reader.magnitudes[-1] - 0.01,
                reader.magnitudes[-1],
                reader.magnitudes[-1] + 0.01,
            ]
        )
        for z in [
            max(0, reader.z[0] - 0.01),
            reader.z[0],
            reader.z[0] + 0.01,
            reader.z[-1] - 0.01,
            reader.z[-1],
            reader.z[-1] + 0.01,
        ]:
            sample = d.sample(z, mags)
            old = t.get_dn_dzdm(np.full(len(mags), z), mags)
            expected = old.copy()
            expected[old < 0] = 1e-20
            compare(t.name + "/density-policy", sample["values"], expected)
            extensions.append(
                dict(
                    name=t.name,
                    query=plain(sample["provenance"]),
                    equality_to_legacy=bool(np.all(old >= 0)),
                )
            )
        if not forest:
            continue
        bad = next(
            x for x in history["rejected_initial_samples"] if x["field"] == t.name
        )
        old = t.get_dn_dzdm(
            np.full(len(bad["magnitudes"]), bad["z_source"]),
            np.array(bad["magnitudes"]),
        )
        assert np.any(old < 0)
        sample = d.sample(bad["z_source"], bad["magnitudes"])
        compare(
            t.name + "/negative-extension",
            sample["values"],
            np.where(old < 0, 1e-20, old),
        )
        extensions.append(
            dict(
                name=t.name,
                negative_reference=old.tolist(),
                query=plain(sample["provenance"]),
                equality_to_legacy=False,
            )
        )
        reference_snr = Spectrograph(c, survey, t)
        paths = sorted(
            (root / "lyaforecast/resources/data" / c[key]["snr-file-dir"]).glob("*.dat")
        )
        s = LegacySNR(SNRReader(paths, smoothing="legacy"))
        r = s.reader
        midz = float(np.mean(r.z))
        midw = float(np.mean(r.wavelength))
        midm = float(np.mean(r.magnitudes))
        queries = [
            (m, midz, midw)
            for m in [
                r.magnitudes[0] - 0.01,
                r.magnitudes[0],
                r.magnitudes[-1],
                r.magnitudes[-1] + 0.01,
            ]
        ]
        queries += [
            (midm, z, midw) for z in [r.z[0] - 0.01, r.z[0], r.z[-1], r.z[-1] + 0.01]
        ]
        queries += [
            (midm, midz, w)
            for w in [
                r.wavelength[0] - 0.01,
                r.wavelength[0],
                r.wavelength[-1],
                r.wavelength[-1] + 0.01,
            ]
        ]
        for m, z, w in queries:
            for pixel, count in [(1.0, 4), (2.0, 8), (1e-30, 4)]:
                sample = s.sample(
                    z_source=z,
                    magnitudes=[m],
                    wavelength=w,
                    pixel_width_angstrom=pixel,
                    exposure_count=count,
                )
                old = (
                    np.asarray(
                        reference_snr.get_pixel_rms_noise(m, z, w, pixel, count)
                    ).reshape(-1)
                    ** 2
                )
                compare(t.name + "/SNR", sample["values"], old)
    return dict(
        source_hashes=hashes,
        comparisons=comparisons,
        extensions=extensions,
        reference="actual Tracer and per-population Spectrograph methods; no NewForecast",
        max_relative=max(x["max_relative"] for x in comparisons),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True)
    p.add_argument("--saved", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    out = Path(a.output)
    if out.exists():
        raise FileExistsError(out)
    result = run(a.reference, a.saved)
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        len(result["comparisons"]), "checks; maximum relative", result["max_relative"]
    )
