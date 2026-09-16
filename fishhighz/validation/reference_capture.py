"""Fresh seven-case reference capture with read-only upstream-array recording.

A validation-owned Fisher subclass observes the actual NewForecast calls. Its
return value is always the unchanged reference result. Captured upstream arrays
are inputs to a separate FishHighz covariance/Fisher calculation, never its F.
"""

import configparser
import contextlib
import importlib.metadata
import json
import sys
import time
from pathlib import Path

import numpy as np

from ..adapters.legacy_compat import plain
from .cases import CASE_IDS, bins, selection, verify_inventory
from .evidence import digest


def imported_reference(root):
    """Require the actually imported reference package to be the inventoried tree."""
    import camb
    import lyaforecast

    root = Path(root).resolve()
    actual = Path(lyaforecast.__file__).resolve().parent
    if actual != (root / "lyaforecast").resolve():
        raise ValueError(f"reference import {actual} differs from requested {root}")
    sources = {str(p.resolve()): digest(p) for p in actual.glob("*.py")}
    return dict(
        reference_origin=str(actual),
        camb_origin=camb.__file__,
        python=sys.executable,
        versions={
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "astropy", "camb", "lyaforecast")
        },
        sources=sources,
    )


def resolved_resources(root, config):
    """Verify resolved files/directories against the explicitly requested checkout."""
    from lyaforecast.utils import get_dir, get_file

    root = Path(root).resolve()
    inputs = [Path(get_file(config["cosmo"]["filename"])).resolve()]
    for key in config.sections():
        if key.startswith("tracer "):
            inputs.append(Path(get_file(config[key]["dn dz"])).resolve())
            if "snr-file-dir" in config[key]:
                inputs.extend(
                    p.resolve()
                    for p in get_dir(config[key]["snr-file-dir"]).glob("*.dat")
                )
    for path in inputs:
        if not path.is_relative_to(root / "lyaforecast/resources"):
            raise ValueError(f"unexpected resolved reference resource {path}")
    return {str(p): digest(p) for p in inputs}


def capture_case(root, original, destination, *, case):
    """Run one actual NewForecast, preserving upstream arrays for every bin."""
    import lyaforecast.forecast_new as module
    from lyaforecast.fisher import Fisher

    out = Path(destination)
    out.mkdir(parents=True, exist_ok=False)
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(original)
    original_bytes = Path(original).read_bytes()
    (out / "original.ini").write_bytes(original_bytes)
    config["output"]["filename"] = str((out / "outputs/forecast").resolve())
    with (out / "effective.ini").open("w") as stream:
        config.write(stream)
    identity = imported_reference(root)
    identity["resources"] = resolved_resources(root, config)
    selected = selection(case)
    fields = [f.id for f in selected.fields]
    all_labels = [f"{fields[i]}_{fields[j]}" for i, j in selected.required_pairs]
    selected_labels = [f"{fields[i]}_{fields[j]}" for i, j in selected.selected_pairs]
    single = {}
    records = []
    forecast = None

    class CaptureFisher(Fisher):
        def compute_fisher(self, models, measurements, spectra_list):
            value = super().compute_fisher(models, measurements, spectra_list)
            index = int(self.zbin_index)
            if len(spectra_list) == 1:
                single.setdefault(index, {})[spectra_list[0]] = value.copy()
            if list(spectra_list) != selected_labels:
                return value
            pk = self._power_spec
            derivative = np.stack(
                [
                    self.compute_derivatives(
                        np.stack([models[n][i] for n in spectra_list]), mu, spectra_list
                    ).T
                    for i, mu in enumerate(pk.mu)
                ]
            )
            k = np.tile(pk.k, len(pk.mu))
            mu = np.repeat(pk.mu, len(pk.k))
            observed_j = (
                derivative.reshape(-1, len(spectra_list))[:, :, None]
                * np.column_stack((mu**2, 1 - mu**2))[:, None, :]
            )
            total = np.column_stack(
                [measurements[label].reshape(-1) for label in all_labels]
            )
            mean = np.column_stack(
                [models[label].reshape(-1) for label in spectra_list]
            )
            arrays = dict(
                k=k,
                mu=mu,
                modes=np.tile(self._num_modes, len(pk.mu)),
                mean=mean,
                total=total,
                observed_j=observed_j,
                selected_pairs=selected.selected_pairs,
                required_pairs=selected.required_pairs,
                legacy_fisher=value,
                legacy_pair_fisher=np.stack(
                    [single[index][n] for n in selected_labels]
                ),
                k_axis=pk.k,
                mu_axis=pk.mu,
                dlogk=pk.dlogk,
            )
            name = f"bin-{index:02d}.npz"
            np.savez_compressed(out / name, **arrays)
            pair_inputs = {}
            for label, cov in forecast._covariance.items():
                row = {
                    key: plain(getattr(cov, key, None))
                    for key in (
                        "_z_mean",
                        "_zq",
                        "_pix_kms",
                        "_res_kms",
                        "forest_length",
                        "_distance_to_velocity",
                        "_angle_to_distance",
                        "_aliasing_weights",
                        "_effective_noise_power",
                        "_w_lya",
                        "_volume",
                    )
                }
                weights = cov._weights
                row["magnitudes"] = weights.maglist.tolist()
                if weights._lya_tracer is not None:
                    row.update(
                        density=weights._get_dn_dkmsdm(
                            weights._zq, weights.maglist, weights._lya_tracer
                        ).tolist(),
                        variance=weights._get_pix_var_m().tolist(),
                        auxiliary_signal=plain(weights._p3d_w),
                        auxiliary_p1d=plain(weights._p1d_w),
                    )
                pair_inputs[label] = row
            lo, hi = bins(case)[index]
            records.append(
                dict(
                    case=case,
                    bin=index,
                    bounds=[lo, hi],
                    fields=fields,
                    selected_labels=selected_labels,
                    required_labels=all_labels,
                    mean_z=float(forecast.survey.z_bin_centres[index]),
                    covariance_z=float(np.sqrt((1 + lo) * (1 + hi)) - 1),
                    array=name,
                    sha256=digest(out / name),
                    pair_inputs=pair_inputs,
                    growth_ratio=float(forecast.cosmo.growth_factor_ratios[index]),
                    growth_rate=float(forecast.cosmo.growth_rate_zbins[index]),
                    k_measure="literal endpoint-inclusive linspace, dk=(kmax-kmin)/(nk-1)",
                    response_ownership="reference mean/J already observed; no additional response",
                    parameter_order=["ap", "at"],
                )
            )
            return value

    module.Fisher = CaptureFisher
    start = time.monotonic()
    try:
        forecast = module.NewForecast(out / "effective.ini")
        result = forecast.new_run_forecast()
    finally:
        module.Fisher = Fisher
    if [r["bin"] for r in records] != list(range(len(bins(case)))):
        raise ValueError("incomplete captured reference bin inventory")
    (out / "result.json").write_text(
        json.dumps(plain(result), indent=2, allow_nan=False) + "\n"
    )
    metadata = dict(
        case=case,
        complete=True,
        seconds=time.monotonic() - start,
        original_sha256=digest(out / "original.ini"),
        effective_sha256=digest(out / "effective.ini"),
        identity=identity,
        records=records,
        result_sha256=digest(out / "result.json"),
    )
    (out / "metadata.json").write_text(
        json.dumps(plain(metadata), indent=2, allow_nan=False) + "\n"
    )
    return metadata


def capture(root, output):
    """Serial fresh seven-case capture; failures remain visible and never resume."""
    root, out = Path(root).resolve(), Path(output).resolve()
    inventory = verify_inventory(root / "examples/desi2")
    out.mkdir(parents=True, exist_ok=False)
    manifest = dict(
        kind="fresh-legacy-upstream",
        schema=2,
        requested=list(CASE_IDS),
        records=[],
        complete=False,
    )

    def save():
        (out / "manifest.json").write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n"
        )

    save()
    for case in CASE_IDS:
        print(f"Capturing actual reference {case}", flush=True)
        try:
            with (
                (out / f"{case}.log").open("w") as log,
                contextlib.redirect_stdout(log),
                contextlib.redirect_stderr(log),
            ):
                result = capture_case(
                    root, inventory[case]["path"], out / case, case=case
                )
            record = dict(
                case=case,
                status="completed",
                metadata=f"{case}/metadata.json",
                sha256=digest(out / case / "metadata.json"),
                seconds=result["seconds"],
            )
        except Exception as error:
            import traceback

            with (out / f"{case}.log").open("a") as log:
                traceback.print_exc(file=log)
            record = dict(
                case=case, status="failed", error=f"{type(error).__name__}: {error}"
            )
        manifest["records"].append(record)
        save()
        print(record, flush=True)
    manifest["complete"] = all(r["status"] == "completed" for r in manifest["records"])
    save()
    return manifest
