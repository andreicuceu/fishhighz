#!/usr/bin/env python3
"""Isolated scientific worker for ``lyaforecast_baseline.py``."""

from __future__ import annotations

import argparse
import configparser
import importlib
import importlib.metadata
import json
import math
import os
import pickle
import platform
import sys
import time
from pathlib import Path

THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def write_json(path: Path, value: object) -> None:
    """Write strict JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def prepare_imports(reference: Path):
    """Import lyaforecast from the requested checkout only."""
    for variable in THREAD_VARIABLES:
        if os.environ.get(variable) != "1":
            raise RuntimeError(f"{variable} must be set to 1 before worker startup")
    sys.dont_write_bytecode = True
    reference = reference.resolve()
    sys.path.insert(0, str(reference))
    forecast_module = importlib.import_module("lyaforecast.forecast_new")
    expected = (reference / "lyaforecast").resolve()
    actual = Path(forecast_module.__file__).resolve()
    if actual != expected / "forecast_new.py":
        raise RuntimeError(
            f"lyaforecast.forecast_new imported from {actual}, expected {expected}"
        )
    return forecast_module


def module_origin(module) -> str | None:
    """Return a resolved module origin when one exists."""
    origin = getattr(module, "__file__", None)
    return str(Path(origin).resolve()) if origin else None


def probe(reference: Path) -> dict[str, object]:
    """Describe the isolated scientific environment and module origins."""
    forecast_module = prepare_imports(reference)
    modules = {
        name: importlib.import_module(name)
        for name in ("lyaforecast", "numpy", "scipy", "camb")
    }
    modules["lyaforecast.forecast_new"] = forecast_module
    distributions = {}
    for name in ("numpy", "scipy", "camb", "lyaforecast"):
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "not-installed"
    return {
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "module_origins": {
            name: module_origin(module) for name, module in modules.items()
        },
        "distribution_versions": distributions,
        "thread_settings": {
            variable: os.environ.get(variable) for variable in THREAD_VARIABLES
        },
        "dont_write_bytecode": sys.dont_write_bytecode,
    }


def read_config(path: Path) -> configparser.ConfigParser:
    """Read a case INI with its original option spellings."""
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    if parser.read(path, encoding="utf-8") != [str(path)]:
        raise RuntimeError(f"could not read {path}")
    return parser


def resolve(reference: Path, config_path: Path) -> dict[str, object]:
    """Resolve every scientific input path through lyaforecast utilities."""
    prepare_imports(reference)
    from lyaforecast.utils import get_dir, get_file

    config = read_config(config_path)
    requests = [("cosmo", "filename", "file")]
    for section in config.sections():
        if not section.startswith("tracer "):
            continue
        requests.append((section, "dn dz", "file"))
        if "snr-file-dir" in config[section]:
            requests.append((section, "snr-file-dir", "directory"))
    paths = []
    for section, option, kind in requests:
        setting = config[section][option]
        resolved = get_file(setting) if kind == "file" else get_dir(setting)
        paths.append(
            {
                "section": section,
                "option": option,
                "kind": kind,
                "original_setting": setting,
                "resolved_path": str(Path(resolved).resolve()),
            }
        )
    return {
        "config": str(config_path.resolve()),
        "resolver_module": module_origin(importlib.import_module("lyaforecast.utils")),
        "paths": paths,
    }


def encode_number(value: object) -> object:
    """Encode finite and non-finite numbers as strict JSON."""
    number = float(value)
    if math.isnan(number):
        return "nan"
    if math.isinf(number):
        return "+inf" if number > 0 else "-inf"
    if isinstance(value, int):
        return int(value)
    return number


def encode_value(value: object) -> dict[str, object]:
    """Encode nested forecast output with order, array shape, and dtype metadata."""
    import numpy as np

    if isinstance(value, dict):
        return {
            "kind": "dict",
            "items": [
                {"key": str(key), "value": encode_value(item)}
                for key, item in value.items()
            ],
        }
    if isinstance(value, np.ndarray):
        return {
            "kind": "ndarray",
            "dtype": str(value.dtype),
            "shape": list(value.shape),
            "data": [encode_number(item) for item in value.reshape(-1)],
        }
    if isinstance(value, np.generic):
        return {
            "kind": "numpy-scalar",
            "dtype": str(value.dtype),
            "python_type": type(value.item()).__name__,
            "value": encode_number(value),
        }
    if value is None:
        return {"kind": "none"}
    if isinstance(value, bool):
        return {"kind": "bool", "value": value}
    if isinstance(value, (int, float)):
        return {
            "kind": "python-scalar",
            "dtype": None,
            "python_type": type(value).__name__,
            "value": encode_number(value),
        }
    if isinstance(value, str):
        return {"kind": "string", "value": value}
    raise TypeError(f"unsupported result value {type(value).__name__}")


def decode_value(node: dict[str, object]) -> object:
    """Decode the typed JSON representation for a worker-side round trip."""
    import numpy as np

    kind = node["kind"]

    def number(value: object) -> float | int:
        if value == "nan":
            return float("nan")
        if value == "+inf":
            return float("inf")
        if value == "-inf":
            return float("-inf")
        return value

    if kind == "dict":
        return {item["key"]: decode_value(item["value"]) for item in node["items"]}
    if kind == "ndarray":
        return np.array(
            [number(value) for value in node["data"]], dtype=node["dtype"]
        ).reshape(node["shape"])
    if kind == "numpy-scalar":
        return np.array(number(node["value"]), dtype=node["dtype"])[()]
    if kind == "python-scalar":
        return {"int": int, "float": float}[node["python_type"]](number(node["value"]))
    if kind == "string":
        return node["value"]
    if kind == "bool":
        return node["value"]
    if kind == "none":
        return None
    raise TypeError(f"unsupported encoded kind {kind!r}")


def equal_results(left: object, right: object) -> bool:
    """Compare a result recursively, including dictionary order and ndarray dtypes."""
    import numpy as np

    if isinstance(left, dict):
        return (
            isinstance(right, dict)
            and list(left) == list(right)
            and all(equal_results(left[key], right[key]) for key in left)
        )
    if isinstance(left, np.ndarray):
        return (
            isinstance(right, np.ndarray)
            and left.dtype == right.dtype
            and left.shape == right.shape
            and np.array_equal(left, right, equal_nan=True)
        )
    if isinstance(left, np.generic):
        return type(left) is type(right) and bool(
            np.array_equal(left, right, equal_nan=True)
        )
    if isinstance(left, float) and math.isnan(left):
        return isinstance(right, float) and math.isnan(right)
    return type(left) is type(right) and left == right


def numeric_array(value) -> dict[str, object]:
    """Encode a metadata array using the same self-describing representation."""
    return encode_value(value)


def result_summary(
    data: dict[str, object], selected_pairs: list[str]
) -> dict[str, object]:
    """Build a concise, unrounded result summary."""
    import numpy as np

    summaries = {}
    for owner in [*selected_pairs, "total"]:
        forecast = data if owner == "total" else data[owner]
        summaries[owner] = {
            metric: {
                "minimum": encode_number(np.min(forecast[metric])),
                "maximum": encode_number(np.max(forecast[metric])),
            }
            for metric in ("sigma_at", "sigma_ap", "corr_coef")
        }
    return {
        "result_key_order": list(data),
        "selected_pairs": selected_pairs,
        "forecasts": summaries,
    }


def run_forecast(reference: Path, request_path: Path) -> dict[str, object]:
    """Run one unchanged full forecast and serialize all returned values."""
    forecast_module = prepare_imports(reference)
    request = json.loads(request_path.read_text(encoding="utf-8"))
    config_path = Path(request["effective_config"])
    started = time.perf_counter()
    forecast = forecast_module.NewForecast(config_path)
    initialized = time.perf_counter()
    metadata = {
        "module_origin": module_origin(forecast_module),
        "tracer_order": list(forecast.tracers),
        "all_pairs": list(forecast.correlations),
        "selected_pairs": list(forecast.correlations_to_compute),
        "redshift_centers": numeric_array(forecast.survey.z_bin_centres),
        "redshift_edges_2d": numeric_array(forecast.survey.z_bin_edges),
        "grid": {
            "k": numeric_array(forecast.power_spectrum.k),
            "mu": numeric_array(forecast.power_spectrum.mu),
        },
    }
    data = forecast.new_run_forecast()
    finished = time.perf_counter()
    result_pickle = Path(request["result_pickle"])
    with result_pickle.open("xb") as stream:
        pickle.dump(data, stream, protocol=5)
    with result_pickle.open("rb") as stream:
        pickle_roundtrip = pickle.load(stream)
    encoded = encode_value(data)
    write_json(Path(request["result_json"]), encoded)
    decoded = decode_value(
        json.loads(Path(request["result_json"]).read_text(encoding="utf-8"))
    )
    roundtrip_equal = equal_results(data, pickle_roundtrip) and equal_results(
        data, decoded
    )
    if not roundtrip_equal:
        raise RuntimeError("pickle/JSON serialization round-trip changed the result")
    write_json(Path(request["metadata"]), metadata)
    write_json(
        Path(request["summary"]), result_summary(data, metadata["selected_pairs"])
    )
    return {
        "roundtrip_equal": roundtrip_equal,
        "initialization_seconds": initialized - started,
        "forecast_seconds": finished - initialized,
        "worker_total_seconds": time.perf_counter() - started,
        "module_origin": module_origin(forecast_module),
        "thread_settings": {
            variable: os.environ.get(variable) for variable in THREAD_VARIABLES
        },
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the private worker command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("probe", "resolve", "run"))
    parser.add_argument("--reference-checkout", required=True, type=Path)
    parser.add_argument("--response", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--request", type=Path)
    return parser


def main() -> int:
    """Run one isolated action."""
    args = build_parser().parse_args()
    if args.action == "probe":
        response = probe(args.reference_checkout)
    elif args.action == "resolve":
        if args.config is None:
            raise RuntimeError("resolve requires --config")
        response = resolve(args.reference_checkout, args.config)
    else:
        if args.request is None:
            raise RuntimeError("run requires --request")
        response = run_forecast(args.reference_checkout, args.request)
    write_json(args.response, response)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
