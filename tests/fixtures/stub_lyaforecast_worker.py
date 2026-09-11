"""Standard-library worker stub for portable baseline-tool tests."""

import argparse
import configparser
import json
import os
import sys
from pathlib import Path

CASE_PAIRS = {
    "lbg_lae_3x2pt.ini": (
        ["lbg_lbg", "lbg_lae", "lae_lae"],
        3,
        ["lbg", "lae"],
    ),
    "lya_lbg_lae_3x2pt.ini": (
        [
            "lya(lbg)_lya(lbg)",
            "lya(lbg)_lbg",
            "lya(lbg)_lae",
            "lbg_lbg",
            "lbg_lae",
            "lae_lae",
        ],
        3,
        ["lya(lbg)", "lbg", "lae"],
    ),
    "lya_lbg_lae_6x2pt.ini": (
        [
            "lya(lbg)_lya(lbg)",
            "lya(lbg)_lbg",
            "lya(lbg)_lae",
            "lbg_lbg",
            "lbg_lae",
            "lae_lae",
        ],
        6,
        ["lya(lbg)", "lbg", "lae"],
    ),
    "lya_qso_2x2pt.ini": (
        ["lya(qso)_lya(qso)", "lya(qso)_qso", "qso_qso"],
        2,
        ["lya(qso)", "qso"],
    ),
}
FIVE_TRACERS = ["lya(qso)", "qso", "lbg", "lae", "lya(lbg)"]
FIVE_PAIRS = [
    f"{left}_{right}"
    for index, left in enumerate(FIVE_TRACERS)
    for right in FIVE_TRACERS[index:]
]
for count in (4, 8, 15):
    CASE_PAIRS[f"lya_qso_lbg_lae_{count}x2pt.ini"] = (
        FIVE_PAIRS,
        count,
        FIVE_TRACERS,
    )


def write_json(path, value):
    """Write strict JSON for the controller."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def array(values, dtype="float64", shape=None):
    """Return a typed-array JSON node."""
    return {
        "kind": "ndarray",
        "dtype": dtype,
        "shape": shape if shape is not None else [len(values)],
        "data": values,
    }


def forecast(selected):
    """Return one selected or zero-filled pair result."""
    values = [0.1, 0.2] if selected else [0.0, 0.0]
    correlations = [-0.2, 0.3] if selected else [0.0, 0.0]
    return {
        "kind": "dict",
        "items": [
            {"key": "sigma_at", "value": array(values)},
            {"key": "sigma_ap", "value": array(values)},
            {"key": "corr_coef", "value": array(correlations)},
        ],
    }


def result(case, role):
    """Return a complete synthetic result for one authoritative case."""
    pairs, count, _ = CASE_PAIRS[case]
    selected = pairs[:count]
    items = [
        {"key": "redshifts", "value": array([2.1, 2.3])},
        {"key": "zedges", "value": array([2.0, 2.2, 2.4])},
        {
            "key": "fiducial redshift",
            "value": {
                "kind": "python-scalar",
                "dtype": None,
                "python_type": "float",
                "value": 2.3,
            },
        },
    ]
    items.extend({"key": pair, "value": forecast(pair in selected)} for pair in pairs)
    total_at = [0.05, 0.06]
    mode = (
        os.environ.get("FISHHIGHZ_STUB_REPEAT_MODE")
        if role == "repeat"
        else os.environ.get("FISHHIGHZ_STUB_QUICK_MODE")
    )
    if mode == "number":
        total_at[0] = 0.5
    items.extend(
        [
            {"key": "sigma_at", "value": array(total_at)},
            {"key": "sigma_ap", "value": array([0.07, 0.08])},
            {"key": "corr_coef", "value": array([-0.1, 0.2])},
        ]
    )
    if mode == "structure":
        items[-1], items[-2] = items[-2], items[-1]
    return {"kind": "dict", "items": items}


def read_config(path):
    """Read a test INI without interpolation."""
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser.read(path)
    return parser


def resolve(args):
    """Resolve fixture paths under the synthetic reference checkout."""
    config = read_config(args.config)
    requests = [("cosmo", "filename", "file")]
    for section in config.sections():
        if section.startswith("tracer "):
            requests.append((section, "dn dz", "file"))
            if "snr-file-dir" in config[section]:
                requests.append((section, "snr-file-dir", "directory"))
    paths = []
    for section, option, kind in requests:
        setting = config[section][option]
        paths.append(
            {
                "section": section,
                "option": option,
                "kind": kind,
                "original_setting": setting,
                "resolved_path": str(
                    (args.reference_checkout / "resources" / setting).resolve()
                ),
            }
        )
    return {
        "config": str(args.config.resolve()),
        "resolver_module": str(
            (args.reference_checkout / "lyaforecast" / "utils.py").resolve()
        ),
        "paths": paths,
    }


def run(args):
    """Create the artifacts normally emitted by the real scientific worker."""
    request = json.loads(args.request.read_text())
    case = request["case"]
    if (
        os.environ.get("FISHHIGHZ_STUB_FAIL_CASE") == case
        or os.environ.get("FISHHIGHZ_STUB_FAIL_ROLE") == request["role"]
    ):
        raise RuntimeError(f"requested stub failure for {case}")
    mutation = os.environ.get("FISHHIGHZ_STUB_MUTATE")
    reference = args.reference_checkout
    targets = {
        "source-change": reference / "lyaforecast" / "forecast_new.py",
        "source-add": reference / "lyaforecast" / "added.py",
        "source-delete": reference / "lyaforecast" / "utils.py",
        "input-change": reference / "resources" / "SNR" / "header-and-data.dat",
        "input-add": reference / "resources" / "SNR" / "added.dat",
        "input-delete": reference / "resources" / "SNR" / "header-and-data.dat",
        "config-change": reference / "examples" / "desi2" / case,
        "config-add": reference / "examples" / "desi2" / "added.ini",
        "config-delete": reference / "examples" / "desi2" / case,
    }
    if mutation in {"source-change", "input-change", "config-change"}:
        targets[mutation].write_text(targets[mutation].read_text() + "# mutated\n")
    elif mutation in {"source-add", "input-add", "config-add"}:
        targets[mutation].write_text("# added during run\n")
    elif mutation in {"source-delete", "input-delete", "config-delete"}:
        targets[mutation].unlink()
    pairs, count, tracers = CASE_PAIRS[case]
    selected = pairs[:count]
    encoded = result(case, request["role"])
    write_json(request["result_json"], encoded)
    Path(request["result_pickle"]).write_bytes(b"trusted-stub-pickle\n")
    metadata = {
        "module_origin": str(
            (args.reference_checkout / "lyaforecast" / "forecast_new.py").resolve()
        ),
        "tracer_order": tracers,
        "all_pairs": pairs,
        "selected_pairs": selected,
        "redshift_centers": array([2.1, 2.3]),
        "redshift_edges_2d": array([2.0, 2.2, 2.2, 2.4], shape=[2, 2]),
        "grid": {
            "k": array([0.01, 0.01 + 0.49 / 3, 0.01 + 2 * 0.49 / 3, 0.5]),
            "mu": array([0.125, 0.375, 0.625, 0.875]),
        },
    }
    write_json(request["metadata"], metadata)
    write_json(
        request["summary"],
        {"result_key_order": [item["key"] for item in encoded["items"]]},
    )
    output = Path(request["effective_config"])
    forecast_path = Path(read_config(output)["output"]["filename"])
    forecast_path.parent.mkdir(parents=True, exist_ok=True)
    (forecast_path.parent / "forecast.log").write_text(f"stub {case}\n")
    return {
        "roundtrip_equal": os.environ.get("FISHHIGHZ_STUB_ROUNDTRIP") != "false",
        "initialization_seconds": 0.01,
        "forecast_seconds": 0.02,
        "worker_total_seconds": 0.03,
        "module_origin": metadata["module_origin"],
        "thread_settings": {
            name: os.environ.get(name)
            for name in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
            )
        },
    }


def main():
    """Run one stub worker action."""
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("probe", "resolve", "run"))
    parser.add_argument("--reference-checkout", required=True, type=Path)
    parser.add_argument("--response", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--request", type=Path)
    args = parser.parse_args()
    action_log = os.environ.get("FISHHIGHZ_STUB_ACTION_LOG")
    if action_log:
        case = None
        if args.request and args.request.is_file():
            case = json.loads(args.request.read_text()).get("case")
        with Path(action_log).open("a") as stream:
            stream.write(json.dumps({"action": args.action, "case": case}) + "\n")
    if args.action == "probe":
        response = {
            "python_version": sys.version,
            "python_executable": sys.executable,
            "platform": sys.platform,
            "module_origins": {
                "lyaforecast.forecast_new": str(
                    (
                        args.reference_checkout / "lyaforecast" / "forecast_new.py"
                    ).resolve()
                )
            },
            "distribution_versions": {
                "numpy": os.environ.get("FISHHIGHZ_STUB_NUMPY_VERSION", "stub"),
                "scipy": "stub",
                "camb": "stub",
                "lyaforecast": "stub",
            },
            "thread_settings": {
                name: os.environ.get(name)
                for name in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                )
            },
            "dont_write_bytecode": True,
        }
    elif args.action == "resolve":
        response = resolve(args)
    else:
        response = run(args)
    write_json(args.response, response)


if __name__ == "__main__":
    main()
