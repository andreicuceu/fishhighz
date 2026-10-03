"""Refinement ladder of the integrated forest-source quadrature on an INI survey.

For an INI that selects ``forest_source_integration = integrated`` (for example
the SRD INIs written to ``$DESI2_SRD_DATA/output_fishhighz_integrated/``) this
script

* runs the BAO forecast of every selected bin for a ladder of quadrature orders
  ``(forest_zq_order, forest_lambda_order, forest_lambda_panels)``,
* reports per bin and forest field the first moment N1 (deg^-2), A (deg^2),
  P_pixel (deg^2 km/s) and z_eff, and per spectrum pair and joint sigma(ap),
  sigma(at), each with the relative change with respect to the previous order,
* times the preparation stages: grid construction, S/N query, density query,
  the rest of the survey preparation, the weight solve, the remaining bin
  preparation and the Fisher derivatives,
* optionally cross-checks N1-N3 at fixed seed weights against an independent
  midpoint rule on a much finer grid (selected bins).

It is a validation script, not part of the quick suite; the default ladder is
heavy (the weight solve scales with n_pixel * n_magnitude) and belongs on a
compute node, with one thread per process::

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    python scripts/check_integrated_forest_refinement.py INI --bins 1,6 \
        --midpoint-bins 1

The exit status is nonzero if the change of any BAO error between the last two
orders exceeds ``--tolerance`` (default 1e-3, the 0.1 per cent target).
"""

import argparse
import configparser
import json
import sys
import tempfile
import time
from contextlib import ExitStack
from pathlib import Path

import numpy as np

import fishhighz
from fishhighz import forecast as forecast_module
from fishhighz import public, survey_config
from fishhighz.geometry import SPEED_LIGHT_KMS

DEFAULT_LADDER = "8,8,2;16,16,4;32,32,8"
TIMING_STAGES = (
    "grid",
    "snr_query",
    "density_query",
    "survey_other",
    "weight_solve",
    "bin_other",
    "fisher",
)


class Timer:
    """Accumulate wall times of named stages.

    Attributes
    ----------
    seconds : dict of str to float
        Accumulated time of each stage in seconds.
    """

    def __init__(self):
        """Create an empty timer."""
        self.seconds = dict.fromkeys(TIMING_STAGES, 0.0)

    def wrap(self, function, stage):
        """Return a function that adds its duration to a stage.

        Parameters
        ----------
        function : callable
            Function to time.
        stage : str
            Key of ``seconds`` receiving the elapsed time.

        Returns
        -------
        wrapped : callable
            Function with the same call signature and return value.
        """

        def timed(*args, **kwargs):
            """Call the wrapped function and record the elapsed time."""
            start = time.perf_counter()
            try:
                return function(*args, **kwargs)
            finally:
                self.seconds[stage] += time.perf_counter() - start

        return timed


class TimedReader:
    """Delegating reader proxy that times the vectorised grid queries.

    Parameters
    ----------
    reader : object
        Density or S/N adapter.
    timer : Timer
        Receives the query times.
    stage : str
        Stage name for ``sample_grid`` / ``query_grid`` / ``variance_grid``.
    """

    def __init__(self, reader, timer, stage):
        """Store the delegate and the timing target."""
        self._reader = reader
        self._timer = timer
        self._stage = stage

    def __getattr__(self, name):
        """Delegate attributes, timing the grid-query methods."""
        attribute = getattr(self._reader, name)
        if name in ("sample_grid", "query_grid", "variance_grid"):
            return self._timer.wrap(attribute, self._stage)
        return attribute


def parse_ladder(text):
    """Parse a ladder of order triples.

    Parameters
    ----------
    text : str
        Triples ``zq,lambda,panels`` separated by semicolons.

    Returns
    -------
    ladder : list of tuple of int
        (forest_zq_order, forest_lambda_order, forest_lambda_panels) per entry.
    """
    ladder = []
    for item in text.split(";"):
        values = tuple(int(value) for value in item.replace(" ", "").split(","))
        if len(values) != 3 or min(values) < 1:
            raise ValueError(f"ladder entry {item!r} must be three positive integers")
        ladder.append(values)
    return ladder


def single_bin_ini(source, bin_index, orders, directory):
    """Write an INI holding one bin of the source survey with given orders.

    Parameters
    ----------
    source : pathlib.Path
        Integrated native INI.
    bin_index : int
        One-based bin of the source INI to keep.
    orders : tuple of int
        (forest_zq_order, forest_lambda_order, forest_lambda_panels).
    directory : pathlib.Path
        Directory receiving the INI.

    Returns
    -------
    path : pathlib.Path
        Written INI. Relative resource paths are made absolute first so the
        copy resolves them like the original.
    """
    parser = configparser.ConfigParser(interpolation=None)
    parser.read(source)
    edges = [float(x) for x in parser["survey"]["z_edges"].replace(",", " ").split()]
    n_bins = len(edges) - 1
    if not 1 <= bin_index <= n_bins:
        raise ValueError(f"bin {bin_index} outside 1..{n_bins}")
    selected = parser[f"pairs bin {bin_index}"]["selected"]
    parser["survey"]["z_edges"] = f"{edges[bin_index - 1]}, {edges[bin_index]}"
    parser["survey"]["num_z_bins"] = "1"
    for index in range(1, n_bins + 1):
        parser.remove_section(f"pairs bin {index}")
    parser.add_section("pairs bin 1")
    parser["pairs bin 1"]["selected"] = selected

    # Resource paths in the source INI are relative to its own directory.
    for section in parser.sections():
        for key in ("density", "snr", "camb_ini", "template"):
            value = parser[section].get(key)
            if (
                value
                and not value.startswith("package:")
                and not Path(value).is_absolute()
            ):
                parser[section][key] = str((source.parent / value).resolve())
    if not parser.has_section("numerical"):
        parser.add_section("numerical")
    for key, value in zip(
        ("forest_zq_order", "forest_lambda_order", "forest_lambda_panels"), orders
    ):
        parser["numerical"][key] = str(value)
    path = directory / f"bin{bin_index}_{'_'.join(map(str, orders))}.ini"
    with path.open("w") as stream:
        parser.write(stream)
    return path


def run_case(path, bin_index, background, template, readers, timer):
    """Prepare and run one single-bin forecast with timed stages.

    Parameters
    ----------
    path : pathlib.Path
        Single-bin integrated INI.
    bin_index : int
        One-based bin index of the source INI, used for reporting.
    background, template : object
        Shared prepared background and template.
    readers : mapping
        Timed density/S-N readers keyed by field ID.
    timer : Timer
        Stage timer, updated in place.

    Returns
    -------
    record : dict
        Per forest field N1, N2, N3, A, P_pixel, z_eff, node counts and
        fallback counts; per pair and joint sigma(ap), sigma(at); the elapsed
        time of each stage; and the prepared forecast (key ``forecast``).
    """
    forecast = fishhighz.Forecast(
        path, background=background, template=template, readers=readers
    )
    start = time.perf_counter()
    prepared = forecast.prepare()
    prepare_total = time.perf_counter() - start
    start = time.perf_counter()
    result = forecast.run()
    run_total = time.perf_counter() - start

    block = prepared.provenance["bins"][0]["forest_source_integration"]
    spec = prepared.survey.bins[0]
    forests = {}
    for field_id, weights in prepared.bins[0].weights.items():
        source = spec.forests[field_id].integrated
        forests[field_id] = dict(
            N1=weights.N1,
            N2=weights.N2,
            N3=weights.N3,
            A=weights.A,
            P_pixel=weights.P_pixel,
            z_eff=weights.z_eff,
            updates=weights.convergence.get("updates"),
            n_zq_nodes=int(len(source.nodes.zq_nodes)),
            n_pixel=int(len(source.nodes.lam_obs)),
            fallback_counts=json.loads(
                json.dumps(block["fallback_counts"][field_id], default=dict)
            ),
        )
    sigma = {}
    for item in (*result.joint, *result.individual):
        name = "joint" if item.pair is None else "x".join(item.pair)
        sigma[name] = [item.sigma_ap, item.sigma_at]
    return dict(
        bin=bin_index,
        forests=forests,
        sigma=sigma,
        prepare_seconds=prepare_total,
        run_seconds=run_total,
        forecast=forecast,
    )


def midpoint_check(spec, field_id, alias, readers, config, *, n_y, n_u):
    """Compare production N1-N3 with a fine midpoint rule at seed weights.

    Parameters
    ----------
    spec : BinSpec
        Prepared integrated bin.
    field_id : str
        Forest field to check.
    alias : float
        Weighting alias B in km/s of the prepared weights.
    readers : mapping
        Untimed density/S-N readers keyed by field ID.
    config : SurveyConfig
        Configuration of the (single) bin.
    n_y : int
        Midpoints in ln(1+z_q) per panel between consecutive breakpoints.
    n_u : int
        Midpoints in ln(lambda) per source over the overlap.

    Returns
    -------
    result : dict
        Production and midpoint ``N1``, ``N2``, ``N3`` (deg^-2) and their
        relative differences.

    Notes
    -----
    The weights are the pixel-local seed of the recurrence, w = (B/Delta_v) /
    ((B/Delta_v) + v) with the pixel variance v of each node, so that the check
    isolates the (z_q, lambda) quadrature from the weight solve. The panel
    boundaries (window ends, overlap kinks, density-table edges and S/N
    source-redshift nodes inside the window) are rebuilt here from the
    configuration, not taken from the production nodes; the midpoint rule is
    second-order accurate on every panel. The magnitude nodes are those of the
    production source.
    """
    source = spec.forests[field_id].integrated
    item = next(f for f in config.fields if f.observed.id == field_id)
    density_reader, snr_reader = readers[field_id]["density"], readers[field_id]["snr"]
    magnitudes, quadrature = source.magnitudes, source.quadrature
    pixel = source.pixel_width_velocity
    reference_power = alias / pixel

    # Logarithmic coordinates of the slice and of the rest-frame forest.
    lya_rest = float(config.survey["lya_rest_angstrom"])
    bin_config = config.bins[0]
    u_low = np.log(lya_rest * (1 + bin_config.z_min))
    u_high = np.log(lya_rest * (1 + bin_config.z_max))
    a_rest, b_rest = np.log(item.min_rest_frame_lya), np.log(item.max_rest_frame_lya)
    y_low, y_high = source.nodes.info["window_y"]

    # Panel boundaries in y: window ends, overlap kinks and table non-smooth points.
    table_z = np.concatenate(
        [
            np.asarray(density_reader.reader.z_edges, dtype=float),
            np.asarray(snr_reader.reader.z, dtype=float),
        ]
    )
    candidates = np.r_[u_low - a_rest, u_high - b_rest, np.log1p(table_z)]
    inside = candidates[(candidates > y_low) & (candidates < y_high)]
    breaks = np.unique(np.r_[y_low, y_high, inside])

    def seed_weights(variance):
        """Return the pixel-local seed weights for a variance array."""
        return reference_power / (reference_power + variance)

    def moments(measure, variance):
        """Return (sum mu w, sum mu w^2, sum mu w^2 v)."""
        weights = seed_weights(variance)
        return np.array(
            [
                np.sum(measure * weights),
                np.sum(measure * weights**2),
                np.sum(measure * weights**2 * variance),
            ]
        )

    production = moments(source.measure, source.variance)

    # Midpoint rule per y panel; the rows of one panel are processed together.
    total = np.zeros(3)
    for lower, upper in zip(breaks[:-1], breaks[1:]):
        step = (upper - lower) / n_y
        y_mid = lower + step * (np.arange(n_y) + 0.5)
        u_start = np.maximum(u_low, y_mid + a_rest)
        overlap = np.minimum(u_high, y_mid + b_rest) - u_start
        u_mid = u_start[:, None] + overlap[:, None] * (np.arange(n_u) + 0.5) / n_u
        zq_mid = np.expm1(y_mid)
        wavelength = np.exp(u_mid).ravel()

        density = density_reader.sample_grid(zq_mid, magnitudes)["values"]
        variance = snr_reader.variance_grid(
            z_source=np.repeat(zq_mid, n_u),
            wavelength=wavelength,
            magnitudes=magnitudes,
            pixel_width_angstrom=pixel * wavelength / SPEED_LIGHT_KMS,
            exposure_count=item.num_exposures,
        )["values"]

        # mu = (1+z_q) c dy du dn/(dz dm) w_m / L_bin with c du / L_bin = du/(u2-u1).
        geometry = (
            ((1 + zq_mid) * step)[:, None]
            * (overlap / n_u)[:, None]
            * np.ones((1, n_u))
            / (u_high - u_low)
        ).ravel()
        measure = (geometry[:, None] * np.repeat(density, n_u, axis=0)) * quadrature[
            None, :
        ]
        total += moments(measure, variance)

    return dict(
        production=production.tolist(),
        midpoint=total.tolist(),
        relative_difference=(production / total - 1).tolist(),
    )


def relative_change(new, old):
    """Return new/old - 1 for scalars, NaN where old is zero or missing.

    Parameters
    ----------
    new, old : float or None
        Values to compare.

    Returns
    -------
    change : float
        Relative change.
    """
    if new is None or old is None or old == 0:
        return float("nan")
    return new / old - 1


def main(argv=None):
    """Run the ladder and print the report.

    Parameters
    ----------
    argv : list of str, optional
        Command line; default ``sys.argv[1:]``.

    Returns
    -------
    status : int
        0 if the last-step change of every BAO error is below ``--tolerance``.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("ini", type=Path, help="integrated native survey INI")
    parser.add_argument(
        "--ladder", default=DEFAULT_LADDER, help="orders 'zq,lam,panels;...'"
    )
    parser.add_argument(
        "--bins", default=None, help="one-based bins, e.g. 1,6 (default all)"
    )
    parser.add_argument(
        "--midpoint-bins", default="", help="bins with the midpoint cross-check"
    )
    parser.add_argument(
        "--midpoint-y", type=int, default=24, help="midpoints per y panel"
    )
    parser.add_argument(
        "--midpoint-u", type=int, default=120, help="midpoints per overlap"
    )
    parser.add_argument("--tolerance", type=float, default=1e-3)
    parser.add_argument(
        "--json", type=Path, default=None, help="write the records here"
    )
    args = parser.parse_args(argv)

    ladder = parse_ladder(args.ladder)
    source_config = fishhighz.parse_survey_ini(args.ini)
    if source_config.forest_integration is None:
        raise SystemExit(f"{args.ini} does not select the integrated mode")
    n_bins = len(source_config.bins)
    bins = (
        list(range(1, n_bins + 1))
        if args.bins is None
        else [int(x) for x in args.bins.split(",")]
    )
    midpoint_bins = [int(x) for x in args.midpoint_bins.split(",") if x]

    records = {}
    with ExitStack() as stack, tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        start = time.perf_counter()
        background = public._default_background(source_config, stack, {})
        template = public._default_template(source_config, background, stack)
        raw_readers = survey_config._normalise_readers(
            source_config, source_config.fields, None, stack
        )
        print(f"background, template and readers: {time.perf_counter() - start:.1f} s")

        for orders in ladder:
            timer = Timer()
            readers = {
                field_id: {
                    "density": TimedReader(item["density"], timer, "density_query"),
                    "snr": None
                    if item["snr"] is None
                    else TimedReader(item["snr"], timer, "snr_query"),
                }
                for field_id, item in raw_readers.items()
            }
            # Stage timers wrap the module-level callables used by the package.
            originals = (
                survey_config.integration_nodes,
                forecast_module.prepare_integrated_forest_weights,
                public.prepare_survey,
                public.prepare_bin,
                public.run_forecast,
            )
            survey_config.integration_nodes = timer.wrap(originals[0], "grid")
            forecast_module.prepare_integrated_forest_weights = timer.wrap(
                originals[1], "weight_solve"
            )
            public.prepare_survey = timer.wrap(originals[2], "survey_other")
            public.prepare_bin = timer.wrap(originals[3], "bin_other")
            public.run_forecast = timer.wrap(originals[4], "fisher")
            try:
                for bin_index in bins:
                    path = single_bin_ini(args.ini, bin_index, orders, directory)
                    case = run_case(
                        path, bin_index, background, template, readers, timer
                    )
                    forecast = case.pop("forecast")
                    if bin_index in midpoint_bins:
                        prepared = forecast.prepare()
                        spec = prepared.survey.bins[0]
                        config = prepared.config
                        untimed = {
                            field_id: {"density": item["density"], "snr": item["snr"]}
                            for field_id, item in raw_readers.items()
                        }
                        case["midpoint"] = {
                            field_id: midpoint_check(
                                spec,
                                field_id,
                                prepared.bins[0].weights[field_id].alias,
                                untimed,
                                config,
                                n_y=args.midpoint_y,
                                n_u=args.midpoint_u,
                            )
                            for field_id in prepared.bins[0].weights
                        }
                    records.setdefault(orders, {})[bin_index] = case
                    print(
                        f"orders {orders} bin {bin_index}: prepare "
                        f"{case['prepare_seconds']:.1f} s, run {case['run_seconds']:.1f} s",
                        flush=True,
                    )
            finally:
                (
                    survey_config.integration_nodes,
                    forecast_module.prepare_integrated_forest_weights,
                    public.prepare_survey,
                    public.prepare_bin,
                    public.run_forecast,
                ) = originals
            # prepare_survey includes the grid and both queries: report the rest.
            timer.seconds["survey_other"] = max(
                0.0,
                timer.seconds["survey_other"]
                - timer.seconds["grid"]
                - timer.seconds["snr_query"]
                - timer.seconds["density_query"],
            )
            timer.seconds["bin_other"] = max(
                0.0, timer.seconds["bin_other"] - timer.seconds["weight_solve"]
            )
            records[orders]["timing"] = dict(timer.seconds)

    print_report(records, ladder, bins)
    failures = check_last_step(records, ladder, bins, args.tolerance)
    if args.json is not None:
        args.json.write_text(
            json.dumps(
                {str(orders): value for orders, value in records.items()},
                indent=1,
                default=float,
            )
        )
    return 1 if failures else 0


def print_report(records, ladder, bins):
    """Print the moments, errors, midpoint comparison and timings.

    Parameters
    ----------
    records : dict
        Ladder records keyed by order triple, then bin index.
    ladder : list of tuple
        Orders in increasing refinement.
    bins : list of int
        Bins reported.
    """
    previous = None
    for orders in ladder:
        print(f"\n=== orders (zq, lambda, panels) = {orders} ===")
        for bin_index in bins:
            case = records[orders][bin_index]
            old = None if previous is None else records[previous][bin_index]
            for field_id, values in case["forests"].items():
                line = f"bin {bin_index} {field_id}: n_zq={values['n_zq_nodes']} n_pix={values['n_pixel']}"
                for name in ("N1", "A", "P_pixel", "z_eff"):
                    change = (
                        float("nan")
                        if old is None
                        else relative_change(
                            values[name], old["forests"][field_id][name]
                        )
                    )
                    line += f" | {name}={values[name]:.8g} ({change:+.2e})"
                print(line)
            for name, sigma in case["sigma"].items():
                changes = (
                    [float("nan")] * 2
                    if old is None
                    else [
                        relative_change(a, b) for a, b in zip(sigma, old["sigma"][name])
                    ]
                )
                print(
                    f"   sigma {name}: ap={sigma[0]:.6g} ({changes[0]:+.2e}) "
                    f"at={sigma[1]:.6g} ({changes[1]:+.2e})"
                )
            for field_id, check in case.get("midpoint", {}).items():
                print(
                    f"   midpoint check {field_id} (N1,N2,N3 at seed weights): "
                    + ", ".join(f"{d:+.2e}" for d in check["relative_difference"])
                )
        timing = records[orders]["timing"]
        print("timing [s]: " + ", ".join(f"{k}={v:.1f}" for k, v in timing.items()))
        previous = orders


def check_last_step(records, ladder, bins, tolerance):
    """List the BAO errors changing by more than the tolerance in the last step.

    Parameters
    ----------
    records : dict
        Ladder records.
    ladder : list of tuple
        Orders in increasing refinement.
    bins : list of int
        Bins checked.
    tolerance : float
        Maximum relative change.

    Returns
    -------
    failures : list of str
        Descriptions of exceeding quantities (empty on success).
    """
    failures = []
    if len(ladder) < 2:
        print("single order: no refinement check")
        return failures
    last, before = records[ladder[-1]], records[ladder[-2]]
    largest = 0.0
    for bin_index in bins:
        for name, sigma in last[bin_index]["sigma"].items():
            for label, new, old in zip(
                ("ap", "at"), sigma, before[bin_index]["sigma"][name]
            ):
                change = abs(relative_change(new, old))
                largest = max(largest, change)
                if not change <= tolerance:
                    failures.append(
                        f"bin {bin_index} {name} sigma_{label}: {change:.2e}"
                    )
    print(
        f"\nlargest last-step BAO-error change: {largest:.3e} (tolerance {tolerance:g})"
    )
    for failure in failures:
        print("EXCEEDS:", failure)
    return failures


if __name__ == "__main__":
    sys.exit(main())
