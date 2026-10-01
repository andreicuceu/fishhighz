"""Run DESI Run-2 with the accepted full-sample compatibility weights.

Install FishHighz and lyaforecast in the same environment, then run
``python examples/desi2_fixed_compatibility.py --reference /path/to/lyaforecast
--output /new/result/directory``.
"""

import argparse

from fishhighz.adapters.desi2_compatibility import run_desi2_compatibility

try:
    from ._desi2_results import create_output, print_results, write_results
except ImportError:  # Direct execution places this examples directory on sys.path.
    from _desi2_results import create_output, print_results, write_results


def run(*, reference, output):
    """Run the six-bin DESI-2 fixed-compatibility forecast.

    Parameters
    ----------
    reference : str or pathlib.Path
        Installed lyaforecast checkout containing the reference inputs.
    output : str or pathlib.Path
        New result directory; existing directories are rejected.

    Returns
    -------
    settings : dict
        Forecast prescription, parameter order, and reference identities.
    records : list of dict
        Individual and joint constraints by redshift bin; AP errors and
        correlations are dimensionless.

    Raises
    ------
    FileExistsError
        If the output directory already exists.

    Notes
    -----
    Creates the output directory, runs the forecast, saves settings.json and
    results.npz, and prints the individual and joint constraints.
    """
    destination = create_output(output)
    settings, records = run_desi2_compatibility(
        reference=reference,
        profile="fixed-compatibility",
        work_directory=destination / "lyaforecast",
    )
    write_results(destination, settings=settings, records=records)
    print_results(records)
    return settings, records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference", required=True, help="installed lyaforecast checkout"
    )
    parser.add_argument("--output", required=True, help="new result directory")
    arguments = parser.parse_args()
    run(reference=arguments.reference, output=arguments.output)
