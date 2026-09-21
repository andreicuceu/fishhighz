"""Reproduce the selected DESI Run-2 lyaforecast calculation with FishHighz.

Install FishHighz and lyaforecast in the same environment, then run
``python examples/desi2_full_compatibility.py --reference /path/to/lyaforecast
--output /new/result/directory``.
"""

import argparse

from fishhighz.adapters.desi2_compatibility import run_desi2_compatibility

try:
    from ._desi2_results import create_output, print_results, write_results
except ImportError:  # Direct execution places this examples directory on sys.path.
    from _desi2_results import create_output, print_results, write_results


def run(*, reference, output):
    """Run the six-bin 15x2pt full-compatibility forecast."""
    destination = create_output(output)
    settings, records = run_desi2_compatibility(
        reference=reference,
        profile="full-compatibility",
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
