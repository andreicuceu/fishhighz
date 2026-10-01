#!/usr/bin/env python3
"""Run the bundled DESI-2 accuracy recipe through the public API."""

from fishhighz import Forecast


def main():
    """Run the bundled accuracy recipe and save its results.

    Returns
    -------
    result : SurveyResult
        Native per-bin and combined forecast results.

    Notes
    -----
    Reads desi2_accuracy.ini through the public configuration resolver and
    saves the forecast to accuracy-desi2 in the current working directory.
    """

    result = Forecast("desi2_accuracy.ini").run()
    result.save("accuracy-desi2")
    return result


if __name__ == "__main__":
    main()
