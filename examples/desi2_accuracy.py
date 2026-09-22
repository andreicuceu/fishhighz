#!/usr/bin/env python3
"""Run the bundled DESI-2 accuracy recipe through the public API."""

from fishhighz import Forecast


def main():
    """Prepare and run the native recipe, then save its result directory."""

    result = Forecast("desi2_accuracy.ini").run()
    result.save("accuracy-desi2")
    return result


if __name__ == "__main__":
    main()
