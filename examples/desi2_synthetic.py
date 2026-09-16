"""Seven tiny synthetic selections; this is not a real DESI-2 suite."""

import json

from fishhighz.validation.synthetic import run

if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
