"""Portable checks of installed-package import behavior."""

import os
import subprocess
import sys


def test_import_is_quiet_and_independent(tmp_path):
    """Import in a fresh process without initializing neighboring packages."""
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-c",
            "import fishhighz; import sys; "
            "assert 'vega' not in sys.modules; "
            "assert 'lyaforecast' not in sys.modules",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
