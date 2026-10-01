"""Focused contract tests for the short public DESI-2 accuracy example."""

import ast
import importlib.util
from pathlib import Path

EXAMPLE_PATH = Path(__file__).resolve().parents[1] / "examples/desi2_accuracy.py"


def _example():
    """Import the accuracy example without invoking its main function.

    Returns
    -------
    module : module
        Loaded example module for the mocked public-API test.
    """
    spec = importlib.util.spec_from_file_location(
        "desi2_accuracy_example", EXAMPLE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_accuracy_example_delegates_public_forecast_and_saves(tmp_path, monkeypatch):
    """The example delegates entirely to Forecast and returns its result.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    example = _example()
    calls = []

    class FakeResult:
        def save(self, path):
            """Record the requested output path without writing forecast data.

            Parameters
            ----------
            path : pathlib.Path
                Path of the temporary test artifact to read or write.

            Notes
            -----
            Appends to the enclosing test call log so provider dispatch can be checked.
            """
            calls.append(("save", path))

    result = FakeResult()

    class FakeForecast:
        def __init__(self, source):
            """Initialize the synthetic FakeForecast fixture.

            Parameters
            ----------
            source : pathlib.Path
                Path to the synthetic survey configuration.

            Notes
            -----
            Sets the instance state used by the enclosing test; no scientific calculation is run.
            Appends to the enclosing test call log so provider dispatch can be checked.
            """
            calls.append(("construct", source))

        def run(self):
            """Return the synthetic public forecast result and record the call where required.

            Returns
            -------
            result : FakeResult
                Placeholder forecast result used to test delegation and saving.

            Notes
            -----
            Appends to the enclosing test call log so provider dispatch can be checked.
            """
            calls.append(("run",))
            return result

    monkeypatch.setattr(example, "Forecast", FakeForecast)
    monkeypatch.chdir(tmp_path)

    assert example.main() is result
    assert calls == [
        ("construct", "desi2_accuracy.ini"),
        ("run",),
        ("save", "accuracy-desi2"),
    ]


def test_accuracy_example_has_no_legacy_imports_or_low_level_execution():
    """The shipped example exposes only the native public-API boundary."""
    tree = ast.parse(EXAMPLE_PATH.read_text(), filename=str(EXAMPLE_PATH))
    imported = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.append(node.module or "")
    assert not any(
        name == "lyaforecast"
        or name == "vega"
        or name.startswith("fishhighz.validation")
        for name in imported
    )
    called_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert not called_names.intersection(
        {"prepare_bin", "run_bin", "run_forecast", "prepare_camb"}
    )
