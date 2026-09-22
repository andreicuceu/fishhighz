# Contributing and building documentation

Use Python 3.11 or newer and an isolated environment from the package root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,docs,templates]'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
./scripts/check.sh
```

The runner invokes pytest, Ruff lint, and Ruff formatting checks. It does not
install tools or rewrite formatting. Tests should be self-contained and use
synthetic inputs; real-survey studies and scheduler actions require separate
explicit authorization. Source uses a direct `fishhighz/` layout, NumPy-style
public docstrings, four-space indentation, and an 88-character Ruff line length.

## Documentation

```bash
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
```

Alternatively, run `make -C docs html`. Open `docs/_build/html/index.html`.
The documentation CI and Read the Docs configuration use the same strict HTML
build. Hosting configuration does not itself publish or connect a project.

Edit narrative pages in MyST Markdown. Include existing examples with
`literalinclude` and download links rather than copying their contents. Select
API members explicitly; never execute forecasts during documentation builds.
Add new pages to a toctree and keep units, array ordering, default values, and
scientific qualifications consistent with source and tests.

Historical records belong in the archives. Preserve scientific claims and dated
commands, distinguish historical proposals from current interfaces, and label
local-only evidence rather than creating links to unavailable files. Update the
[migration inventory](migration.md) when relocating a document.

## Distributions and import checks

```bash
python -m build
```

Documentation sources and examples ship in the source distribution. Generated
HTML and research/development archives are not runtime wheel contents. The base
package requires only NumPy; documentation dependencies are optional.

Install the built wheel in an isolated environment and import it from outside
the checkout when checking packaging. A successful import from the checkout
alone does not establish that the wheel is complete. Read package versions with
`importlib.metadata.version("fishhighz")`; `pyproject.toml` is authoritative.

Earlier distribution-validation commands and machine-local helper paths are
preserved in the [development archive](../archive/index.md), not required for
routine contributor checks.
