# Contributing and building documentation

Use Python 3.11 or newer and an isolated environment from the package root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,docs,templates,camb,cosmology]'
./scripts/check.sh
```

`scripts/check.sh` sets single-threaded BLAS/OpenMP and runs pytest, Ruff lint
and the Ruff format check; it neither installs tools nor rewrites files. Tests
are self-contained and use synthetic inputs. Source uses a direct `fishhighz/`
layout, NumPy-style docstrings, four-space indentation and an 88-character line
length. Numba kernels under `fishhighz/kernels/` define the reference
arithmetic of the forecast; changes to them need an accompanying comparison of
saved results (see `AGENTS.md`).

The environment variables `FISHHIGHZ_CAMB_CACHE` (set to `0` to disable the CAMB
background cache) and `FISHHIGHZ_CACHE_DIR` (cache location) control caching.

## Documentation

```bash
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
python scripts/check_docs.py docs/_build/html
```

`make -C docs html` runs both. Open `docs/_build/html/index.html`.

Narrative pages are MyST Markdown. Include existing examples with
`literalinclude` and download links instead of copying them. Select API members
explicitly and wrap autodoc and autosummary directives in `eval-rst` blocks;
`scripts/check_docs.py` requires a rendered signature and anchor for every
documented object and member. Never execute forecasts during a documentation
build. Add new pages to a toctree and keep units, array ordering and defaults
consistent with the source.

## Distributions

```bash
python -m build
```

Documentation sources and examples ship in the source distribution; generated
HTML does not. Verify packaging by installing the wheel in a clean environment
and importing it from outside the checkout. `pyproject.toml` is authoritative
for dependencies and extras.
