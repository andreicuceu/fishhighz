# Step 1 — DESI Run-2 input packaging handoff

Status: implementation and bounded checks complete; ready for independent and
user review.  This is not a scientific-acceptance claim and does not authorize
Step 2 or later work.

## Implemented

* Added byte-for-byte copies of the three density tables, 12 DESI-2 QSO SNR
  tables, 12 DESI-2 LBG SNR tables, `Planck18.ini`, and the Vega
  `Planck18_z_2.406.fits` K/PK/PKSB template under `fishhighz/data/`.
* Added `fishhighz.resources.bundled_resource`, `bundled_path`, and
  `resolve_input_path`.  The first two use `importlib.resources`; the latter
  anchors a relative user path at the survey INI directory and never at the
  process CWD.  No native INI parser was added.
* Added exact inventory, units/semantics, source commits, sizes, and SHA-256
  records to `fishhighz/data/README.md`.  Copied the upstream data notes and
  the GPL notices from both read-only source checkouts into
  `fishhighz/data/licenses/`.
* Declared package-data and source-distribution inclusion in `pyproject.toml`
  and `MANIFEST.in`, with focused deterministic tests in
  `tests/test_resources.py`.

## Checks

All 29 scientific assets were compared with their stated upstream paths using
`cmp`: byte-identical.  The focused checks passed:

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest tests/test_resources.py tests/test_import.py -q
5 passed in 0.15s
.venv/bin/ruff check fishhighz/resources.py fishhighz/__init__.py tests/test_resources.py
All checks passed!
```

A local no-isolation build (using the already installed NERSC setuptools because
the isolated build could not resolve PyPI DNS) produced:

```text
/tmp/fishhighz-step1-final-dist.fvJBdM/
  fishhighz-0.1.0.dev0-py3-none-any.whl
  fishhighz-0.1.0.dev0.tar.gz
wheel: 33 data/README/license/input files, byte-identical
sdist: 33 data/README/license/input files
wheel sha256: 58c44915ed032f54cb4f384407e91270de1e74f690a8f5fbaf8a484b7912619f
sdist sha256: 58aef7e754157ac1d9a6ab2920723511b673e4e3c35b65a84ce2b1934d2a0cec
zip-import resource lookup: passed
```

The isolated `python -m build --wheel --sdist` attempt was an infrastructure
failure only: temporary build-environment creation could not resolve
`pypi.org`.  No forecast, Slurm action, commit, push, or planning-document edit
was made.

## Provenance and redistribution limitation

The density-table origins and their limitations are stated by the upstream
`lyaforecast` data README.  The SNR headers provide generation command fragments
and numerical axes, while the FITS header provides Vega table and cosmological
metadata.  Neither source checkout supplies a separate per-data-product
copyright or redistribution statement for those SNR/FITS products.  Exact
GPL-3.0-or-later repository notices and the upstream data README are included,
but that does not establish independent data rights.  Maintainer/legal
confirmation is therefore required before publishing a release containing the
copied SNR and FITS assets; no missing provenance has been invented.

Stop for independent/user review.
