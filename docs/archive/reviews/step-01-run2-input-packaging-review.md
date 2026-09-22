# Step 1 — DESI Run-2 input packaging review

Status: review passed with no consequential correction requested. This review
supports progression to Step 2 under the user's 2026-09-21 standalone-forecast
assignment; it is not scientific acceptance or public-release authorization.

## Findings

The implementation satisfies the bounded Step-1 contract. All 29 scientific
assets are byte-identical to the stated files in the read-only `lyaforecast`
and Vega checkouts. The package inventory distinguishes the three density
tables, 24 population-specific SNR tables, CAMB configuration and signed
K/PK/PKSB template, and records their units, checksums and available source
history without inventing missing generation metadata.

`bundled_resource()` and `bundled_path()` use `importlib.resources`, including
when FishHighz is imported directly from a newly built wheel. Relative external
paths are anchored at the survey-INI directory and are independent of the
process working directory. The data are declared for both wheel and source
distribution inclusion. No cosmology, bias, survey-preparation or forecast
logic was introduced prematurely.

The implementation also states the material release limitation correctly:
repository GPL notices and the upstream density-table notes are preserved, but
the source checkouts do not separately establish redistribution rights for the
SNR and FITS products. Maintainer confirmation remains necessary before a
public distribution containing those assets.

## Independent checks

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest tests/test_resources.py tests/test_import.py -q
5 passed in 0.17s

.venv/bin/ruff check fishhighz/resources.py fishhighz/__init__.py \
  tests/test_resources.py
All checks passed!

git diff --check
passed

cmp against all 29 stated upstream scientific inputs
29 scientific assets byte-identical
```

A fresh Python-3.14 `pip wheel --no-build-isolation --no-deps` build succeeded.
An isolated interpreter imported FishHighz from that wheel as a zip archive,
read `Planck18.ini` as an `importlib.resources` traversable, and materialized
the 28,800-byte FITS template through `bundled_path()`.

The project `.venv` could not independently run `python -m build` because that
environment lacks `setuptools.build_meta`; this is an environment limitation,
not a failure of the package metadata. The implementation handoff records a
successful wheel and source-distribution build using an available local
setuptools backend.

No real forecast, Slurm action, commit or push was performed.
