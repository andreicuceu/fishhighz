# Step 2 — native cosmology and bias review

Status: implementation review passed with no further code correction requested.
Progression to Step 3 is supported under the user's standalone-forecast
assignment. This is not scientific acceptance of forecast equivalence.

## Scientific and numerical assessment

The native preparation boundary preserves the intended quantities and keeps
CAMB optional. A single bulk CAMB solve requests the complete unique redshift
set in strictly decreasing-redshift order. The returned CAMB redshift metadata
must match that request before `sigma8` and `fsigma8` are associated with
redshift and reordered to caller order. Deterministic controls reject altered
metadata, unprepared growth redshifts and ambiguous template-redshift aliases.

The accuracy example now uses the explicit native `H(z)` and transverse
comoving `D_M(z)` callables for geometry. Template-growth and damping-reference
sigma8 normalizations remain separate. The analytic density-bias constants,
fixed Ly-alpha beta, galaxy `f/b` definition and linear tabulated extrapolation
match the inspected `lyaforecast` formulae. No Kaiser, damping, AP derivative,
weight, covariance or Fisher implementation changed.

## Independent checks

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q tests/test_import.py \
  tests/test_forecast_imports.py tests/test_cosmology.py tests/test_biases.py \
  tests/test_desi2_accuracy_example.py
15 passed in 1.13s

.venv/bin/ruff check <changed Python files>
All checks passed!

.venv/bin/ruff format --check <changed Python files>
7 files already formatted

git diff --check
passed
```

The three broader-suite failures reported in the implementation handoff were
reproduced independently: two stale synthetic fixtures construct
`prepared.p3d=None` in `test_step12_performance.py`, and the plotting test lacks
the optional Matplotlib dependency. They do not exercise the Step-2 code.

## Remaining evidence boundary

CAMB 2.0.1 successfully completed the native eight-redshift bulk preparation
in 59.4 s using the bundled Planck18 INI. A native-versus-legacy CAMB comparison
was stopped when the combined solves exceeded the bounded login-node runtime;
there is therefore no certified live numerical equivalence result at this
stage. The synthetic CAMB-surface controls validate ordering and extraction,
not CAMB's physical output. Live intermediate-quantity and forecast-level
comparison remains part of Step 5 and may require an explicitly authorized
compute allocation.

No real forecast, Slurm action, commit or push was performed.
