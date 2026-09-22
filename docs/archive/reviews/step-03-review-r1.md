# Step 03 revision 1 review

Outcome: review passed; no findings requiring changes. Ready to proceed when the
user requests the next step. Step 04 has not been drafted or authorized.

## Reviewed scope and provenance

Reviewed the current assignment, `reviews/step-03.md`, all six new implementation
files, four new test files, README changes, and dependency metadata. The private
`_arrays.py` helper keeps shared validation outside the numerical loops without
adding a backend or plugin framework. Scientific providers remain ordinary
callables, and NumPy is the only new runtime dependency.

All 12 inline source/test/document hashes in the handoff match the current files.
The final wheel hash matches
`373bce8b93bb7af97213348ca04c2bac14aaae28bada9377dd36e1c6d08427f7`.
All seven Python modules in that wheel are byte-identical to the current source.
Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228`; this step's
implementation is uncommitted. No implementation code or tests were changed by
the reviewer.

## Independent checks

From the FishHighz root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step03-review-r1/verify.py
git diff --check
```

Results:

- **146 tests passed in 23.78 seconds**, including the 52 retained scaffold and
  synthetic reference-tool checks. Ruff lint and format checks passed
  (28 files already formatted).
- Thirty deterministic randomized selections across one to six fields matched
  independently indexed dense field-matrix covariance products. This checks
  required-pair ordering and all four lookup maps without implementing a new
  covariance kernel.
- Nine Jacobian mappings, including empty local dependencies and repeated global
  indices, agreed with independent matrix contractions.
- Quadrature orders 1, 2, 4, and 6 on unequal k bins reproduced the corresponding
  analytic k/mu polynomial moments.
- Source/wheel hashes and all module bytes passed independent comparison.
  The verifier and machine-readable results are retained under
  `.validation/step03-review-r1/`.

The reviewer also reran `.validation/step03/probe.py` using the reported
`.validation/step03/wheel-env/bin/python -I`, from `/tmp`, with the exact final
wheel and README example as arguments. The three thread settings were one.
It confirmed quiet package import without eagerly importing NumPy or neighboring
packages, installed-module bytes matching the wheel, execution of the README
integration example, and custom-grid/P1D behavior. Imports resolved to the
isolated environment's site-packages. Python was 3.13.15 and NumPy 2.5.3.

This reran the existing installed-wheel probe and verified its source/artifact
identity; it did not claim a second fresh installation by the reviewer. The
original handoff records the fresh installation and its retained logs.

## Contract assessment

- Field IDs stay opaque; selection order is preserved while pair identities and
  covariance-required spectra are canonical. Unused fields add no dependencies.
- Parameter sharing is explicit. Duplicate local contributions sum correctly;
  fixed settings remain outside the free vector, and bin-local nuisance IDs
  remain independent.
- Tensor flattening, dk/dmu weights, positive-mu mode normalization, and fixed
  reference h follow the reviewed conventions. Grid preparation owns and protects
  its stored arrays, with no parameter-dependent cuts or weights.
- Callable output validation enforces shapes and real finite numerical types,
  preserves signed powers/derivatives, and keeps P1D independent of P3D.
- Tests cover array ownership, invalid inputs, integration of the contracts,
  and polynomial normalization. No covariance/Fisher engine, scientific model,
  or external-package dependency was introduced prematurely.

## Limits and progression

No real lyaforecast quick capture or full suite was run. Neither is needed for
this contracts-only assignment; full execution remains restricted to explicit
user requests. Validation exercised Python 3.13, not a Python-version matrix.
Polynomial exactness does not establish convergence of eventual survey forecasts,
and this review does not claim a working forecast engine or JIT performance.

The implementation and report satisfy Step 03 revision 1. No repair revision is
requested. Await the user's instruction to prepare Step 04.
