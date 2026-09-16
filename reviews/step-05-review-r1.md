# Step 05 revision 1 review

Outcome: review passed; no findings requiring changes. Ready to proceed when the
user requests the next step. Step 06 has not been drafted or implemented.

## Reviewed scope and provenance

Reviewed revision 1 of `IMPLEMENTATION_STEP.md`, the handoff, all four new
implementation files, both test files, and the README example. The private
`_information.py` helper shares boundary validation between assembly and results
without changing the accepted Step 04 covariance contract. No implementation
source or tests were changed by the reviewer.

All seven inline handoff hashes match current files. The exact wheel hash matches
`72cdba93af0a2c7d84c4c73fcc393f2137b3be7f4cc5308bcd16718d20e732e8`,
and all 14 Python modules match the reviewed source byte for byte. Independently
checked all 37 entries in the implementer's pre-dispatch snapshot: README.md was
the only changed preexisting file. This comparison preceded the reviewer's
documentation/status edits. Earlier source, tests, reports, and tools remain
preserved. The base commit remains `0d69786`; implementation is uncommitted.

## Independent validation

From the package root:

```bash
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step05-review-r1/verify.py
git diff --check
```

- **261 tests passed in 25.95 seconds**, including all 75 new tests and 186
  retained tests. Ruff lint and formatting passed (42 files already formatted).
- Thirty-six deterministic assembly/result cases covered 1, 2, 5, and 15
  observables, 1, 2, and 5 parameters, and three independent nodes. Constructing
  `J=L B` and `C=L L.T` gave the independent expected Fisher `sum(B.T B)`.
  Signed observable rescalings spanning 120 orders of magnitude, factor reuse,
  and read-only inputs preserved the answers.
- Those cases also checked data combination with a correlated prior, the full
  covariance identity, ordered marginalized subsets against an independent
  Schur complement, explicit fixing against principal information submatrices,
  and general non-orthogonal linear parameter transformations.
- Eight additional rank-deficient cases across 2, 3, and 6 parameters verified
  numerical rank, null-direction residuals in the documented normalized basis,
  and rejection of marginalized errors. Parameter scales spanned 120 orders of
  magnitude.
- Source/report hashes, wheel/source bytes, and preexisting-file preservation
  passed. The independent verifier and JSON results are retained under
  `.validation/step05-review-r1/`.

The reviewer reran `.validation/step05/probe.py` using the reported
`.validation/step05/wheel-env/bin/python -I` from `/tmp`, passing the exact wheel,
saved README example, and source root, with all three thread limits set to one.
It passed quiet import, all 14 source/wheel/installed module comparisons, the
README example, factor reuse, named uncertainty checks, singular-result errors,
and prior-once combination. Module origins were within the isolated environment's
site-packages. Python was 3.13.15 and NumPy 2.5.3.

This reran the existing installed-wheel probe and verified its artifact identity;
it does not claim a second fresh installation. The implementation handoff records
the original fresh installation and retained logs.

## Contract assessment

- Cholesky factors and forward substitution implement the fixed-covariance mean
  Fisher without extra mode weights or covariance-derivative information.
  Factors are reusable across supplied Jacobians. Numerical functions remain
  separate from metadata, with node-bounded workspace and no covariance inverse.
- Covariance factorization judges the selected covariance, correctly allowing
  invertible selections from singular parent field power. Normalized rank and
  symmetry checks diagnose unresolved systems without jitter or projection.
- Named results preserve data and prior information separately, enforce the
  common registry, combine before marginalizing, and apply a shared prior once.
  Conditional, marginalized, and explicitly fixed uncertainties have distinct,
  correctly tested meanings.
- Singular information remains inspectable. Null vectors have an explicit
  normalized coordinate convention and conversion scales. Joint marginalized
  errors fail with diagnostics instead of returning pseudoinverse errors.
- Input ownership, finite arithmetic checks, and the README's supplied-array
  scope agree with the report. No derivative engine, survey work, or compilation
  was introduced prematurely.

## Limits and progression

No real lyaforecast quick capture or full suite ran; neither is required for this
assignment. Full execution still requires an explicit user request. Validation
uses synthetic arrays on Python 3.13, without a Python-version matrix or claims
of survey convergence, external-model derivative accuracy, or JIT performance.
Directly supplied factors receive structural checks; callers own their provenance
and conditioning. Use `factor_covariance` for the normalized solvability check.
Partially identifiable subsets of a singular joint Fisher remain unsupported as
specified; explicit fixing or a finite prior is required before marginalization.

The implementation and handoff satisfy Step 05 revision 1. No repair revision is
requested. Await the user's instruction to prepare Step 06.
