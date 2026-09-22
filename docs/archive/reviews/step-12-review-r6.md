# Step 12 revision 6 independent review

**The bounded revision-6 evidence repair passes independent review. R2 is
resolved for the reviewed contract. No further code correction was found.
R3 remains unresolved; Step 12 is not scientifically accepted.**

Reviewed the [handoff](step-12.md), live source/tests, exact installed wheel and
saved results against the [revision-6 instructions](../notes/IMPLEMENTATION_STEP.md).
The exact reviewed instructions, handoff, initial source hashes and governance
are snapshotted in .validation/step12-review-r6 (local-only path: `../.validation/step12-review-r6/`).
No real forecast, reference recapture, new installation, scientific prescription
change, production edit, agent dispatch, Slurm action, commit or push occurred
during review. No subsequent revision or Step 13 is planned.

## R2 correction and independent checks

The changes are confined to validation numerics, schema and trial binding, the
README, handoff/archive and new tests. Pair Fisher matrices and derived summaries
now have individual consistency checks. Primary/final trials, metric operands
and replayed operands use each matrix's own scale. Scaling before subtraction
and norm calculation removes the earlier false-zero comparison for weak
information. The controller replay reconstructs metrics and the convergence
verdict from saved trials and requires agreement even under `require_pass=False`.
Failed attempts remain explicit; truthful unconverged reports remain inspectable.

The original r5 independent reproducer was rerun with its full synthetic
five-field, 15-spectrum shape against the installed r6 wheel. Its true
individual-spectrum error change is **0.01941932430908011 (1.9419%)**. Copying
the weak spectrum's upper operand into the lower operand and setting the claimed
change to zero still leaves every actual study array unchanged. It now rejects
through all three routes: semantic validator, writer, and consistently rehashed
offline bundle. All reject at:

```
metric/trial operand weights 0 metric_pair_fisher block (0,): inconsistent numerical content
```

original-probe.json (local-only path: `../.validation/step12-review-r6/original-probe.json`) records
the result. This is a synthetic per-record corruption test with fixture provenance,
not a real forecast. The offline test promotes only outer completion flags so
that rejection reaches numerical semantics rather than stopping at an incomplete
bundle. The complete synthetic 12-primary/72-diagnostic inventory test is a
separate regression.

The new tests cover information ratios 1, 1e-8 and 1e-16; weak auto/cross spectra
at positions 0, 1, 8 and 14; isolated, combined and complete substitutions;
primary/final/replay bindings; pair summaries; inconsistent metrics/verdicts;
and valid weak, roundoff, null and partially constrained controls. Historical R1,
r2 repeated-operand, r5 inventory and performance regressions remain present.
The before/after table correctly distinguishes old-wheel test failures from
actual falsely certified forecasts.

A separate Decimal oracle (local-only path: `../.validation/step12-review-r6/range_review.py`)
checks 21 matrix comparisons at scales from **1e-320 to 1e308**, including harmless
roundoff and 4% discrepancies. Four mixed-scale stack mutations verify that
neither another trial, side nor spectrum can conceal the changed block. Exact-zero
and nonfinite controls pass. This establishes range behavior for these cases;
it is not a claim of arbitrary-range covariance inversion.

The compatibility comparison's change from exact dictionary equality to the
existing metric-consistency rule is appropriate: it admits roundoff from the
range-safe norm while rejecting a false metric. Evidence consistency remains
5e-12; scientific thresholds remain 1e-3 for Fisher/error changes and 1e-6 for
volume. No physical weighting, general Fisher/covariance, model or response
kernel changed. The reviewed performance implementation and NumPy default remain.

## Fresh validation and source identity

- Ordinary suite: **1243 passed, 25 optional-compiler skips, 181.52 s**; Ruff lint
  and format passed. Command: `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`;
  the script sets OMP/OPENBLAS/MKL thread counts to one.
- Installed affected regressions: **251 passed, 155.38 s**
  (`installed-tests.log`). Executed `-I -m pytest -q tests/test_step12_revision6.py
  tests/test_step12_revision5.py tests/test_step12_revision2.py
  tests/test_step12_performance.py` with one-thread settings from a fresh copy of
  fixtures at `/tmp/fishhighz-step12-review-r6-tests`. The performance tests select
  both NumPy and Numba; the process default is NumPy.
- Independent saved-array Wick covariance, direct C solves and analytic 2x2
  inverses reproduce **12 joint forecasts and 180 individual-spectrum results**.
  Maximum discrepancies: C 0; joint F 1.34e-15; per-spectrum F 2.20e-15;
  joint errors 1.33e-15; per-spectrum errors 2.11e-15. All six accuracy reports
  retain correct individual trial operands and truthful failed verdicts.
- All six examples pass outside the checkout with isolated imports. All **52**
  source/wheel/installed modules and installed METADATA/WHEEL bytes match the
  final wheel, SHA256
  `7286a2c34014807c2d3309f9b9f140594aeb64fd8fe415932ba234cae399697c`.
- All **155** implementation snapshot entries match their preserved copies.
  All **35** non-validation modules are unchanged; all **144** saved r5
  primary/diagnostic array and report hashes match. The r5 handoff archive and
  original r3 instructions are unchanged.

New commands, scripts, outputs and identities are saved separately from the
implementation evidence. The existing r6 environment uses NumPy 2.5.3, SciPy
1.18.1, Numba 0.67.0, llvmlite 0.49.0 and Astropy 8.0.1. No fresh installation
was needed. The identity script's first attempt omitted the documented README
change from its expected changed-file list; it was corrected after inspection,
and both logs remain. Nonfatal site MUNGE messages in the ordinary log do not
change the passing test result.

The handoff's 251 explicitly compiled and 435 NumPy installed checks remain
implementation evidence, distinct from this review's installed run. Its disclosed
explicit-Numba/SciPy-blocked subprocess failure is inherited from r4 and outside
this repair; that combination was not rerun or counted as a pass here. The
NumPy-only matrix performance target remains unmet. These limitations are
preserved, not newly introduced or waived.

## Scientific status and user control

The unchanged saved results retain six passing compatibility bins, six
unconverged accuracy bins, 60 completed diagnostics and 12 unavailable diagnostics
out of 72 requests. The scoped scientific gate still exits 1 with
`partial or scientifically failed validation`, as expected. This does not fail
the bounded evidence repair and does not establish scientific acceptance.

The fixed-weight magnitude integral and cumulative-weight diagnosis remain the
historical r5 results linked in the [handoff](step-12.md). No figure or scientific
comparison was regenerated. The existence and quadrature independence of a
normalized asymptotic noise limit remain unestablished. The three alternatives
remain unselected: a specified finite-iteration convention with reviewed
acceptance criteria; investigation of a normalized limit of the same rule;
or an independently justified estimator. R2 closure selects none of them.

The active step status, design and roadmap record this bounded review pass.
Scientific weighting decisions, acceptance, any further real run and progression
remain with the user. Wait for the user's next request.
