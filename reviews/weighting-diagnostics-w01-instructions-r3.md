# Forest weighting diagnostic W01, revision 3

Status: revision 2 requires corrections after independent review. Revision 3
is proposed for user approval and separate implementation-agent dispatch.
W01-R1/R2/R3 and the numerical result pass independent review, but non-finite
JSON handling and the execution/check inventory require W01-R4/R5 repairs. This
remains W01; W02 and production changes are not authorized. Preserve the exact
revision-1 instructions at
`.validation/forest-weight-diagnostics/w01-r1-20260915T203132Z/WEIGHTING_DIAGNOSTIC_STEP.md`,
the [revision-1 handoff](reviews/weighting-diagnostics-w01-r1.md), its evidence,
and the [revision-1 review](reviews/weighting-diagnostics-w01-review-r1.md).
Also preserve the exact revision-2 instructions at
`.validation/forest-weight-diagnostics/w01-r2-20260915T212017Z/WEIGHTING_DIAGNOSTIC_STEP.md`,
the [revision-2 handoff](reviews/weighting-diagnostics-w01-r2.md), its evidence,
and the [revision-2 review](reviews/weighting-diagnostics-w01-review-r2.md).

Created: 2026-09-15. Last revised: 2026-09-15. This is a diagnostic assignment,
not package Step 14 or a revision of package Step 13.

## Revision-3 required corrections

The revision-2 source/provenance, recurrence, Decimal, projection and
failed-baseline binding remains required and already passes review for finite
values. Do not repeat counts 12/24 or change scientific inputs. Repair and
recertify only final evidence parsing and execution/check provenance.

### W01-R4: fail closed on non-finite JSON numbers

1. Reject `NaN`, `Infinity` and `-Infinity` while reading every JSON artifact,
   including summary, manifest and check inventory. JSON unavailable quantities
   remain `null`; do not replace them with non-finite numbers.
2. Require every actual and expected JSON numerical scalar to be finite before
   applying zero-safe relative comparison. Preserve intentional NPZ NaNs only
   where the reconstructed array has an unavailable error at the same position.
3. Add full `--check-only` mutations with refreshed dependent hashes for a
   non-finite coefficient-trigger operand, forecast comparison, finite table
   value and projection check. Add non-finite created-time/check timing or other
   applicable metadata controls for every JSON reader. Each must reject before
   it can support a scientific verdict.

### W01-R5: validate chronological execution and required checks

1. Define the exact pre-finalization check inventory: focused pytest, Ruff lint,
   Ruff format check and the bounded real-source numerical run. Require each
   category exactly once, with its expected command scope, zero exit code,
   nonempty result and finite nonnegative wall/numerical timing where recorded.
   Reject missing, duplicate, unknown or fabricated check rows.
2. Validate the summary execution mapping rather than excluding it from semantic
   reconstruction. Cross-check its numerical command against the corresponding
   check row and required input/output/instruction paths; require the recorded
   interpreter and version fields to be nonempty and all three thread settings
   to equal `1`. Require a finite nonnegative numerical duration.
3. Do not place claims that finalization or final `--check-only` already passed
   inside the checks file that finalization takes as an input and hashes. Seal
   only operations completed before the manifest. Run finalization next, then
   run authenticated `--check-only`; record these later observations in a
   clearly ordered external handoff/review record, or use an outer record that
   does not claim to authenticate itself.
4. Add full-route controls replacing all checks by one fabricated pass, omitting
   and duplicating each required category, altering execution metadata, changing
   thread limits, and inserting self-referential post-finalization rows. Each
   must reject after dependent hashes are refreshed.

Implement only W01-R4/R5 in the existing standalone script and focused tests.
Write a revision-3 handoff and unique evidence bundle, run the focused tests,
Ruff and the same bounded one-thread real-source replay, then stop for review.
Do not edit package source, change the supported numerical result, execute new
scientific counts or start another diagnostic step.

## Revision-2 corrections retained as reviewed baseline

The original question, fixed scientific inputs, stopping rule and sections
below remain applicable. W01-R1/R2/R3 pass revision-2 review and remain
regression requirements. Do not recompute counts 12 or 24: the independently
verified t=3 to t=6 coefficient trigger retains the revision-1 early stop.

### W01-R1: bind the complete evidence and conclusion

1. Derive authoritative identity, provenance and scientific inputs from
   `load_source`; do not treat paths, hashes, columns or labels in the submitted
   summary as their own oracle. Compare the recorded source manifest, report,
   arrays, producer metadata/arrays, live reference sources/resources and source
   fingerprint with the resolved authoritative values.
2. From authenticated source arrays, reconstruct every completed recurrence
   state and all coefficient changes relative to t=3 and the preceding attempted
   count. Reconstruct the first trigger or failure, later stop statuses, terminal
   classification, verdict list, signs, table and projected-count selection.
   The terminal reason must equal the first actual terminal event; failure,
   coefficient sensitivity and finite-cap no-sensitivity are not interchangeable.
3. Recompute and compare both 80- and 160-digit Decimal states and their
   float64/refinement discrepancies for the baseline and first comparison. A
   string copied from the report is not validated evidence.
4. Reconstruct P1D, response, projected noise, total, covariance, Fisher matrix,
   rank/constrained mask, errors, scalar-oracle discrepancies, changes and
   `forecast_sensitive` verdict. Bind every named NPZ projection array. Require
   t=3 to match both authenticated `reference_pair_fisher` and `pair_fisher`
   within `5e-12` independently.
5. Preserve an immutable copy of the exact diagnostic script used for the final
   calculation and authenticate it, the instruction snapshot, arrays, summary,
   final check inventory and revised handoff without a circular hash. The full
   `--check-only` route must reject altered provenance or artifact hashes.
6. Add full-route corruption tests that independently alter: terminal/verdict;
   trigger operands; table coefficients/errors; forecast comparison; projected
   errors/totals/covariances/P1D/response; Decimal coefficients; script identity;
   and source report path/hash. Each must reject. Retain the revision-1 tests for
   altered source scalar/density, population, pair column, count/trajectory and
   omitted attempts.

### W01-R2: enforce exact-zero scalar baseline checks

Use a zero-safe scalar comparison for captured A and P_pixel. A zero reference
matches only exact zero; otherwise apply the established per-quantity `5e-12`
relative tolerance. Add separate zero/nonzero mismatch tests in both directions
for A and P_pixel, plus exact-zero matches. Preserve the existing array-weight
zero handling.

### W01-R3: retain bounded evidence for failed baselines

Initialization failure, failure at the three-update evaluation, Decimal/float64
disagreement, and a finite but mismatched captured baseline must produce a
serializable terminal result and a new exclusive evidence bundle. Retain all
completed attempts and the classified failure/source discrepancy; mark later
counts not attempted, perform no Fisher projection, and do not emit
`reproduced_legacy`. Add full CLI write-and-`--check-only` tests for each route,
including the Decimal distinction between float64 arithmetic failure and an
algebraic signed singularity. A successful real-source rerun must still preserve
the t=6 early stop and the independently verified revision-1 numbers.

Revision 2 implemented these repairs in the standalone script and focused test
file. Its source/provenance, finite-value semantic, exact-zero and failed-baseline
controls remain mandatory regressions for revision 3. Do not edit package source
or execute another diagnostic step.

## Question and stopping point

Does increasing the update count change the original lyaforecast forest-noise
coefficients or a single-spectrum BAO forecast, with every other legacy input
and the original magnitude grid held fixed?

Use one population in one saved bin: `lya(qso)` auto-power, bin 0,
`2.000 <= z < 2.235`, from the saved `lya_qso_lbg_lae_15x2pt` case. The case name
identifies the source archive; it does not authorize a 15-spectrum calculation.
Do not evaluate the LBG population, other bins or cross spectra in W01.

The first independently verified iteration sensitivity answers whether the
effect predates the profile changes for this example. Stop there for review.
It does not identify a universal cause, invalidate a deliberately specified
three-update estimator, establish optimality, or resolve the continuum limit.

## Read order and relation to the main assignment

Read workspace/package AGENTS.md, design section 0, the main roadmap handover,
and [the diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md).
Read package `IMPLEMENTATION_STEP.md` only for context; do not execute or edit it.
Then inspect:

- `reviews/step-12-r5.md`, `reviews/step-12-review-r5.md`, and
  `reviews/step-12-review-r6.md` for the scientific diagnosis and evidence repair;
- `reviews/step-13.md`, `reviews/step-13-review-r1.md` and
  `reviews/step-13-weight-limit-r2.md` for limits on existing diagnostic claims;
- `../lyaforecast/lyaforecast/weights.py`, `covariance.py`, `power_spectrum.py`
  and `fisher.py`, read-only;
- `fishhighz/validation/reference_capture.py`, `compatibility.py`,
  `numerics.py`, `schema.py`, and the literal replay in
  `scripts/attribute_desi2.py:negative_floor_noise`;
- `fishhighz/kernels/weights.py` to understand why its positive-support handling
  and arithmetic guards are not a literal signed-legacy replay.

The Step 13 revision-2 diagnostic is implemented but not independently reviewed.
It is not the numerical oracle for W01; its nonnegative-mass assumptions do not
cover the signed legacy inputs below.

## 1. Verify and extract the existing inputs

Start from
`.validation/step12-r5-20260914T191855Z/profiles-checked/`:

- `manifest.json`;
- `records-000.report.json` and `records-000.npz` (compatibility, bin 0).

Resolve these through the manifest and verify task, report and array hashes;
do not trust the record number alone. Verify recorded source revisions/hashes,
including any referenced producer metadata needed to establish provenance.
If live reference code differs, inspect the recorded source snapshot. Missing
evidence does not authorize a new NewForecast run or a different case.

Planning inspected metadata and array shapes only. The selected saved row
`settings.pair_inputs['lya(qso)_lya(qso)']` has 107 magnitude nodes from 16.1 to
26.75, nine negative density interpolants and no zero density entries. Its
reported auxiliary values are S=1.7189455896517238 and B=16.333147249645222;
the saved three-update A and P_pixel are 0.022118395206768268 and
0.5480269062497443. Re-read authoritative bytes and validate these observations;
do not use rounded numbers in this document as calculation inputs.

Extract only this field's magnitudes, signed density per source velocity,
variance, auxiliary S/B, forest length L, pixel width Delta_v, response width,
conversion factors, evaluation/source redshifts, and saved weights/coefficient
arrays. Density in this row is already per source velocity: do not convert or
normalize it again. Set q_i=m[1]-m[0] for every node, including both endpoints.

For the optional scalar-covariance projection below, extract k, mu, mode counts,
the selected auto-power column of `total`, the matching `observed_j` column,
and the corresponding saved `reference_pair_fisher`/`pair_fisher` result.
Resolve field and pair indices by identity in the manifest, separately for the
required powers and selected Jacobians. Never use the joint Fisher matrix as
the one-spectrum reference. There are 5000 saved Fourier nodes; using these
arrays avoids generating a new model or a smaller, unvalidated Fourier grid.

## 2. Reproduce the exact finite prescription before changing its count

Implement a narrowly scoped offline diagnostic, preferably
`scripts/diagnose_legacy_forest_iterations.py`, with pure array functions
separated from CLI/file handling. Add self-contained tests in
`tests/test_legacy_forest_iterations.py`. These are suggested new paths, not
existing commands. No new runtime dependency, packaged API, generic validation
framework or production method is required.

Preserve literal legacy operation order for the baseline:

    b = B / Delta_v
    w(0) = b / (b + v)
    I1_prefix = cumsum(rho * w * dm)
    n_eff = I1_prefix * (L / Delta_v)
    w(next) = S / (S + v / n_eff)
    I1 = cumsum(rho * w * dm)[-1]
    I2 = cumsum(rho * w**2 * dm)[-1]
    I3 = cumsum(rho * w**2 * v * dm)[-1]
    A = I2 / (L * I1**2)
    P_pixel = Delta_v * I3 / (L * I1**2).

Signed densities, signed intermediate sums and any resulting signed weights
remain literal diagnostic inputs. Record them; do not floor, clip, drop cells,
renormalize weights, insert a positive-support mask, or infer a physical
positive-density interpretation. A division singularity or nonphysical result
must be recorded distinctly from iteration sensitivity and floating-point loss.
Do not apply positive-matrix spectral-radius arguments to this signed array.

Independently reproduce the recurrence and moment sums using scalar Decimal
arithmetic, starting from the exact saved float64 values. Use 80 and 160 decimal
digits for the baseline and first comparison requiring confirmation. Separate
the effects of representation/operation order from changes with iteration count.
Retain explicit exception/availability records; local floating-point monitoring
must not modify global settings or silently ignore lost contributions.

Baseline gate: t=3 weights, A and P_pixel must reproduce the captured values,
with each nonzero quantity checked at its own scale (relative tolerance 5e-12;
exact zeros separately). Decimal precision refinement must support the result.
If the baseline does not reproduce, stop with a source/arithmetic finding.
Do not continue to later counts using a silently adjusted input or tolerance.

## 3. Bounded update-count comparison

Record t=0 as initialization, t=3 as the reference, then test t=6, 12, 24 in
that order. Each count starts from the same w(0), or from a verified shared
trajectory; counts are absolute updates, not additional batches of updates.
All other inputs, including the rectangular magnitude masses, remain identical.

At each attempted count save weights, I1/I2/I3, A/P_pixel, signs, nonfinite or
range failures, and signed changes relative to t=3 and the previous checkpoint.
Compare A and P_pixel independently; never use a combined norm that can hide
one coefficient. Record absolute-weight changes for context only. Common
amplitude decay is not a coefficient-convergence test.

Use 1e-3 relative change in either coefficient as a diagnostic trigger for a
controlled forecast projection, not as a new package acceptance threshold.
Confirm a trigger against the independent Decimal calculation at both
precisions. Then perform section 4 for t=3 and that count only, and stop for
review regardless of whether its forecast effect exceeds 1e-3. Later counts
are `not_attempted_after_identification`, not missing or successful tests.

If an arithmetic failure appears first, use the same bounded Decimal control to
distinguish float64 range failure from an algebraic singularity or cancellation.
Report the distinction and stop. A recoverable high-precision coefficient does
not authorize changing a production arithmetic guard.

If no coefficient trigger occurs through t=24, project only t=3 and t=24 and
report sensitivity over that finite range. Do not certify an asymptotic limit
from a finite plateau. Maximum real scope is one 107-sample trajectory through
24 updates and two single-spectrum 2x2 Fisher sums.

## 4. One-spectrum forecast consequence using saved arrays

Hold the saved observed mean Jacobian J, Fourier nodes, modes, geometry,
instrument response and intrinsic covariance signal fixed. At each Fourier
node n use the legacy expression

    N_n(t) = [A(t) * P1D(k_parallel,n) * W_n**2 + P_pixel(t)] * d_deg**2 / a_v
    T_n(t) = T_n(3) + N_n(t) - N_n(3)
    C_n(t) = 2 * T_n(t)**2 / modes_n
    F_ab(t) = sum_n J_na * J_nb / C_n(t).

Use the captured legacy z, a_v, d_deg, pixel/resolution widths and the unchanged
PD2013 P1D/floor formula. A small analytic P1D/response evaluation on these saved
nodes is allowed; no external model provider or CAMB/background preparation is
needed. Do not subtract the saved mean to infer the covariance signal: the
legacy mean and covariance were evaluated at different redshifts.

This changes covariance between diagnostic runs because the observing weights
change; covariance and weights stay fixed within each parameter derivative.
There is no covariance-derivative Fisher term. No extra factor of modes or two
belongs in the Fisher sum. Do not allocate or contract a 15x15 covariance.

Require positive finite physical total auto-power and a well-resolved 2x2
Fisher matrix before reporting marginalized ap/at uncertainties. Use the
existing rank-aware conventions; no jitter or pseudoinverse error bars. Report
unavailable errors explicitly if those conditions fail.

Validate F(3) against the selected saved single-spectrum reference within 5e-12
relative discrepancy. Independently cross-check direct summation and 2x2
inversion against a separate small scalar/analytic implementation. Report each
ap/at error change and the Fisher change; 1e-3 is the existing forecast diagnostic
budget. Small forecast changes can coexist with significant coefficient changes
when that noise contribution carries little weight.

## 5. Concrete tests and independent review requirements

Run these synthetic checks before the saved-array calculation:

1. One-cell analytic control below, at and above d=LSr/(Delta_v*v)=1;
   A=1/(Lr), P_pixel=Delta_v*v/(Lr) when positive and finite. Verify that
   changing amplitude alone does not generate a noise-change finding.
2. Two/three-cell unequal signed and positive examples: compare literal
   vectorized operation order with independent scalar Decimal updates and
   integrals at 0/3/6 updates. Do not claim a positive-domain theorem for signed
   controls. Include a cancellation/zero-denominator case with an explicit
   singular or arithmetic outcome rather than invented finite coefficients.
3. A synthetic coefficient change that is hidden by small absolute weights;
   ensure the coefficient comparison detects it. Conversely, stable ratios
   with decaying amplitude must not be called nonconvergent noise.
4. One-spectrum analytic Fisher normalization, two parameter directions,
   field/pair permutation, and a null-direction unavailable-error control.
5. Source/attempt binding: altered scalar or density values with old provenance,
   swapped population identity, incorrect pair column, changed count with
   unchanged trajectory, and omitted failed/capped attempt must reject.
   Truthful early-stop results must remain readable without claiming full scope.
6. Exercise the first-discrepancy stop with a stub evaluator: no later count,
   second population or broad forecast is called. CLI input/output paths are
   explicit; the output directory must not pre-exist.

Test only this small diagnostic and directly affected helpers, then Ruff on the
new/changed files. A focused run is sufficient for a new standalone script;
this user-requested small-diagnostic scope does not require another whole-package
suite, installation or wheel build. If package source unexpectedly needs repair,
stop and report it for a revised assignment instead of expanding the edits.
Use the existing `.venv` with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1.
Target <=30 seconds for each diagnostic invocation; retain partial outcomes
if this bound is reached. No Slurm, parallel workers or new environment installs.

The reviewer independently recomputes t=3 and the first reported differing count
(or t=24 if no difference), checks exact inputs and count semantics, validates
the single-spectrum covariance/Fisher normalization, reruns relevant synthetic
tests, and challenges the conclusion. Review must distinguish a valid diagnosis
from scientific acceptance of any weighting prescription.

## 6. Deliverables and permitted edits

Write new evidence exclusively under `.validation/forest-weight-diagnostics/w01-r3-<UTC>/`.
Keep a small JSON summary, named NPZ arrays with allow_pickle=False, input/source
hashes, exact instructions snapshot, actual commands, versions, timings and
attempt outcomes. Distinguish `reproduced_legacy`, `coefficient_sensitive`,
`forecast_sensitive`, `arithmetic_failure`, `singular_signed_recurrence`, and
`no_sensitivity_detected_within_cap`; they are not interchangeable verdicts.
Supply a simple table of t, A, P_pixel, available ap/at errors and discrepancies.
Only add plots if the table fails to explain the finding.

Write `reviews/weighting-diagnostics-w01-r3.md`. Reviewer writes a separate
`reviews/weighting-diagnostics-w01-review-r3.md`. Preserve the revision-1 and
revision-2 handoffs, reviews and evidence byte-for-byte.
Include the exact changed-file list and evidence that production code,
`IMPLEMENTATION_STEP.md`, previous reports and historical input bundles were
preserved. Archive this instruction revision before any later repair revision.

Permitted implementation edits: the existing standalone script, its focused
tests, the new revision-3 handoff, and its unique evidence directory only. No automatic rewrite of
roadmaps, AGENTS.md, the package README, source kernels, old validators or reports.
No NewForecast/controller run, raw-reader interpolation, external P3D call,
full 15x2pt assembly, additional population/bin, commit, push or agent dispatch.

Stop for independent/user review after W01. Recommend the smallest next
comparison if necessary; do not plan or execute it without the user's request.
