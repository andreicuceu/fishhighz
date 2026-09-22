# Forest weighting diagnostic W04, revision 1: aliasing in the weighting signal

Status: prepared at the user's request after W03's passing scientific review;
awaiting user approval and separate implementation-agent dispatch. This proposes
an explicitly defined diagnostic variant, not a production weighting choice or
a reconstruction of the historical C++ algorithm. Package Step 13 is unchanged.
Created: 2026-09-15.

## 1. Scientific question and scope

For the full-sample update on the saved W01/W03 QSO inputs, how does including
sampling aliasing in the weighting signal affect the coefficients and their
approach to finite-count stability, compared with the existing fixed intrinsic
signal? Can the observed three-to-six-update drift be reduced simply by more
iterations of the existing full-sample rule?

W03 passes review: replacing prefixes alone does not remove the t=3/6 drift.
A changes by -0.9196% and P_pixel by +2.0171% in its full-sample branch. W03
stopped at six updates, so it does not distinguish a short transient from
persistent drift. W04 compares just two full-sample trajectories on the same
107-node signed inputs; extending both to at most 24 updates is inexpensive and
prevents attributing an ordinary count effect to the aliasing term.

The proposed variant recomputes S_t=P0+A[w_t]B before every update, starting from
the common legacy seed. W02 established that the historical signal-refresh and
initialization convention is unresolved. This explicit-refresh variant is a
conditional experiment; it is not asserted to be the uniquely intended ME07/FR14
algorithm. Approval of this assignment selects only this diagnostic experiment.
An alternative frozen-aliasing signal is not silently substituted or also run.

Only prepared coefficient calculations are assigned. No prefix extension,
Fisher sum, magnitude refinement, new physical inputs or new forecast is needed.

## 2. Read only relevant evidence and live source

Read workspace/package AGENTS.md, design section 0 and main roadmap handover/
progress register, then the diagnostic roadmap (local-only path: `../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md`),
W03 report (local-only path: `reviews/weighting-diagnostics-w03-r1.md`),
W03 review (local-only path: `reviews/weighting-diagnostics-w03-review-r1.md`), and W02's signal
assessment in its review (local-only path: `reviews/weighting-diagnostics-w02-review-r1.md`).
The reviewed W03 instructions (local-only path: `reviews/weighting-diagnostics-w03-instructions-reviewed-r1.md`)
are archived. `IMPLEMENTATION_STEP.md` belongs to package Step 13, not W04.

Inspect `scripts/compare_forest_weight_integrals.py` for initialization, moment
products, full-sample update and scalar independent calculation; follow live
legacy/kernel source only where needed to confirm these operations. Preserve
all dirty/untracked source and reports. Do not reopen historical code/literature
searches or repair W01 evidence tooling without a concrete scientific reason.

The authoritative prepared inputs remain:

```text
.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json
SHA-256: 3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6
```

Read `settings.pair_inputs["lya(qso)_lya(qso)"]`. Check compatibility bin 0,
z=[2.0,2.235], field identity, finite arrays, 107 magnitude nodes from 16.1 to
26.75 and nine negative densities. Retain all signed inputs literally. W03's
checkpoint vectors are available for replay comparison at:

```text
.validation/forest-weight-diagnostics/w03-r1-20260915T224027Z/arrays.npz
.validation/forest-weight-diagnostics/w03-r1-20260915T224027Z/summary.json
```

Direct report loading and the replay below suffice. Do not run W03/W01 controllers
or recreate their artifact-validation machinery. Missing or changed authoritative
inputs require a bounded handoff, not reference recapture.

## 3. Two explicit recurrences, all other quantities fixed

Keep the grid, rectangular dm, signed rho_i, variance v_i, L, l_p, response,
P0=`auxiliary_signal`, B=`auxiliary_p1d`, and seed from W03. The names P0 and S_t
separate the fixed intrinsic power from the signal actually used at each update.
P0 and B already contain the relevant instrumental smoothing. Apply no additional
response, cosmological conversion, density floor or inter-update normalization.

For each branch, independently evolve from

```text
w_i^(0) = (B/l_p)/(B/l_p + v_i)
I1_t = sum_i rho_i w_i^(t) dm
I2_t = sum_i rho_i (w_i^(t))^2 dm
A_t = I2_t/(L I1_t^2)
N_t = (L/l_p) I1_t

fixed intrinsic signal:  S_t = P0
refreshed aliasing:      S_t = P0 + A_t B

w_i^(t+1) = S_t/(S_t + v_i/N_t)
```

Use one common N_t and S_t across all magnitudes within each branch. Compute A_t
from that branch's current w_t, before the simultaneous weight update. Do not
use the other branch's weights, a fixed W03 A, an updated w_(t+1), a prefix A_i,
or P_pixel in S_t. B is evaluated at the saved weighting mode and is fixed.

Preserve W03's moment product order and `cumsum(...)[-1]` reduction in float64.
Final coefficients in both branches are the same functions of their own weights:

```text
I3_t = sum_i rho_i (w_i^(t))^2 v_i dm
P_pixel,t = l_p I3_t/(L I1_t^2).
```

rho_i has units deg^-2 (km/s)^-1 mag^-1 and the integrated moments have units
deg^-2 (km/s)^-1. v_i is dimensionless; L, l_p and B have units km/s, A is
deg^2, and P0, A B and
P_pixel have units deg^2 km/s. Explicitly check that the added aliasing term
has the same units and response convention as P0. The variant changes weight
preparation only; it does not add a second aliasing contribution to final noise.

## 4. Minimal bounded numerical experiment

1. Replay fixed-intrinsic full-sample t=3 and t=6 from the seed. Compare its
   weights and moments with W03 and its coefficients with:

   | t | A [deg^2] | P_pixel [deg^2 km/s] |
   | --- | --- | --- |
   | 3 | 0.022021379679306282 | 0.5298546380249203 |
   | 6 | 0.021818865253617570 | 0.5405425825641871 |

   Require rtol=5e-12, atol=0, finite operands and exact reference-zero handling.
   Stop on a material replay discrepancy. No new prefix trajectory is required.
2. Evaluate the refreshed-aliasing branch at the same counts. Report matched-count
   coefficient shifts, within-branch t=3/6 changes, and A_t B/P0 to show the size
   of the added signal. Save S_t for each executed update, including t=0.
3. Continue both branches to t=12, and to t=24 only if the stopping condition
   below has not been satisfied. Evolve from previous checkpoints rather than
   reseeding them. Both branches have at most 24 updates of 107 weights each.
   W03's sensitivity-at-six early stop is deliberately replaced here: larger
   counts are needed to separate finite-count drift from the signal contrast.
4. After three checkpoints a<b<c, define bounded coefficient stability for a
   branch by requiring both A and P_pixel to differ by at most 1e-3 in all three
   comparisons b/a, c/b and c/a, using signed fractional changes and their
   absolute magnitudes. Test (3,6,12), then (6,12,24) if needed. Stop as soon as
   both branches satisfy this condition, otherwise stop at 24. Classify each
   branch separately and retain all attempted results; a sensitivity finding
   at t=3/6 alone is not an early stop in W04.

A zero reference has no fractional change: report values and the absolute
change, treating exact zero-to-zero as unchanged and zero-to-nonzero as changed.
Negative coefficients, nonpositive effective density or nonpositive S_t must
be reported as leaving the physical diagnostic domain, not as convergence.
Stop on these, undefined moment normalization, zero update denominator, or
nonfinite arithmetic; retain the last finite states. Confirm an obstruction
with the independent scalar arithmetic before attributing it to the recurrence.
Do not floor, regularize or use a different solver to continue past it.

The 1e-3 criterion describes stability over the stated checkpoints only. It
neither proves an asymptotic limit nor accepts a physical signed-density measure.
Report max(abs(w)), I1, signs and output-only w/max(abs(w)) where they help explain
the coefficients; never feed a normalized shape into the next update.

## 5. Independent checks sufficient for this comparison

- One-bin positive analytic control: rho*dm=L=l_p=P0=B=v=1, w0=1/2. A=1 for
  nonzero w, so the aliasing branch has S=2 and leaves w=1/2 fixed. The intrinsic
  branch maps w to w/(1+w). Derive this by hand; confirm the two update functions.
- Two-bin control: rho*dm=(1,1), v=(1,4), L=l_p=P0=B=1, w0=(1/2,1/5).
  Initially A=29/49, S_alias=78/49, and its first full-sample update is
  (39/74,39/179). Derive the next update with exact rational arithmetic and
  recomputed A to check refresh timing, rather than a frozen aliasing term.
  Compare the float implementation with this independent calculation.
- Independently reproduce the executed real-input checkpoint weights, moments,
  signals and coefficients with scalar 80-digit Decimal arithmetic, converting
  the exact saved binary floats. Compute sums and S_t in the independent loop;
  do not call the vectorized implementation's moments or update functions.
  Compare with the replay tolerance. Increase precision once only if a concrete
  mismatch, cancellation or singularity affects the conclusion.

These small controls and the real-input contrast are sufficient. Prefer one
standalone script with assertions; a separate focused test file is optional.
No broad test suite, new wheel/environment, schema, replay framework or general
provenance-hardening work is required. Run only these checks and Ruff on new
Python. Use the existing interpreter, one process, and OMP_NUM_THREADS,
OPENBLAS_NUM_THREADS and MKL_NUM_THREADS set to 1. Cap each numerical invocation
at 30 seconds and retain partial results if reached; no automatic extension.

## 6. Handoff, interpretation and review

Allowed additions: `scripts/compare_forest_weight_signals.py`, optional focused
tests, `reviews/weighting-diagnostics-w04-r1.md`, and small outputs in a new
`.validation/forest-weight-diagnostics/w04-r1-<UTC>/` directory. Preserve earlier
scripts, reports, artifacts and production code. Do not modify planning files
during implementation. Save the input identity, checkpoint vectors, moments,
S_t history, exact attempts/stop, command and independent-check results in
transparent small arrays/tables; no generalized validator is needed.

The report must lead with the scientific answer and its limits. Include:

- The explicit-refresh convention and the historical uncertainty it does not
  resolve; list all held quantities and the single signal change.
- Matched-count coefficient differences, each branch's iteration changes, signal
  contribution and bounded-stability classification, with actual counts shown.
- Whether additional iterations alone reduce the original full-sample drift,
  and whether adding aliasing alters that behavior or mainly the coefficient
  levels. If both stabilize, do not claim aliasing was necessary for stability.
  If neither stabilizes within the cap, do not claim either has no finite limit.
- Signs, any arithmetic/physical-domain obstruction, independent agreement, and
  the distinction between historical W03 results and new t=12/24 measurements.
- A concise remaining scientific question, not a future numbered assignment.

The reviewer checks the same scientific contrast with independent small
arithmetic, verifies the replay and signal-refresh timing, and writes
`reviews/weighting-diagnostics-w04-review-r1.md`. Review only the source/tests
needed for this conclusion. Propose corrections only when likely to change the
scientific conclusion or its justified scope. Begin with **Scientific conclusion
and impact of proposed changes**; for each proposed correction give its possible
scientific consequence and the smallest resolving check. State explicitly when
none are needed. No optional cleanup or hypothetical hardening list.

Stop for user review. No production prescription, aliasing-refresh convention for
production, magnitude/instrument/input-policy change, additional population/bin,
McQuinn–White adoption, Fisher sum, new survey/model evaluation, full forecast,
package Step 13 work, Slurm, agent dispatch, commit or push is authorized.
