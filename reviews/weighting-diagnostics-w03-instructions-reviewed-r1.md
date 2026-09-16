# Forest weighting diagnostic W03, revision 1: prefix versus full-sample update

Status: revision 1 is implemented and passes independent scientific review;
awaiting user review and acceptance. No scientifically consequential changes
are needed. The exact reviewed assignment is archived in
[W03 r1 instructions](reviews/weighting-diagnostics-w03-instructions-r1.md);
see [the independent review](reviews/weighting-diagnostics-w03-review-r1.md).
No progression or production prescription is authorized. W01 evidence-only
findings remain deferred. Package Step 13 is unchanged.
Created: 2026-09-15.

## 1. Scientific question and controlled comparison

On the exact saved W01 QSO-forest example, does replacing each magnitude-prefix
I1 by a common full-sample I1 change the prepared noise coefficients and their
sensitivity to iteration count, with every other input held fixed?

W02 established two distinct differences from the published discussion: prefixes
enter the live update, and its signal contains no aliasing contribution. W03
isolates the first difference. Retain the existing fixed, smoothed intrinsic
P3D signal S. This is a controlled equation substitution, not a complete
implementation of ME07/FR14 or a decision about their aliasing-signal convention.

Evaluate coefficients only, for one population and one saved bin. No physical
input preparation, magnitude refinement, Fisher calculation or new forecast is
needed. A changed finite result, persistent sensitivity, bounded stability, or a
confirmed arithmetic obstruction can each answer this bounded assignment.

## 2. Context, source and minimal input identification

Read workspace/package AGENTS.md, design section 0 and main roadmap handover/
progress register, then the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W02 report](reviews/weighting-diagnostics-w02-r1.md),
[W02 review](reviews/weighting-diagnostics-w02-review-r1.md), and the scientific
result in [W01's review](reviews/weighting-diagnostics-w01-review-r2.md).
The [reviewed W02 instructions](reviews/weighting-diagnostics-w02-instructions-reviewed-r1.md)
are archived. `IMPLEMENTATION_STEP.md` remains package Step 13 context only.

Inspect the relevant functions in `../lyaforecast/lyaforecast/weights.py` and
`scripts/diagnose_legacy_forest_iterations.py`: initialization, moment products,
prefix density, update, and coefficient extraction. Inspect only dependencies
needed for the comparison. Preserve dirty/untracked code and historical evidence.
Do not reopen W02's literature or C++ search without a specific contradiction.

The authoritative prepared inputs are already in:

```text
.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json
```

This report has SHA-256
`3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6`.
Use `settings.pair_inputs["lya(qso)_lya(qso)"]`; verify the record context is
compatibility, bin 0, z=[2.0,2.235], and the selected population is lya(qso).
Do not read a different forest's row because it has matching shapes.
The historical W01 result is in:

```text
.validation/forest-weight-diagnostics/w01-r2-20260915T212017Z/summary.json
.validation/forest-weight-diagnostics/w01-r2-20260915T212017Z/arrays.npz
```

Directly read the required report fields or reuse a suitable read-only loader.
Do not invoke W01's diagnostic/forecast controller or build another evidence
validator. The report hash, identity, finite input checks and baseline replay
below are sufficient for this question. If the report is missing or its identity
has changed, report the specific obstacle; do not recapture survey inputs.

## 3. Held quantities and equations

Keep all 107 ordered magnitude nodes (16.1 to 26.75), endpoint-inclusive
rectangular dm, signed density rho_i, pixel variance v_i, forest length L,
pixel width l_p, auxiliary S and B, and initialization exactly as in W01.
Verify the nine negative-density nodes; retain them literally. Do not floor,
mask, alter support, normalize weights between updates, or pass this example
through FishHighz's nonnegative-only weighting API.

Use the saved `_pix_kms`, `forest_length`, `auxiliary_signal`, `auxiliary_p1d`,
`magnitudes`, `density`, `variance`, `_w_lya`, `_aliasing_weights`, and
`_effective_noise_power`. Record the actual scalar values and units used.
Initialization and alternatives are:

```text
w_i^(0) = (B/l_p) / (B/l_p + v_i)
I1,i^(t) = sum_(j<=i) rho_j w_j^(t) dm

prefix: N_i^(t) = (L/l_p) I1,i^(t)
full:   N_i^(t) = (L/l_p) I1,last^(t) for every i

w_i^(t+1) = S / (S + v_i/N_i^(t))
```

Retain W01's floating-point product order and cumulative summation in both
branches. The only changed operation is selecting/broadcasting the last I1
rather than using each prefix. In particular, do not simultaneously replace
`cumsum(...)[-1]` with a differently ordered reduction. Each branch evolves its
own weights from the same seed; do not borrow the other branch's later I1.

After each evaluated checkpoint use full-sample moments in both branches:

```text
I1 = sum_i rho_i w_i dm
I2 = sum_i rho_i w_i^2 dm
I3 = sum_i rho_i v_i w_i^2 dm
A = I2/(L I1^2)
P_pixel = l_p I3/(L I1^2)
```

Density and moments have units deg^-2 (km/s)^-1; v is dimensionless, L/l_p is
dimensionless, B is km/s, S and P_pixel are deg^2 km/s, and A is deg^2.
Keep the same response already represented in S, B and v. No aliasing is added
to S; the final aliasing coefficient A remains part of the comparison.

## 4. Bounded execution and stopping rules

1. **Replay control.** Reproduce prefix t=3 and t=6, including the t=3 saved
   weight vector. Check the historical coefficients independently of any W01
   verdict label:

   | t | A [deg^2] | P_pixel [deg^2 km/s] |
   | --- | --- | --- |
   | 3 | 0.022118395206768268 | 0.5480269062497443 |
   | 6 | 0.022308973469159228 | 0.5544833486755246 |

   Use per-quantity rtol=5e-12 with exact handling of reference zeros; verify
   finite values before comparison. Stop on a material baseline discrepancy.
2. **Single substitution.** Evaluate the full-sample branch at t=3 and t=6.
   Save the weight vectors and moments. Report full/prefix differences at each
   common count, and the signed fractional t=3 to t=6 change within each branch,
   separately for A and P_pixel. This distinguishes a change in coefficient
   level from a change in iteration sensitivity.
3. **Conditional confirmation only.** If either full-sample coefficient changes
   by more than 1e-3 between t=3 and t=6, confirm the result independently and
   stop: the substitution has not removed that finite-count sensitivity.
   This does not establish its asymptotic behavior. If both changes are at most
   1e-3, evaluate only the full-sample branch at t=12. Compare t=12 against both
   t=6 and t=3; stop with either bounded stability over those checkpoints or
   renewed sensitivity. Do not extend either branch beyond these bounds.
4. Stop an affected branch on nonfinite arithmetic, zero update denominators,
   or undefined coefficient normalization. Preserve the last finite state and
   distinguish mathematical singularity from float64 failure using the small
   independent check below. Do not regularize, floor or substitute a solver.

The 1e-3 threshold is the existing diagnostic sensitivity scale, not a new
production acceptance tolerance. For a zero reference coefficient, report both
values and the absolute difference; do not invent a fractional change. Exact
zero-to-zero is unchanged; zero-to-nonzero is a change requiring interpretation.
A negative coefficient or singular signed measure must be reported explicitly
and must not be certified as physically converged even if finite changes are small.

Report max(abs(w)) and the relative shape w/max(abs(w)) at each checkpoint,
with signs retained and all-zero weights treated as undefined shape. This is
output-only normalization. Inspect the evolution of I1 and the update
denominators if needed to explain amplitude decay or signed cancellation; no
asymptotic or continuum derivation is required. The one-bin W02 result already
shows why a common integral does not guarantee nonzero limiting weights.

## 5. Minimal independent checks

- Retain W02's two-bin analytic control: r=(1,1), v=(1,4), L=l_p=S=B=1,
  w0=(1/2,1/5). The first prefix update is (1/3,7/47), and the full update is
  (7/17,7/47). Check the common last-bin first update; later equality is not
  expected. One focused assertion is sufficient.
- Independently reconstruct the decisive real-input coefficients with a short
  scalar Decimal calculation at 80 digits, using the exact saved float inputs
  and the specified recurrence, not the vectorized implementation's update
  function. Check prefix t=3/6 and full t=3/6; include full t=12 only if reached.
  Compare finite weights and each coefficient at 5e-12 relative tolerance,
  treating exact zeros explicitly. Increase precision once only if a mismatch
  or near-singularity affects the conclusion; otherwise no precision ladder.
- Use a direct loop/full sum in the independent calculation. Its independent
  arithmetic is intended to detect a wrong prefix, shared state, missing dm,
  or cancellation error; it must not simply echo the vectorized implementation.
- Verify both branches read identical input arrays/scalars and share only their
  initialization, not their evolving states. A test of artifact formats is not
  a substitute for checking the scientific contrast.

One small standalone script, with assertions and optionally one focused test
file, is enough. No reusable framework, schema, replay engine, broad mutation
suite, full pytest suite, wheel or fresh environment is required. Run the small
checks and Ruff on new Python files only. Use the existing interpreter, one
process and OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1. Target under
30 seconds per numerical invocation; stop with partial evidence at that cap.

## 6. Handoff and independent review

Add only a standalone `scripts/compare_forest_weight_integrals.py`, optional
focused tests, `reviews/weighting-diagnostics-w03-r1.md`, and small outputs in a
new `.validation/forest-weight-diagnostics/w03-r1-<UTC>/` directory. Save inputs,
checkpoint weights and coefficients in a small NPZ/CSV or similarly transparent
form. Record source/report identity, the instruction revision, actual commands,
tolerances, attempted counts and early stop. Preserve W01/W02 scripts and reports.
Do not modify production modules or planning files during implementation.

The report must lead with the answer and its scope, give the two branches'
coefficient levels and within-branch changes, and explain which quantities were
held fixed. State whether the prefix substitution removes, reduces, increases
or retains the measured sensitivity. Keep that finite result separate from any
claim of convergence, physical admissibility or fidelity to the full paper
prescription. Explain weight amplitude/shape only where useful to that conclusion.
State the smallest remaining scientific question; do not plan another step.

The reviewer independently checks the small recurrence and decisive numerical
contrast against the identified inputs, and writes
`reviews/weighting-diagnostics-w03-review-r1.md`. Reuse W02's established source
assessment unless new evidence contradicts it. Propose repairs only if likely to
affect the scientific conclusion or its justified scope; begin with a short
**Scientific conclusion and impact of proposed changes** summary. For each
proposed correction, explain its possible scientific consequence and the
smallest check that resolves it. No unrelated hardening or optional cleanup list.

Stop after handoff/review. User approval, dispatch, acceptance and progression
remain separate. No aliasing-signal change, magnitude refinement, input-policy
change, McQuinn–White adoption, other population/bin, Fisher calculation, new
survey run, full 15x2pt forecast, package Step 13 work, Slurm, agent dispatch,
commit or push is authorized by this assignment.
