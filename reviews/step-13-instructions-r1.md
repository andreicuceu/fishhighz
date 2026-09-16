# Step 13, revision 1: normalized limit of the cumulative forest weights

Status: proposed investigation following the user's selection of option 2.
Await user plan approval and implementation-agent dispatch. The user explicitly
authorizes progression to this investigation; Step 12's R2 repair passed review,
but its unconverged accuracy results are not scientifically accepted. R3 is now
carried by Step 13. Do not interpret this assignment as selection of another
estimator, changed acceptance criteria or permission for a real forecast run.

Read [AGENTS.md](AGENTS.md), [design section 0 and weighting conventions](../../FISHHIGHZ_DESIGN.md),
the [roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md), the
[planning assessment](reviews/step-13-planning-r1.md),
[r6 handoff](reviews/step-12.md), [r6 review](reviews/step-12-review-r6.md), and the
[r5 scientific diagnosis](reviews/step-12-r5.md). The previous active assignment
is preserved in [step-12-instructions-r6.md](reviews/step-12-instructions-r6.md).
Implement only this bounded investigation when dispatched; write its handoff
and stop for independent/user review. Do not dispatch agents or advance steps.

## 1. Scientific question and scope

Determine whether the existing cumulative forest-weight update has a finite,
physically usable limiting pair of noise coefficients A and P_pixel that is
independent of magnitude discretization. Distinguish this question from whether
unnormalized weights vanish, whether a numerical representation underflows,
and whether the final two-parameter Fisher inversion is well conditioned.

A supported negative result is a valid outcome of this investigation. Do not
assume that normalization cures R3, optimize an alternative estimator, or tune
settings to reproduce the legacy errors. A bounded unresolved result must state
what is missing; it is not a successful convergence demonstration.

Authorized work on dispatch: mathematical derivation, small deterministic
synthetic tests, a validation-only equivalent numerical prototype, and bounded
read-only studies using the existing saved 15x2pt weight inputs. No new real-survey CAMB or
external P3D/P1D evaluation, raw-reader resampling, full 15x2pt forecast, other six
benchmark cases, NewForecast capture, or real sensitivity suite. Ordinary synthetic model examples remain in the quick checks. No changes to
production weighting, general covariance/Fisher/model/response APIs, strict
reader defaults, input floors/support or overlap assumptions. Do not relax the
current underflow guards or silently use a diagnostic result in a primary run.

Expected new code is narrowly confined to a validation module, a small explicit
analysis script and self-contained tests. Suggested locations are
`fishhighz/validation/weight_limit.py`, `scripts/diagnose_weight_limit.py`, and
`tests/test_weight_limit.py`; choose similarly clear names if the actual layout
requires it. Array calculations must remain separate from file access and
plotting. Use NumPy and the standard library; no new runtime dependency or
compiler/backend redesign. An equivalent log/amplitude representation is a
numerical investigation, not a new weighting prescription.

## 2. Baseline, provenance and preserved evidence

Verify the live dirty/untracked checkout before editing. The reviewed r6 baseline
has 1268 parametrized cases: 1243 passes and 25 optional-compiler skips. Its exact
52-module wheel is in `.validation/step12-r6-20260914T222237Z/dist/`, SHA256
`7286a2c34014807c2d3309f9b9f140594aeb64fd8fe415932ba234cae399697c`.
The r6 independent review also passed 251 installed regressions and six examples.
These are historical checks, not new Step 13 passes. Preserve all existing tests,
including R1/R2 corruption, range, optional-import and r4 performance regressions.

Read-only scientific inputs are under
`.validation/step12-r5-20260914T191855Z/`:

- `weight-diagnosis/diagnosis.json` and `bin-0.npz` through `bin-5.npz` contain
  both forest populations, magnitude orders 4/8/16/32/64, sampled masses,
  magnitudes, variances, finite weights and fixed/coupled coefficient evidence.
- `profiles-checked/` contains 12 primary records (six bins, two profiles), 180
  individual-spectrum results and 72 diagnostic outcomes. The report settings
  supply L, geometry, response conventions and sample metadata. Verify the source
  relationship to the older inputs used by `diagnose_desi2_weights.py`; record
  any quantities reconstructed from unchanged scalar settings explicitly.
- The r5 joint, pair, attribution, sensitivity and weight figures are historical.
  Six compatibility bins pass; all six accuracy bins fail; 12 diagnostics are
  unavailable. Preserve these verdicts and error reasons unchanged.

Use hashes, report contexts, pair/field identities, units and actual controls,
not file names alone. Require exactly six bins and both `lya(qso)`/`lya(lbg)`
populations for a claimed complete real-input diagnosis. Missing files do not
authorize new reader/model calls; finish independent synthetic work and report
which real-input checks are unavailable. Preserve r2/r4/r5/r6 reports, arrays,
archives and review directories. New outputs go to a unique
`.validation/step13-r1-.../` directory with an exclusive output-directory policy.

## 3. Phase A — equations, units and exact finite-step equivalence

Start with `fishhighz/kernels/weights.py`, `fishhighz/weights.py`,
`fishhighz/noise.py`, `validation/weight_diagnosis.py`, and
`scripts/diagnose_desi2_weights.py`. Inspect live caller semantics and preserve
magnitude ordering, inclusive prefix sums and positive quadrature masses.

Use notation r_i=rho_i*q_i, v_i=variance_i, Delta_v=pixel width, L=forest length,
S=auxiliary signal, B=auxiliary alias, alpha=L/Delta_v and a_i=v_i/S. On supported
samples with positive variance the accepted recurrence and coefficients are

    w_i(0) = 1/(1 + Delta_v*v_i/B)
    N_i(t) = alpha * sum(j<=i, r_j*w_j(t))
    w_i(t+1) = N_i(t)/(N_i(t)+a_i)
    I1 = sum(r*w), I2 = sum(r*w*w), I3 = sum(r*w*w*v)
    A = I2/(L*I1^2)
    P_pixel = Delta_v*I3/(L*I1^2).

r has units deg^-2 (km/s)^-1; v and w are dimensionless; S is deg^2 km/s;
B and L,Delta_v are in km/s. N and a have matching units. A is deg^2 and
P_pixel is deg^2 km/s. Retain the same auxiliary signal/alias, variance and
population definitions. P1D remains independent of P3D.

Derive an exact amplitude/shape representation before coding its evaluator.
One useful parameterization is w=s*u, max(u)=1. For positive denominators,

    h_i = alpha * sum(j<=i, r_j*u_j)
    b_i = h_i/(a_i+s*h_i)
    u_next = b/max(b)
    s_next = s*max(b).

Retain the amplitude, e.g. log(s), throughout the nonlinear update. Simply
renormalizing w and feeding it back into the original update is inequivalent.
Derive a log-domain or equivalently range-safe calculation; do not reconstruct
tiny s merely to feed it into an unstable product. Track relative components as
well as the global scale: u itself may develop extreme dynamic range. Exact
zero support and zero variance require explicit branches, not tiny floors.
Positive supported weights must not silently become exact zeros.

Derive A and P_pixel directly from the shape with cancellation of common scale
before arithmetic. Retain log coefficients if their mathematical values exceed
float64, with explicit availability flags for ordinary values. A finite log
representation is not a finite converged physical coefficient. If a contribution
must be dropped because of relative range, justify an upper bound on its effect
on I1/I2/I3 and the final ratios; record the bound and fail when it is insufficient.
Do not declare equivalence by globally ignoring floating-point exceptions.

## 4. Phase B — discrete and continuum asymptotics

Analyze before undertaking large iteration scans. For strictly positive a,
the zero-weight linearization is

    M_ij = alpha*r_j/a_i for j<=i, and zero otherwise.

Establish the domains of any bound, fixed-point result or asymptotic reduction.
Analyze spectral radius below, equal to and above one, positive fixed points,
unique and repeated dominant diagonal entries, near-degeneracy, initial-support
accessibility and nonnormal transients. The real-input maximum spectral radius
0.4072476439 is historical evidence for amplitude decay; it is not a bound on
noise convergence or on the rate of normalized-shape convergence.

At fixed mesh, determine what controls the limiting shape and coefficient
moments. Under mesh refinement, a quadrature cell's mass changes. Examine whether
the shape concentrates on cells whose mass tends to zero, and whether A or
P_pixel diverges despite convergence on each fixed mesh. Respect the physical
magnitude coordinate: changing node order is not an invariant of a prefix rule.
Subdivision means resampling the same continuous density/variance and conserving
the underlying measure, not inventing more sources.

Write the continuum prefix operator for a fixed support and density measure.
Distinguish limit(t -> infinity) at fixed grid, limit(grid -> continuum) at fixed
t, and joint limits. Establish any assumptions needed to exchange these limits;
do not claim interchangeability from a single diagonal refinement sequence.
For the full nonlinear rule, identify when a small-amplitude linearization is
controlled and estimate its accumulated error. Counterexamples from a linearized
operator alone must not be presented as proofs for the full real-data recurrence.

Mandatory analytic controls include:

- One positive cell, d=alpha*r/a: w_next=d*w/(1+d*w). For d!=1,
  w_t=d^t/[1/w_0+d*(1-d^t)/(1-d)]; for d=1,
  w_t=w_0/(1+t*w_0). For this cell, A=1/(L*r) and
  P_pixel=Delta_v*v/(L*r), even when w_t -> 0. This is a finite-coefficient
  positive control against conflating vanishing weights with divergent noise.
- Constant-coefficient linearized continuum operator
  Kf(x)=c*integral(0..x,f(y)dy), x in [0,1], starting from f=1.
  The normalized shape is x^t and, for total measure R,
  A_t=(t+1)^2/[(2*t+1)*L*R], P_pixel=Delta_v*v*A_t.
  This diverges with t. It is a required negative control against certifying
  coefficient convergence from bounded normalized weights or small amplitudes.
- Constant and nonconstant finite grids, including repeated/near-equal diagonal
  eigenvalues and subdivision of the same smooth measure. Derive expected
  behavior explicitly rather than treating a dense numerical eigensolver as
  an independent proof. Use tiny high-precision controls for the nonlinear map.

The planning note and its tiny checks are starting evidence; independently
verify them and extend the analysis. A general theorem need not cover arbitrary
unphysical inputs, but its stated assumptions must include each real sample to
which a conclusion is applied. Record exceptions separately.

## 5. Phase C — prototype tests before real-input application

Keep the new evaluator diagnostic-only, preserving `_iterate`, `_integrals`,
`prepare_forest_weights`, current method names and their error behavior.
No production method switch or automatic fallback is authorized in revision 1.

Concrete acceptance tests:

1. Compare direct float64 recurrence, the amplitude/shape or log evaluator and
   an independent Decimal recurrence on small positive arrays. Include 0/1/3/6/
   12/24 updates, multiple initial amplitudes, unequal masses/variances and
   scalar d below/equal/above one. Wherever representable, match weights or
   reconstructed log weights and A/P_pixel within 5e-12 relative discrepancy.
   Scale each tested quantity independently. Test coefficient scale invariance
   separately from the non-invariance of the nonlinear update.
2. Continue small Decimal controls past the original underflow boundary. Compare
   log amplitude, relative shape and both noise ratios; do not require float64
   reconstruction of an unrepresentable w. Include the Step 10 mixed-product
   case and the saved r5 approximately 8.12e-566 product. The old production
   path must still reject; a mathematically justified diagnostic result must
   retain that historical failure and its new representation's provenance.
3. Verify exact-zero density, supported zero variance, all-zero/invalid support,
   nonfinite/negative inputs and ordered quadrature validation. Either reproduce
   the existing allowed behavior or explicitly declare a narrower mathematical
   domain for the diagnostic; do not silently alter caller semantics. Preserve
   zero coefficients as exact zeros where physically required.
4. Independently verify the scalar and continuum controls in Phase B. Use direct
   analytic moments or sufficiently exact quadrature, plus precision refinement
   (e.g. Decimal 80/160 digits) for selected small nonlinear cases. Include a
   deliberately false normalized-only recurrence and ensure an equivalence test
   distinguishes it from the exact amplitude-retaining evolution.
5. Expose nonconvergence with repeated/near-degenerate spectra, slow drift and
   concentration. A tiny change between adjacent late iterations or vanishing
   unnormalized updates cannot by itself produce a convergence verdict.
6. Bind reported iteration/refinement comparisons to their actual saved inputs
   and trajectory snapshots. Test swapped field/order, mismatched coefficient
   operands, missing capped/failed attempts, and a divergent control mislabeled
   converged. Preserve the r6 per-spectrum evidence semantics; do not use a
   stack-wide norm that hides a weak population or small coefficient.

Use tests that exercise the derived mathematics and scientific verdicts rather
than duplicating implementation expressions. Do not repeatedly allocate a real
forecast grid for unit tests. No ordinary test may depend on local DESI data.

## 6. Phase D — bounded diagnosis of the saved 15x2pt populations

After finite-step equivalence and counterexample tests pass, use only saved
sample arrays for both populations in all six bins. Recover S/B/L/Delta_v and
sample provenance from verified saved metadata and the frozen recipe; document
any scalar unit conversions. No model calls or new reader interpolation.

First reproduce all available saved 3/6/12/24-update coefficients and classify
all recorded arithmetic failures. Then extend the same discrete recurrence in
the diagnostic representation, recording actual iterations and failures. Start
with checkpoints 3, 6, 12, 24, 48, 96, 192, 384, 768. Use an initial cap of
1024 updates per population/order; this is a work bound, not a convergence
criterion. Include all existing magnitude orders 4/8/16/32/64. Use analysis of
the spectrum/concentration to stop a disproved branch early with a documented
reason rather than exhausting every combination mechanically.

Run one population/order or similarly small batch per invocation, with one-thread
settings, bounded memory and a practical wall-time cap (target <=30 s for each
batch). Preserve partial progress and report a cap honestly. If the initial cap
cannot distinguish finite convergence from slow drift, estimate the needed
additional work analytically and report it for user review; do not automatically
launch a long job, increase grid orders, or request/use Slurm.

Save log amplitude, a measure-normalized weight distribution
p_i=r_i*u_i/sum(r*u), A/P_pixel (and logs/availability where needed), moments of
magnitude, concentration widths/effective support and appropriate residuals.
A useful continuous-measure concentration statistic is
sum(p_i^2/r_i)=L*A on supported cells; it detects shrinking effective sightline
measure without confusing amplitude loss with information change. Include
spectral gap/degeneracy diagnostics where analytically applicable. Plot population
results separately; do not allow the dominant QSO forest to hide LBG behavior.

For a candidate finite limit, require agreement across at least three successive
iteration doublings and three nested existing grid orders, with relative changes
in both positive coefficients <=1e-4 as a screening target. Handle an exact zero
P_pixel analytically. This is a numerical screening target for the investigation,
not a replacement for the original 1e-3 forecast-error/Fisher requirements.
Support the candidate with a derivation or controlled remainder/concentration
analysis that rules out a transient plateau at the claimed accuracy. If existing
data do not supply sufficient refinements, classify it as unresolved. Fixed-mesh
agreement alone cannot establish a continuum limit.

A negative conclusion requires an explicit divergence/grid-dependence argument
and numerical controls under its assumptions, identifying precisely which real
populations it covers. Report any remaining populations separately. An unresolved
case is neither proof of nonexistence nor a numerical pass.

Do not reassemble a new full 15x2pt forecast in this phase. Explain the implication
for forest noise using the existing relation
N_F(k,mu)=[A*P1D(k_parallel)*W(k_parallel)^2+P_pixel]*d_deg^2/a_v.
Hold all other quantities fixed. A future user-authorized forecast remains needed
to establish joint and individual-spectrum Fisher convergence under any adopted
formulation; coefficient screening alone cannot scientifically accept Step 12.

## 7. Phase E — conclusion and subsequent decision

The handoff must distinguish three outcomes, by population and globally:

- **Finite equivalent limit supported:** state assumptions, limiting coefficients,
  error estimates and grid/iteration tests. Supply the validated diagnostic
  method and a concrete proposal for later production integration. Do not
  switch the primary accuracy profile in this revision. User review decides
  adoption and any required API/arithmetic-policy change or real validation run.
- **No suitable finite grid-independent limit:** provide the derivation and
  counterexample/real-input evidence. State which condition fails and why more
  precision or more iterations cannot remedy it. Recommend the scientific
  question the user must settle next, without selecting a new estimator or
  finite-iteration acceptance convention.
- **Unresolved within the stated bounds:** identify the missing mathematical or
  numerical evidence, the limiting resource/precision issue and a bounded next
  investigation. Do not label R3 resolved or request an automatic forecast rerun.

The first two can complete this investigation after independent review; only a
validated adopted prescription and its forecast checks can establish scientific
acceptance of the accuracy profile. Preserve R2 closure and all earlier evidence.
The user has chosen to investigate the same rule, not approved a different
estimator, density floor, support change, convergence waiver or production policy.

## 8. Quick validation and handoff requirements

On implementation dispatch, run focused mathematical/diagnostic tests first,
then the ordinary suite and Ruff with one-thread settings:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
```

Retain all 1268 baseline cases. Run affected installed tests and the six examples
outside the checkout with isolated imports. If validation code is packaged,
build the final exact wheel and verify all source/wheel/installed module hashes
and METADATA/WHEEL bytes. Record interpreter, dependencies and actual backend;
NumPy remains the default. Reuse appropriate local environments without modifying
siblings. Preserve the disclosed explicit-Numba/SciPy-blocked subprocess failure
and unmet NumPy matrix speed target; neither is assigned for repair here.

Write `reviews/step-13.md` and a technical derivation at
`reviews/step-13-weight-limit-r1.md`, with immutable figures/tables under the new
artifact directory. Include equations and assumptions, exact recurrence/domain,
source/input hashes, before/after finite-step agreement, analytic/high-precision
controls, all 12 population diagnoses across the existing grid orders, capped
and failed attempts, concentration and A/P_pixel plots, and a decision table.
Record commands, timings of targeted diagnostics only, test outcomes, wheel
identity and a changed-file manifest. Do not present a forecast runtime estimate
as measured performance or historical r5 figures as new results.

Complete the authorized mathematical and synthetic work even if some historical
inputs are missing; document the incomplete real-input coverage. Do not alter
planning/governance documents during implementation without user instruction;
report any needed revision. Stop for independent/user review. No production
adoption, new real forecast, agent dispatch, Step 14, commit or push.
