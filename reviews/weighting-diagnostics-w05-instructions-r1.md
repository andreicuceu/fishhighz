# Forest weighting diagnostic W05, revision 1: explicit inverse-variance reference

Status: proposed for user approval and implementation-agent dispatch. W04 passes
independent scientific review; the user requested progression to this diagnostic.
This assignment compares a literature-motivated reference, without selecting a
production prescription. Package Step 13 and forecast acceptance are unchanged.
Created: 2026-09-15.

## 1. Scientific question and scope

Do the original three-update weights and the two finite-count-stable W04
full-sample weights agree with the explicit inverse-variance reference at the
same weighting mode? How do their combined aliasing and pixel-noise powers
compare? Stability of A and P_pixel separately does not establish optimality.

W04 found that both full-sample branches satisfy 0.1% coefficient stability over
t=6/12/24. Aliasing is not necessary for that bounded stability, but including
it changes the coefficient levels: at t=24, A decreases by 5.08747% and P_pixel
increases by 12.60455% relative to fixed P0. Neither change alone determines
whether the combined noise is lower.

W05 has two explicitly separate scopes:

1. Derive the inverse-variance reference and its minimum-noise property for a
   nonnegative density measure under the stated covariance assumptions; check
   this on a tiny positive synthetic example.
2. Evaluate that same algebraic prescription on the unchanged signed legacy
   grid and compare four weight vectors. This is an algebraic comparison, not
   a claim of physical optimality for a population with negative densities.

Only fixed-input, fixed-mode coefficient calculations are assigned. No new
iterations, magnitude refinement, density policy, Fisher calculation or forecast
is necessary. Stop after this comparison; do not select a replacement estimator.

## 2. Relevant sources and historical evidence

Read workspace/package AGENTS.md, design section 0, the main roadmap handover/
progress register, the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W04 report](reviews/weighting-diagnostics-w04-r1.md) and
[W04 review](reviews/weighting-diagnostics-w04-review-r1.md). The
[reviewed W04 instructions](reviews/weighting-diagnostics-w04-instructions-reviewed-r1.md)
are preserved. `IMPLEMENTATION_STEP.md` is package Step 13, not this assignment.

Inspect the seed and moment calculation in `scripts/compare_forest_weight_signals.py`
and the relevant live legacy/kernel definitions. Use W02's
[source audit](reviews/weighting-diagnostics-w02-r1.md) for response ownership.
Do not modify or rerun W01–W04 scripts or reopen the historical C++ search.

Use [McQuinn & White (2011), section 2.2, equations 11–13](https://arxiv.org/pdf/1102.1752)
for the reference: sightline weights proportional to the inverse sum of 1D
forest and instrumental noise powers. Its variance argument neglects
subdominant off-diagonal covariance and assumes factorized sightline weights.
State those approximations. Our B denotes the saved 1D forest power, not the
paper's arbitrary normalization B in equation 11. Distinguish the published
formula from the derivation below in FishHighz units. This is not an assertion
that the ME07/FR14 recurrence is equivalent to this estimator.

Authoritative inputs:

```text
.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json
SHA-256: 3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6
row: settings.pair_inputs["lya(qso)_lya(qso)"]

W04 checkpoints:
.validation/forest-weight-diagnostics/w04-r1-20260915T225726Z/arrays.npz
.validation/forest-weight-diagnostics/w04-r1-20260915T225726Z/summary.json
```

Verify compatibility bin 0, z=[2.0,2.235], 107 magnitude nodes from 16.1 to
26.75, and nine negative densities. Retain the signed rho, variances v,
rectangular dm=mag[1]-mag[0], forest length L, pixel length l_p, P0 and B exactly.
Match the saved W04 inputs/scalars to the authoritative row before using its
vectors. Record file hashes and the decisive source identity, without creating
an artifact-validation framework. If inputs are unavailable or materially
inconsistent, hand off the obstruction; do not recapture them.

## 3. Scientific conventions and derivation

For this diagnostic, all forests share the supplied L, l_p and smoothed B at
the saved weighting mode. The instrumental response is already included in
P0 and B; v is the saved dimensionless pixel variance. Apply no extra smoothing,
noise deconvolution, cosmological conversion or alteration of response ownership.
The individual 1D noise power is l_p v_i; P_pixel below is the effective 3D noise.

Define r_i=rho_i dm, D_i=B+l_p v_i and the reference weights

```text
nu_i = B / D_i.
```

Show that nu has exactly the relative shape of the legacy initialization
(B/l_p)/(B/l_p+v_i); its amplitude is also the same algebraically. Compare the
float expressions numerically, allowing rounding. This reference has no
iteration count or dependence on P0. Do not insert nu into a nonlinear update.

For any supplied weight vector w, retain the existing coefficient definitions:

```text
I1 = sum_i r_i w_i
I2 = sum_i r_i w_i^2
I3 = sum_i r_i w_i^2 v_i
A = I2 / (L I1^2)
P_pixel = l_p I3 / (L I1^2)
Q[w] = A B + P_pixel = sum_i r_i D_i w_i^2 / (L I1^2).
```

Q is the additive noise at this one weighting mode. The intrinsic P0 is common
to every vector. Do not add aliasing a second time or include P0 in D. Lower Q
means lower diagonal auto-power variance under these assumptions for a physical
nonnegative measure; this is not a multi-mode BAO or joint-forecast optimum.

For r_i>=0, at least one r_i>0, B>0, L>0, l_p>0, v_i>=0 and I1!=0, derive using
Cauchy–Schwarz (or completing a square):

```text
K = sum_i r_i / D_i
Q[w] >= 1/(L K), with equality for w_i proportional to 1/D_i on r_i>0.
n_eff = L sum_i r_i nu_i = L B K
Q[nu] = 1/(L K) = B/n_eff.
```

For clarity, rho has units deg^-2 (km/s)^-1 mag^-1; I1/I2/I3 have units
deg^-2 (km/s)^-1; L, l_p, B and D have units km/s; A has units deg^2;
n_eff has units deg^-2; Q and P_pixel have units deg^2 km/s.

With signed r, the identity Q[nu]=1/(L K) remains algebraic when defined, but
the inequality and physical interpretation of n_eff do not follow. Label the
real-input quantity a formal n_eff and nu an algebraic reference. Keep the nine
negative nodes; do not floor, discard, interpolate, or replace them. A result
above or below this reference on the signed grid is not evidence of an optimum.
No new physically admissible density prescription is selected in W05.

## 4. Minimal real-input comparison

Compute coefficients for just these four vectors, without evolving them:

| Vector | Source |
| --- | --- |
| Explicit reference nu | B/(B+l_p v), checked against W04 `seed` |
| Original prefix rule at t=3 | Authoritative row `_w_lya` |
| Full-sample, fixed P0 at t=24 | W04 `fixed_24_weights` |
| Full-sample, refreshed aliasing at t=24 | W04 `refreshed_24_weights` |

Use the historical float64 moment product order and `cumsum(...)[-1]` reduction
for matching old coefficients. Confirm the prefix coefficients against the last
elements of the row's `_aliasing_weights` and `_effective_noise_power`; confirm
the two W04 coefficient pairs against its summary and saved moments. This
recomputes moments of saved vectors; it does not rerun any trajectory.

For each vector report I1/I2/I3, A, A B, P_pixel and Q. Report signed fractional
changes of A, P_pixel and Q relative to nu, plus the refreshed-versus-fixed Q
contrast. Show whether A B or P_pixel dominates at this mode. Compare weight
shapes w/max(abs(w)) with the reference and report the maximum absolute shape
difference. This normalization is for display only and changes no recurrence.

For nu also compare Q from moments with 1/(L K)=B/n_eff. Record finite values,
signs and any cancellation affecting the conclusion. If D is nonpositive,
normalization is undefined, or arithmetic is nonfinite, verify the obstruction
independently and stop. If signed inputs yield nonpositive K, A, P_pixel or Q,
record the values and stop interpretation at the algebraic/domain finding; do
not repair inputs or claim improved forecast sensitivity. For zero references
report absolute differences instead of undefined fractional changes.

There is no requirement that nu give a particular ranking on the signed inputs.
The scientific result may be agreement, a difference, or a restricted domain
conclusion. Stop after answering the comparison; do not enlarge the grid or
add counts, modes, populations or survey bins to seek a preferred result.

## 5. Concrete acceptance checks and execution bounds

- Derive the positive-measure inequality, equality condition, and dimensions.
  Explain explicitly why it cannot prove optimality for the saved signed rho.
- One exact two-bin control suffices: r=(1,1), v=(1,4), L=l_p=B=1.
  Derive nu=(1/2,1/5), K=7/10, n_eff=7/10, A=29/49,
  P_pixel=41/49 and Q=10/7. Uniform weights give Q=7/4, which is larger.
  Check the moment and effective-density expressions independently with rational
  arithmetic. Multiplying nu by 3 must preserve A, P_pixel and Q.
- Compare nu to the legacy seed and independently reconstruct coefficients for
  all four saved-input vectors. Use an independent scalar 80-digit Decimal
  calculation from the exact binary inputs and saved weights, with direct sums;
  do not call the float moment helper. For the reference also compute nu in
  Decimal from B, l_p and v. Verify Q via the direct numerator sum(r D w^2)
  as well as the A B+P_pixel decomposition and effective-density identity.
- Require rtol=5e-12, atol=0 for nonzero weights/moments/coefficients and exact
  reference-zero handling. For fractional contrasts compare the underlying
  coefficients and require absolute agreement <=5e-12 in the dimensionless
  contrast (avoid relative tests on nearly zero differences). A scientifically
  relevant mismatch requires explanation or bounded handoff, not relaxed limits.

Prefer one short standalone script with assertions and a derivation in its
report; a separate focused test file is optional. No optimizer, generalized
validator, full pytest, wheel or new environment is needed. Run the tiny checks
and Ruff on new Python only, using the existing interpreter, one process and
OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1. Each numerical invocation
has a 30-second cap; retain partial evidence if reached and do not auto-extend.

## 6. Handoff and independent review

Allowed additions: `scripts/compare_inverse_variance_forest_weights.py`, optional
focused tests, `reviews/weighting-diagnostics-w05-r1.md`, and small arrays/tables
under a new `.validation/forest-weight-diagnostics/w05-r1-<UTC>/` directory.
Preserve production code, earlier scripts, artifacts and reports. Do not edit
planning documents during implementation.

The report starts with the scientific answer and includes the derivation,
assumptions, four-vector comparison, exact commands, actual checks, source/input
identity, and enough small numerical output to reproduce the comparison. Mark
W04 stability and its checkpoint vectors as historical evidence; W05's new
measurements are their combined-noise/reference comparisons. Explain separately:

- Whether the explicit reference is already the legacy seed.
- Whether the stable iterative shapes/coefficient levels agree with that seed,
  and how their Q compares; a stable recurrence need not minimize Q.
- What the positive synthetic control establishes and what the signed example
  cannot establish. No continuum or physical-forecast conclusion follows here.

The reviewer inspects only the decisive source/equations and independently
checks the tiny control and four-vector comparison. Write
`reviews/weighting-diagnostics-w05-review-r1.md`, beginning with **Scientific
conclusion and impact of proposed changes**. Propose corrections only if likely
to affect the scientific conclusion or its justified scope; give the possible
impact and smallest resolving check for each. Explicitly state when none are
needed. No optional cleanup or hypothetical provenance-hardening list.

Stop for user review. No production adoption, density/response change, new
iteration or magnitude refinement, Fisher sum, real model/survey evaluation,
full forecast, package Step 13 repair, subsequent step, Slurm action, agent
dispatch, commit or push is authorized by this assignment.
