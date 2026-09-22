# FishHighz compatibility weighting tests

Created: 2026-09-17. Status updated 2026-09-18: stages 1–5 are complete.
The stage-4 BAO calculation and additional magnitude-grid check passed independent
scientific review. The user selected early lyaforecast weights as the baseline,
with McDonald weights retained as an explicit alternative. See the
[weighting decision](../../research/FOREST_WEIGHTING_DECISION.md) for results and
limits, and the [roadmap](FISHHIGHZ_IMPLEMENTATION_PLAN.md) for the next proposed
three-profile comparison. Step 13 and the W-series remain closed.
The specification below records the completed study, including its historical
all-15 bin-1 results. All future forecasts exclude bin-1 correlations involving
LBG, LAE or lya(lbg), as specified in the revised roadmap. Source:
[LYAFORECAST_WEIGHTING_FINDINGS.md](../notes/LYAFORECAST_WEIGHTING_FINDINGS.md).
The user confirmed five variants, individual and joint BAO results, and omission
of poorly constrained bins from plots (the precise cut is specified in stage 4).

Historical execution decision, 2026-09-17: the user authorized the coordinating agent to
dispatch Astra (Light) implementers and Sol (High) scientific reviewers, check
each stage and progress through stages 1–3. Review corrections must be likely
to affect the scientific conclusion. This supersedes earlier requirements for
separate user dispatch/progression within these three stages. It does not
authorize stage 4 forecasts, Slurm actions, commits or pushes.

## Execution record

| Stage | State and evidence |
| --- | --- |
| 1 | Implemented; Sol (High) scientific review passes, with no consequential correction. All 54 saved weight/coefficient preparations reproduce exactly; six-bin baseline Fisher discrepancies are at roundoff. Coordinator reran all 18 focused tests and checked coverage. [Handoff](../reviews/compatibility-weighting-stage1.md); [review](../reviews/compatibility-weighting-stage1-review.md). |
| 2 | Implemented; Sol (High) scientific review passes with no consequential correction. All 810 trajectories and 12 plots saved; coordinator reran 21 tests, checked 849 converged classifications against their saved tails, and verified input hashes. [Handoff](../reviews/compatibility-weighting-stage2.md); [review](../reviews/compatibility-weighting-stage2-review.md). |
| 3 | Implemented; Sol (High) scientific review passes with no consequential correction. Opt-in full-sum positive-signal forest-auto stopping only: 192 converged, 24 capped and 1,404 ineligible validation results; 810 exactly-three-update checks unchanged. Coordinator reran all 30 focused tests and checked Ruff lint/format for all ten new Python files. [Handoff](../reviews/compatibility-weighting-stage3.md); [review](../reviews/compatibility-weighting-stage3-review.md). |
| 4 | Complete; independent scientific review passes. Six-bin individual and joint BAO comparison at three updates and demonstrated convergence. [Handoff](../reviews/compatibility-weighting-stage4.md); [review](../reviews/compatibility-weighting-stage4-review.md). |
| 5 | Assessment and accepted baseline recorded in the [weighting decision](../../research/FOREST_WEIGHTING_DECISION.md). Additional grid-to-BAO check covers bins 2–6 and passes [independent review](../reviews/early-lyaforecast-grid-bao-review.md); bin 1 is outside that refinement result. Accuracy integration remains future work under the revised roadmap. |

The implemented stopping defaults are rtol=1e-4, minimum 3 updates, three stable
transitions, doubled-count confirmation and cap 96. Full-sum intrinsic and
moment-aliasing weights converge in 10/12 forest-auto cases on each tested grid;
the earliest Python formula converges in 12/12. LBG bins 0–1 remain nonconverged
for the first two variants. Both prefixes and non-auto preparations are ineligible
for this initial adaptive API. These are finite-grid results: the historical
formula retains 0.194% aliasing and 0.117% pixel-noise sensitivity at the final
magnitude-spacing refinement. Those stage-2/3 coefficient results alone make no
BAO-error convergence claim; subsequent bounded BAO and refinement results are
summarized in the weighting decision.

## 1. Add five weighting variants and reproduce compatibility

Retain each compatibility forest input: magnitude nodes and measure `r_i=n_i*q_i`,
pixel variance `v_i`, length `L`, pixel width `lp`, intrinsic signal `P`,
response-smoothed 1D power `B`, and seed `w_i(0)=B/(B+lp*v_i)`.
Three iterations means three updates after this seed. At each update form
`J1=sum(r*w)` and `J2=sum(r*w*w)` over the full sample (broadcast scalars) or
each prefix (arrays), then update simultaneously using
`PN_i=lp*v_i/(L*J1_i)` and `w_i(t+1)=S_i/(S_i+PN_i)`:

| Variant | Moments entering update | Iterative signal S |
| --- | --- | --- |
| `prefix_intrinsic` | cumsum | P; unchanged compatibility control |
| `sum_intrinsic` | full-sample sum | P |
| `prefix_aliasing` | cumsum for both J1 and J2 | P + B*J2/(L*J1²) |
| `sum_aliasing` | full-sample sums | P + B*J2/(L*J1²) |
| `sum_historical` | full-sample sum | P + B/(L*J1); earliest Python approximation |

Recompute aliasing from the current iterate before every update, starting from
the common seed. Full-sample moment aliasing is the paper-supported interpretation;
refreshing it is an explicit implementation choice, not a uniquely specified
published algorithm. The earliest Python formula is a separate historical
comparison. Exclude the dimensionally inconsistent intermediate historical formula.

Final coefficients always use full-sample moments:
`A=I2/(L*I1²)`, `P_pixel=lp*I3/(L*I1²)`, `I3=sum(r*v*w²)`.
Aliasing stays in the final covariance for all five variants. Apply response once.
Preserve signed compatibility inputs and literal arithmetic in this validation
path; do not silently mask them through the nonnegative public preparation API.

**Checks:** reproduce saved three-update weights, coefficients and forecasts;
use one/two-cell calculations, equal-noise positive populations, the homogeneous
fixed-point examples in the findings, and full-sample permutation invariance.
Keep reduction order common where possible (`cumsum(...)[-1]` is a full-sample
sum), distinguishing roundoff effects from the normalization change.

The compatibility runner currently reads saved total powers. Add only the
weight preparation and noise/covariance reassembly needed for these comparisons;
relabelled cached totals or substituted accuracy inputs cannot test the variants.

## 2. Measure weight convergence

Study both forest populations in all six 15×2pt bins, including any distinct
preparation contexts used by compatibility pair calculations. At fixed magnitude
grids, save trajectories through 3, 6, 12, 24, 48 and 96 updates. This is a
proposed finite cap, not an assumed sufficient count. Stop invalid branches with
their last finite state and reason; do not regularize weights or moment ratios.

Track amplitude, normalized weight shape, fixed-point residual, moments, A and
P_pixel. Plot weights versus magnitude and convergence measures versus iteration.
Separate finite nonzero convergence, amplitude decay with stable normalized
quantities, continued drift/oscillation, and invalid arithmetic. Signed inputs
support algebraic comparisons but do not inherit positive-density proofs.

Propose relative stability of amplitude, shape, A and P_pixel at 1e-3 for three
successive updates, confirmed against a doubled count within the cap. Measure
vector changes relative to their own norm, not a fixed absolute floor that
mistakes decay for convergence; report exact zeros separately. Repeat with 1e-4.
Only finite nonzero convergence qualifies for the initial adaptive prescription;
stable coefficients alone do not establish converged weights. A candidate too
late for doubled-count confirmation remains unconfirmed within this study.
Execution clarification: a doubled-count candidate must also remain stable
against every later saved state through the cap. The stage-2 calculation found
prefix-aliasing plateaus that subsequently drift; such plateaus do not establish
eligibility for adaptive stopping. This is bounded trajectory evidence, not a
proof of asymptotic convergence.

Then make a small, separate magnitude-spacing refinement with the same
compatibility interpolation, endpoints and input policies. Do not substitute
accuracy quadrature. This checks discretization sensitivity without reopening
Step 13 or attempting a continuum proof. Report any unresolved result within
the cap for review instead of extending the calculation indefinitely.

## 3. Add a conditional stopping prescription

If stage 2 demonstrates reliable convergence in a modest number of updates,
add opt-in adaptive stopping for eligible variants: a minimum count, consecutive
stable updates, tolerances and a maximum count. Return actual counts, residuals
and explicit status. Preserve fixed-count legacy behavior and the exactly-three-
update option for all variants.

Choose defaults from the measured trajectories; the stage-2 settings are proposed
study controls, not adopted guarantees. Test early stopping against the longer
trajectories. A cap or arithmetic failure means reported nonconvergence, never
a silent fallback. Do not adopt a stopping rule for a variant whose weights
have not demonstrated convergence.

## 4. Repeat the 15×2pt BAO comparison

Redo all six bins for lya(qso), lya(lbg), QSO, LBG and LAE: 15 auto/cross spectra.
For every variant calculate (a) exactly three updates and (b) converged weights
where established. Mark unavailable converged results explicitly; keep valid
three-update results. Do not substitute the last attempted iterate or W12 weights.

Hold compatibility signal, BAO Jacobian, Fourier modes, volume, redshift
conventions, response, density/SNR treatment, galaxy noise, parameters and priors
fixed. Recompute forest noise and all affected covariance entries consistently;
keep weights fixed in mean derivatives. Identify and preserve pair-specific
legacy preparation where present. Verify that saved inputs supply the required
signal/noise decomposition before reusing them.

Produce two matched sets of 15-panel plots, at three updates and convergence:
marginalized sigma(alpha_parallel) and sigma(alpha_perp) versus redshift for
each spectrum alone. Show fractional changes relative to the three-update
compatibility control and converged/three-update changes within each variant.
Add joint 15×2pt curves with the full inter-spectrum covariance; summing individual
Fisher matrices does not give the joint forecast.

**Plot cut selected by the user:** omit a variant's correlation/bin point from
both BAO-error curves when either marginalized uncertainty exceeds 0.2, or is
unavailable/nonfinite. Apply this to joint curves using the joint uncertainties.
Mask derived ratios when either operand is omitted; leave gaps rather than
connecting across excluded bins. Keep all results and reasons in tables, with
the cut stated in captions. This is only a plotting cut: retain all valid spectra
and bins in the joint calculation. Do not regularize singular individual fits.

Save tables of variant, bin, spectrum, count/status, coefficients and BAO errors.
Check baseline reproduction in every bin, unchanged galaxy-only individual
forecasts, and independent covariance/Fisher sums in a small subset. Tightening
the stopping tolerance must change each available BAO uncertainty by less than
0.1%, including valid points hidden by the plot cut; otherwise revise the
stopping choice before interpreting convergence. This assesses iteration
accuracy on compatibility grids, not total physical forecast accuracy.

## 5. Assess the results and replan

Separate normalization, iterative-aliasing and stopping effects, including their
interactions. Summarize unresolved convergence/grid sensitivity and the changes
in individual and joint BAO constraints. Then reassess the accuracy profile and
other FishHighz changes with the user and write the package-completion plan.

Each stage receives its own bounded implementation/review assignment; use minimal
analytic and saved-input checks. The current execution request covers stages
1–3, including saved-input baseline reproduction and weight diagnostics.
The stage-4 15×2pt execution is a later explicit assignment, using compute resources for
intensive work and explicit approval for Slurm actions. Preserve historical
evidence; other reference cases and fresh upstream captures are outside scope.
