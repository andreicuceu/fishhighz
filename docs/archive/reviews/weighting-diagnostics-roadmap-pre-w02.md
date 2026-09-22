# FishHighz forest-weighting diagnostic roadmap

Status: W01 revisions 1 and 2 are implemented and independently reviewed. The
numerical result and W01-R1/R2/R3 repairs pass, but revision 2 retains
non-finite-JSON and execution/check-provenance defects. W01 revision 3 is written
for user approval; no acceptance or progression to W02.
Created: 2026-09-15.

## Purpose and relation to the package plan

Identify whether the observed finite-update and magnitude-resolution sensitivity
comes from the inherited cumulative rule, changed physical inputs, discretization,
or an arithmetic/implementation defect. Separately quantify forecast changes
conditional on a specified finite weighting prescription. A change in forecast
amplitude is not itself a convergence failure.

This investigation runs alongside the package roadmap (local-only path: `FISHHIGHZ_IMPLEMENTATION_PLAN.md`).
It does not replace package Step 13, approve its pending revision-2 review,
adopt a weighting limit, relax Step 12 convergence requirements or authorize
Step 14. Parallel means separate documents and review sequence, not concurrent
agents or automatic execution. User control is unchanged.

- This file owns the diagnostic questions, sequence and progress register.
- WEIGHTING_DIAGNOSTIC_STEP.md (local-only path: `lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md`)
  contains exactly one detailed diagnostic assignment: W01 revision 3.
- Agent prompts (local-only path: `FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PROMPTS.md`) explicitly select
  this investigation instead of package `IMPLEMENTATION_STEP.md`.
- New handoffs/reviews are `lib/fishhighz/reviews/weighting-diagnostics-wNN-rR.md`
  and `weighting-diagnostics-wNN-review-rR.md`.
- New numerical evidence is exclusive to
  `lib/fishhighz/.validation/forest-weight-diagnostics/wNN-rR-<UTC>/`.

## Questions and scientific distinctions

1. Was the original three-update calculation already sensitive to more updates,
   on its exact original grid and inputs?
2. How much do changes outside weight preparation affect the forecast with the
   original weighting prescription retained?
3. How much do signal/response changes affect the forecast if their feedback into
   weight preparation is deliberately suppressed?
4. Do updated weighting inputs change the behaviour even with the original
   rectangular magnitude sum and three updates?
5. What happens when iteration count and magnitude sampling are changed separately,
   and is there an interaction requiring a joint test?

Three questions remain distinct: faithful evaluation of a specified three-update
prescription; stability over a stated finite range of updates/grid resolutions;
and existence of an iteration/refinement limit. This diagnostic sequence first
addresses the first two. It does not require solving Step 13's continuum problem
before learning which profile changes matter. Nor does it prove optimality.

There are two independently adjustable direct controls: update count and
magnitude discretization. The third comparison requested in the discussion is
their interaction, not a third independent parameter.

## Minimal comparison design

Start with one forest auto-power, `lya(qso)_lya(qso)`, bin 0 of the saved 15x2pt
case. A forest auto-power is sufficient to measure the weighting contribution
to aliasing/pixel noise and its BAO error. Reusing its column from a 15x2pt archive
is not running a 15x2pt forecast. Preserve k/mu nodes and use a scalar covariance
with a 2x2 ap/at Fisher sum. No within-bin evolution or configuration-space
correlation calculation is introduced.

Evaluate coefficients first. At most two one-spectrum Fisher sums are needed to
establish whether the first coefficient discrepancy matters for that observable.
Keep the saved 5000 Fourier nodes for W01 rather than introduce a second,
unvalidated Fourier-grid change. Smaller memory/operation count comes from
one selected spectrum, not necessarily fewer nodes.

### Operational definitions for controlled comparisons

| Comparison type | Held fixed | Allowed changes and interpretation |
| --- | --- | --- |
| Strictly indirect | Every recurrence input, original magnitude grid/count, weight vector and A/P_pixel; the underlying spectrum/template and response inputs remain fixed where changing them would affect weighting. | Operations genuinely outside preparation, such as volume/mode integration or a derivative-operator comparison on an explicitly identical supplied signal. Identify exactly which operators are well-defined; do not import an updated template under an indirect label. |
| Frozen-weight signal/response | The same legacy recurrence inputs, weight vector and A/P_pixel. | Deliberately change the selected mean/signal or response outside preparation. Report the missing feedback into S/B/weights. This is a controlled partial-effect calculation, not a self-consistent updated survey. |
| Updated inputs with legacy direct controls | Original 107-point rectangular magnitude rule and three updates. | Recompute all intended inputs, S/B, weights and noise consistently under the accuracy conventions. This is the intermediate profile requested by the user. |
| Direct-control tests | One fully specified physical/input convention. | Change count alone, magnitude sampling alone, then both only if needed. |

The user explicitly selected both the strictly indirect and frozen-weight
signal/response comparisons on 2026-09-15. Keep them separately labelled.
Changing response while freezing only the weights need not freeze the observed
noise, because aliasing still carries W squared. Each step must list the actual
held/changed arrays and conversions; labels are insufficient.

The earlier profile comparison (local-only path: `lib/fishhighz/reviews/cumulative-forest-weighting-profile-comparison-2026-09-15.md`)
is a starting inventory, not proof that its categories are causally independent.
Before each later assignment, derive the dependency list from live code. Changes
to density/variance/S/B/L/pixel scales are coupled; their effects on stability
must be measured rather than assumed small or independent.

## Sequence: small assignments with a review after each

Only W01 has detailed implementation instructions. Later rows define questions,
bounded candidate experiments and decision points; they are not dispatchable
assignments. The reviewer recommends the next smallest useful row; the user
decides whether to stop, skip attribution-only work, or request its plan.

| Step | Question and smallest comparison | Evidence and early exit |
| --- | --- | --- |
| **W01 — exact legacy iteration sensitivity** | Saved QSO forest auto, bin 0, original signed density/SNR arrays, S/B, units and 107 rectangular samples. Reproduce three updates; compare 6, then 12/24 only if earlier counts do not identify sensitivity. Project the first confirmed difference with at most two scalar-covariance Fisher sums. No new model or reader calls. | If a verified change already exists, report that profile changes are unnecessary for sensitivity in this example and stop. If the exact replay fails, stop on the identity/arithmetic problem. A finite plateau is only bounded evidence. |
| **W02 — strictly indirect forecast effects** | Retain the W01 three-update preparation; select only truly independent operator/volume/Fourier changes whose effects remain of interest. One spectrum/bin, baseline plus one composite endpoint initially; use exact native evaluations or verified saved common-node arrays. | Quantifies the conditional forecast shift without changing the recurrence. It cannot explain a change in A/P_pixel when those are identical. Omit changes with zero effect on this auto-power; do not add a cross spectrum merely to fill an inventory. Split a composite effect only if the user needs attribution. |
| **W03 — frozen-weight signal/response effects** | Retain the same legacy weight vector and coefficients. Add the requested signal/template/response changes only outside preparation, with an explicit list of feedback paths suppressed. One baseline/endpoint pair on the same spectrum. | Isolates the partial effect the user requested. It is not the full effect of a coupled change and is not a converged accuracy profile. Stop if the desired forecast attribution is already answered. |
| **W04 — all updated inputs, legacy direct controls** | Recompute the updated physical/model/input conventions at the original 107 magnitude nodes and three updates. Compare to the legacy three-update baseline and W03; separate observable changes from changed preparation. Compare 3 versus 6 updates in coefficient space under both input conventions if not already available. | If sensitivity appears or changes here, the new magnitude method is not required. Resolve a consequential input group with one forward substitution and its reverse, in a revision of W04; do not launch all combinations. Differences between W03/W04 are preparation feedback only where every other downstream convention is identical. |
| **W05 — iteration count alone** | On the input convention selected from W01/W04, retain the rectangular grid and extend only the count: use existing 3/6 results, then 12 and at most 24 if still informative. Coefficients first, at most one new one-spectrum consequence. | Distinguishes a finite-update change from a grid change. Reuse W01/W04 results; skip this step if they already answer it. Stop on the first confirmed sensitivity or arithmetic explanation. |
| **W06 — magnitude sampling alone** | Fix three updates and the same physical density/variance functions. First reproduce the rectangular rule with explicit weights; then compare one partitioned Gauss–Legendre sampling. If needed, refine the rectangular grid and GL rule separately, retaining the same support/policies. | Determine whether fixed-count results approach a common value, differ due to prefix discretization, or fail because of interpolation/support/range handling. Include a fixed-weight control. Do not interpret nominal total-integral GL order as prefix accuracy. Stop once the origin of a discrepancy is isolated. |
| **W07 — interaction only if unresolved** | A minimal 2x2 experiment: rectangular/GL crossed with 3/6 updates, same physical inputs and support. Reuse all available corners and calculate only missing ones. | Quantify whether the count effect depends on the magnitude rule. If necessary, inspect a tiny synthetic prefix example before more real sampling. Do not escalate to infinite-iteration/continuum claims. |
| **W08 — one targeted confirmation, optional** | Only after a mechanism is identified, repeat the single decisive coefficient contrast for `lya(lbg)` in the same bin, or one other bin if a specific hypothesis predicts different behaviour there. Choose one, not both. | Extends or restricts the inference. A new spectrum forecast is justified only if the coefficient result does not answer the stated question. No automatic all-population validation. |

W02/W03 answer attribution questions; they are not prerequisites for diagnosing
the recurrence. If W01 already identifies the inherited sensitivity, the user
may stop or request a targeted mechanism test instead. If W04 identifies a
coupled input as the decisive factor, revise W04 to isolate it before spending
time on W05–W07. The aim is a minimal supported explanation, not completion of
every row. One example cannot establish prevalence across the survey.

### Conditional refinement bounds for future plans

Use at most three levels per magnitude method initially. A suitable proposed
rectangular family is 107/213/425 endpoint-inclusive nodes, halving the original
spacing while preserving its endpoint convention. A GL family can reuse the
fixed partition with orders 4/8/16; extend only in a later approved revision.
Evaluate the same frozen interpolants at new nodes; do not interpolate the
107 sampled weights and call that the same evolving rule. New real interpolation
and any template/background evaluation must be explicitly listed in that step's
approved assignment. Missing saved arrays do not authorize broad regeneration.

For W07, report both within-rule count differences and their interaction, e.g.
for a positive quantity X,

    [log X(GL,6) - log X(GL,3)] - [log X(rect,6) - log X(rect,3)].

Use signed/zero-safe alternatives with stated denominators when X is not
positive. Do not add sequential forecast percentage changes as independent
physical errors. Agreement at one count/grid does not prove convergence or
exclude compensating effects.

## Common scientific and numerical requirements

- Preserve the same cumulative recurrence. No replacement optimizer,
  renormalization inside its nonlinear update, support removal, automatic
  flooring or relaxed production guards is a repair authorized by this plan.
- Legacy signed inputs require literal signed arithmetic and a separate physical
  interpretation. Positive-mass proofs and the Step 13 log evaluator do not
  automatically apply. W01 planning found nine negative QSO and twelve negative
  LBG density interpolants in the saved bin-0 legacy rows; counts are metadata
  observations to verify, not a diagnosis of their effect.
- Compare A and P_pixel separately, using each quantity's own scale. Establish
  arithmetic fidelity with independent scalar/high-precision controls before
  interpreting changes as scientific sensitivity. Underflow, singularity,
  negative/unphysical coefficients and finite coefficient drift are different
  outcomes. Keep finite failed results inspectable.
- Preserve P3D/P1D independence, signed cross powers, fixed h units, fixed
  observed cuts and response ownership. Changing weights between experiments
  changes their covariance; each Fisher derivative still uses fixed weights and
  covariance, with no covariance-derivative information or added priors.
- For a single auto-power, C=2*T^2/modes. A cross-power, if later justified,
  requires both autos and their noise for covariance closure. Never isolate a
  cross spectrum by dropping the powers required for its covariance.
- Use 5e-12 for baseline/replay consistency, with per-quantity scaling and exact
  zero handling. The inherited forecast sensitivity budget is 1e-3 for Fisher
  and individual ap/at errors. The W01 1e-3 coefficient trigger is only a
  diagnostic stopping rule, not a new scientific acceptance criterion.
- A supported sensitivity needs an independent arithmetic reproduction and
  unchanged controls apart from the declared contrast. A stable finite range
  must be reported as such; no plateau or vanishing amplitude proves a limit.

## Cost and execution boundaries

Default to saved-array or synthetic work. One selected spectrum, one redshift
bin and one process; OMP/OPENBLAS/MKL thread counts are one. Target <=30 seconds
per numerical invocation, retain partial evidence at a cap, and do not increase
scope automatically. Scalar-coefficient calculations precede Fourier work.
Use existing local interpreters; do not modify shared environments.

Do not invoke `NewForecast`, `desi2_profiles.run`, the full study controller or
another helper that unconditionally iterates over all fields/bins. Inspect helper
scope before reuse. Later targeted new evaluations must be explicitly included
in the approved step and implementation prompt. No full seven-case suite,
joint 15x2pt forecast, reference recapture or Slurm work is authorized by this
planning request or by general approval of this roadmap.

Only affected self-contained tests and independent small numerical checks are
routine completion requirements for these standalone diagnostics. No mandatory
new wheel, environment, full ordinary suite or installed examples for every
diagnostic script. If production edits become necessary, stop and request a
revised scope; their appropriate regressions belong to that separate assignment.

## Plan / implementation / independent review cycle

1. Planning agent reads live source, existing evidence and user feedback; writes
   only the requested W-step/revision in WEIGHTING_DIAGNOSTIC_STEP.md. Archive
   the previous instructions before replacement. Ask before consequential new
   estimator/input choices; do not reopen confirmed choices.
2. User reviews/approves the concrete instructions and dispatches an implementer.
3. Implementer executes only that assignment, its bounded diagnostics and focused
   tests; saves an immutable handoff and stops. An early scientific finding is
   a valid completion if the assigned checks supporting it pass.
4. Independent reviewer inspects exact source/input identity, reruns focused
   tests and a separate numerical oracle, then writes a separate review. If
   repairs are needed, revise the same W-step with tests closing each finding;
   do not silently fix production or diagnostic code during review.
5. User controls acceptance, continuation, skipping and termination. A review
   pass neither accepts the weighting science nor advances either roadmap.

Reports must state: question; held/changed quantities; selected population/bin;
source and instruction hashes; actual attempts and early stops; independent
checks; coefficient and available forecast results; strongest supported
conclusion and its limits; smallest remaining question. Failures or omitted
attempts are not zero effects. Preserve historical reports and source artifacts.

## Progress register

| Step/revision | State | Evidence / user decision |
| --- | --- | --- |
| W01 r1 | Implemented; independent review requires repair | The exact t=3 replay and first t=6 difference are independently reproduced: A changes by 0.8616%, P_pixel by 1.1781%, and both one-spectrum BAO errors exceed the 0.1% diagnostic budget. Review (local-only path: `lib/fishhighz/reviews/weighting-diagnostics-w01-review-r1.md`) finds incomplete semantic/provenance binding, incorrect exact-zero scalar handling and no serializable failed-baseline handoff. The unmodified numerical evidence is truthful; no weighting prescription is accepted. |
| W01 r2 | Implemented; independent review requires finalization repair | W01-R1/R2/R3 pass and the numerical result remains independently reproducible. Review (local-only path: `lib/fishhighz/reviews/weighting-diagnostics-w01-review-r2.md`) finds that NaN can bypass finite JSON comparisons and that execution/required-check claims are not semantically or chronologically validated. |
| W01 r3 | Proposed; awaiting user approval/dispatch | Same-step W01-R4/R5 repair: reject non-finite JSON and bind a non-self-referential required-check/execution record. Preserve the t=6 early stop; no new scientific evaluation or production change. |
| W02–W08 | Conditional roadmap only | Plan only the requested next comparison after review and user direction. |

## Planning evidence and unresolved choices

Live source, package Git status, main instructions and Step 12/13 reports were
inspected on 2026-09-15. The legacy bin-0 archive is present with 107 magnitude
nodes, 5000 Fourier nodes, saved per-pair Fisher results and per-field legacy
coefficient inputs. W01 revision 1 verified the source hashes and reproduced the
selected baseline before changing the iteration count; its independent review
reproduces the numerical result. Revision 2 is restricted to recertifying that
result with complete evidence binding and failure-path coverage. Revision 2
closes the original three findings for finite payloads; revision 3 is limited to
the remaining JSON and execution/check provenance boundary.

The user chose to include both strictly indirect and frozen-weight signal/
response comparisons, separately labelled. No unresolved scientific choice is
needed for the W01 repair. Choosing an alternative estimator, physical
density-cell widths/support, accepting three updates as definitive, or adopting
a continuum limit remains outside this investigation and requires the user's
decision.
