# FishHighz weighting diagnostic agent prompts

Current assignment: **W10 revision 1**, proposed for user approval and dispatch.
W09 passed independent scientific review; the user requested progression and an
outline of the remaining stages. W11–W13 are high-level prospects only. They do
not authorize another assignment, profile adoption or a full forecast.

## Implementation agent

Use this prompt only when approving and dispatching the current assignment.

```text
Act as the implementation agent for FishHighz weighting diagnostic W10 r1 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. I approve W10 revision 1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and authorize only its bounded
saved-input calculations and quick checks.

Read workspace/package AGENTS.md, design section 0 and the main-roadmap handover/
progress register; IMPLEMENTATION_STEP.md is separate package context. Then read
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md, the exact W10 assignment, relevant W06,
W07 and W09 reviews, and only the required source. Check live Git state and
preserve uncommitted/untracked work and original evidence.

Answer W10's question for exactly two accuracy samples, in order: lya(lbg) bin 0
and lya(qso) bin 5. For each, first prepare fixed-reference coefficients at
magnitude orders 16/32/64 and check independent scalar sums and 0.1% refinement
stability. Only then invert its small saved auto-spectrum Fisher matrices for
reference order-32/64 stability and historical t=3/t=6 contrasts at order 32.
Anchor inputs at the reports' actual native orders (32 for bin 0, 64 for bin 5).
Do not regenerate the unavailable bin-0 order-64/t=6 diagnostic. Stop on a
numerical or reference-stability failure and report the supported finding.
A reversal of the historical improvement is a result, not a reason to change
samples or a test that the implementation must force to pass.

Use one short script, one process/thread, and the assigned 60-second cap. No
new mode sums, model or derivative evaluations, production/profile/test/plan/
manuscript changes, full forecasts, additional samples/orders, broad suites,
installations, Slurm, agent dispatch, commits or pushes. Report a production
obstruction rather than silently fixing it.

Write the specified W10 r1 handoff and small outputs with source/input identities,
policies, held scalars, independent arithmetic, residuals, refinement ratios,
BAO errors/contrast denominators if reached, and limitations/early stopping.
Distinguish newly prepared coefficients from historical Fisher matrices.
Stop for my review and acceptance. Do not proceed to W11.
```

## Independent review agent

Use this Action after the W10 implementation and handoff exist. This planning
update does not claim that W10 has been implemented.

```text
Act as the planning and independent review agent for FishHighz weighting W10 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. Follow
FISHHIGHZ_PLANNING_REVIEW_PROMPT.md and workspace/package AGENTS.md. Read design
section 0 and main-roadmap handover/progress register; IMPLEMENTATION_STEP.md is
separate package context. The active assignment is W10 r1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md. Read its handoff and relevant W06/W07/W09
reviews. Verify live source and preserve existing work/evidence.

Action: W10 r1 has been implemented. Review the bounded two-sample diagnostic.
Check identity, native-order anchoring, nonnegative measure, fixed B_star and
physical inputs, coefficient/refinement arithmetic, and the selected individual
Fisher matrices and error-ratio denominator. Rerun only the small diagnostic
and independent checks in the assignment. No new W09 mode sum or W08 campaign.

Lead the review with the scientific finding and how any proposed corrections
could affect it. Propose changes only if scientifically consequential, with the
smallest closing check; explicitly state when none are needed. A supported
instability or improvement reversal can be a valid finding, not an implementation
defect. Save a separate W10 review and revise this same assignment only if needed.
Do not silently fix code, change profiles, accept forecasts, dispatch, commit,
push or prepare W11. I control approval, dispatch, acceptance and progression.
```

Previous W09 prompts are preserved in
[their archive](lib/fishhighz/reviews/weighting-diagnostics-w09-prompts-r1.md).
