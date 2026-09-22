# FishHighz weighting diagnostic agent prompts

Current assignment: **W11 revision 1**, proposed for approval and dispatch.
W10 passed independent scientific review and the user requested progression.
W12–W13 remain prospective stages; no reference convention/default or full run
has been approved by this planning update.

## Implementation agent

Use this prompt only when approving and dispatching the current assignment.

```text
Act as the implementation agent for FishHighz weighting diagnostic W11 r1 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. I approve W11 revision 1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and authorize only its bounded
saved-input reference-mode sensitivity calculation and quick checks.

Read workspace/package AGENTS.md, design section 0 and the main-roadmap handover/
progress register; IMPLEMENTATION_STEP.md is separate package context. Then read
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md, the exact W11 assignment, relevant W09/
W10 reviews and only the required source. Check live Git state and preserve all
existing uncommitted/untracked work and original evidence.

Use only the original accuracy bin-0 lya(qso) auto-spectrum, magnitude order 32.
Reproduce the q_star=0.00035 s/km baseline, then test 0.001 s/km; test 0.003 only
if neither BAO error changes by more than 0.1% at the first alternative. Verify
that the alternatives probe beyond the retained P1D floor and inside the saved
mode domain. Compute B_star=P1D(q_star)*W(q_star)^2 and change only the fixed
weights it prepares. Hold saved signal, mean Jacobian, mode counts, geometry,
physical P1D and response fixed. Covariance changes between tests through noise,
but is fixed during differentiation. Never use constant B_star as the physical
mode-dependent P1D, and never introduce mode-dependent weights.

Use the public chain and independent scalar checks specified in W11. Compare
Q_0 at the common baseline physical mode; use baseline errors in all sensitivity
ratios. Stop on a numerical failure or a detected change above 0.1% of either
sign. A sensitivity finding is a valid outcome, not a production defect or an
instruction to select a new reference mode. No optimization or expanded scan.

Add only the short diagnostic script, handoff and small outputs assigned. One
process/thread, 60-second cap, at most three preparations and one-spectrum mode
sums. No production/profile/test/plan/manuscript changes, new P3D or derivative,
full forecasts, other samples or refinements, broad suites, installation, Slurm,
agent dispatch, commit or push. Report an obstruction instead of silently fixing it.

Write the W11 r1 handoff with identities, held/varied quantities, scalar checks,
reference/P1D/response/coefficients, 2x2 Fisher matrices/errors, Q_0 and sensitivity
ratios, actual stopping point, omitted trials and limitations. Preserve historical
baseline evidence. Stop for my review and acceptance; do not proceed to W12.
```

## Independent review agent

Use this Action only after the W11 implementation and handoff exist.

```text
Act as the planning and independent review agent for FishHighz weighting W11 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. Follow
FISHHIGHZ_PLANNING_REVIEW_PROMPT.md and workspace/package AGENTS.md. Read design
section 0 and the main-roadmap handover/progress register; IMPLEMENTATION_STEP.md
is separate context. The active assignment is W11 r1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md. Read its handoff and W09/W10 reviews.
Verify live source and preserve existing work/evidence.

Action: W11 r1 has been implemented. Review the exact bounded reference-mode
sensitivity test. Check the baseline, sample/parameter identity, reference P1D
floor/domain, response ownership, changed B_star versus held physical mode power,
fixed weights within each run, recomputed covariance between runs, common-mode
Q_0 comparison and baseline-denominator BAO ratios. Rerun only the assigned small
calculation and independent checks. Verify the 0.1% classification and early stop.

Begin the review with the scientific finding and how any proposed corrections
could affect it. Propose changes only if scientifically consequential, with the
smallest closing check; explicitly state when none are needed. Sensitivity of
either sign is a valid finding, not a reason to change the scientific convention
automatically. Save a separate W11 review and revise the same assignment only
if needed. No production fix, extra scan, new population, profile adoption, full
forecast, dispatch, commit, push or W12 planning. I control scientific choices,
approval, dispatch, acceptance and progression.
```

Previous W10 prompts are preserved in
their archive (local-only path: `lib/fishhighz/reviews/weighting-diagnostics-w10-prompts-r1.md`).
