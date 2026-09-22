# FishHighz weighting diagnostic agent prompts

Current assignment: **W09 revision 1**, proposed for user approval and dispatch.
W08 revision 2 passed independent scientific review; the user requested
progression. W09 replays one saved QSO-forest auto-spectrum/bin through the
implemented fixed-reference weights, noise, covariance and Fisher calculation.
It does not adopt a profile or reopen package Step 13.

## Implementation agent

Use this prompt only when approving and dispatching the current assignment.

```text
Act as the implementation agent for FishHighz weighting diagnostic W09 r1 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. I approve W09 revision 1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and authorize its bounded saved-input
calculation and quick checks.

Read the workspace/package AGENTS.md, design section 0 and main-roadmap handover
and progress register. IMPLEMENTATION_STEP.md is separate package context.
Then read FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md, the exact W09 assignment,
W07 and W08 r2 reviews, and only the relevant source. Verify live Git state and
preserve all uncommitted/untracked work and original evidence.

Answer only W09's scientific question: does the new public fixed-reference
calculation reproduce the saved W07 BAO result for accuracy bin 0, lya(qso)
auto-spectrum, magnitude order 32? Use the specified saved signal, derivatives
and mode counts. Execute the historical control first, then the reference
noise and independent scalar covariance/Fisher checks. Stop on a consequential
failure; do not broaden the calculation. Explain any diagnostic geometry adapter
and never use its artificial volume. Distinguish newly summed reference results
from historical cumulative comparators.

Use the minimal script and checks in the assignment, one process/thread and
the 60-second cap. Do not change production code, profiles, tests, other packages,
planning documents or the scientific note. Do not run full forecasts, other bins
or spectra, new integrations or derivatives, a broad suite, installations,
Slurm, commits, pushes or other agents. Report a production obstruction rather
than silently repairing it.

Write the specified W09 r1 handoff and small numerical outputs with source/input
identities, independent formulas, residuals, errors, contrast denominators and
limitations. Stop for my review and acceptance; do not advance the sequence.
```

## Independent review agent

The following Action is for use after the implementation/handoff exists; this
planning update makes no claim that W09 has already been implemented.

```text
Act as the planning and independent review agent for FishHighz weighting W09 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. Follow
FISHHIGHZ_PLANNING_REVIEW_PROMPT.md and the workspace/package AGENTS.md.
Read design section 0, the main-roadmap handover/progress register and separate
package IMPLEMENTATION_STEP.md as context; the active assignment is
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md and
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md, W09 r1. Read its handoff and relevant
W07/W08 reviews. Verify live source and preserve existing work/evidence.

Action: W09 r1 has been implemented. Review the assigned one-spectrum/bin replay,
its historical control and independent numerical checks. Check response/unit
ownership, fixed weights/covariance, sample and parameter identity, scalar auto
covariance and Fisher normalization, and the percentage denominator. Use only
minimal source and checks needed to decide this scientific question. No full
forecast, repeated W08 test campaign, production repair or agent dispatch.

Begin the review with the scientific finding and a short explanation of how
any proposed changes could affect that conclusion. Propose changes only if
scientifically consequential; give the smallest check closing each finding.
Explicitly say when none are needed. Save a separate W09 review and update the
same step's instructions only if a consequential revision is needed. Do not
silently fix code, change profiles, accept forecasts, commit/push, or plan W10.
I control approval, dispatch, acceptance and progression.
```

Historical W08 prompts are preserved in
their archive (local-only path: `lib/fishhighz/reviews/weighting-diagnostics-w08-prompts-reviewed-r2.md`).
