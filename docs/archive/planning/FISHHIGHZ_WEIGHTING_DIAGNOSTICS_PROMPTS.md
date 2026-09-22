# Retired weighting diagnostic prompts

Status: inactive as of 2026-09-17. W12 and the diagnostic plan are done and
closed by user decision; do not dispatch these historical prompts. Use the
[compatibility weighting test plan](FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md)
and [implementation roadmap](FISHHIGHZ_IMPLEMENTATION_PLAN.md) for current planning.

The original prompts are preserved below as historical text only.

---

# FishHighz weighting diagnostic agent prompts

Current assignment: **W12 revision 1**, proposed for approval and dispatch.
W11 passed review and detected reference-mode sensitivity. The user explicitly
selected retention of q_star=0.00035 s/km and planning accuracy-profile adoption.
That convention decision does not itself dispatch implementation or authorize
W13/full forecasts.

## Implementation agent

Use this prompt only when approving and dispatching the current assignment.

```text
Act as the implementation agent for FishHighz weighting diagnostic W12 r1 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. I approve W12 revision 1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and authorize only its bounded profile
integration, focused tests and saved-input local calculation. I have selected
fixed inverse-variance weighting at q_star=0.00035 s/km for accuracy adoption.

Read workspace/package AGENTS.md, design section 0 and main-roadmap handover/
progress register. IMPLEMENTATION_STEP.md is separate package context. Then read
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md, the exact W12 assignment, relevant W09/
W11 reviews, and only the required source. Check live Git state and preserve all
existing uncommitted/untracked work and historical evidence.

Connect the existing fixed-reference method to the accuracy recipe, record the
reference convention/B_star per forest, and retain explicit historical legacy
routing and compatibility behavior. Use direct reference P1D and field response;
no P3D weighting signal or auxiliary coordinates. The new fixed-reference
refinement schedule excludes iterations while retaining all other numerical
families, combined checks, thresholds and per-spectrum trial/operand binding.
Preserve version-1 legacy evidence semantics and explicitly identify the narrow
version-2 fixed-reference contract. Do not reinterpret historical failures or
reuse a cumulative cached primary as a fixed-reference calculation.

Validate the actual profile forest-input construction with focused self-contained
fixtures and a one-bin saved-input calculation selecting only FF, FG and GG
for lya(qso) and qso. Use forest magnitude orders 32/64; hold saved signal,
Jacobian, modes and galaxy noise fixed. Check independent cross-covariance,
joint/individual Fisher results and the W09 FF reference. Synthetic controller
checks and this local comparison do not establish full-survey convergence.

Follow the exact edit boundary and minimal checks. Cap focused tests and the
saved script at 60 seconds each, one process/thread. No core kernel, compatibility
recipe, physical policy, tolerance, INI/CLI, general framework, plan/manuscript or
sibling changes. No real-model controller, full forecast, broad suite, installation,
Slurm, dispatch, commit or push. Report an obstruction rather than expanding scope.

Write the W12 r1 handoff with selected convention, exact identities, changed files,
controller/legacy checks, small matrices/errors/refinement ratios, commands,
limitations and stopping point. Preserve old evidence. Stop for my review and
acceptance; do not proceed to W13.
```

## Independent review agent

Use this Action only after W12 implementation and its handoff exist.

```text
Act as the planning and independent review agent for FishHighz weighting W12 in
/global/cfs/cdirs/desicollab/users/acuceu/vega_dev. Follow
FISHHIGHZ_PLANNING_REVIEW_PROMPT.md and workspace/package AGENTS.md. Read design
section 0 and main-roadmap handover/progress register; IMPLEMENTATION_STEP.md is
separate context. The active assignment is W12 r1 in
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md and
FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md. The user has selected retention of the
baseline q_star=0.00035 s/km for planned accuracy adoption; do not reopen that
choice merely because W11 detected sensitivity. Read the handoff and W09/W11
reviews, check live source and preserve existing work/evidence.

Action: W12 r1 has been implemented. Review the exact profile adoption, metadata/
cache routing, appropriate refinement families, historical contract semantics
and one-bin cross-spectrum calculation. Verify the weight-iteration family is
inapplicable rather than fabricated as passing, while other convergence tests
and legacy failed verdicts remain intact. Rerun only focused checks and the
assigned saved-input calculation; independently reconstruct selected covariance,
one joint Fisher sum and error/refinement ratios. No full survey run.

Begin the review with the scientific finding and how proposed corrections could
affect it. Require changes only if scientifically consequential, with the smallest
closing test; explicitly state when none are needed. Local/synthetic success
must not be presented as survey convergence or validation of inherited physical
policies. Save a separate W12 review and revise the same assignment if needed.
No silent production repair, new reference scan, full forecast, dispatch, commit,
push or W13 planning. I retain implementation acceptance and execution decisions.
```

Previous W11 prompts are preserved in
[their archive](../reviews/weighting-diagnostics-w11-prompts-r1.md).
