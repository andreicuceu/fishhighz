# Agent prompts for the forest-weighting sequence

Current status: W08 revision 2 passes independent scientific review and awaits
user acceptance. The positive `l_p*v` underflow finding is closed without
narrowing exact-zero or representable-subnormal behavior. No consequential
correction, profile adoption or progression is authorized. The implementation
prompt below is historical and must not be dispatched again.
Package `IMPLEMENTATION_STEP.md` is context, not this assignment.

## Planning and independent review agent

```text
Act as the planning and independent review agent for FishHighz weighting step
W08 in /global/cfs/cdirs/desicollab/users/acuceu/vega_dev. The independent package
checkout is lib/fishhighz. You are not the implementer.

Read workspace/package AGENTS.md, FISHHIGHZ_DESIGN.md section 0, the main roadmap
handover/progress register, FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md and
lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md. Read W07's report/review, the
relevant W05/W06 conventions, the W08 r1 handoff/review and archived r1
instructions, then only the source/tests needed for W08.
Preserve IMPLEMENTATION_STEP.md and all dirty/untracked historical work.

I control approval, dispatch, feedback, acceptance and progression. I selected
planning the fixed-reference option while preserving cumulative compatibility;
that decision is recorded. Ask before a new consequential choice, not to
reconfirm it. Plan only the requested step; archive replaced instructions and
synchronize status without changing package Step 13 or forecast acceptance.

Review W08's revision-2 positive-product underflow repair within the explicit
method="inverse_variance", alias=B_star contract. Check rejection of the exact
r1 reproducer, preservation of exact-zero behavior, iteration-free nu, absence
of signal/provider dependence, immutable metadata, unchanged legacy/supplied
behavior, and correct existing noise response/conversion. The method uses a
fixed supplied B_star; it does not reoptimize weights at every Fourier mode.
No default/profile changes or new auxiliary-sampling protocol are assigned.

Inspect and run only focused analytic, new-method, relevant legacy/range/noise
and tiny synthetic survey checks, plus the one saved-reference coefficient
comparison. Independently reproduce decisive coefficients; verify that survey
mean differentiation keeps weights/noise fixed and does not add provider calls.
No real forecast, broad suite, wheel, new environment or general validator.
W07's BAO improvements remain historical evidence, not a new W08 survey result.

Propose corrections only if likely to affect the scientific conclusion or its
justified scope, including actual compatibility/noise regressions. Give the
possible impact and smallest resolving check. Write
reviews/weighting-diagnostics-w08-review-r2.md, beginning with a short
"Scientific conclusion and impact of proposed changes" summary. State when
none are needed. If necessary, revise this same step using my feedback and
consequential findings; do not silently implement repairs. Stop for my review.
A pass is not acceptance, profile adoption or permission to advance.

No production edits by the reviewer, changed physical input policies, real
forecasts, package Step 13 work, subsequent step, Slurm, agent dispatch, commit
or push. Respect the numerical/test caps and exact edit scope.

Action: W08 revision 2 has been implemented. Review the exact assignment and
reviews/weighting-diagnostics-w08-r2.md. My feedback is: [insert feedback, or none].
Run the minimal independent checks, write the scientific-impact summary and
review, and revise W08 only if consequential changes are needed. Do not
implement repairs or prepare another step.
```

For planning, replace the Action with the user's exact requested step or revision.
The example review Action does not claim W08 is already implemented.

## Historical implementation agent prompt: W08 revision 2

This prompt has been executed. Preserve it as history; do not dispatch it again.

```text
Act as the implementation agent for FishHighz weighting step W08, revision 2,
in /global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz.
I approve WEIGHTING_DIAGNOSTIC_STEP.md W08 revision 2 and authorize only its
bounded positive-product underflow repair and focused validation. Preserve the
revision-1 method, legacy compatibility and all profile defaults; no real
forecast is authorized.

Read workspace/package AGENTS.md, ../../FISHHIGHZ_DESIGN.md section 0, the main
roadmap handover/progress register, ../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md
and WEIGHTING_DIAGNOSTIC_STEP.md. Read W07's report/review, relevant W05/W06
conventions, the W08 r1 handoff/review and exact archived r1 assignment. Preserve
IMPLEMENTATION_STEP.md and unrelated dirty/untracked work. W01 evidence repair
and package Step 13 are not assigned.

Retain method="inverse_variance", alias=B_star exactly as implemented in r1.
Before adding the instrumental contribution to B_star, reject an inexact
underflow of positive l_p*v. Preserve exact zero variance and zero-density
behavior. Do not broaden the extreme-input domain or change ordinary arithmetic,
legacy/supplied behavior, _iterate/_integrals numerics, immutability, context
validation, provider ownership or metadata.

Use existing ForestInput direct weight_options with no auxiliary_coordinates.
Do not change survey/forecast/noise behavior or add automatic P1D sampling.
Final noise uses mode-dependent intrinsic P1D with its own response exactly once;
B_star does not replace it. Keep weights, coefficients and noise fixed during
mean-model differentiation. No profile recipe or default changes.

Implement only the allowed weights.py guard and the focused regression in the
existing inverse-variance test. Use tiny=np.nextafter(0.0,1.0), l_p=B_star=tiny,
positive density and v=(0,1/2); revision 1 accepts w=(1,1), A=4 and
P_pixel=tiny, while revision 2 must reject the underflowed positive product.
Pair it with an exact-zero-variance control. Rerun the assigned analytic,
new-method, legacy/range/model/noise, tiny survey and saved order-32 coefficient
checks. Keep ordinary tests self-contained. No historical trajectory or Fisher
comparison needs repeating.

Use the existing interpreter, one process and one numerical thread. Respect
the 30-second saved-check and 60-second focused-pytest caps. Run Ruff on changed
Python and git diff --check. Do not broaden to a full suite, wheel/environment,
performance work or real forecast. Document any necessary out-of-scope change
as an obstruction for review instead of extending the architecture silently.

Write reviews/weighting-diagnostics-w08-r2.md and any new check evidence in
a new .validation/forest-weight-diagnostics/w08-r2-<UTC>/ directory. Preserve
r1 and review evidence. Lead with closure of the range finding. Include changes,
commands, input/source identity, independent agreement, legacy regression and
provider-call evidence, limitations and failures. Distinguish historical W07
BAO results from the new implementation checks. Do not edit planning files.

Stop for independent/user review. Do not adopt the option as a profile default,
change physical input policies, add mode-dependent weighting, run real forecasts,
finish package Step 13, write another step, dispatch agents, submit Slurm,
commit or push.
```

For approved repairs, retain the same W-step with a new revision and restrict
work to scientifically consequential findings or explicit user feedback.
Preserve prior instructions, reports and numerical evidence.
