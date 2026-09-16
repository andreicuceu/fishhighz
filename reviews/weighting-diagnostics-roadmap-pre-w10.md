# FishHighz forest-weighting diagnostic roadmap

Status: W08 revision 2 passed independent scientific review; the user requested
progression to W09. The positive `l_p*v` underflow finding is closed, with
ordinary coefficients unchanged. **W09 revision 1 is proposed for approval and
dispatch:** replay one saved QSO-forest auto-spectrum/bin through the explicit
fixed-reference implementation, noise, covariance and Fisher calculation.
Profile defaults, W01 evidence repair, package Step 13 and forecast acceptance
remain unchanged.
Created: 2026-09-15. Last updated: 2026-09-16.

## Purpose and current scope

Understand the inherited forest-weight sensitivity through individually requested,
controlled comparisons. W01 established finite-update sensitivity. W02 identified
prefix feedback and absence of aliasing in the weighting signal. W03 found larger
absolute t=3/6 drift with full-sample feedback; W04 showed both full-sample branches
satisfy 0.1% coefficient stability over t=6/12/24. W05 identifies the explicit
inverse-variance reference with the legacy seed. On its signed grid, prefix t=3
and full-sample fixed/refreshed t=24 have 14.3618%, 12.8358% and 19.1309% higher
combined noise Q than the reference. This signed-grid ordering establishes no
physical optimum or forecast improvement.

Step 12 r5 already independently verified magnitude stability of the fixed seed.
W06 connects that evidence to the explicit reference in one nonnegative accuracy
QSO-forest bin. The reference is stable to roundoff over orders 16/32/64;
cumulative t=3 passes the 0.1% coefficient criterion, while t=6 fails only
P_pixel at 64/16 (+0.129453%). At fixed order, increasing the count from 3 to 6
raises Q by 1.7039–1.7533%. All cumulative Q values exceed the conditional
reference minimum. This establishes neither nonconvergence nor a continuum
limit or a forecast improvement.

W07 connects the coefficient diagnosis to one QSO-forest auto-spectrum/bin:
relative to the fixed reference, three cumulative updates raise radial/transverse
BAO errors by 7.63057%/10.32563%, and six by 8.61956%/11.70862%. The reference
order-32/64 check passes to roundoff. These are analyses of historical Fisher
matrices, not new forecasts or a general multi-mode optimum.

The user chose the fixed-reference option for W08 planning after being offered
that option, a mode-dependent comparison, or further cumulative-rule diagnosis.
This selects an explicit opt-in implementation scope; it does not adopt a new
profile default, validate the retained physical input policies or accept the
unconverged package forecasts.

This investigation is separate from the [package roadmap](FISHHIGHZ_IMPLEMENTATION_PLAN.md).
Package Step 13, its pending review, and Step 12 scientific acceptance remain
unchanged. The user controls approval, dispatch, acceptance and progression.
Progression to W09 does not retrospectively close W01's outstanding findings
or adopt a profile default.

Only completed work and the current requested step are listed here. Add another
step only when the user requests it; no future sequence is prescribed.

## Active assignment

**W09 — fixed-reference replay of one saved BAO forecast, revision 1.**
Reconstruct only accuracy bin 0, `lya(qso)` auto-spectrum, magnitude order 32.
W07 inverted saved Fisher matrices; W08 verified coefficients and synthetic
integration. W09 checks the implemented reference through a new mode sum on
saved signal, derivatives and mode counts. First reconstruct the historical
six-update noise/Fisher result as a convention control, then compare the public
reference calculation with independent scalar formulas and the W07 matrix.
Reuse historical three/six-update matrices for the quoted BAO contrasts.

The detailed contract, minimal acceptance tests and stopping conditions are in
[WEIGHTING_DIAGNOSTIC_STEP.md](lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md).
[Agent prompts](FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PROMPTS.md) select W09. Approval
and dispatch remain with the user. One small saved-input script suffices; no
production edit, profile adoption, other spectrum/bin, full survey run or
repeated magnitude-refinement experiment is assigned. This checks numerical
consistency, not physical input validity or a general multi-mode optimum.

## Scientific conventions and minimum necessary work

- Preserve units, response ownership, independent P3D/P1D inputs and the
  difference between 1D pixel-noise power and effective 3D noise.
- Separate published equations, executable behavior, code comments and inference.
  Do not claim that a cumulative implementation is required by a paper merely
  because the code cites it, or that ancestry proves historical behavior.
- Separate finite-update sensitivity, weight-amplitude decay, normalized noise
  limits and continuum convergence. A tiny example supports its stated scope.
- Preserve literal signed legacy inputs as historical evidence. Positive-measure
  arguments do not apply automatically; flooring or changing the scientific
  prescription requires a user decision.
- Implementation and review use only the minimum code/tests needed to answer
  the active scientific question. Prefer a derivation and one short script over
  new infrastructure. No routine full suite, wheel, fresh environment, evidence
  schema or provenance-validator repair is required for these diagnostics.
- Reviewers propose further changes only if they are likely to affect the
  scientific conclusion or its justified scope. Each proposed change must state
  the scientific consequence and the smallest check that resolves it. Lead the
  review with a short summary of how the proposed changes could affect the
  conclusion; explicitly state when none are needed.
- Sufficient source/input identification and independent arithmetic remain
  necessary to assess a scientific claim. Resolve concrete identity or numerical
  ambiguities; do not add hypothetical hardening requirements unrelated to it.

The user's prior interest in both strictly indirect forecast effects and
separately labelled frozen-weight signal/response effects is retained. Neither
is assigned a future step or authorized for execution by this document.

## Plan / implementation / independent review cycle

1. The planning agent writes only the requested step/revision, including its
   question, scientific conventions, minimal independent checks and stopping
   condition. Archive replaced instructions and synchronize current status.
2. The user approves the concrete assignment and separately dispatches an
   implementer. The implementer performs only that assignment and writes a
   concise report with the equations, relevant source identity and actual checks.
3. The independent reviewer checks the decisive source and calculations and
   writes a separate report. Propose repairs only under the scientific-impact
   criterion above; a restricted conclusion can be sufficient when uncertainty
   is documented. Do not silently fix code or generate optional repair lists.
4. If consequential repairs are needed, revise the same step with the minimal
   check for each finding, preserving earlier reports. The user controls
   acceptance and any request for another step. A passing review does not
   automatically select an estimator or advance the investigation.

Use reports `lib/fishhighz/reviews/weighting-diagnostics-wNN-rR.md` and
`weighting-diagnostics-wNN-review-rR.md`. Any small numerical outputs go in a
new `lib/fishhighz/.validation/forest-weight-diagnostics/wNN-rR-<UTC>/` directory.
Preserve historical evidence and distinguish it from newly performed checks.

No production changes, NewForecast/reference recapture, full profile controller,
joint 15x2pt forecast, Slurm action, agent dispatch, commit or push is authorized
by planning or reviewing this diagnostic. Future real evaluations require
explicit scope in a user-approved assignment; full forecasts require an explicit
user request.

## Progress register

| Step/revision | State | Scientific result and remaining scope |
| --- | --- | --- |
| W01 r1/r2 | Implemented and independently reviewed; scientific result supported | The original 107-node signed legacy example changes between 3 and 6 updates: A by 0.8616%, P_pixel by 1.1781%, and both one-spectrum BAO errors by more than 0.1%. This establishes finite-update sensitivity for that example. [Review](lib/fishhighz/reviews/weighting-diagnostics-w01-review-r2.md). |
| W01 r3 | Evidence-only repair deferred, not implemented or accepted | W01-R4/R5 remain recorded in the historical review. Their proposed repair is archived; it is not a prerequisite for the user-requested W02 audit. |
| W02 r1 | Independent review passes; user requested progression to W03 | Prefix versus full-sample updates differ; live weighting S has no aliasing contribution. Tiny independent checks pass. Historical C++ and aliasing-signal iteration details remain unresolved; no survey conclusion or prescription adopted. [Review](lib/fishhighz/reviews/weighting-diagnostics-w02-review-r1.md). |
| W03 r1 | Independent scientific review passes; user requested progression to W04 | At fixed saved inputs and S/B, full-sample feedback increases absolute t=3/6 sensitivity: A changes by -0.9196% and P_pixel by +2.0171%, versus +0.8616% and +1.1781% for prefixes. Independent 80-digit checks pass; conditional stop at t=6. No asymptotic or production conclusion. [Review](lib/fishhighz/reviews/weighting-diagnostics-w03-review-r1.md). |
| W04 r1 | Independent scientific review passes; user requested progression to W05 | Both full-sample branches satisfy 0.1% coefficient stability over t=6/12/24; aliasing is not necessary for this bounded stability. At t=24 refreshed aliasing lowers A by 5.08747% and raises P_pixel by 12.60455%. Independent 80-digit replay/contrast checks pass; no consequential correction required. Fixed signed grid only, no asymptotic or production conclusion. [Review](lib/fishhighz/reviews/weighting-diagnostics-w04-review-r1.md). |
| W05 r1 | Independent scientific review passes; user requested progression to W06 | Reference equals the legacy seed; positive-measure minimum verified. Signed-grid Q is higher by 14.3618%, 12.8358%, 19.1309% for prefix t=3 and fixed/refreshed t=24. Refreshed Q exceeds fixed by 5.57897%. No consequential correction; no physical optimum or production adoption. [Review](lib/fishhighz/reviews/weighting-diagnostics-w05-review-r1.md). |
| W06 r1 | Independent scientific review passes; user requested progression to W07 | nu and t=3 meet the 0.1% coefficient refinement criterion; t=6 fails only 64/16 P_pixel (+0.129453%). Fixed-order t=3/6 Q increases 1.703884–1.753327%; all cumulative Q exceed the conditional reference minimum. No consequential correction or continuum claim. [Review](lib/fishhighz/reviews/weighting-diagnostics-w06-review-r1.md). No production adoption. |
| W07 r1 | Independent scientific review passes; user requested progression to W08 | Relative to the reference, order-32 t=3 raises radial/transverse BAO errors by 7.63057%/10.32563%, and t=6 by 8.61956%/11.70862%. The reference order-32/64 check passes to roundoff. No consequential correction; one spectrum/bin only, no new mode sum or production adoption. [Review](lib/fishhighz/reviews/weighting-diagnostics-w07-review-r1.md). |
| W08 r1 | Implemented; one range correction required after independent review | The explicit fixed-B_star method reproduces exact rational, supplied/legacy-seed, noise, survey and saved W06 coefficients. One positive subnormal `l_p*v` product underflows to zero without rejection while final coefficients remain representable, so the claimed range policy is incomplete. No ordinary/saved coefficient or compatibility/noise regression was found. [Review](lib/fishhighz/reviews/weighting-diagnostics-w08-review-r1.md). |
| W08 r2 | Independent scientific review passes; user requested progression to W09 | The W08-R1 reproducer now rejects with field context. Exact-zero and exactly representable subnormal controls pass; 111 focused tests and independent rational/saved-reference coefficient reconstructions pass. No consequential correction, compatibility/noise regression, new forecast, profile adoption or expanded scientific claim. [Review](lib/fishhighz/reviews/weighting-diagnostics-w08-review-r2.md). |
| W09 r1 | Proposed for user approval and dispatch | One saved accuracy QSO-forest auto-spectrum/bin: public fixed-reference weights → noise → covariance → Fisher, checked against scalar formulas and historical W07 results. No calculation performed during planning, production change or default adoption. |

## Historical documents

- [W01 revision-3 instructions](lib/fishhighz/reviews/weighting-diagnostics-w01-instructions-r3.md)
  preserve the former current-step file byte for byte.
- [Pre-W02 roadmap](lib/fishhighz/reviews/weighting-diagnostics-roadmap-pre-w02.md)
  preserves the previous conditional sequence as historical evidence. Its future
  rows are superseded and are not the active plan.

- [Original W02 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w02-instructions-r1.md)
  and [reviewed W02 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w02-instructions-reviewed-r1.md)
  preserve its original assignment and the status after independent review.

- [Original W03 r1 assignment](lib/fishhighz/reviews/weighting-diagnostics-w03-instructions-r1.md)
  and [reviewed W03 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w03-instructions-reviewed-r1.md)
  preserve the scientific assignment and the subsequent review status.

- [Original W04 r1 assignment](lib/fishhighz/reviews/weighting-diagnostics-w04-instructions-r1.md)
  and [reviewed W04 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w04-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review.

- [Original W05 r1 assignment](lib/fishhighz/reviews/weighting-diagnostics-w05-instructions-r1.md)
  and [reviewed W05 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w05-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review.

- [Original W06 r1 assignment](lib/fishhighz/reviews/weighting-diagnostics-w06-instructions-r1.md)
  and [reviewed W06 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w06-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review. No repair
  revision was required.

- [Original W07 r1 assignment](lib/fishhighz/reviews/weighting-diagnostics-w07-instructions-r1.md)
  and [reviewed W07 r1 instructions](lib/fishhighz/reviews/weighting-diagnostics-w07-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review. No repair
  revision was required.

- [Original W08 r1 assignment](lib/fishhighz/reviews/weighting-diagnostics-w08-instructions-r1.md)
  preserves the exact implemented assignment. Its independent review requires
  only the revision-2 positive-product underflow repair.

- [Exact W08 r2 assignment](lib/fishhighz/reviews/weighting-diagnostics-w08-instructions-r2.md)
  preserves the implemented repair instructions. The
  [revision-2 review](lib/fishhighz/reviews/weighting-diagnostics-w08-review-r2.md)
  closes W08-R1 without requesting another revision.

- [Reviewed W08 r2 instructions](lib/fishhighz/reviews/weighting-diagnostics-w08-instructions-reviewed-r2.md)
  and [W08 prompts](lib/fishhighz/reviews/weighting-diagnostics-w08-prompts-reviewed-r2.md)
  preserve the pre-W09 assignment/status and prompts byte for byte. The prompts
  retain links relative to their former workspace-root location.

Exact instruction archives retain relative links from their original location
`lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md`; the pre-W02 roadmap retains links
from the workspace root. Use current documents for navigation and instructions.
Historical wait-for-progression statements describe their earlier scope. The
user has now requested progression after W08 review. Only W09 planning is
complete; approval/dispatch, profile adoption and forecast acceptance remain
separate user decisions.
