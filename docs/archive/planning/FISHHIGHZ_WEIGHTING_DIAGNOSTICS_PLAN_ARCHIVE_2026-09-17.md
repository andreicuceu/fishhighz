# FishHighz forest-weighting diagnostic roadmap

Status: W11 revision 1 passes independent scientific review and detects
reference-mode sensitivity. The user explicitly chose **retain the baseline
q_star=0.00035 s/km and plan accuracy-profile adoption**. W12 revision 1 is
proposed for approval/dispatch: connect fixed-reference weighting, make the
iteration control inapplicable, and validate one saved Lyalpha–QSO three-spectrum
bin. W13 remains a separately requested broader-validation stage. No profile
code has changed during planning; W01 repairs, package Step 13 and historical
forecast acceptance remain unchanged.
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
Progression to W12 does not retrospectively close W01's outstanding findings
or accept historical forecasts. The user has now selected the baseline reference
convention for planned accuracy-profile adoption; its implementation remains
subject to the W12 approval/dispatch and review cycle.

The user requested the remaining sequence at high level. Only W12 has current
detailed instructions; later assignments depend on the preceding findings and
separate user requests. Listing a stage is not approval or dispatch.

## Active assignment

**W12 — adopt the declared fixed reference in the accuracy profile, revision 1.**
Following W11, the user selected “Retain baseline and plan adoption”: use fixed
inverse-variance weights with q_star=0.00035 s/km explicitly recorded, retaining
cumulative weighting for compatibility. This is a modelling convention, not an
optimized reference or a claim of reference-mode independence.

W11 found radial/transverse error changes of -0.00147679%/-0.0000622378% at
0.001 s/km and +0.116214%/+0.220774% at 0.003 s/km relative to the baseline,
with all other inputs fixed. Its [review](../reviews/weighting-diagnostics-w11-review-r1.md)
passes without a consequential correction. Both larger-mode changes exceed the
0.1% screen; the user's explicit decision resolves the prerequisite for planning
adoption. No further reference-mode scan is assigned.

The live accuracy recipe and refinement controller still implement cumulative
iterations. W12 changes only their necessary method routing and applicable
refinement families: new fixed-reference studies omit iteration trials, while
legacy studies and historical contracts retain their original semantics. Bind
the new method/convention explicitly so an old cumulative result cannot be
relabelled or reused as the new calculation. Other convergence criteria remain.

Local numerical validation uses only saved bin-0 QSO-forest/QSO inputs and
FF, FG, GG spectra at forest magnitude orders 32 and 64, with fixed signal,
Jacobian, mode counts and quasar noise. Check the new profile's forest preparation,
cross-spectrum covariance, independent joint/individual Fisher sums and bounded
magnitude refinement. No real-model/full-profile controller run is assigned.

The [current instructions](../notes/WEIGHTING_DIAGNOSTIC_STEP.md) specify the
narrow edit scope, scalar/matrix tests, historical compatibility checks, and
handoff. [Agent prompts](FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PROMPTS.md) select W12.
Approval, dispatch, implementation acceptance and full-run authorization remain
separate user actions. W10's two-sample findings and weak-LBG limitations remain
in its [review](../reviews/weighting-diagnostics-w10-review-r1.md).

## Remaining sequence — user-requested outline

These stages implement the proposed continuation. Plan and review each separately;
stop or revise the sequence if a scientific finding changes what is needed.
Later entries deliberately omit detailed implementation contracts.

| Stage | Scientific question and smallest intended scope | Decision before proceeding |
| --- | --- | --- |
| **W10: population/redshift extension — reviewed** | Two additional samples pass bounded reference refinement checks and retain smaller reference BAO errors than the historical cumulative comparison. | User requested progression; no broader physical/optimality claim. |
| **W11: reference-mode sensitivity — reviewed** | At 0.003 s/km, errors rise by 0.116214%/0.220774%; the 0.001 point is below the screen. No numerical failure or optimum claim. | User explicitly chose to retain q_star=0.00035 s/km and plan adoption. |
| **W12: accuracy adoption/local validation — current proposal** | Connect the selected fixed reference, remove the inapplicable iteration-convergence family with historical semantics preserved, and check one saved FF/FG/GG bin at two magnitude orders. Detailed assignment above. | Convention selected; plan approval/dispatch and acceptance remain with the user. Physical policies and full forecast authorization are unchanged. |
| **W13: return to package validation** | After local validation, reassess the Step 12 convergence criteria and perform the explicitly requested broader forecast comparison. Distinguish new fixed-reference results from historical cumulative forecasts; record which acceptance criteria are met. | A real full suite requires an explicit execution request, followed by independent review and user scientific acceptance. This outline does not authorize that run or close package Step 13. |

W13 here is a **weighting-sequence label**, distinct from package Step 13's
cumulative-limit investigation. The two records must not be conflated.

### Separate or conditional work

- Indirect-only profile comparisons and separately labelled frozen-weight
  signal/response comparisons remain unperformed. They can attribute remaining
  compatibility/accuracy differences if needed; they are not prerequisites for
  W12 or automatic additions to its execution.
- W01's deferred evidence repairs and package Step 13's pending review/unresolved
  continuum question remain recorded. A useful replacement does not resolve
  the mathematical limit of the inherited recurrence.
- Physical density/SNR policies remain an independent modelling question. Any
  change needs an explicit decision; numerical stability cannot validate them.
- After diagnostic review, documentation may be updated in a separately requested
  task. Package performance/distribution work remains in the main roadmap.

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
separately labelled frozen-weight signal/response effects is retained. These
remain conditional attribution work, outside the active assignment.

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
| W01 r1/r2 | Implemented and independently reviewed; scientific result supported | The original 107-node signed legacy example changes between 3 and 6 updates: A by 0.8616%, P_pixel by 1.1781%, and both one-spectrum BAO errors by more than 0.1%. This establishes finite-update sensitivity for that example. [Review](../reviews/weighting-diagnostics-w01-review-r2.md). |
| W01 r3 | Evidence-only repair deferred, not implemented or accepted | W01-R4/R5 remain recorded in the historical review. Their proposed repair is archived; it is not a prerequisite for the user-requested W02 audit. |
| W02 r1 | Independent review passes; user requested progression to W03 | Prefix versus full-sample updates differ; live weighting S has no aliasing contribution. Tiny independent checks pass. Historical C++ and aliasing-signal iteration details remain unresolved; no survey conclusion or prescription adopted. [Review](../reviews/weighting-diagnostics-w02-review-r1.md). |
| W03 r1 | Independent scientific review passes; user requested progression to W04 | At fixed saved inputs and S/B, full-sample feedback increases absolute t=3/6 sensitivity: A changes by -0.9196% and P_pixel by +2.0171%, versus +0.8616% and +1.1781% for prefixes. Independent 80-digit checks pass; conditional stop at t=6. No asymptotic or production conclusion. [Review](../reviews/weighting-diagnostics-w03-review-r1.md). |
| W04 r1 | Independent scientific review passes; user requested progression to W05 | Both full-sample branches satisfy 0.1% coefficient stability over t=6/12/24; aliasing is not necessary for this bounded stability. At t=24 refreshed aliasing lowers A by 5.08747% and raises P_pixel by 12.60455%. Independent 80-digit replay/contrast checks pass; no consequential correction required. Fixed signed grid only, no asymptotic or production conclusion. [Review](../reviews/weighting-diagnostics-w04-review-r1.md). |
| W05 r1 | Independent scientific review passes; user requested progression to W06 | Reference equals the legacy seed; positive-measure minimum verified. Signed-grid Q is higher by 14.3618%, 12.8358%, 19.1309% for prefix t=3 and fixed/refreshed t=24. Refreshed Q exceeds fixed by 5.57897%. No consequential correction; no physical optimum or production adoption. [Review](../reviews/weighting-diagnostics-w05-review-r1.md). |
| W06 r1 | Independent scientific review passes; user requested progression to W07 | nu and t=3 meet the 0.1% coefficient refinement criterion; t=6 fails only 64/16 P_pixel (+0.129453%). Fixed-order t=3/6 Q increases 1.703884–1.753327%; all cumulative Q exceed the conditional reference minimum. No consequential correction or continuum claim. [Review](../reviews/weighting-diagnostics-w06-review-r1.md). No production adoption. |
| W07 r1 | Independent scientific review passes; user requested progression to W08 | Relative to the reference, order-32 t=3 raises radial/transverse BAO errors by 7.63057%/10.32563%, and t=6 by 8.61956%/11.70862%. The reference order-32/64 check passes to roundoff. No consequential correction; one spectrum/bin only, no new mode sum or production adoption. [Review](../reviews/weighting-diagnostics-w07-review-r1.md). |
| W08 r1 | Implemented; one range correction required after independent review | The explicit fixed-B_star method reproduces exact rational, supplied/legacy-seed, noise, survey and saved W06 coefficients. One positive subnormal `l_p*v` product underflows to zero without rejection while final coefficients remain representable, so the claimed range policy is incomplete. No ordinary/saved coefficient or compatibility/noise regression was found. [Review](../reviews/weighting-diagnostics-w08-review-r1.md). |
| W08 r2 | Independent scientific review passes; user requested progression to W09 | The W08-R1 reproducer now rejects with field context. Exact-zero and exactly representable subnormal controls pass; 111 focused tests and independent rational/saved-reference coefficient reconstructions pass. No consequential correction, compatibility/noise regression, new forecast, profile adoption or expanded scientific claim. [Review](../reviews/weighting-diagnostics-w08-review-r2.md). |
| W09 r1 | Independent scientific review passes; user requested progression to W10 | The new public fixed-reference mode sum reproduces W07 for the original QSO bin. Independent primary residuals are below 1e-14; t=3/t=6 radial/transverse error contrasts remain +7.63057%/+10.32563% and +8.61956%/+11.70862%. No consequential correction or default adoption. [Review](../reviews/weighting-diagnostics-w09-review-r1.md). |
| W10 r1 | Independent scientific review passes; user requested progression to W11 | Both added samples have roundoff-level reference refinement stability. Cumulative t=3/t=6 error increases are about 116–118%/250–254% for bin-0 LBG and 6.1–6.7%/8.1–8.9% for bin-5 QSO; weak-LBG Fisher limitation retained. No consequential correction. [Review](../reviews/weighting-diagnostics-w10-review-r1.md). |
| W11 r1 | Independent scientific review passes; user selected baseline convention and requested progression | Sensitivity detected at 0.003 s/km: +0.116214% radial/+0.220774% transverse. At 0.001, changes are -0.00147679%/-0.0000622378%. Independent reconstruction finds no consequential correction. [Review](../reviews/weighting-diagnostics-w11-review-r1.md). |
| W12 r1 | Proposed for user approval/dispatch | Fixed-reference accuracy-profile adoption at the selected q_star=0.00035 s/km, appropriate numerical controls, and one saved three-spectrum local check. No implementation or diagnostic run during planning. |

## Historical documents

- [W11 r1 assignment](../reviews/weighting-diagnostics-w11-instructions-r1.md)
  and [W11 prompts](../reviews/weighting-diagnostics-w11-prompts-r1.md)
  preserve the replaced documents byte for byte. Their pending-choice statements
  describe the earlier stage; the user has now selected retention of the baseline.

- [W10 r1 assignment](../reviews/weighting-diagnostics-w10-instructions-r1.md)
  and [W10 prompts](../reviews/weighting-diagnostics-w10-prompts-r1.md)
  preserve the replaced documents byte for byte. Their proposed status is
  superseded by the W10 review/current register; prompts retain root-relative links.

- [W09 r1 assignment](../reviews/weighting-diagnostics-w09-instructions-r1.md),
  [W09 prompts](../reviews/weighting-diagnostics-w09-prompts-r1.md), and
  [pre-W10 roadmap](../reviews/weighting-diagnostics-roadmap-pre-w10.md)
  preserve the prior documents byte for byte. W09's original assignment still
  says proposed; its review and the current register give its later disposition.
  Prompts/roadmaps retain relative links from the workspace root.

- [W01 revision-3 instructions](../reviews/weighting-diagnostics-w01-instructions-r3.md)
  preserve the former current-step file byte for byte.
- [Pre-W02 roadmap](../reviews/weighting-diagnostics-roadmap-pre-w02.md)
  preserves the previous conditional sequence as historical evidence. Its future
  rows are superseded and are not the active plan.

- [Original W02 r1 instructions](../reviews/weighting-diagnostics-w02-instructions-r1.md)
  and [reviewed W02 r1 instructions](../reviews/weighting-diagnostics-w02-instructions-reviewed-r1.md)
  preserve its original assignment and the status after independent review.

- [Original W03 r1 assignment](../reviews/weighting-diagnostics-w03-instructions-r1.md)
  and [reviewed W03 r1 instructions](../reviews/weighting-diagnostics-w03-instructions-reviewed-r1.md)
  preserve the scientific assignment and the subsequent review status.

- [Original W04 r1 assignment](../reviews/weighting-diagnostics-w04-instructions-r1.md)
  and [reviewed W04 r1 instructions](../reviews/weighting-diagnostics-w04-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review.

- [Original W05 r1 assignment](../reviews/weighting-diagnostics-w05-instructions-r1.md)
  and [reviewed W05 r1 instructions](../reviews/weighting-diagnostics-w05-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review.

- [Original W06 r1 assignment](../reviews/weighting-diagnostics-w06-instructions-r1.md)
  and [reviewed W06 r1 instructions](../reviews/weighting-diagnostics-w06-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review. No repair
  revision was required.

- [Original W07 r1 assignment](../reviews/weighting-diagnostics-w07-instructions-r1.md)
  and [reviewed W07 r1 instructions](../reviews/weighting-diagnostics-w07-instructions-reviewed-r1.md)
  preserve the assignment and its status after independent review. No repair
  revision was required.

- [Original W08 r1 assignment](../reviews/weighting-diagnostics-w08-instructions-r1.md)
  preserves the exact implemented assignment. Its independent review requires
  only the revision-2 positive-product underflow repair.

- [Exact W08 r2 assignment](../reviews/weighting-diagnostics-w08-instructions-r2.md)
  preserves the implemented repair instructions. The
  [revision-2 review](../reviews/weighting-diagnostics-w08-review-r2.md)
  closes W08-R1 without requesting another revision.

- [Reviewed W08 r2 instructions](../reviews/weighting-diagnostics-w08-instructions-reviewed-r2.md)
  and [W08 prompts](../reviews/weighting-diagnostics-w08-prompts-reviewed-r2.md)
  preserve the pre-W09 assignment/status and prompts byte for byte. The prompts
  retain links relative to their former workspace-root location.

Exact instruction archives retain relative links from their original location
`lib/fishhighz/WEIGHTING_DIAGNOSTIC_STEP.md`; the pre-W02 roadmap retains links
from the workspace root. Use current documents for navigation and instructions.
Historical wait-for-progression statements describe their earlier scope. The
user has now requested progression after W11 review and selected the baseline
convention for planned adoption. W12 planning is complete; approval/dispatch,
implementation acceptance, later assignments and full-run/forecast acceptance
remain separate user decisions.
