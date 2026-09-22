# FishHighz implementation roadmap

Status: Steps 01–11 accepted for progression; Step 12 revision 6 evidence repair passed review (R2 resolved). Step 13 revision 2 is implemented with a handoff awaiting independent/user review; closure of the revision-1 findings is not yet independently verified. Its handoff classifies the continuum limit as unresolved for both forest populations in all six saved bins. Step 12 R3 and scientific acceptance remain open; no new real forecast is authorized.
Created: 2026-09-10. Last updated: 2026-09-16.

This document describes the implementation sequence and review process at a high
level. Scientific requirements and conventions live in
[FISHHIGHZ_DESIGN.md](FISHHIGHZ_DESIGN.md). Only the currently proposed step is
specified in detail in
[lib/fishhighz/IMPLEMENTATION_STEP.md](../notes/IMPLEMENTATION_STEP.md).
Standing guidance for builders lives in
lib/fishhighz/AGENTS.md (local-only path: `lib/fishhighz/AGENTS.md`).

## Fresh-agent handover

Parallel weighting update, 2026-09-16: W07 revision 1
[passes independent scientific review](../reviews/weighting-diagnostics-w07-review-r1.md).
For one accuracy QSO-forest auto-spectrum/bin, three cumulative updates raise
radial/transverse BAO errors by 7.63057%/10.32563% relative to the fixed reference;
six updates raise them by 8.61956%/11.70862%. Reference order-32/64 refinement
is stable to roundoff. No consequential correction is required; these are
historical forecast comparisons, not a general multi-mode optimum.
The user requested progression and selected the fixed-reference option. W08
revision 1 implements explicit `method="inverse_variance", alias=B_star`
preparation with no signal or iterations. Independent review reproduces the
ordinary analytic, fixed-noise/survey and saved-reference results, but finds
one range-safety gap: positive `l_p*v` can underflow to zero without rejection
while final coefficients remain representable. W08 revision 2 rejects only that
product underflow and passes independent scientific review. Exact zeros and an
exactly representable subnormal product remain valid; ordinary analytic,
legacy/supplied, fixed-noise/survey and saved-reference checks pass. No
consequential correction is needed. Direct scalar B_star is already
response-smoothed; automatic auxiliary sampling and per-mode weights remain
outside this bounded change.

W09 now passes independent scientific review: public fixed-reference noise,
covariance and Fisher calculations reproduce W07 for the original QSO bin,
with primary independent residuals below 1e-14. W10 also passes independent
scientific review: bin-0 LBG and bin-5 QSO reference coefficients/errors are
stable to roundoff under the assigned refinements, and historical cumulative
errors are larger. The weak LBG constraints retain the local-Fisher limitation.
W11 passes independent scientific review and detects reference-mode sensitivity:
changing q_star from 0.00035 to 0.003 s/km raises radial/transverse errors by
0.116214%/0.220774% in the original QSO bin; the 0.001 test is below the screen.
No numerical correction is required. The user explicitly selected retention of
q_star=0.00035 s/km and planning accuracy-profile adoption. W12 revision 1 is
proposed for approval/dispatch: fixed-reference recipe routing, an applicable
refinement schedule without iteration trials, preserved legacy evidence, and one
saved bin with FF/FG/GG spectra at two forest magnitude orders. No optimization,
physical-policy change or accuracy-tolerance change is selected.
The [diagnostic roadmap](FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md) retains W13
as separately requested broader validation. Only W12 has current detailed
instructions. No code or numerical result changed during this planning update;
the chosen convention awaits implementation and review. [Agent prompts](FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PROMPTS.md)
retain minimal work and scientific-impact-based review. W01 evidence repairs,
package Step 13 and scientific acceptance of the failed forecasts remain unchanged.

Overview update, 2026-09-15: read the current [Step 13 revision-2 handoff](../reviews/step-13.md)
alongside the r1 independent review and r2 instructions. Its eight saved
source/report hashes match the live files. The implementation reports 1318 tests
passed and 25 optional skips; no tests or forecasts were rerun for this overview.
Revision-2 independent review remains pending. This update supersedes the
proposal/dispatch status in the dated snapshot below, without accepting the
implementation or changing its scientific scope.

This section supplies working context without requiring the preceding chat.
Snapshot date: 2026-09-14. Recheck live files/Git state at the start of each
assignment; historical test counts and artifact paths are evidence pointers,
not a claim that a changed checkout has passed.

### Read order and document ownership

1. Read workspace AGENTS.md and lib/fishhighz/AGENTS.md.
2. Read the status/progress register here and design section 0 for the actual
   implemented API map and remaining design choices. Read the scientific design
   sections relevant to the requested step before drafting its instructions.
3. Read lib/fishhighz/IMPLEMENTATION_STEP.md for Step 13 revision 2, then
   reviews/step-13-review-r1.md and the preserved r1 handoff/derivation. Revision 1
   requires scientific-verdict, source/attempt-binding and scalar-range repairs.
   The original proposal is in reviews/step-13-instructions-r1.md. The user selected
   investigation of a normalized
   limit of the existing cumulative weighting rule and authorized this new step's
   planning. No implementation dispatch is implied. Step 12 r6 handoff/review
   (reviews/step-12.md; reviews/step-12-review-r6.md) closes R2. Its exact active
   instructions are archived in reviews/step-12-instructions-r6.md. R3 is carried
   into Step 13: all six accuracy bins remain unconverged, with 12 unavailable
   diagnostic variants. Scientific acceptance and real runs remain user-controlled.
   Read reviews/step-12-r5.md for the weight diagnosis and the r4 performance
   handoff/review for the preserved optimization. NumPy's matrix target remains unmet.
   Step 11 revision 2 is accepted for progression. Earlier
   accepted prerequisites include reviews/step-10.md and
   reviews/step-10-review-r2.md; Step 10 revision 2 closes R1. Read Step 09
   reports for geometry/response and earlier reports as dependencies require;
   Step 06's latest accepted repair is revision 2, not revision 1.
4. Inspect README.md, pyproject.toml, relevant source/tests, and the six examples
   to verify interfaces. Do not infer implemented behavior from the target module
   diagram or copy an older plan's API guesses into a new assignment.

The design owns scientific conventions and the current API/dependency map; this
roadmap owns sequencing, status, and the review workflow. IMPLEMENTATION_STEP.md
contains exactly one detailed assignment. AGENTS.md owns standing working rules.
Implementation reports describe what was built; independent review reports
record findings and actual verification. Preserve old review reports rather
than rewriting them to imply a later state. The
[draft prompt](FISHHIGHZ_PLANNING_REVIEW_PROMPT.md) can start a new planning/review
agent. Keep these documents synchronized when decisions or status change.

### Checkout and environment

- Workspace: /global/cfs/cdirs/desicollab/users/acuceu/vega_dev. Package checkout:
  lib/fishhighz; source: lib/fishhighz/fishhighz, with no src directory. This is
  an independent child repository. Origin is configured as
  https://github.com/andreicuceu/fishhighz; this handover did not fetch or verify
  remote contents.
- At this snapshot local branch main is at
  0d69786a06d5d676564a51fad14a7156951c2228, the Step 02 commit. Implemented Steps
  03–12 include uncommitted and untracked files. Inspect both, not just git diff;
  do not reset, clean, stash away, or recreate them. A fresh remote clone is not
  assumed to contain this reviewed implementation.
- The root design, roadmap, and prompt sit OUTSIDE the FishHighz Git checkout.
  Preserve/carry those files separately when transferring workspaces or handing
  off a clone. Ignored .validation artifacts and virtual environments are local
  evidence, not guaranteed to be present elsewhere. Missing evidence is not a
  reason to invent results or launch an unrequested full capture.
- Python >=3.11; NumPy is the sole unconditional runtime dependency. The optional
  templates and cosmology extras have Astropy/SciPy; survey has SciPy;
  dev has pytest/Ruff/build.
  Imports remain lazy. Current tested version is 0.1.0.dev0; declarations are in
  pyproject.toml. Do not interpret its old scaffold description as absent code.
- Use the package-local .venv when present and appropriate. The umbrella shell
  may have no python executable; explicitly select the interpreter. Do not
  modify shared environments or assume a previous shell activation persists.

From the workspace root, the usual quick check is:

```bash
cd lib/fishhighz
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
```

The runner sets OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, and MKL_NUM_THREADS to 1,
then runs pytest, Ruff lint, and Ruff format checking. For targeted numerical
scripts outside it, set the same three limits explicitly. If a new development
environment is needed, follow README's isolated install of .[dev,templates,cosmology,survey].
Do not reinstall solely because an agent/session changed. Lightweight quick
checks run on the login node; long/intensive resource use follows the applicable
NERSC guidance, and Slurm actions always require explicit user authorization.

### Evidence and reference navigation

Step 12 revision-2 review passes 1059 ordinary tests, existing installed probes,
all 78 independent primary C/F reconstructions (maximum F discrepancy 3.161e-15)
and the 720-series plot check. All 139 implementation-manifest files, 48 wheel
modules and 75 input hashes match. The five original R1 probes reject; a new R2
probe accepts false convergence because metric operands are not bound to saved
trials. Independent binding checks on all 39 untouched accuracy records pass.
The handoff correctly reports 39/39 compatibility passes, 5/39 accuracy passes,
and 46/468 unavailable diagnostics. R3 is the unresolved forest accuracy result.
No new real forecast or installation was executed during review. The user has
archived the later 15x2pt-only revision-3 repair for future resumption. Revision 4
implements performance optimization, subsequent synthetic profiling and a
conditional next-run estimate. Independent review passes 1163 ordinary tests
(25 optional skips), 358 installed regressions, 129 performance regressions
outside the checkout, independent numerical/study checks and a short benchmark.
All 101 manifest entries and 49 wheel modules match. Compiled execution meets
both targets; the NumPy-only 5x matrix target remains unmet. Historical 1167.37 s
accuracy and 12.61 s compatibility timings cover unequal workloads. No real
forecast ran in r4. Revision 5 completed the restored bounded 15x2pt assignment. Independent review
passes 1183 quick tests (25 skips), 12 C/F reconstructions, 180 individual-spectrum
Fisher/error checks and the 240-series plot check. The saved six accuracy records
truthfully fail, but a synthetic weak-spectrum mutation defeats the new binding
check and falsely certifies a 1.9419% error change. See reviews/step-12-review-r5.md
and .validation/step12-review-r5/. Revision 6 now implements per-spectrum binding and
saved-trial metric/verdict checks. Its independent review passes 1243 quick tests
(25 skips), 251 installed regressions, six examples, Decimal/range checks and
all 12 joint plus 180 individual-spectrum reconstructions. The original false-pass
probe now rejects through all three routes. R2 is resolved; R3 remains open.
See reviews/step-12-review-r6.md and .validation/step12-review-r6/. No new forecast,
scientific choice, repair revision or next-step planning occurred during review.

Step 11 revision-2 review passed 833 tests in 33.77 seconds, lint/format, new
physical-width and legacy-divisor numerical oracles, all 20 original dtype
rejection cases, and the existing installed-wheel probe. All 34 reference hashes
match and 33 newly evaluated bounded comparisons pass with none blocked. These
are historical Step 11 review results; source changes require fresh verification.
The Step 12 request accepts Step 11 for progression. See reviews/step-11-review-r2.md.

Step 10 revision-2 review passed 638 ordinary tests in 36.52 seconds, lint/format,
the original R1 source/wheel replays, 12 independent boundary rejection cases,
six zero controls, one exact subnormal control, 96 Decimal weighting comparisons,
and the existing installed-wheel probe. R1 is resolved, and the user requested
Step 11, accepting Step 10 for progression. These are historical results, not new
Step 11 planning checks or claims about later source changes.
The 81 earlier bounded reference comparisons remain historical; all 34 reference
hashes match. Versions were Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, and
Astropy 8.0.1. No real legacy forecast ran in Steps 03–10.

Paths below are relative to lib/fishhighz and should be checked before use:

| Evidence | Location and use |
| --- | --- |
| Latest accepted implementation and review | reviews/step-11.md; reviews/step-11-review-r2.md. Revision 2 closes R1/R2 and is accepted for progression after the Step 12 request. Historical reports remain in reviews/step-11-r1.md and reviews/step-11-review-r1.md. |
| Earlier accepted prerequisite | reviews/step-10.md; reviews/step-10-review-r2.md. Step 10 numerical-range R1 is closed; do not confuse it with Step 11 R1. |
| Historical Step 10 revision 1 | reviews/step-10-r1.md and reviews/step-10-review-r1.md retain the original implementation and failed range finding; revision 2 is the accepted repair. |
| Latest completed evidence repair | reviews/step-12.md and reviews/step-12-review-r6.md: r6 passes review and resolves R2. The user authorizes progression to Step 13 investigation; six unconverged accuracy bins and 12 unavailable diagnostics remain historical scientific limitations, not accepted results. |
| Current Step 13 review/repair | reviews/step-13-review-r1.md records verified finite-step/fixed-grid numerics and required S13-R1/R2/R3 corrections. IMPLEMENTATION_STEP.md is revision 2, proposed for approval/dispatch. The r1 instructions, handoff, derivation and .validation/step13-r1-20260914T235633Z/ are preserved. New review checks and corruption copies are in .validation/step13-review-r1/. No production adoption or real run. |
| Step 12 revision 6 review | .validation/step12-review-r6/: exact instruction/handoff/governance snapshots, 155-entry and 52-module identity, 1243 quick tests/25 skips, 251 installed regressions, six examples, original false-pass rejection, Decimal/range controls and 12 joint/180 individual-spectrum saved-array checks. No new forecast or installation. |
| Step 12 revision 5 review | .validation/step12-review-r5/: exact instructions/handoff/governance snapshots, identity checks, quick/installed logs, 12 C/F and 180 individual-spectrum checks, 40 coupled-weight trial comparisons, 56 Decimal products, 240-series plot check and a false-pass weak-spectrum reproducer. reviews/step-12-instructions-r5.md preserves the reviewed assignment. |
| Step 12 revision 5 planning | .validation/step12-plan-r5/ preserves the preceding active instructions/governance, source identities and the r3-to-r5 restoration diff. The original r3 archive and r4 handoff/review remain unchanged. No forecast or implementation during planning. |
| Step 12 performance planning | reviews/step-12-instructions-r3.md is the exact archived scientific assignment. reviews/step-12-performance-assessment-r4.md and .validation/step12-performance-plan-r4/ contain short synthetic measurements, historical timing extraction and source-preservation evidence. No real forecast or production edits. |
| Step 12 r2 implementation evidence | .validation/step12-r2-implementation/: seven fresh actual references, profiles-streamed with 78 primary and 468 diagnostic requests, exact streamed wheel, plots-final, attribution-final, sensitivities-final and failed scientific gate. Preserve all interrupted attempts and cached-input provenance. |
| Step 12 r2 review evidence | .validation/step12-review-r2/: revision-2 assignment/governance snapshots, 139-file/48-module identity, 1059-test log, all 78 C/F oracles, 39 trial-binding checks, false-convergence reproducer, installed probes, plot and seven attribution checks. No new full forecast or installation. |
| Step 12 review evidence | .validation/step12-review-r1/: assignment/governance snapshots, 118-entry identity, 1017-test quick log, existing-wheel probes, independent saved-array F and five invalid-bundle reproductions, 104 new raw-reference checks. No full forecast or fresh installation. |
| Step 12 r1 implementation evidence | .validation/step12-r1-implementation/: exact 40-module wheel, source manifests, 104/112 reference-policy checks, one-bin real-02 and convergence evidence. Real results are one bin, not a full seven-case comparison. |
| Step 12 planning preservation | .validation/step12-plan-r1-20260913T161841Z/: previous Step 11 assignment, governance and source snapshots; 32-module wheel identity; seven original INI inventories; offline full/quick baseline checks with schema-v1 limitation. No new forecast or test-suite execution. |
| Step 11 implementation artifacts | .validation/step11-r2-implementation/: source manifest, exact wheel, installed environment and passing repair/reference reports. Preserve .validation/step11-r1-implementation/ with its historical blocked comparisons. |
| Step 11 independent review | .validation/step11-review-r2/: 104-entry source identity, quick log, independent density/measure/legacy and original dtype rejection checks, rerun 33 reference comparisons and existing-wheel probe, original assignment/governance snapshots. No new installation or full forecast. Preserve original failing evidence in .validation/step11-review-r1/. |
| Step 11 planning preservation | .validation/step11-plan-r1-20260913T145817Z/: pre-planning hashes, previous assignment and governance snapshots. Verified 28 reviewed production modules against the exact Step 10 wheel; no numerical-suite rerun. |
| Step 10 planning preservation | .validation/step10-plan-r1-20260913T021406Z/: pre-planning hashes, previous Step 09 assignment and design/roadmap/AGENTS copies. Planning verified all 13 Step 09 manifest entries; no numerical suite was rerun. |
| Step 09 planning preservation | .validation/step09-plan-r1-20260913T002414Z/: pre-planning hashes, previous Step 08 assignment, design/roadmap/AGENTS snapshots, and planning checks. No new forecast/test-suite evidence is implied. |
| Accepted Step 10 implementation artifacts | .validation/step10-r2-implementation/: before/after manifests, exact dist/ wheel, wheel-env/, probe and R1 closure results. Revision-1 artifacts remain in .validation/step10-r1-implementation/, including bounded reference evidence. |
| Accepted Step 09 artifacts | .validation/step09-r1-20260913-implementation/: final wheel and original evidence; see its reports for the historical reinstall disclosure. |
| Accepted Step 10 review artifacts | .validation/step10-review-r2/: snapshots, identity.json, production-r1-r2.diff, passing range replay, verify.py and numerical results, copied examples and wheel probe. Original failing evidence remains in .validation/step10-review-r1/. Manifest checks preceded governance status edits. |
| Accepted Step 09 review evidence | .validation/step09-review-r1/: preserved independent geometry/P1D/response and installed checks. |
| Template format/interpolation evidence | reviews/step-07.md and reviews/step-07-review-r1.md; .validation/step07-r1-20260912T194343/ and .validation/step07-review-r1/. Includes two read-only Vega FITS hashes and source/wheel identity. |
| Historical full legacy bundle | .validation/baseline/20260911T184946Z-7a6a795d/: seven primary DESI-2 cases plus one independent 15x2pt repeat, schema v1. See the Step 02 reports for its explicit inventory limitations. |
| Historical quick legacy bundle | .validation/baseline/20260911T205607Z-81e76c36/: one full-resolution 15x2pt case, schema v2, matching its full baseline. See reviews/step-02.md and reviews/step-02-review-r3.md. |

Offline legacy evidence checking, from the package root, is read-only with
respect to the saved bundles and does not run a forecast:

```bash
.venv/bin/python scripts/lyaforecast_baseline.py check \
  .validation/baseline/20260911T184946Z-7a6a795d --require-suite full
.venv/bin/python scripts/lyaforecast_baseline.py check \
  .validation/baseline/20260911T205607Z-81e76c36
```

Ordinary pytest uses synthetic workers and does not run those real cases. A real
quick capture is separate work (historically about two minutes), requiring a
compatible saved full baseline and the explicitly selected legacy interpreter.
Only include it where it validates the requested step; re-running unchanged
legacy code alone cannot test a new FishHighz implementation. Full capture is
seven primaries plus a repeat and remains explicit-user-request only, including
at Step 12. The user's revision-2 request was used for fresh legacy reference,
both profiles and bounded studies during implementation. The user subsequently
narrowed revision 3 to 15x2pt only, then archived it during performance-only
revision 4. Revision 5 used the restored assignment for bounded 15x2pt work, preserving
saved reference arrays. Revision 6 completed synthetic/offline evidence
repair; any further real run requires an explicit user request. No other real
case or new reference capture is authorized.

Read-only scientific reference entry points for subsequent planning:

| Topic | Files under lib/lyaforecast/lyaforecast |
| --- | --- |
| Default P1D | analytic_p1d_PD2013.py:P1D_z_kms_PD2013; power_spectrum.py:compute_p1d_kms and compute_p1d_hmpc. Preserve the redshift-dependent low-k floor; smoothing is a separate squared field-transfer factor. |
| Geometry and mode conventions | cosmoCAMB.py:velocity_from_distance, distance_from_wavelength, distance_from_degrees; covariance.py:_get_survey_volume and _get_num_modes. These are references, not a requirement to depend on CAMB. |
| Survey/source/instrument input | survey.py, tracer.py, spectrograph.py, and the seven original examples/desi2 INIs one directory above the Python package. |
| Weighting/noise | weights.py:compute_weights, _compute_weights_lya, compute_int_1/2/3, get_np_eff_lya; covariance.py:compute_eff_density_and_noise and _compute_total_power_lya. Inspect density normalization and auxiliary model queries before porting. |

Vega template/scaling references and their conventions are in design section 6.
Read sibling AGENTS.md before running sibling tests or modifying anything there;
changes are outside the default FishHighz scope. Check source/provenance and
licenses before copying formulas/code/assets. Baseline end results may not
contain every intermediate needed by a later port; identify missing evidence
explicitly and propose bounded new checks instead of fabricating a comparison.

### Planning and review procedure for fresh agents

When asked to plan, inspect prerequisites and unresolved scientific choices;
ask about consequential alternatives before fixing the dependent design. Write
only the requested next step, with revision, objective, current API inputs,
deliverables/exclusions, equations/units/ordering, independent acceptance tests,
exact quick-validation expectations, evidence requirements, and a stop-for-review
instruction. Routine implementation choices can use judgment. Do not implement
the step or dispatch an implementer. The user requested Step 12 after Step 11
revision 2 passed review, accepting it for progression. Step 12 revision 2 is
implemented and reviewed with R2/R3 requiring further work. Revision 3 is
preserved as an archive; revision 4 supplies reviewed performance optimization.
Revision 6 implements and passes review of the weak-spectrum evidence correction
using synthetic/offline checks. No further code correction was found. The user
now selects investigation of the existing rule's normalized/asymptotic limit and
requests Step 13 revision 1. This explicitly authorizes progression to the
investigation, while preserving R3 and withholding scientific acceptance of the
unconverged forecasts. Plan approval, dispatch, production adoption and any new
real run remain with the user. Further performance work and distribution move
to Steps 14 and 15; no detailed scope is added to them.

When asked to review, read the exact assignment/revision and implementation
handoff, verify the changed-file manifest against live files, and inspect code,
tests, and examples. Account for earlier dirty/untracked work; a HEAD diff alone
does not isolate the step. Run its applicable quick checks and focused independent
analytic/convergence checks. Verify packaging evidence when relevant; do not
expand every review into the original wheel/sdist/version matrix.

For a new build, use one exact artifact in a fresh output directory and a fresh
isolated install; verify module/metadata bytes and imports outside the checkout
with -I. An unchanged reported environment's probe may be rerun instead when
that is sufficient for review, but label it as a rerun, not a new fresh install.
Read evidence scripts before executing: some write output files or assume their
original run directory. Use a new reviewer directory and preserve prior fixtures,
logs, snapshots, and reports. Never select an arbitrary same-version wheel glob.

Report findings with severity, exact source locations, a concrete failure or
missing requirement, and the check that would close it. Separate required fixes
from optional improvements and reviewer-test mistakes; do not weaken accepted
scientific tolerances to obtain a pass. If repair is needed, increment the same
step's revision, retain still-applicable scope/tests, and add precise repair
instructions incorporating both the user's feedback and review findings.
Do not fix production code during a review-only assignment. Save a separate
reviews/step-NN-review-rN.md and update the status/register consistently.

When no required findings remain, say the step is ready for user progression;
do not replace the current plan until the user requests the next step. Do not
commit, push, delegate, or advance automatically. Report actual commands/results,
failed/blocked checks, skipped unrequested full execution, artifact identity,
and limits. Historical successes are not current reruns. NERSC MUNGE messages
have accompanied successful checks here; inspect exit status and test results
rather than treating those messages alone as either proof of failure or a reason
to ignore a failing command.

## User-controlled review cycle

1. The planning/review agent writes one detailed step, including scope, expected
   deliverables, and concrete acceptance checks. Its initial state is proposed.
2. The user reviews or edits the proposal and chooses an implementation agent.
   Asking that agent to implement the identified step authorizes its in-scope
   implementation and quick tests. Real full-suite execution requires a separate
   explicit user request; dispatch or plan approval alone does not authorize it.
   No additional approval is needed for the authorized quick checks.
3. The implementation agent completes only that step and supplies a reproducible
   handoff report. Passing tests means ready for review, not accepted.
4. The user reviews the work and supplies feedback, then asks the reviewing agent
   to assess the implementation and evidence against the same acceptance criteria.
5. If either reviewer requests changes, the planning/review agent rewrites the
   current step as a revision of the same numbered step, incorporating both sets
   of feedback. Resolved requirements remain satisfied; new tests target the
   remaining gaps. Disagreements affecting scope or scientific behavior are
   presented to the user rather than resolved by silently weakening criteria.
6. The user decides whether the step is accepted and explicitly requests the
   next detailed step. Only then is the step document replaced with the next
   step. No agent advances, dispatches another agent, or marks user acceptance
   on the user's behalf.

The user can pause, reorder, split, or revise any step. The roadmap is a planning
framework, not blanket permission to implement it.

## Step sequence

Steps are ordered to expose interoperability and numerical issues early. The
acceptance summaries below describe the intended outcome; exact tests belong in
the detailed step when it is proposed. A step that grows too large for a focused
review should be split by the user and planning agent before implementation.

| Step | Deliverable | Assessment before acceptance |
| --- | --- | --- |
| 01 | Minimal installable package using a direct fishhighz/ source directory and development checks | Source and built-wheel imports work independently of neighboring repositories; repeated validation installs the exact new artifacts into fresh isolated environments |
| 02 | Reproducible scientific reference capture and quick validation | A one-case quick mode compares against saved full evidence; portable checks validate full/quick coverage and integrity; new real seven-case captures plus repeat run only at the user's explicit request |
| 03 | Core array contracts, field/pair identity, parameter mapping, and fixed 2D grids | Synthetic cases validate ordering, covariance dependencies, units, cuts, shared parameters, callable contracts, and Gauss–Legendre/custom quadrature; fixed fiducial covariance/weights is the confirmed Fisher convention |
| 04 | Generic Gaussian covariance kernel from supplied observed signal/noise arrays | Analytic one/two-field limits, mode normalization, pair selection/permutation, signed cross power, PSD/singular cases, and bounded-batch equivalence pass; validation stays separate from the numerical kernel |
| 05 | General Fisher assembly and named results from supplied Jacobians | Analytic normalization, reusable Cholesky factors, independent-bin combination, priors added once, conditional/marginalized errors, and unit-aware degeneracy diagnostics pass |
| 06 | External callable integration and derivative engine | Explicit multi-provider pair routing, independent P1D invocation, mixed analytic/numerical global derivatives, tied-parameter perturbations, bounds-aware explicit steps, and opt-in convergence checks pass; a synthetic forecast reuses fixed covariance |
| 07 | Vega-format template input and interpolation | Generated FITS fixtures and read-only reference checks verify metadata, signed smooth/wiggle components, h-unit conversion, explicit redshift amplitude, cubic interpolation/derivative accuracy, and strict domains; confirmed optional Astropy/SciPy preparation feeds NumPy evaluation kernels |
| 08 | Built-in Kaiser, separate component scaling, and fixed BAO broadening | Analytic limits, three scaling bases, tracer-specific RSD, cross-width rules, domain/derivative checks, individual Vega conventions, and a synthetic fixed-covariance forecast exercise the existing provider interface |
| 09 | Survey geometry, optional Astropy adapter, instrument response, and default P1D | Integrated-volume convergence and analytic limits, explicit units/sigma/FWHM and legacy compatibility, unchanged P1D floor, fixed-covariance composition, bounded reference checks, and optional installed-wheel paths pass |
| 10 | Fixed per-field forest weighting/noise and galaxy sampling noise from prepared arrays | Explicit legacy cumulative iterations or supplied weights, normalized density/unit contracts, independent P1D, declared independence or full noise PSD, scalar/matrix tests, bounded DESI-2 input comparisons, and installed checks pass |
| 11 | Python survey-to-forecast orchestration and explicit raw density/SNR adapters | One built-in/external pipeline; explicit legacy-compatible optional SciPy preparation, independent bins/shared registry, fixed noise/factors, node batching, priors once, analytic and bounded-reader tests. INI/CLI/serialization deferred; revision 2 accepted for progression |
| 12 | DESI-2 validation and performance optimization | Revision 2 delivered the seven-case comparison, with R2/R3 still open. Scientific repair r3 is archived; r4 implements targeted optimization, subsequent profiling and a conditional runtime report using synthetic checks. Compiled targets pass independent review; NumPy-only matrix target remains unmet. Revision 5 supplies 15x2pt evidence and a weighting diagnosis; review finds weak-spectrum R2 binding incomplete and R3 unresolved. Revision 6 passes the bounded offline evidence repair review, resolving R2 while preserving r4 and strict defaults. R3 transfers to Step 13 investigation; the failed accuracy results are not scientifically accepted. |
| 13 | Investigate the normalized limit of the existing cumulative forest weights | User-selected option 2. Derive exact amplitude/shape evolution, distinguish fixed-grid from continuum limits, validate a diagnostic numerical representation, and assess both forest populations in all six saved 15x2pt bins. A supported finite limit or demonstrated absence is reviewable; no production adoption, alternative estimator, convergence waiver or real forecast in revision 1. |
| 14 | Further performance work if needed after the weighting investigation | Former prospective Step 13, renumbered after insertion of the weighting investigation. Preserve r4 optimizations; reassess remaining need after scientific review. No detailed plan or dispatch; any JAX prototype remains separately scoped. |
| 15 | Distribution and user-facing documentation | Clean installation, packaged resources, standalone examples, public API documentation, and agreed regression checks pass. Scope deferred INI translation, CLI and serialization explicitly when this step is planned |

The sequence deliberately verifies covariance/Fisher machinery with analytic
arrays before connecting realistic survey noise. It also exercises an external
model before the built-in model becomes the default path, so interoperability is
tested as a core requirement.

NumPy with dedicated Numba-friendly kernels is the confirmed starting backend
(user decision, 2026-09-11). Targeted Numba compilation follows validated reference
calculations and profiling; a JAX backend is optional future work. The user also
confirmed 2D P(k,mu) on 2026-09-11, as recorded in design section 5.1. No multipole
forecast implementation is required. Step 03 introduces NumPy as the first
runtime scientific dependency. On 2026-09-11 the user also confirmed
Gauss–Legendre integration in mu and within k-bin edges with custom-grid support,
and fixed fiducial covariance and survey weights. These two Step 03 design gates
are resolved. Other proposed contracts remain subject to the user's plan review.

## Handoff and review evidence

Use a stable step ID and a revision number, such as `01 / revision 1`. Each
implementation handoff records the plan revision, changed files, exact commands
and outcomes, artifact locations, limitations, and feedback addressed. Identify
the reviewed code with its Git revision and dirty diff when available, or a
file-hash manifest when the package is not yet a Git repository.

Keep concise, completed-step reports under `lib/fishhighz/reviews/` as they are
produced; keep generated build/test/baseline outputs in ignored validation
directories. Reports should refer to logs without committing large generated
outputs. The single-step file must not accumulate old plans or future detailed
steps. Accepted reports preserve the evidence and decision history.

Revisions should list unresolved findings and the specific checks that close
them. Preserve the accepted behavior of earlier steps. A change to an accepted
public contract or scientific convention requires the user's explicit decision
and a corresponding design-document update.

## Validation cadence

Implementation and review agents run quick checks by default: self-contained
ordinary tests, relevant full-resolution 15x2pt comparisons against compatible
saved evidence, and offline bundle checks. The capture CLI should default to
quick mode. Synthetic tests may exercise full-mode orchestration without running
real forecasts.

Execute the real full suite only when the user explicitly asks for that run.
This applies at review checkpoints and step completion, and after failures or
relevant changes. General implementation/testing requests and approval of a plan
are not requests to execute the full suite. Never fall back to it automatically;
report missing or incompatible baseline evidence. One explicit full-run request
does not authorize subsequent full reruns.

Do not routinely ask for a full run or make an unrequested one a prerequisite for
a review handoff. Report it as not run under the user's validation policy, with
historical full evidence clearly distinguished from current quick evidence.
Preserve accepted reference bundles for later FishHighz comparisons; running
unchanged lyaforecast does not itself validate new FishHighz code. Exact quick
checks remain part of each user-reviewed detailed assignment. This policy
governs all roadmap steps, including future complete survey validation.

Historical explicit exception to the default quick cadence: the user's Step 12
revision-2 request authorized seven fresh legacy cases, both FishHighz profiles
and bounded studies. Implementation has attempted and preserved this work.
The user subsequently restricted revision 3 to 15x2pt, then archived it during
performance-only revision 4. Revision 5 used the restored bounded scientific assignment for
15x2pt only: all six bins and 15 individual spectra plus the joint forecast.
Review does not authorize another real run. Revision 6 completed its synthetic
and offline evidence checks; the six unconverged accuracy bins and 12 unavailable
diagnostics remain explicit. No other real case, fresh NewForecast capture or
Slurm action is authorized. Scientific choices, dispatch and acceptance remain
with the user.

Step 13 revision 2 proposes, on user dispatch, derivations, synthetic controls
and bounded diagnostics on saved weight-input arrays. It does not authorize
regeneration of raw samples, external model/background calls or a full 15x2pt
forecast. A normalized numerical representation must retain the amplitude in the
nonlinear update and establish finite-step equivalence. No outcome is assumed:
fixed-grid convergence, continuum convergence and possible divergent limits are
separate questions. Production adoption or a changed estimator requires later
user review. Preserve all earlier scientific failures and source evidence.

## Progress register

The parallel weighting register is maintained in the
[diagnostic roadmap](FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md). W08 revision 2
passes scientific review, as do W09–W11. The user selected the baseline reference
for accuracy-profile adoption. W12 revision 1 is proposed for approval/dispatch;
W13 remains separately requested broader validation.
This separate weighting sequence does not alter the package-step states below
or accept the failed Step 12 forecasts.

| Step | Current state | Evidence and acceptance |
| --- | --- | --- |
| 01, revision 4 | Accepted for progression | 2026-09-11: user reported no comments; independent review passed with no outstanding findings, and the user subsequently requested Step 02. See [revision-4 review](../reviews/step-01-review-r4.md). |
| 02, revision 3 | Accepted for progression | 2026-09-11: user reported no comments; [independent revision-3 review](../reviews/step-02-review-r3.md) passed with 52 tests and saved-evidence checks. Changes were committed/pushed as `0d69786`; the user subsequently requested Step 03. No new real full run was required. |
| 03, revision 1 | Accepted for progression | 2026-09-11: user reported no comments; [independent review](../reviews/step-03-review-r1.md) found no required changes. All 146 tests, lint, source/wheel identity, installed-wheel probe, and independent array checks passed. No real reference forecasts were run. The user subsequently requested Step 04. |
| 04, revision 1 | Accepted for progression | 2026-09-11: user reported no comments; [independent review](../reviews/step-04-review-r1.md) found no required changes. All 186 tests, lint/format, 63 independent quadratic-form cases, source/wheel identity, and installed-wheel probe passed. No real reference forecasts ran. The user subsequently requested Step 05. |
| 05, revision 1 | Accepted for progression | 2026-09-11: user reported no comments; [independent review](../reviews/step-05-review-r1.md) found no required changes. All 261 quick tests, lint/format, 36 independent assembly/result cases, 8 singular-diagnostic cases, source/wheel identity, and installed-wheel probe passed. No real reference forecasts ran. The user subsequently requested Step 06. |
| 06, revision 2 | Accepted for progression | 2026-09-12: [revision-2 review](../reviews/step-06-review-r2.md) confirms R1 resolved with no further changes. All 332 quick tests, lint/format, 16 independent zero-call rejection cases, 6 valid one-point refinements, source/wheel identity, and installed-wheel probe passed. No real reference forecasts ran. The user subsequently requested Step 07. |
| 07, revision 1 | Accepted for progression | 2026-09-12: user reported no comments; [independent review](../reviews/step-07-review-r1.md) found no required changes. All 408 quick tests, lint/format, 24 independent analytic cases, two read-only reference-file checks, source/wheel identity, and installed-wheel probe passed. No real reference forecasts ran. The user subsequently requested Step 08. |
| 08, revision 1 | Accepted for progression | 2026-09-12: user reported no comments; [independent review](../reviews/step-08-review-r1.md) found no required changes. All 468 quick tests, lint/format, 18 independent five-field model cases and 216 derivative columns, source/wheel identity, and installed-wheel probe passed. Saved Vega convention checks retain matching source hashes. No real reference forecasts ran. The user subsequently requested Step 09 (2026-09-13). |
| 09, revision 1 | Accepted for progression | 2026-09-13: user reported no comments; [independent review](../reviews/step-09-review-r1.md) found no required changes. All 555 quick tests, lint/format, 48 independent volume cases, 272 high-precision P1D points, signed-response covariance/Fisher checks, source hashes, and the installed-wheel probe passed. Saved 33 legacy-method comparisons retain matching source hashes; no real forecast ran. The user subsequently requested Step 10 (2026-09-13). |
| 10, revision 2 | Accepted for progression | 2026-09-13: user reported no comments. [Revision-2 independent review](../reviews/step-10-review-r2.md) closes R1 with no further required changes. All 638 quick tests, source/wheel R1 replays, 12 boundary rejection cases, six zero controls, one exact subnormal control, 96 Decimal comparisons, and installed probe passed. All 34 reference hashes match; 81 earlier bounded comparisons remain historical. No real forecast ran. The user subsequently requested Step 11 on 2026-09-13. |
| 11, revision 2 | Accepted for progression | 2026-09-13: [revision-2 review](../reviews/step-11-review-r2.md) closes R1/R2 with no further required changes. All 833 quick tests, 24 independent physical-density queries and count measures, 12 legacy-divisor tables, 20 original dtype rejection cases, three metadata controls, and existing 32-module installed-wheel probe passed. All 34 reference hashes match; 33 new comparisons pass and none are blocked. Historical r1 findings remain preserved. No real full forecast, production fixes or Step 12 planning. User subsequently requested Step 12 on 2026-09-13. |
| 12, revision 4 | Optional compiled performance path passes independent review; user acceptance pending | 2026-09-14: [r4 handoff](../reviews/step-12-performance-r4.md) and [r4 review](../reviews/step-12-review-r4.md): 1163 quick tests/25 skips, 358 installed regressions, 129 performance regressions outside checkout, 32 C/F combinations, 28 boundary cases, four layouts, five baseline-matched studies and installed examples pass. All 101 manifest entries/49 modules match. Composed speedups 5.77x NumPy / 19.21x compiled; matrix 2.77x / 17.56x. Default NumPy matrix target remains unmet; no required code corrections found for compiled path. R2/R3 remain in [archived r3](../reviews/step-12-instructions-r3.md). No real forecast, scientific acceptance, repair revision or Step 13 planning. |
| 12, revision 5 planning | Historical restoration plan; subsequently implemented and reviewed below | 2026-09-14: user requested [revision 5](../reviews/step-12-instructions-r5.md) using the [archived r3 instructions](../reviews/step-12-instructions-r3.md). Restores R2 trial binding and R3 convergence diagnosis/repair for 15x2pt only: six bins, 12 primaries, 180 individual-spectrum results and 72 diagnostics; preserves reviewed r4 optimizations and 1188 current tests. Scientific prescription changes require user approval. No forecast, production edit or test-suite rerun during planning; Step 12 is not scientifically accepted. |
| 12, revision 5 review | Changes required; scientific acceptance withheld | 2026-09-14: [r5 review](../reviews/step-12-review-r5.md) passes 1183 quick tests/25 skips, 191 installed regressions, 12 saved C/F and 180 individual-spectrum checks, 40 coupled-trial comparisons and 240 plotted series. Weak-spectrum corruption still produces false convergence in a one-record writer/offline bundle (R2); no such inconsistency was found in the submitted arrays. Six accuracy bins remain unconverged and 12 diagnostics unavailable (R3). The fixed-weight quadrature and zero-weight Jacobian diagnosis are verified. No real forecast or production fix during review. |
| 12, revision 6 | Bounded evidence repair passes independent review; scientific acceptance open | 2026-09-14: [r6 review](../reviews/step-12-review-r6.md) resolves R2. The original weak-spectrum corruption rejects in validator, writer and offline reader; 1243 quick tests/25 skips, 251 installed regressions, six examples, 21 Decimal comparisons, four mixed-scale mutations and 12 joint/180 individual-spectrum numerical checks pass. All 52 source/wheel/installed modules match; 35 non-validation modules and historical r5 evidence are unchanged. Six accuracy bins remain unconverged, 12 diagnostics unavailable (R3). No scientific choice, real run, new revision or Step 13; await the user. |
| 13, revision 1 | Requires corrections after independent review | [r1 review](../reviews/step-13-review-r1.md): 1303 quick passes/25 skips, 121 installed regressions, six examples, all 53 module identities and all-sixty independent coefficient/source checks pass. The refinement trend is verified; continuum overclaim (S13-R1), source/attempt binding (S13-R2) and scalar range (S13-R3) require repair. Original reports/evidence are preserved; Step 12 R2 remains closed and R3 remains open. |
| 13, revision 2 | Implemented; awaiting independent/user review | 2026-09-15 overview: [handoff](../reviews/step-13.md) and live diagnostic source are present; all eight saved source/report hashes match. The handoff reports 1318 quick passes/25 skips and classifies the continuum limit as unresolved for both forest populations in all six bins. Closure of S13-R1/R2/R3 remains unreviewed; Step 12 R3 remains open. No tests or forecasts rerun for this overview; no production adoption, scientific acceptance or progression. |

The planning/review agent updates this register after user decisions. Add an
acceptance date and report link only after the user accepts the identified work.
