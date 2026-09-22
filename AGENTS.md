# FishHighz agent guidance

## Project and authoritative documents

FishHighz is a standalone Fourier-space covariance/Fisher forecasting package.
Steps 01–11 are accepted for progression; Step 12 evidence repair passed review,
while historical accuracy forecasts remain scientifically unaccepted.

On 2026-09-17 the user marked **Step 13 and W12 done**, closed the entire
previous weighting diagnostic plan, and directed that these assignments not be
resumed. W12's bounded fixed-reference accuracy-profile review passed; preserve
its implementation. Step 13's continuum conclusion remains unresolved in the
historical evidence. Deferred W01 repairs and the old W13 proposal are retired.

Current planning update, 2026-09-18: the five-stage compatibility weighting
study and additional grid-to-BAO check are complete. `docs/research/FOREST_WEIGHTING_DECISION.md`
records the accepted early-lyaforecast baseline and retained McDonald alternative.
S1–S3 now implement full-compatibility, fixed-compatibility and updated accuracy,
qualify all six revised-accuracy bins within the stated finite-refinement tests,
and complete the selected six-bin DESI-2 comparison. Astra (light) implementations
and Sol (high) independent scientific reviews are recorded in `docs/archive/reviews/s1.md`,
`docs/archive/reviews/s2.md`, `docs/archive/reviews/s3.md` and their `-review.md` reports; all three reviews
pass with no consequential correction requested. S4 attribution is also complete with a passing Sol scientific review:
`docs/archive/reviews/s4.md`, `docs/archive/reviews/s4-impact-table.md`, `docs/archive/reviews/s4-profile-inventory.md`
and `docs/archive/reviews/s4-review.md`. Scientific acceptance remains with the user;
S5 reassessment began on 2026-09-21 as a joint planning discussion. The user
requested an evidence summary, interpretation and proposed next work; that
initial request authorized planning only. The later S5 execution authorization
below now governs implementation and review. The user
confirmed a research-ready Python API as the first S5 milestone and the tested
S2–S4 accuracy recipe as its research baseline, retaining mixed-pair damping and
full wiggle AP derivatives. The reference-mode comparison is now a separate
later study in `docs/research/DEFERRED_SCIENTIFIC_TESTS.md`; S5 covers only prescription
documentation, Python API documentation/examples, and validation. No further
scientific decision is pending before implementation within that scope. The
legacy weighting reference remains in use, with a 1% BAO-error attribution
trigger and 0.1% numerical-refinement target.

For every future forecast, exclude bin-1 correlations involving LBG, LAE or
lya(lbg). Retain only lya(qso) auto, QSO auto and their cross in bin 1, including
its joint covariance. Bins 2–6 retain all 15 spectra. This is a forecast selection,
not a plot cut; preserve historical files and their original selection.

The completed validation study implements five fixed-count variants and opt-in
full-sum forest-auto adaptive stopping. Its trajectory, BAO and magnitude-grid
results and independent reviews are linked from the weighting decision. Iteration
convergence does not establish magnitude-grid convergence or physical accuracy.
The historical extra grid-to-BAO check covers bins 2–6 only; S2 adds the retained
three-spectrum bin-1 check. W12's explicit implementation is retained. Revised
accuracy uses its own fiducial per-field P/B, response and magnitude quadrature,
with frozen weights in BAO derivatives. Its maximum tested individual/joint
BAO-error refinement changes are 0.0011411%/0.0012685%. S3 finds material
accuracy/fixed differences, triggering the completed S4 study. Its 364 controlled calculations identify
mixed forest–galaxy damping as the dominant isolated reduction, partly opposed
by the grouped AP derivative change. Maximum endpoint mismatch is 0.000129%;
maximum tested error refinement is 0.010783%. These quantify model assumptions
and numerical behavior, not scientific acceptance of the reconstruction physics.
Execution authorization, 2026-09-18: the user subsequently authorized the
coordinator to dispatch Astra (light) implementation agents and Sol (high)
scientific reviewers, check results and progress through S1–S3 inclusive.
Reviews should request changes only when likely to affect scientific conclusions.
The user subsequently authorized S4 attribution with the same agent roles and
scientific review criterion, including a report table of all scientifically
meaningful profile changes and their joint BAO effects. This includes the scoped
S2/S3 and S4 real forecasts; it did not itself authorize S5 calculations, commits,
pushes or Slurm actions. S5 now has its separate authorization below. The user explicitly directed login-node
execution for these short forecasts. Use one numerical thread and bounded runs.
The remaining standing review controls apply; scientific acceptance stays with
the user.

- [Archived scientific design](docs/archive/planning/FISHHIGHZ_DESIGN.md): historical scientific conventions and API
  context; dated former assignments are historical.
- [Archived implementation roadmap](docs/archive/planning/FISHHIGHZ_IMPLEMENTATION_PLAN.md): historical
  sequence and status; it does not authorize current work.
- [Archived compatibility weighting test plan](docs/archive/planning/FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md):
  completed study and evidence, including the five variants, individual/joint
  plots and historical plotting cut for uncertainties above 0.2.
- [Deferred scientific tests](docs/research/DEFERRED_SCIENTIFIC_TESTS.md): later, non-critical
  studies with agreed and proposed statuses kept distinct.
- [Archived closed diagnostic plan](docs/archive/planning/FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md): historical closure
  record and link to its full archive.
- `docs/archive/notes/IMPLEMENTATION_STEP.md` and `docs/archive/notes/WEIGHTING_DIAGNOSTIC_STEP.md` retain retired
  assignments; neither authorizes implementation or dispatch. New detailed
  assignments require a user request. Do not follow their historical pending-work
  instructions as current tasks.

This file supplies standing context and working rules. Keep step-specific tasks,
acceptance criteria, and feedback in the current-step file, not here. Follow the
user's latest explicit instructions if they change a document's requirements;
report inconsistencies so the planning agent can reconcile the documents.

### S5 execution authorization, 2026-09-21

The user has now explicitly requested S5 implementation, with **Sol (medium)**
implementation agents and **Sol (high)** independent review agents. The coordinator
controls dispatch and checks results. S5 reviews cover scientific correctness,
software and documentation. This supersedes the earlier planning-only status
and narrower review scope for this assignment. Implement the roadmap's three
S5 parts, including bounded synthetic/example, focused-test and installed-package
checks. Keep the adopted numerical recipe and historical evidence unchanged;
report any consequential correction before expanding numerical scope. D1/D2,
new real-survey studies, commits, pushes and Slurm actions remain outside scope.

### S5 completion, 2026-09-21

S5 is complete with a passing Sol (high) scientific, software and documentation
review. `docs/research/RESEARCH_BASELINE.md` defines the adopted S2–S4 recipe;
`examples/research_bao_forecast.py` demonstrates the standalone Python API with
synthetic inputs. `docs/archive/reviews/s5.md` and `docs/archive/reviews/s5-review.md` record 101 focused
tests, two installed checks, direct covariance reconstruction and unchanged
numerical code/evidence. Only package docstrings and three validation error-message
occurrences changed in existing source. The guide and examples ship in the source
distribution. D1/D2 remain deferred; no next step is authorized by this completion.

## Roles and review control

The user chooses and prompts implementation agents and retains control over plan
review, implementation review, acceptance, and progression. An instruction to
implement the current step includes its quick tests. The real full validation
suite requires an explicit user request; a general implementation, testing,
review, or step-completion request does not authorize it. Complete the bounded
assignment and quick checks, report evidence, and stop for review; do not ask for
extra permission to run already authorized quick tests.

Do not advance the roadmap, write future detailed steps, or treat passing tests
as user acceptance. Do not spawn or dispatch other agents unless the user
explicitly requests delegation. Do not edit planning/governance documents during
an implementation assignment unless the user authorizes that revision. Report
necessary scope or design changes rather than silently expanding the step.

A review assignment assesses code and test evidence against the identified plan
revision. Report findings with concrete file references and missing or failing
checks. Do not silently implement fixes during a review-only assignment. The
planning/review agent incorporates feedback into the same step when asked; only
the user's acceptance and request for the next step permit advancement.

## Minimal diagnostic implementation and scientific review

For forest-weighting diagnostics, use only the minimum source, code and tests
needed to answer the assigned scientific question. A report, derivation and tiny
standalone check can be sufficient. Do not require full suites, fresh wheels or
environments, generalized validators or hypothetical provenance hardening unless
they are needed to resolve a concrete scientific ambiguity in the assignment.
These diagnostic-specific instructions supersede broader routine validation
requirements below; they do not alter historical package-step evidence.

During review, propose further changes only if they are likely to affect the
scientific conclusion or its justified scope. For each proposed correction,
identify the evidence, smallest necessary check and possible scientific impact.
Prefer a justified narrower conclusion or explicit uncertainty when sufficient.
Do not generate optional cleanup/hardening lists. Begin the review with a short
scientific summary explaining how proposed changes could affect the conclusion,
or state that no scientifically consequential changes are needed. Preserve user
feedback, control of approval/dispatch/acceptance and the stop after each step.

## Scientific invariants

- Forecast 2D P(k,mu), as confirmed on 2026-09-11, with model/noise evaluation
  on a 2D grid and independent Fourier-cell covariance blocks. Multipole
  compression is not required. Correlation functions, sampling likelihoods,
  and samplers are outside scope.
- Treat redshift bins as independent and use a common volume across fields in
  each bin. Preserve covariance between different spectra within a bin.
- Distinguish physical tracers from observed samples. Forests measured through
  QSO and LBG spectra can share physics but have different noise and response.
- Preserve the spectra required to construct covariance even when those spectra
  are not selected as forecast observables.
- Accept arbitrary external callable models; automatic differentiation is
  optional. Never require an external model to depend on FishHighz or JAX.
- The built-in template loader reads Vega-format K/PK/PKSB: no-wiggle is PKSB,
  wiggle is PK-PKSB. Do not regenerate a decomposition in the forecasting path.
- Built-in forest Kaiser factors use b and nuisance beta; galaxy factors use
  b+f*mu², with a shared growth-rate target in each bin.
- BAO broadening affects wiggles only. Per-tracer parallel/transverse widths are
  fixed inputs; cross widths squared are the mean of the two squared auto widths.
  Varying f never changes these input widths.
- Wiggle and no-wiggle scaling can be separate or tied, using the three confirmed
  parameterizations. Follow the design's Fourier mapping and normalization.
- k_min/k_max define fixed observed-coordinate cuts for a run. Model-domain
  padding for rescaling must not alter those cuts or mode counts.
- Preserve lyaforecast's default P1D formula and low-k treatment. External P3D
  does not imply replacement P1D; model response/noise ownership must be explicit
  to avoid counting smoothing or noise twice.

The user confirmed NumPy with dedicated Numba-friendly kernels as the starting
backend on 2026-09-11. Validate NumPy reference calculations before targeted
compilation; a JAX backend is not an initial requirement. On 2026-09-11 the user
confirmed Gauss–Legendre quadrature in mu on [0,1] and within user-supplied k-bin
edges, with custom nodes/weights supported, and fixed fiducial covariance and
survey weights for the initial Fisher calculation. Differentiate the predicted
mean; do not add covariance-derivative information. Other proposed contracts
remain subject to the user's plan review. Read the design
rather than inferring scientific choices from similar code in a neighboring
package.

## Structure and implementation practice

Keep file/configuration handling, validation, names, units, adapter setup, and
caching separate from intensive numerical functions. Kernels accept explicit
arrays/scalars and integer indices, with no filesystem access, logging, dynamic
plugin lookup, or configuration objects in numerical loops.

Prefer small, functional interfaces and batched evaluations. Preserve useful
NumPy/vectorized operations; optimize measured bottlenecks rather than rewriting
for style. Use explicit float64 for scientific calculations and inspect numerical
conditioning. Do not silently change tolerances, regularize singular results, or
enable aggressive fastmath. Keep dependencies minimal and optional integrations
out of package import. Do not add an abstraction solely for possible future use.

The user confirmed optional Astropy/SciPy template preparation on 2026-09-12.
Astropy reads FITS; SciPy prepares spline coefficients. Keep both in the
`templates` extra with lazy imports. Repeated template evaluation uses dedicated
NumPy kernels, and the external-model base installation requires only NumPy.
On 2026-09-13 the user also confirmed an optional Astropy cosmology adapter for
Step 09, integrated bin volumes, and explicit sigma/FWHM instrument conversions
with separately labeled legacy compatibility. The cosmology extra
keeps Astropy/SciPy optional; supplied-background preparation remains NumPy-only.
On 2026-09-13 the user confirmed prepared-array inputs for weighting/noise,
explicit legacy cumulative iterations plus caller-supplied weights, and explicit
independent sampling or supplied full noise. Keep per-field weights consistent
across pairs; overlap-derived noise remains deferred. For Step 11 the user
confirmed a Python API and standalone raw readers, deferring INI/CLI/serialization.
A labeled optional SciPy survey adapter preserves explicit legacy density/SNR
interpolation and normalization settings with strict domain/value errors; no
legacy floors or clamps. Keep the normalized-array forecast path NumPy-only.
The user also confirmed nonuniform redshift support with explicit raw cell widths
and a separately labeled legacy first-spacing conversion. Forecast-bin bounds
and raw source-table cells are separate; do not infer physical widths from centres.
For Step 12 the user confirmed wiggle-only BAO ap/at validation, a separately
labeled legacy input floor/clamp adapter, and a real lyaforecast intrinsic-P3D
example. A negative interpolated-density floor is an explicitly labeled extension
beyond exact legacy behavior. Strict reader defaults stay unchanged; retain
fallback provenance and keep real full-suite execution explicitly opt-in.
The user-authorized Step 12 revision-2 full comparison has been attempted and
preserved: seven reference cases, both profiles and bounded diagnostic studies.
The subsequent 15x2pt-only revision-3 repair was archived during performance-only
revision 4. Revision 5 used the restored scientific assignment on the reviewed r4
baseline. Preserve its saved 15x2pt primaries, diagnostics and weighting diagnosis.
Revision 6 completes the reviewed synthetic/offline evidence repair. Step 13's
normalized-limit investigation is now closed; preserve its diagnostic and saved
results without resuming its former assignment. In the new compatibility study,
retain amplitude in nonlinear updates and distinguish iteration from grid
sensitivity. The planned 15×2pt execution needs its own explicit assignment;
do not run the other six cases or recapture NewForecast. Preserve r4 optimizations
and regression tests. Further weighting/input-policy adoption requires the user's
scientific decision.
The literal legacy estimator path
is validation-only. The accuracy profile preserves existing-model scope and
retained input policies with sensitivity tests; neither changes strict defaults.

Use four-space indentation, snake_case functions/modules, PascalCase classes,
NumPy-style public docstrings, pytest, and the package's configured Ruff rules.
Use a direct `fishhighz/` source directory under this repository root, as requested
by the user; do not add an intervening `src/`. Restrict package discovery to the
intended import package and retain isolated installed-wheel checks outside the
checkout so working-directory imports cannot mask packaging defects. Read
development commands from the current README/pyproject rather than assuming
planned scripts already exist.

## Reference packages and validation

`../lyaforecast` is the scientific survey/noise reference; `../vega` supplies
template format and selected model conventions. Treat both as read-only during
FishHighz steps unless the user explicitly expands scope. Read a reference
package's own AGENTS.md before any work that runs its tests or changes it.
Do not alter authoritative example INIs, FITS templates, or existing outputs.
Retain provenance and applicable license notices for any copied scientific code;
do not invent a package license or redistribute assets without checking rights.

Capture the seven authoritative DESI-2 baselines before numerical implementation,
as scheduled in the roadmap. Record inputs, versions, source state, hashes, logs,
outputs, and timings. A scaffold-only step does not require a forecast baseline.
Existing baselines are reference evidence, not assumed-current results. Match
the environment for numerical comparisons and explain expected differences from
the legacy BAO recipe rather than forcing agreement by changing physics.

Run only quick checks by default: fast ordinary tests, relevant one-case reference
checks, and offline validation of saved bundles. Execute the real full suite only
when the user explicitly asks for that run, including at final step/review
checkpoints. Failures, relevant changes, or a missing/incompatible baseline do not
authorize automatic full execution or fallback. A request for one full run does
not authorize later full reruns. Synthetic tests of full-mode orchestration are
allowed in the quick suite because they do not run real forecasts.

Do not routinely ask for a full run or hold an otherwise complete review handoff
pending one. Report an unrequested full run as not run under the user's policy;
do not call it a failed or outstanding acceptance check. Report any actual quick
validation failure or blocked comparison separately.
Reuse accepted immutable reference evidence for FishHighz comparisons when its
provenance remains applicable. Rerunning unchanged reference code alone does not
test newly implemented FishHighz behavior; future steps need their own focused
numerical tests and comparisons.

Tests should be self-contained, deterministic, and independent of NERSC paths or
neighboring installed packages. Use synthetic arrays and small generated
fixtures for ordinary tests; keep local cross-package comparisons as explicitly
identified validation. Prefer analytic limits, finite-difference convergence,
normalization checks, and observable behavior over tests that repeat the code.

Implement the current step's tests, run its required checks, and report actual
commands/results, failures, and limitations. Keep generated artifacts separate
from concise review reports; tie evidence to the code and plan revision reviewed.
Do not describe planned or skipped tests as passed. After a revision, rerun the
affected tests and agreed regression checks without unrelated test expansion.

## Workspace operations

Follow the parent NERSC guidance. Use bounded searches and keep writes inside
this package and explicitly authorized artifact locations. Keep source and test
changes local; do not modify shared environments or sibling repositories.

Use isolated environments where practical. Lightweight development checks run
on the login node with OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, and
MKL_NUM_THREADS=1. Follow lyaforecast's own guidance for its small DESI-2 reference
runs. Long, memory-intensive, MPI, or accelerator jobs require appropriate
resources and explicit user authorization for any Slurm action; do not launch
them merely because a performance step exists in the roadmap.

Never overwrite unrelated user changes, delete unverified paths, or expose
credentials. Do not commit, push, publish, or create remote resources unless the
user requests it. Reference code, data, and command outputs must not be sent to
third-party services independently.
