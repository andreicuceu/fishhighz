# FishHighz agent guidance

## Project and authoritative documents

FishHighz is a new, standalone forecasting package in this directory. It computes
Fourier-space Gaussian covariances and Fisher forecasts for combinations of
galaxy/quasar and Lyα forest fields, with interchangeable external P3D and P1D
providers. The scaffold, reference capture, core array contracts, and covariance
from supplied observed powers, Fisher assembly, and named results from supplied
Jacobians, external model integration, and derivative generation are accepted for
progression, as are built-in template input and interpolation. The intrinsic
Kaiser/BAO model is also accepted for progression. Step 09 geometry, instrument
response, and default P1D are accepted for progression. Step 10 fixed forest
weighting/noise and galaxy sampling noise are accepted for progression after
revision-2 review resolved the original numerical finding. Step 11 revision 2
implements Python survey orchestration and explicit raw density/SNR readers;
independent review passed, closing nonuniform redshift and input dtype findings.
The user accepted Step 11 for progression by requesting Step 12. Step 12
revision 1 implements DESI-2 validation and opt-in legacy input/P3D adapters.
Revision 2 implements the full seven-case comparison and plots. Independent
review verifies all 39 compatibility bins and closes the five original R1
probes; R2 requires binding convergence evidence to its recorded trials.
Only five galaxy-only accuracy bins converge; all 34 forest bins remain
unresolved (R3). The user archived the exact revision-3 scientific repair in
reviews/step-12-instructions-r3.md for later use. Revision 4 implements performance
optimization and passes independent review for the optional compiled path;
the default NumPy matrix target remains unmet. Revision 5 restores that
assignment and supplies a 15x2pt-only handoff.
Revision-5 review found weak-spectrum trial binding incomplete. Revision 6
implements per-spectrum checks and saved-trial metric/verdict reconstruction;
independent review passes, resolving R2. All six accuracy bins remain unconverged
(R3), with a verified weight diagnosis and 12 unavailable diagnostics. The user
selected investigation of a normalized/asymptotic limit of the same
cumulative rule and requested Step 13 revision 1. This authorizes planning and
progression of that investigation, not scientific acceptance of the failed
forecasts. Step 13 revision 1 is implemented and independently reviewed. Finite-step and
fixed-grid numerics are verified, but continuum-verdict, source/attempt-binding
and diagnostic scalar-range findings require revision 2. The current assignment
remains diagnostic-only on synthetic and saved inputs; production adoption,
scientific choices, new real runs, plan approval and dispatch remain with the user. Later performance/distribution
steps are numbered 14/15. INI translation, production CLI and serialization
remain future work.

The parallel forest-weighting diagnostic W01 supports inherited finite-update
sensitivity; W01-R4/R5 remain deferred. W02–W07 pass scientific review. W05
identifies the explicit inverse-variance reference with the legacy seed, W06
verifies its bounded magnitude stability, and W07 finds smaller BAO errors for
that reference in one accuracy QSO-forest auto-spectrum/bin; see
[W07 review](reviews/weighting-diagnostics-w07-review-r1.md). These results
establish no general multi-mode optimum or physical validation of input policies.
On 2026-09-16 the user selected the fixed-reference option. W08 revision 1
implements explicit `method="inverse_variance"` with direct scalar B_star.
Independent review verifies its ordinary analytic, noise, survey and saved
coefficients but finds one range-safety gap: positive `l_p*v` can underflow to
zero without rejection while final coefficients remain representable. W08
revision 2 rejects only that product underflow and passes independent scientific
review. Exact zeros and an exactly representable subnormal product remain valid;
ordinary analytic, legacy/supplied, noise, survey and saved-reference checks
pass. No consequential correction is needed.

W09 also passes independent scientific review: the public fixed-reference
noise/covariance/Fisher chain reproduces W07 for the original QSO bin. W10 also
passes independent scientific review: two additional population/redshift samples
retain roundoff-level reference refinement stability and smaller reference errors
than the historical cumulative results. The weak-LBG local-Fisher limitation
remains explicit. W11 passes independent review and detects reference-mode
sensitivity above 0.1% at 0.003 s/km. The user explicitly selected retention of
q_star=0.00035 s/km and planning accuracy-profile adoption. W12 revision 1 is
proposed for approval/dispatch: host-side fixed-reference routing, applicable
refinement controls with historical legacy semantics retained, and one saved
FF/FG/GG bin at two forest magnitude orders. This choice is not optimization or
physical-policy validation. Only W12 has current detailed instructions; W13
broader validation requires a separate request. No production change or numerical
run occurred during this planning update. Package Step 13 and acceptance of
historical failed forecasts remain unchanged.

- [Scientific design](../../FISHHIGHZ_DESIGN.md): conventions, architecture, and
  confirmed versus proposed decisions.
- [High-level roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md): sequence, review
  process, and acceptance register.
- [Current package step](IMPLEMENTATION_STEP.md): the detailed scope for the
  main package assignment.
- [Parallel forest-weighting diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md)
  and [current diagnostic step](WEIGHTING_DIAGNOSTIC_STEP.md): a separate
  user-requested plan/implementation/review sequence, created 2026-09-15.
  Select the assignment explicitly from the user's prompt; diagnostic work must
  not execute or replace the package Step 13 assignment. Both retain user control
  of approval, dispatch, acceptance and progression. Diagnostic-specific scope,
  checks and execution bounds belong in its current-step file. The new roadmap
  authorizes no automatic agent dispatch or scientific acceptance.

This file supplies standing context and working rules. Keep step-specific tasks,
acceptance criteria, and feedback in the current-step file, not here. Follow the
user's latest explicit instructions if they change a document's requirements;
report inconsistencies so the planning agent can reconcile the documents.

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
Revision 6 completes the reviewed synthetic/offline evidence repair. Step 13
investigates the same cumulative rule's normalized limit using a diagnostic
prototype and saved samples. It must retain amplitude in the nonlinear update
and test iteration and grid limits separately; it must not assume a finite limit.
Production weighting/guards remain unchanged in this assignment. A further real
forecast requires an explicit user request. Do not run the other six cases or
recapture NewForecast. Preserve r4 optimizations and regression tests. Weighting/input
policy changes require the user's scientific decision.
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
