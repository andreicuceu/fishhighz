# FishHighz agent guidance

## Project and authoritative documents

FishHighz is a new, standalone forecasting package in this directory. It computes
Fourier-space Gaussian covariances and Fisher forecasts for combinations of
galaxy/quasar and Lyα forest fields, with interchangeable external P3D and P1D
providers. The initial scaffold is implemented and under review; numerical
forecasting functionality has not been implemented yet.

- [Scientific design](../../FISHHIGHZ_DESIGN.md): conventions, architecture, and
  confirmed versus proposed decisions.
- [High-level roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md): sequence, review
  process, and acceptance register.
- [Current step](IMPLEMENTATION_STEP.md): the only detailed implementation scope
  for the current assignment.

This file supplies standing context and working rules. Keep step-specific tasks,
acceptance criteria, and feedback in the current-step file, not here. Follow the
user's latest explicit instructions if they change a document's requirements;
report inconsistencies so the planning agent can reconcile the documents.

## Roles and review control

The user chooses and prompts implementation agents and retains control over plan
review, implementation review, acceptance, and progression. An instruction to
implement the current step includes its tests. Complete that bounded assignment,
report evidence, and stop for review; do not ask for an extra permission between
implementation and its already authorized tests.

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
compilation; a JAX backend is not an initial requirement. Other recommendations,
such as the specific angular quadrature and fixed-covariance policy, remain to
be resolved before the relevant numerical step. Do not turn recommendations
into accepted requirements without the user resolving them. Read the design
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
