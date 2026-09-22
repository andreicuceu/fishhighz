# FishHighz design proposal

Status: Steps 01–11 accepted for progression; Step 12 evidence repair retained.
Step 13, W12 and the previous weighting diagnostic plan are done and closed by
user decision on 2026-09-17. The five-stage compatibility weighting study is
complete and early lyaforecast weights are the accepted baseline; McDonald
weights remain an explicit alternative. S1–S3 now implement the three profiles,
qualify revised accuracy and complete the comparison, with passing independent
scientific reviews. S4 attribution is complete with a passing independent
scientific review. Historical
failed forecasts and unresolved continuum conclusions remain unchanged.
Created: 2026-09-10. Status updated: 2026-09-21.

This is the working design document for FishHighz. Update the relevant sections
and decision register when feedback changes the design; keep the main proposal
consistent with the latest decisions. Section 0 and sections marked implemented
or reviewed describe existing APIs; proposed interfaces and the target module
diagram describe future architecture.

Implementation sequencing and user-controlled review gates are maintained in
[FISHHIGHZ_IMPLEMENTATION_PLAN.md](FISHHIGHZ_IMPLEMENTATION_PLAN.md). The package
directory is `lib/fishhighz`; its
[former current-step plan](../notes/IMPLEMENTATION_STEP.md) is retired.
The [compatibility weighting test plan](FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md)
records the completed weighting study; AGENTS.md (local-only path: `lib/fishhighz/AGENTS.md`) supplies
standing guidance. These planning documents do not constitute a package implementation.

## 0. Start here: current implementation and handover

### Current planning decision, 2026-09-18

The [weighting decision](../../research/FOREST_WEIGHTING_DECISION.md) records the
completed study and accepted early-lyaforecast baseline. The revised
[roadmap](FISHHIGHZ_IMPLEMENTATION_PLAN.md) records full-compatibility,
fixed-compatibility and updated accuracy profiles, completed S1–S3 tests and
comparison, and the completed S4 attribution.
The user confirmed retaining the legacy representative weighting mode initially,
with a separate fixed-comoving-mode test, and a 1% BAO-error attribution trigger
with a 0.1% numerical-refinement target. Updated accuracy uses the early recurrence
with its own fiducial per-field model, response and magnitude quadrature.
All future forecasts exclude bin-1 correlations involving LBG, LAE or lya(lbg):
only lya(qso) auto, QSO auto and their cross remain in that bin's forecast and
joint covariance. Bins 2–6 retain all 15 spectra. Historical results are unchanged.
S1–S3 implementations and their independent scientific reviews are complete:
[S1](../reviews/s1-review.md),
[S2](../reviews/s2-review.md), and
[S3](../reviews/s3-review.md) all pass. All six revised-accuracy bins
meet the finite-refinement criteria; maximum individual/joint error changes are
0.0011411%/0.0012685%. The [comparison](../reviews/s3.md) finds joint
accuracy/fixed radial reductions of 2.6693–3.8610% and transverse reductions of
0.6665–1.6089%, exceeding the attribution trigger in every joint radial result.
S4 attribution is complete using Astra (light) implementation and a
[passing Sol (high) scientific review](../reviews/s4-review.md).
The [controlled impact table](../reviews/s4-impact-table.md) identifies
mixed forest–galaxy damping as the dominant isolated reduction (joint radial
2.376–4.105%, transverse 0.954–2.383%), partly opposed by the grouped AP
derivative change. The 364 direct controls recover the endpoints within
0.000129% in BAO errors and have a largest tested refinement change of 0.010783%.
Effects are conditional and nonadditive; the study does not establish the
physical preference of the mixed-pair reconstruction prescription. Scientific
acceptance remains with the user. S5 began on 2026-09-21 as a joint planning
assessment of remaining sensitivity, scientific interpretation and package
completion; its proposed priorities and open decisions are in the roadmap.
The user confirmed a research-ready Python API as the first milestone and the
tested S2–S4 accuracy recipe as its research baseline, retaining mixed-pair
damping and full wiggle AP derivatives. S5 now consists of documenting that
prescription, aligning Python API documentation/examples, and validating the
changes. The reference-mode comparison is deferred to a separate later study in
[DEFERRED_SCIENTIFIC_TESTS.md](../../research/DEFERRED_SCIENTIFIC_TESTS.md).
No further scientific decision is pending within this S5 scope. S5 is complete
with Sol (medium) implementation and a
[passing Sol (high) scientific, software and documentation review](../reviews/s5-review.md).
The [research baseline guide](../../research/RESEARCH_BASELINE.md) and
standalone joint-BAO example (local-only path: `lib/fishhighz/examples/research_bao_forecast.py`)
are available; the [handoff](../reviews/s5.md) records focused tests,
installed-distribution checks and preservation of the numerical implementation.
Deferred studies remain outside this completed assignment.
W12's explicit method and historical evidence remain intact.
Older profile prescriptions below describe their dated implementation and are
superseded where the revised roadmap specifies a new proposed behavior.

### Historical planning decision, 2026-09-17

The user closed Step 13, W12 and the entire previous weighting diagnostic plan;
none will be resumed. W12's fixed-reference implementation and bounded passing
review are retained. The next study holds compatibility inputs and forecast
conventions fixed while comparing sum/cumsum, iterative aliasing and the earliest
Python formula, first at three updates and then at demonstrated convergence.
The user selected five variants and individual plus joint 15×2pt plots, omitting
poorly constrained points according to the new test plan's 0.2 uncertainty cut.
Accuracy-profile reassessment and a new completion plan follow those results.

Subsequently the user authorized coordinator dispatch and progression through
stages 1–3 using Astra (Light) implementers and Sol (High) scientific reviewers.
All three reviews pass. The [test plan execution record](FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md)
links the implementations, convergence evidence and bounded adaptive prescription.
These additions are validation-only; production defaults and the accuracy profile
remain unchanged. The stage-4 BAO comparison has not been executed.

The current roadmap and separate test plan supersede older scheduling, dispatch
and pending-repair statements throughout this document. Scientific/API history
below remains reference material; closing an investigation does not establish
its unresolved scientific claims. The initial planning revision changed no code
or evidence; the subsequently authorized stage-1–3 work is recorded above.

### Historical implementation and planning record

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

Overview update, 2026-09-15: the [revision-2 handoff](../reviews/step-13.md)
and diagnostic source are present; all eight entries in its saved source/report
hash manifest match the live files. The handoff reports 1318 tests passed and
25 optional skips; these are implementation evidence, not new review checks.
Independent revision-2 review remains pending. The dated account below preserves
the preceding review history; its proposal/dispatch status is superseded by this
handoff. Production weighting and the failed Step 12 forecast verdicts remain
unchanged. This overview did not rerun tests or forecasts.

This section is the entry point for a fresh planning/review agent. The status at
the top and the roadmap's progress register distinguish implementation review
from the user's acceptance. As of 2026-09-14, Steps 01–11 are accepted for
progression following the user's Step 12 request. Step 10 implements
fixed forest weights/noise and galaxy sampling noise from prepared arrays. Its
[revision-1 review](../reviews/step-10-review-r1.md) found mixed
intermediate underflow could silently lose pixel-noise contributions. Revision 2
now diagnoses those operations explicitly and [passed independent review](../reviews/step-10-review-r2.md).
Step 11 revision 2 implements Python orchestration and standalone raw density/SNR
readers and [passed independent review](../reviews/step-11-review-r2.md).
R1 is resolved by explicit nonuniform raw redshift-cell widths and a separately
labeled legacy first-spacing conversion, with unequal forecast-bin coverage.
R2 is resolved by preserving input dtypes until scientific validation. INI
translation, CLI and serialization remain deferred. Step 12 revision 1 now
implements seven-case synthetic validation, one-bin real BAO validation, opt-in
legacy input policies and a real lyaforecast intrinsic-P3D adapter/example.
[Revision-1 review](../reviews/step-12-review-r1.md) identified R1
semantic evidence gaps. Revision 2 implements the user's full seven-case,
two-profile comparison with plots and attribution. Its
[independent review](../reviews/step-12-review-r2.md) passes 1059 quick
tests and all 78 saved primary C/F reconstructions. All 39 compatibility bins
match legacy; only five galaxy-only accuracy bins converge. The original five
R1 probes reject, but R2 requires convergence operands to be bound to recorded
trials; R3 retains the 34 unconverged forest bins and 46 unavailable diagnostics.
The user subsequently suspended the 15x2pt-only scientific repair and archived
its exact [revision-3 instructions](../reviews/step-12-instructions-r3.md).
[Revision 4](../reviews/step-12-performance-r4.md) implements batched matrix
validation/factorization, optional lazy Numba Fisher contraction and reuse of
prepared factors/study payloads at unchanged scientific settings. Its
[handoff](../reviews/step-12-performance-r4.md) and
[independent review](../reviews/step-12-review-r4.md) establish that
the compiled path meets both performance targets; NumPy meets the composed
3x target but misses the 5x matrix target. Review found no required code
corrections; the NumPy-only matrix target remains unmet. Revision 5 restored the
archived 15x2pt assignment and is implemented in
[the preserved r5 handoff](../reviews/step-12-r5.md). Its
[independent review](../reviews/step-12-review-r5.md) passes 1183 quick
tests (25 optional skips), all 12 saved C/F reconstructions and 180 individual
Fisher/error checks. Compatibility passes all six bins; accuracy passes none,
and 12 of 72 diagnostic variants remain unavailable. The fixed-weight magnitude
integral converges, but the cumulative update remains grid dependent. Its tested
zero-weight Jacobian spectral radius is below one; a scale-normalized noise limit
has not been established. Scientific weighting/acceptance choices remain open.
The r5 review found that aggregate pair-array norms allowed a weak spectrum's
1.9419% error change to be replaced by zero and falsely certified. The original
r2 mutation was rejected, and the submitted r5 arrays remained truthful.
[Revision 6](../reviews/step-12.md) implements per-spectrum evidence
binding and saved-trial reconstruction of metrics/verdicts. Its
[independent review](../reviews/step-12-review-r6.md) passes 1243 quick
tests (25 optional skips), 251 installed regressions and all six installed
examples. The original false-pass probe now rejects through validator, writer
and offline reader. Independent Decimal/range controls and all 12 saved joint
forecasts plus 180 individual-spectrum checks pass. R2 is resolved for this
contract; no further code correction was found. All historical forecasts, the
r3/r5 archives and r4 performance changes are preserved. No real forecast or
production fix occurred during that review. The user now chooses investigation
of a normalized/asymptotic limit of the same cumulative rule and requests
Step 13. Its diagnostic revision 1 is implemented and
[independently reviewed](../reviews/step-13-review-r1.md): 1303 ordinary
tests/25 skips, 121 installed regressions, six examples and independent all-sixty
fixed-grid/finite-step checks pass. The large refinement trend is verified, but
the report overstates the continuum result; source/attempt binding and scalar
range defects require S13-R1/R2/R3 repairs. All 52 r6 modules remain unchanged.
[Revision 2](../notes/IMPLEMENTATION_STEP.md) proposes corrections to this
same diagnostic step. Its approval and dispatch remain with the user. No new
estimator, production adoption, convergence waiver or real forecast is authorized.
Step 12 R2 remains closed; R3 and scientific acceptance remain open.

Read the workspace and package AGENTS.md, this section, the
[roadmap handover](FISHHIGHZ_IMPLEMENTATION_PLAN_ARCHIVE_2026-09-17.md#fresh-agent-handover), the
single current assignment, and the latest implementation/review reports before
working. The [reusable agent prompt](FISHHIGHZ_PLANNING_REVIEW_PROMPT.md) describes
the user-controlled role. Current source and tests must be inspected when
assessing an implementation; a report or older conversation is not evidence
that a new checkout contains the same code.

### Implemented API map

All paths below are relative to lib/fishhighz/fishhighz. Import from the named
modules; these are not promised top-level package exports. README examples and
the corresponding step handoffs document full signatures and validation.

| Module | Implemented entry points | Boundary for subsequent work |
| --- | --- | --- |
| fields.py | ObservedField, PairSelection | Fields currently hold identity/kind/physical label/background, not survey data. Preserve original field order, selected order, required covariance closure, and selected_to_required gather. |
| parameters.py | Parameter, ParameterRegistry, ParameterBinding, gather_local, map_jacobian | Global IDs are opaque; roles are metadata; bounds are not priors. Equality ties gather local values and sum Jacobian contributions. Fixed settings are outside the free vector. |
| grids.py | IntegrationGrid, gauss_legendre_grid | Fixed 2D tensor integration grid, C-order flattened paired nodes, h_fid, cuts, explicit weights, q_mode. Arbitrary paired provider queries need not be tensor grids. |
| models/protocols.py | P3D/P1D/Jacobian protocols and output validators | Plain callable contracts, not required inheritance. Intrinsic P3D and independent velocity-space P1D. |
| models/external.py | BoundParameters, P3DProvider, PreparedP3D, evaluate_p3d, evaluate_p1d | Explicit pair ownership with complete closure. One identical registry object per prepared collection. No implicit AP, noise, response, or default P1D. |
| derivatives.py | evaluate_derivatives, check_convergence | Returns required-pair powers/Jacobians, stencils and call counts. Mixed supplied/numerical derivatives; explicit steps; tied global perturbations; unused columns zero; convergence is opt-in. |
| models/templates.py | PowerTemplate, prepare_template, load_template | Signed smooth/wiggle templates, provenance, h conversion, explicit G, log-k cubic coefficients, dP/dk, strict domains. Astropy/SciPy are optional preparation dependencies only. |
| models/kaiser.py | Scaling, KaiserModel | One-bin built-in intrinsic callable. Fixed/free slots, explicit ties, galaxy f versus forest beta, separate component scaling, fixed widths/G. Production derivatives use the existing numerical engine. |
| geometry.py | BinGeometry, prepare_geometry, prepare_astropy_geometry, mode_counts, explicit unit converters | Integrated common-bin volume and fixed h_fid/units; immutable arrays; optional cosmology extra; no background state in derivative loops. |
| response.py | InstrumentResponse, prepare_response, velocity_response, pair_response, explicit width converters | Signed observed-coordinate field transfers; required-pair products; explicit sigma/FWHM and labeled legacy conventions. Noise is separate. |
| models/p1d.py | default_p1d, p1d_floor | Intrinsic velocity-space PD2013 with stationary low-k floor and zero local parameters; response/conversion remain explicit consumer operations. |
| weights.py | ForestWeights, AuxiliarySamples, prepare_forest_weights, sample_auxiliary, density_per_velocity | Fixed prepared arrays and per-field auto/P1D sampling; explicit cumulative iteration, supplied weights, or W08 fixed inverse-variance weights. W08 r2 passes independent review and closes the positive `l_p*v` underflow finding while preserving exact zeros and representable subnormals. Integral-product underflow remains explicitly rejected. |
| noise.py | ForestNoise, forest_noise, galaxy_noise, local_galaxy_density, prepare_noise | Explicit units and independence or full replacement; independent noise PSD; pixel/Poisson noise unsmoothed. No raw readers or overlap inference. |
| survey.py | BinSpec, ForestInput, PreparedBin | Step 11 revision 2 accepted: input dtypes preserved until validation; immutable fixed state and explicit ownership; object-array snapshots rejected. |
| forecast.py | prepare_bin, run_bin, run_forecast | Step 11 revision 2 reviewed: independent unequal bins, one registry, fixed noise/factors, node batching and one combined prior; accepted for progression. |
| adapters/legacy_inputs.py | DensityReader, SNRReader, sample_forest_readers | Optional survey/SciPy preparation. Explicit nonuniform redshift-cell widths or labeled legacy_first_spacing; uniform magnitude axis, strict domains, no floors. Revision 2 reviewed. |
| adapters/legacy_compat.py | LegacyDensity, LegacySNR, sample_legacy_forest | Step 12 r1 implemented: explicit floor/clamp policies and provenance; negative-density floor labeled as an extension. Strict readers unchanged. |
| adapters/lyaforecast.py | IntrinsicP3D | Caller-prepared intrinsic model; explicit pair routing, k/z domain and h conversion; independent P1D. |
| validation/ | cases/recipes, RealRecipe, accuracy/compatibility profiles, execute/check, schema, plots, trial replay | Step 12 r6 implements schema-3 per-spectrum evidence checks and trial-derived convergence metrics/verdicts; R2 passes independent review. Historical r2 has 5/39 accuracy passes; r5 retains six failed 15x2pt accuracy bins. Step 13 r1 implements validation/weight_limit.py; review verifies finite-step/fixed-grid numerics but requires scientific-verdict, source/attempt-binding and scalar-range corrections. Revision 2 is proposed; production weighting is unchanged. |
| covariance.py | combine_observed_power, gaussian_covariance | Accepts already observed signal/noise and explicit mode counts, in required-pair order. Validates total field PSD and returns selected-pair covariance. It does not prepare W, N, or survey volume. |
| fisher.py | factor_covariance, fisher_from_factors, fisher_matrix | Reusable Cholesky factors and mean-Jacobian contraction. No extra mode factor or covariance derivatives. |
| results.py | FisherResult, diagonal_prior, combine_results | Data/prior separation, rank/null diagnostics, fixing versus marginalization, independent data combination and one shared prior. No pseudoinverse-based finite errors for unconstrained directions. |
| kernels/ | covariance, fisher, derivatives, templates, kaiser, response, weights | Existing dedicated numeric functions. Private helpers in _arrays.py and _information.py support host validation; do not add I/O or configuration work to kernels. |

Implemented reference examples are
external_forecast.py (local-only path: `lib/fishhighz/examples/external_forecast.py`),
builtin_forecast.py (local-only path: `lib/fishhighz/examples/builtin_forecast.py`),
survey_primitives.py (local-only path: `lib/fishhighz/examples/survey_primitives.py`),
weighted_noise.py (local-only path: `lib/fishhighz/examples/weighted_noise.py`), and
survey_forecast.py (local-only path: `lib/fishhighz/examples/survey_forecast.py`). They demonstrate
model-to-Fisher and geometry/response/P1D composition plus generated-file survey
orchestration. Step 12 adds
desi2_synthetic.py (local-only path: `lib/fishhighz/examples/desi2_synthetic.py`) and a bounded real
validation worker. Revision 2 adds the full two-profile comparison and offline
plots. Revision 4 implements batched matrix checks and optional lazy Numba
contraction selected with FISHHIGHZ_FISHER_BACKEND=numba; the compiled extra
requires Numba/SciPy. NumPy stays the default. Revision 5 restores the archived
revision-3 scientific repair while retaining these performance changes.
Later module names in the architecture diagram remain planned where they go
beyond these reviewed primitives.

### Contracts to preserve when adding survey preparation

- P3D uses paired k in h_fid/Mpc and mu in [0,1], with power shaped
  (node,required_pair) in (Mpc/h_fid)^3. P1D separately uses k_parallel in s/km
  and power in km/s. Keep h_fid fixed even when external cosmology varies.
- Prepare response W(node,field), known sampling noise N(node,required_pair),
  and common bin volume V once. Form observed total powers
  T_ij=W_i*W_j*P_ij+N_ij and pass mode counts V*grid.q_mode to covariance.
  q_mode=k^2*w_k*w_mu/(2*pi^2). Preserve signed cross terms.
- For fixed W and known noise subtraction, observed mean derivatives are
  W_i*W_j*dP_ij/dtheta. Gather selected_to_required before Fisher assembly;
  the derivative engine returns the required-pair Jacobian, not that selection.
  If a future model parameter describes a residual stochastic/noise mean term,
  declare its mean derivative explicitly. Varying covariance alone supplies
  no information under this convention.
- Covariance blocks already contain mode normalization. Factor once at the
  fiducial model and reuse while varying means/parameters. Combine independent
  bin data matrices in the same global metadata basis before applying one prior
  and marginalizing. Preserve null directions and fail clearly when the requested
  marginalized result is unresolved.
- Auxiliary P3D/P1D nodes for forest weights are preparation requirements, not
  additional forecast modes. Freeze weights, response, covariance, cuts, and
  volumes for derivative runs. User-supplied P3D does not select or derive P1D.
- Models own their intrinsic AP/RSD/damping physics; survey preparation owns
  instrumental response and sampling noise. Apply each effect once. In
  particular distinguish a field transfer W from its power multiplier W^2.

### Remaining scientific decisions for upcoming plans

These are questions to resolve when the relevant step is requested, not new
approved defaults or a second detailed assignment:

| Upcoming work | Existing constraints | Choices still requiring concrete design |
| --- | --- | --- |
| Step 09: geometry, response, default P1D | Resolved by user, 2026-09-13: optional Astropy adapter alongside supplied backgrounds; integrated bin volume; velocity sigma/pixel widths with explicit wavelength/FWHM and labeled legacy conversions. Fixed h units and exact PD2013 floor remain unchanged. | Revision 1 is accepted for progression following the Step 10 request. See step-09.md and step-09-review-r1.md for implementation and review evidence. |
| Step 10: weighting and noise | Resolved by user, 2026-09-13: prepared normalized density/variance arrays; legacy cumulative iterations with explicit count and per-field auto-power, or supplied fixed weights; explicit independent sampling or full supplied noise. | Revision 2 is accepted for progression after the Step 11 request. Scientific choices are unchanged. Raw file readers are scoped in Step 11; overlap-derived noise remains deferred. |
| Step 11: orchestration | User confirmed Python API plus standalone raw readers, with INI/CLI/serialization deferred; optional SciPy adapter with explicit legacy normalization/interpolation and strict value/domain errors. | Revision 2 accepted for progression after the Step 12 request: R1/R2 resolved. Explicit raw cell widths plus labeled legacy first-spacing mode; fixed preparation, shared bindings and one combined prior. |
| Step 12: validation and compatibility | User confirmed wiggle-only ap/at, opt-in legacy policies/P3D adapter, validation-only exact-legacy algebra, and a converged existing-model accuracy profile with retained-input sensitivity tests. | R2 metric/trial binding is resolved by r6 review; R3 forest convergence remains open. Full r2 evidence and the exact r3 instructions are preserved. Revision 4 implements performance-only optimization and a conditional runtime estimate. Independent review passes for the compiled path; NumPy-only matrix target remains unmet. Revision 5 supplies the 15x2pt-only comparison and weight diagnosis; independent review finds weak-spectrum R2 binding incomplete. Revision 6 passes independent review of per-spectrum evidence binding and trial-derived metrics/verdicts. The user selected same-rule limit investigation in Step 13; no replacement estimator or scientific acceptance. Retain r4 optimizations. |
| Step 13 | User selected investigation of a normalized/asymptotic formulation of the same cumulative weighting rule. | Derive finite-step equivalence, investigate discrete and continuum limits and diagnose saved inputs. A finite usable limit is not assumed; production adoption, alternate estimator and new real runs remain user decisions. |
| Steps 14–15 | Further performance work if needed; standalone packaging | Former prospective Steps 13–14, renumbered after insertion of the weighting investigation. Preserve r4 optimizations; no detailed plan or dispatch. Final distribution/documentation scope includes deferred INI/CLI/serialization. |

Ask the user before fixing a consequential unresolved choice that changes the
scientific result or architecture. Continue independent source inspection and
plan preparation while awaiting an answer. Do not re-ask confirmed decisions in
section 1 or treat deferred choices as reasons to reopen the accepted cores.

## 1. Purpose and decisions

FishHighz will compute Gaussian covariances and Fisher forecasts for arbitrary
combinations of galaxy/quasar and Lyα forest measurements. External packages can
supply P3D models and their parameters, and optionally P1D models for forest
noise. Built-in Kaiser/BAO P3D and the existing lyaforecast P1D prescription will
make the package usable without an external tracer model.

The same Fisher engine should handle BAO dilations, template full-shape AP/RSD
parameters, and cosmological parameters evaluated by an external theory package,
with anything from a few biases to many EFT nuisance parameters. These are model
and parameter choices, not separate forecasting engines.

| Decision | Status | Consequence |
| --- | --- | --- |
| Arbitrary external models take priority; autodiff is optional | Confirmed by user, 2026-09-10 | No requirement that external packages use JAX or be traceable |
| NumPy reference implementation with dedicated Numba-friendly kernels | Confirmed by user, 2026-09-11 | Start with NumPy; add targeted compilation after validation/profiling; JAX remains optional future work |
| Optional Astropy and SciPy for template preparation | Confirmed by user, 2026-09-12 | Astropy reads FITS; SciPy prepares spline coefficients; repeated evaluation uses NumPy kernels, and the external-model base remains NumPy-only |
| Optional Astropy cosmology adapter alongside supplied H(z)/D_M(z) | Confirmed by user, 2026-09-13 | Step 09 provides lazy optional preparation through a cosmology extra; the core remains NumPy-only with explicit h_fid and no chosen default cosmology |
| Integrate the comoving volume element over each redshift bin | Confirmed by user, 2026-09-13 | One common volume from transverse distance and H(z); explicit evaluation redshift and quadrature order; legacy centre approximation is comparison evidence only |
| Explicit sigma/pixel inputs and wavelength/FWHM conversions, with labeled legacy resolving-power compatibility | Confirmed by user, 2026-09-13 | Apply sinc times Gaussian as a field transfer; distinguish ordinary FWHM resolving power from legacy c/R interpreted as sigma |
| Prepared normalized source-density and pixel-variance arrays for Step 10 | Confirmed by user, 2026-09-13 | No raw density/SNR readers or implicit target renormalization in the weighting core; bounded real-input checks stay separate |
| Explicit legacy cumulative forest updates or caller-supplied fixed weights | Confirmed by user, 2026-09-13 | Use forest auto-power once per field and freeze across pairs/derivatives; three updates are compatibility behavior, not a convergence claim |
| Explicit independent sampling or supplied full known noise | Confirmed by user, 2026-09-13 | No assumed independence or inferred overlaps; generated diagonal noise and full replacement are distinct, with physical PSD validation |
| Use 2D P(k, mu) as the forecast observable | Confirmed by user, 2026-09-11 | Evaluate models/noise on a 2D grid and accumulate independent Fourier-cell Fisher contributions; multipole compression is not an implementation requirement |
| Independent redshift bins and a common survey volume within each bin | Confirmed by user, 2026-09-10 | Covariance is independent by Fourier cell in the initial idealized survey treatment |
| Fourier-space forecasts only | Confirmed by user, 2026-09-10 | No correlation-function support now or planned for the future |
| Cosmological parameters are varied for Fisher forecasts | Confirmed by user, 2026-09-10 | No sampling likelihood or sampler integration |
| Gauss–Legendre integration in mu on [0,1] and within supplied k-bin edges; custom nodes/weights supported | Confirmed by user, 2026-09-11 | Prepare configurable quadrature outside kernels; keep integration distinct from template sampling |
| Fixed fiducial covariance and fixed measurement weights for the standard Fisher calculation | Confirmed by user, 2026-09-11 | Differentiate the predicted mean spectra, with no additional covariance-derivative information term |
| Built-in model reads a pre-decomposed linear matter template from a Vega-format FITS file | Confirmed by user, 2026-09-10 | Read K, PK, PKSB; derive the wiggle residual as PK minus PKSB |
| Gaussian BAO widths specified separately for each tracer and fixed throughout a run | Confirmed by user, 2026-09-10 | Store widths outside the forecast parameter vector; never differentiate them as nuisance parameters |
| Independent or shared wiggle/no-wiggle coordinate rescaling | Confirmed by user, 2026-09-10 | Support ap/at, alpha/phi, and alpha_iso/epsilon for both components |
| Lyα beta is a nuisance parameter; galaxy growth rate is a full-shape target | Confirmed by user, 2026-09-10 | Use different Kaiser factors for forest and galaxy fields |
| User sets k_min and k_max before a run | Confirmed by user, 2026-09-10 | Fixed cuts in fiducial Fourier coordinates; template evaluation may need a wider domain |
| Cross-spectrum BAO width uses the mean of the two squared auto widths | Confirmed by user, 2026-09-10 | Recovers each tracer's auto damping when the two indices coincide |
| User dispatches implementers and controls review and advancement one step at a time | Confirmed by user, 2026-09-10 | Keep a separate roadmap, one current detailed step, and standing agent guidance |
| Direct fishhighz/ source directory under lib/fishhighz, without src/ | Confirmed by user, 2026-09-11 | Restrict package discovery and verify installed artifacts independently of the checkout |

The previous `FULL_SHAPE_FORECAST_PLAN.md` concerns a Vega-based configuration-space
package. FishHighz has a different, confirmed scope: Fourier-space Fisher
forecasts with customizable models. The older configuration-space architecture
does not apply to this package.

## 2. Lessons from the current lyaforecast implementation

Inspected local `lib/lyaforecast` at commit `5abe8bc`; its working tree was clean.
This original assessment read code and configurations, without rerunning
forecasts or measuring backend performance. Step 02 subsequently captured the
legacy reference forecasts; the roadmap locates that historical evidence.

| Existing component | Keep scientifically | Change structurally |
| --- | --- | --- |
| `survey.py`, `tracer.py`, `spectrograph.py` | Density distributions, forest lengths, pixel/resolution and S/N inputs | Separate input loading and interpolation from numerical arrays; replace named-tracer restrictions with field descriptions |
| `weights.py` | Forest weighting and magnitude integrals | Accept sampled P3D/P1D and survey arrays; extract dedicated numerical functions |
| `covariance.py` | Volume/mode counting, shot noise, forest aliasing and effective pixel noise | Separate geometry, noise preparation, and observed field spectra |
| `fisher.py:gaussian_covariance_array_func` | Gaussian covariance between arbitrary selected pairs | Move to an independent array kernel with integer pair indices |
| `fisher.py:compute_fisher` | Derivative–covariance contraction | Replace the hard-coded two-parameter BAO construction with a general Jacobian and linear solves |
| `analytic_p1d_PD2013.py` | Existing analytic formula, including the low-k flattening | Expose it as the default standalone P1D provider |
| `power_spectrum.py` | Kaiser physics and optional observation smoothing | Separate intrinsic theory from observing response and from grid ownership |
| `forecast_new.py` | Survey orchestration and availability of covariance-required spectra | Replace pair-specific stateful objects and BAO-only outputs with a prepared forecast and general results |

Important details:

- lyaforecast already assumes independent redshift bins and a common volume
  across tracers within each bin. `Survey` reads one survey area;
  `Covariance._get_survey_volume` combines it with the shared bin depth and
  geometry, and `NewForecast` uses one mode-count array for all correlations in
  that bin. Volumes can differ between redshift bins. Its combined result joins
  correlations within each bin, without cross-bin covariance.
- The 15x2pt example has five observed fields: Lyα from QSO spectra, QSO, LBG,
  LAE, and Lyα from LBG spectra. Five fields give fifteen unique auto/cross pairs.
  The two forest samples may share physical model parameters while having
  different observing noise and response.
- A requested cross spectrum still needs both auto spectra for its variance.
  Covariance between two requested pairs can require additional cross spectra.
  Pair selection must never remove these dependencies.
- Current forest weights call both P3D and P1D at a reference mode and use three
  iterations. Custom models must feed this preparation step as well as the final
  forecast grids. Weighting choices should be explicit and configurable.
- Current density normalization differs for continuous and discrete tracers.
  Compatibility requires preserving those differences, not unifying them blindly.
- The current mean and noise paths use different redshift-centre conventions;
  weights and spectrographs are also pair-specific. A joint field covariance
  requires a consistent definition of each measured field. Resolve differences
  explicitly in validation rather than promising exact equality after redesign.
- Existing BAO derivatives subtract an eighth-order fit in log-amplitude and
  apply damping before a finite grid difference. This is a particular BAO recipe,
  not a general parameter derivative engine.

The package repository location is `lib/fishhighz`, with Python import `fishhighz`
and source files directly in `lib/fishhighz/fishhighz/` (no intervening `src/`).
Port small scientific functions with provenance and
license attribution, and use lyaforecast as a validation reference. Avoid making
its full stateful object graph a permanent runtime dependency.

## 3. Architecture: preparation, models, and numerical kernels

Historical architecture diagram (Mermaid source):

```text
flowchart TD
    A[Survey inputs and field definitions] --> B[Preparation and validation]
    P[Parameter definitions] --> B
    B --> G[Geometry, grids, pair indices, survey arrays]
    E[External or built-in P3D and P1D] --> M[Model evaluation and derivatives]
    G --> M
    M --> N[Response and noise kernels]
    G --> N
    N --> C[Gaussian covariance]
    M --> D[Observable Jacobian]
    C --> F[Fisher solves and accumulation]
    D --> F
    F --> R[Named results, priors, marginalization, diagnostics]
```

**Preparation layer:** readable Python objects, file/config parsing, units,
parameter names, validation, adapter construction, interpolation setup, and cache
management. Prepare numerical quadrature nodes and tabulated density/S/N inputs
here. No file I/O or object construction inside numerical loops.

**Model layer:** small contracts for external or built-in models. Ordinary Python
adapters can manage package-specific objects; JAX-compatible adapters can expose
pure numerical functions. FishHighz owns the requested evaluation coordinates.

**Kernel layer:** array/scalar functions for responses, weighting iterations,
integrals, noise, covariance construction, derivatives where possible, and Fisher
accumulation. No configuration parsers, tracer-name parsing, logging, or dynamic
plugin lookup inside these functions.

Target module organization (mixes implemented and future modules; section 0 is
the authoritative implemented API map):

```text
fishhighz/
    fields.py          # observed-field identity, pair selection, covariance dependencies
    survey.py          # survey/source inputs and validation
    geometry.py        # fiducial distances, volumes, unit conversions
    grids.py           # evaluation nodes, integration weights, domain planning
    parameters.py      # names, transforms, bounds, sharing, local/global mapping
    forecast.py        # prepare and run orchestration
    models/
        protocols.py   # P3D/P1D and optional derivative contracts
        templates.py   # read K/PK/PKSB, metadata, interpolation preparation
        kaiser.py      # Lyα b/beta and galaxy b/f factors
        bao.py         # fixed per-tracer widths and pair damping kernels
        scaling.py     # component-specific Fourier rescaling, three bases
        p1d.py         # lyaforecast PD2013 prescription
    adapters/          # optional external-package and legacy-input adapters
    kernels/
        response.py
        weights.py
        noise.py
        covariance.py
        fisher.py
    derivatives.py     # finite differences, supplied Jacobians, optional autodiff
    results.py         # named Fisher matrices and diagnostics
    io.py              # optional configuration and result serialization
```

No correlation-function projection or likelihood modules are needed. Backend
implementation files should be introduced when needed; do not build a general
array-backend framework before measuring kernels.

## 4. Observed fields, spectra, and external contracts

An `ObservedField` describes the measured sample: a stable ID, galaxy or forest
kind, physical-model label, selection, noise inputs, and response. A Lyα field
also identifies its background-source population. Physical parameters can be
shared across fields without merging their noise or counting the same data twice.

Use integer field indices and ordered pair tuples internally, never split strings
such as `lya(qso)_qso` inside kernels. With T fields there are T(T+1)/2 possible
pairs; the user selects any subset. Preparation computes its covariance closure.

Illustrative minimal contracts, evaluated one redshift block at a time:

```python
p3d(theta, z, k, mu, pairs) -> power       # shape (n_node, n_pair)
p1d(theta_1d, z, k_parallel) -> power_1d  # shape (n_node,)
jacobian(theta, z, k, mu, pairs) -> dp   # optional: (n_node, n_pair, n_parameter)
```

Here `k` and `mu` are equal-length coordinate arrays (a tensor grid may be
flattened into them), and `pairs` is a prepared integer array. A simple callable
for one pair can be adapted to the batched contract. Multiple model providers
may serve different pairs; routing happens outside compiled kernels. Each
provider declares its parameter dependencies and supported domain.

The preferred P3D output is intrinsic redshift-space clustering power before
instrument response and sampling noise. It may include physical EFT stochastic
terms if their ownership is declared. A model that already includes AP, smoothing,
or shot noise must declare that through its adapter so FishHighz applies each
effect exactly once. Providers must supply all required cross spectra explicitly;
geometric means of auto spectra are not a general substitute for an EFT model.

Implemented canonical Fourier convention: k in h_fid/Mpc and P3D in
(Mpc/h_fid)^3, with the reference h fixed for a run. Adapters convert external
units, including parameter-dependent h, explicitly. The data grid is always in
fixed fiducial coordinates. P1D's native contract uses k_parallel in s/km and
P1D in km/s, matching lyaforecast; conversion to comoving units belongs to geometry.
For a=H(z)/[(1+z)h_fid], k_velocity=k_comoving/a and
P1D_comoving=P1D_velocity/a. Instrument response is applied separately.

P1D providers are assigned per forest population, with a shared default. P1D is
an independent input to aliasing/noise; do not integrate a finite-range custom
P3D to construct it implicitly. Overlapping forest samples may require cross-P1D
or a directly supplied cross-noise model; diagonal P1D alone cannot specify that.

## 5. Coordinate grids, interpolation, and derivatives

Use two concepts: the **model evaluation grid**, chosen for accurate sampling,
and the **forecast integration grid**, chosen for converged covariance/Fisher
quadrature. They can be identical initially. Direct evaluation on requested nodes
is simplest; use interpolation only when it reduces expensive external calls.

The user sets run-wide `k_min` and `k_max` at preparation, with
0 < k_min < k_max, in the fiducial Fourier convention. Build integration nodes
and weights over that fixed interval. Parameter perturbations never move the
cuts, alter the observed data vector, or change its mode counts. Both rescaled
components may require model evaluations outside this interval; validate the
template's coverage of those coordinates and derivative stencils in advance.
Auxiliary evaluations used to prepare forest weights are also model-domain
requirements, not additional modes included in the Fisher integral.

For a smooth model, a reusable table over k and μ at each z can reduce repeated
work. Request enough domain padding for AP transformations and derivative
perturbations. Refine until derivatives and marginalized errors converge, not
just until model values look smooth. Preserve signed cross powers; never
interpolate their logarithms blindly. Reject out-of-domain evaluation instead of
silently clipping or extrapolating. Grid smoothness alone says nothing about
smoothness or numerical precision in cosmological parameter directions.

Separate two interpolation uses:

- Spatial interpolation can provide derivatives with respect to k and μ for AP
  or template dilation, subject to convergence.
- It cannot provide derivatives with respect to cosmology or EFT parameters
  from a single fiducial table. Those require new parameter evaluations, supplied
  derivatives, or an explicitly parameterized differentiable emulator.

Derivative strategies share one output contract:

1. Supplied analytic derivatives/Jacobians, with declared parameter ordering.
2. Central finite differences for arbitrary external models, with parameter
   scaling, bounds-aware stencils, and checks at multiple step sizes.
3. Optional JAX autodiff for compatible model components.

Allow mixed strategies: for example, external cosmology finite differences plus
analytic nuisance derivatives. For p independent parameters, basic central
differences cost about 2p perturbed evaluations plus the fiducial evaluation;
step convergence checks cost additional calls. Dependency-aware caching can
reuse a matter spectrum while varying only tracer nuisance parameters. Cache
keys must cover every relevant parameter, grid, redshift, and model setting.

Parameter definitions include fiducial value, derivative scale/step, bounds or
transforms, prior, and scope. Support global cosmological parameters and
field/pair/redshift-specific nuisance parameters through an explicit mapping.
Combine redshift contributions in the shared parameter basis before marginalizing
shared nuisance parameters. Add a shared prior once. Retain rank diagnostics;
do not silently turn a singular Fisher matrix into finite errors with a pseudoinverse.

### 5.1 Two-dimensional spectra versus multipoles

The user confirmed P(k,mu) as the forecast observable on 2026-09-11. Evaluate
model, noise, and covariance on a common two-dimensional integration grid.
Multipole compression is not part of the implementation requirements. The
comparison below records the rationale for this choice; multipoles P_ell(k)
are also Fourier-space observables and do not imply correlation-function support.

Under the assumed homogeneous volume and independent Fourier modes, the 2D
calculation solves one selected-pair covariance block per (k,mu) node and sums
the Fisher contributions. It does not require one dense covariance spanning the
whole grid. This directly accommodates the k_parallel=k*mu dependence of forest
P1D aliasing and instrument response, along with anisotropic BAO damping and AP
rescaling. Angular cuts, if later requested, would also be simple in this basis.

Multipoles instead integrate angular information. For the real, even spectra
currently in scope,

\[
P_\ell(k)=(2\ell+1)\int_0^1 d\mu\,
\mathcal L_\ell(\mu)P(k,\mu),\qquad \ell=0,2,4,\ldots.
\]

A complete, invertible angular basis change preserves Fisher information when
the same data and full transformed covariance are used. Retaining only a few
multipoles is a compression and need not preserve it. Low galaxy multipoles
can perform well for particular models/surveys, but their success cannot be
assumed for the forest and mixed-tracer forecasts here; see
[Taruya, Saito & Nishimichi](https://arxiv.org/abs/1101.4723).

An undamped Kaiser spectrum at identity AP scaling has a mean that terminates
at ell=4. That alone does not prove that conventional ell=0,2,4 estimators
retain all Fisher information: the optimal angular weights involve the inverse
covariance, which is generally anisotropic, and AP derivatives can contain higher
angular orders. Damping, response, and AP rescaling also produce angular
structure beyond a finite Kaiser polynomial. Higher multipoles can have zero
mean but nonzero covariance with informative multipoles; see the discussion in
[Inoue et al.](https://arxiv.org/html/2406.19669v2). The appropriate truncation
must be tested on marginalized forecast errors, not just the fiducial mean.

| Consideration | 2D P(k,mu) | Selected multipoles P_ell(k) |
| --- | --- | --- |
| Angular information | Retained to numerical quadrature accuracy | Potential loss from finite ell_max |
| Covariance structure | Small pair blocks per k and mu | Coupled (pair,ell) blocks per k; cross-multipole terms are required |
| Lyα noise/response and AP | Evaluated directly at the relevant coordinates | Usually still require angular quadrature before projection |
| External packages | Natural for models that return arbitrary requested coordinates | Convenient when a package or intended analysis supplies/uses particular multipoles |
| Computational cost | More angular evaluations but independent small solves | Fewer observables but larger coupled solves and projection work; not automatically faster |

With s selected spectra and L retained multipoles, the matrix dimension is s
for each 2D angular node versus s*L for each multipole k block. For the 15x2pt
case, these are 15x15 blocks versus 45x45 blocks for ell=0,2,4. The winning
runtime depends on angular resolution, model cost, and the implementation.

Even a multipole forecast generally still needs a 2D model/noise evaluation grid
to compute its covariance, unless the angular integrals are available analytically.
The observable basis and the model's sampling grid are separate design choices.
If multipoles are requested, form a fixed Legendre/quadrature projection A and
transform the mean, Jacobian, and covariance consistently:
m_ell=A m, J_ell=A J, C_ell=A C A^T. Do not project the inverse covariance, drop
cross-multipole covariance, or assume an external model's truncated multipoles
specify all angular information needed for covariance. Validate convergence with
both quadrature resolution and ell_max against the 2D forecast. No multipole
implementation is added to the current roadmap merely by documenting this option.

The user confirmed Gauss–Legendre quadrature in mu on [0,1] and in k within
user-supplied bin edges on 2026-09-11, with explicit weights and configurable
orders. These default nodes lie inside their integration intervals. Also accept
custom nodes/weights; midpoint rules remain available through this interface.
Choose resolution through derivative/Fisher convergence, especially near mu=1;
do not inherit lyaforecast's ten angular samples as a universal accuracy
requirement. Keep integration nodes distinct from any finer template interpolation
grid and prepare all nodes/weights outside the kernels.

The inspected lyaforecast implementation already uses 2D arrays of shape
(n_mu,n_k), not multipoles: `power_spectrum.py` creates linearly spaced k values
and the midpoints of uniformly spaced mu bins; `fisher.py:compute_fisher` loops
over mu and sums contributions over k. All seven current DESI-2 INIs specify
500 k values from 0.01 to 0.5 h/Mpc and ten mu bins over [0,1], giving mu=0.05,
0.15,...,0.95. These are 5,000 evaluation nodes per pair per redshift bin.

### 5.2 Step 03 numerical contracts

Step 03 implements the following conventions and has passed independent review.
Public class/function names are recorded in its
[implementation handoff](../reviews/step-03.md); these contracts do not
imply an implemented forecast engine.

- Preserve input field order and use opaque field IDs plus canonical integer
  pairs. Preserve selected-pair order; prepare the required covariance spectra
  and integer product/gather mappings separately.
- Bind ordered local parameter names to explicit global IDs. Repeated local
  slots tied to one global ID sum their Jacobian contributions. Fixed model
  settings remain outside the free vector; parameter transformations and priors
  are later work.
- Store separate tensor axes and dk/dmu weights; flatten `(n_k,n_mu)` in C order,
  with mu varying fastest. Model outputs have shape `(n_node,n_pair)` and local
  Jacobians `(n_node,n_pair,n_local)`. Custom integration rules initially use
  tensor axes; providers themselves accept arbitrary paired evaluation points.
- Store `q_mode=k²*w_k*w_mu/(2*pi²)`, giving `N_modes=V_fid*q_mode` for the
  positive-mu, real/even-spectrum convention. This quantity contains neither
  model dependence nor a second conjugate-mode factor.
- Copy and protect stored preparation arrays; normalize numerical boundary
  values to contiguous float64 and index maps to explicit integer arrays.
  Keep name/configuration processing outside standalone array functions.

Step 03 validates these conventions with synthetic identities, polynomial
quadrature checks, model stubs, and isolated packaging checks. No new legacy
forecast capture is needed for this contracts-only step.

### 5.3 Step 06 external evaluation and derivatives

Prepare explicit pair-to-provider routes covering every covariance-required pair
exactly once, with parameter bindings in a declared common registry. Plain
callables evaluate fixed paired k/mu nodes directly; spatial interpolation and
automatic unit/response transformations are deferred. Independent P1D invocation
uses its own binding and velocity coordinates and is never inferred from P3D.

Supplied local Jacobians map into global columns through existing equality
bindings. Numerical derivatives perturb a global parameter before gathering
local values, so tied entries always move together. A provider may declare some
global columns analytic and others numerical; every tied local contribution to
one global column uses a single strategy within that provider. Only affected
providers are evaluated for each numerical perturbation.

Numerical derivatives require explicit absolute steps from parameter metadata or
validated overrides. Prefer second-order central stencils, with second-order
one-sided stencils near bounds. Reject requests that cannot fit a complete
stencil or produce distinct float64 points; do not silently clip or shrink steps.
Use actual offsets and record methods/call counts. Step-size convergence studies
are explicit opt-in operations, with absolute and relative diagnostics, so normal
evaluations do not incur repeated studies automatically. Refinement preflight
must also reject consecutive requested steps that round to the same actual
perturbation points, before any model call; repeated stencils cannot substantiate
convergence. This guard was implemented and passed review in Step 06 revision 2.

The [Step 06 review](../reviews/step-06-review-r2.md) records its passing
synthetic checks and forecast example with supplied fixed response/noise and the
reviewed covariance/Fisher APIs. It does not add production survey orchestration
or recompute covariance and weights during derivative evaluation.

## 6. Built-in P3D and P1D

### 6.1 File-based linear matter template

The built-in P3D reads a linear matter spectrum already decomposed by
`lib/vega/bin/make_template.py` (local-only path: `lib/vega/bin/make_template.py`). Its FITS binary
table, extension `PK`, has columns:

| Column | Meaning |
| --- | --- |
| `K` | Wavenumber in the template's h/Mpc convention |
| `PK` | Full linear matter power spectrum |
| `PKSB` | Smooth/no-wiggle matter power spectrum |

Prepare P_nw=PKSB and P_w=PK-PKSB once. PK is not the wiggle-only column.
Do not fit a new broadband polynomial or run a wiggle/no-wiggle decomposition
inside FishHighz. Read the selected file once, retain its hash and metadata,
and prepare interpolation coefficients outside numerical kernels. Interpolate
the signed wiggle residual without a logarithm of its amplitude.

Both inspected reference files also contain negative smooth-component samples.
Preserve those values; the template loader does not enforce component positivity
or take logarithms of smooth amplitudes. The full observed-total covariance
boundary remains responsible for physical validity in the chosen forecast domain.

Both inspected examples in
`lib/vega/vega/models/Planck18` (local-only path: `lib/vega/vega/models/Planck18`) contain 814 rows
and these columns: `Planck18_z_2.406.fits` and `DESI-2024_z_2.33.fits`.
Their `ZREF` values are 2.406 and 2.33 respectively. Read reference redshift,
cosmology, and any available normalization/growth metadata explicitly; missing
metadata needed by a chosen operation must be supplied, never silently invented.
The generator does not write explicit column-unit cards, so the loader must
document the Vega unit convention and any conversion to h_fid units.

Provide a file-path input and a documented example template, with attribution
and distribution rights checked before bundling. Using this built-in model
requires no runtime Vega or CAMB installation: generation is an external step.
FITS reading belongs to the host preparation layer, not JIT execution.

At a forecast redshift z, use a declared amplitude factor
G(z)=[D(z)/D(z_ref)]², supplied or prepared from the fiducial background, for both
components. G is fixed for a template-based AP/RSD run. A template at the target
redshift needs no growth rescaling. A template header's F_ZREF is a possible
fiducial value at z_ref, not a rule that fixes all forecast growth-rate parameters.
Cosmology-dependent shapes and evolution remain the external model's responsibility.

#### Step 07 reviewed preparation and interpolation

Require explicit h_fid, with h_template read from H0/100 or supplied if missing.
Convert k_fid=k_file*h_template/h_fid and
P_fid=P_file*(h_fid/h_template)^3 once, preserving source provenance. Missing
required reference-redshift/unit metadata must be supplied; contradictory inputs
fail. This conversion is not AP scaling or cosmological evolution.

The implemented interpolator is a not-a-knot cubic spline in ln(k_fid), acting on
linear signed smooth and wiggle amplitudes. Prepare coefficients outside kernels;
array-only evaluation provides powers and dP/dk at arbitrary paired-query k
values, preserves order, and rejects all extrapolation. Closed-domain endpoints
are accepted. Explicit G scales both components and derivatives once; G defaults
to unity only at z_ref. No growth solver is added in Step 07.

The user confirmed an optional `templates` extra with Astropy for FITS and SciPy
for coefficient preparation on 2026-09-12. Retain a NumPy-only external-model base
path and array-only repeated evaluation; neither preparation library is imported
eagerly. The [Step 07 handoff](../reviews/step-07.md) documents
`prepare_template`, `load_template`, and `PowerTemplate.evaluate`, including
component ordering, metadata precedence, and validation conventions. The
[independent review](../reviews/step-07-review-r1.md) passed with no
required changes. No cosmological FITS asset is bundled. The Kaiser/BAO model
is implemented in Step 08, as described below.

### 6.2 Tracer-specific Kaiser factors

Define the factor for an individual field i as

\[
B_i(\mu,z)=
\begin{cases}
b_i(z)[1+\beta_i(z)\mu^2], & \text{Ly}\alpha\text{ forest},\\
b_i(z)+f(z)\mu^2, & \text{galaxy/quasar}.
\end{cases}
\]

Use B_i B_j for each auto/cross spectrum. Biases are nuisance parameters for
both tracer kinds, with negative Lyα bias allowed. Each forest beta is an
independent nuisance parameter unless explicitly shared between forest samples.
Galaxy fields use a common physical f(z) within a redshift bin, a main full-shape
target parameter, instead of independent galaxy beta nuisance parameters.
The effective galaxy beta=f/b is derived. Forest-only spectra have no direct
Kaiser sensitivity to f in this parameterization; a resulting unconstrained f
must be reported, not regularized into a finite error.

Keep template normalization fixed for forecasts targeting f directly. Optional
b sigma8 / f sigma8 parameter transforms must state that normalization; freeing
an additional matter amplitude changes the degeneracies and interpretation.
Varying f in the template RSD model does not implicitly change G, BAO widths,
or the background distances. In direct cosmological forecasts an external
provider can consistently link these physical quantities instead.

### 6.3 Fixed Gaussian BAO broadening

Each tracer has nonnegative input widths Sigma_parallel,i and Sigma_perp,i in
Mpc/h_fid. Values may be specified per redshift bin, but are fixed for that run,
kept outside the differentiable parameter vector, and never marginalized.
Distinct forest samples may share or override these fixed settings explicitly.

Follow the Gaussian convention in
`PowerSpectrum.compute_peak_nl` (local-only path: `lib/vega/vega/power_spectrum.py`):

\[
D_{ij}(k_\parallel,k_\perp)=
\exp\left[-\tfrac12\left(k_\parallel^2\Sigma_{\parallel,ij}^2
                         +k_\perp^2\Sigma_{\perp,ij}^2\right)\right].
\]

Confirmed cross-pair convention:
Sigma_axis,ij²=(Sigma_axis,i²+Sigma_axis,j²)/2. This gives each specified
auto-spectrum width for i=j and the geometric mean of auto damping factors
for a cross spectrum. Vega's cited routine accepts one width pair for a
correlation; it does not itself define this per-tracer combination rule.

Apply D only to the wiggle component. Zero widths recover the undamped model.
Require both widths explicitly in the first implementation. In particular,
do not copy Vega's optional inference Sigma_parallel=(1+f)Sigma_perp into the
evaluation path: varying f must not change the input widths. Any future helper
that estimates widths or incorporates reconstruction must run once using fixed
fiducial inputs and store the resulting widths before forecasting.

### 6.4 Separate wiggle and no-wiggle Fourier rescaling

Provide independent scaling blocks for `wiggle` and `no_wiggle`, with explicit
parameter ties when common scaling is wanted. Each block supports the three
requested parameterizations from
`lib/vega/vega/scale_parameters.py` (local-only path: `lib/vega/vega/scale_parameters.py`):

| Inputs for one component | a_parallel | a_perp |
| --- | --- | --- |
| ap, at | ap | at |
| alpha, phi | alpha / sqrt(phi) | alpha * sqrt(phi) |
| alpha_iso, epsilon | alpha_iso * (1+epsilon)² | alpha_iso / (1+epsilon) |

Thus phi=a_perp/a_parallel, alpha=sqrt(a_parallel*a_perp), and
alpha_iso=(a_parallel*a_perp²)^(1/3). All scale factors and phi must be positive,
and epsilon must exceed -1. Identity uses unit alpha/phi/ap/at and zero epsilon.
The inspected Vega module also offers `aiso_aap`; it is outside the three
requested options and is not part of this design.

Use Vega's distance-dilation definitions for parameter values and apply the
inverse mapping directly in Fourier space, separately for each component c:

\[
k_{\parallel,c}=k\mu/a_{\parallel,c},\qquad
k_{\perp,c}=k\sqrt{1-\mu^2}/a_{\perp,c},\qquad
k_c=\sqrt{k_{\parallel,c}^2+k_{\perp,c}^2},\qquad
\mu_c=k_{\parallel,c}/k_c.
\]

Use the Fourier volume prefactor
Q_c=1/(a_parallel,c*a_perp,c²) for each rescaled component. Independent component
scalings are a phenomenological template model; tied geometric scaling recovers
the usual AP transformation. This is a Fourier implementation of the parameter
definitions, not a call to Vega's correlation-function scaling code.

The complete built-in intrinsic signal is

\[
P_{ij}(k,\mu,z)=G(z)\left[
Q_{\rm nw}B_i(\mu_{\rm nw},z)B_j(\mu_{\rm nw},z)P_{\rm nw}(k_{\rm nw})
+Q_{\rm w}B_i(\mu_{\rm w},z)B_j(\mu_{\rm w},z)
D_{ij}(k_{\parallel,\rm w},k_{\perp,\rm w})P_{\rm w}(k_{\rm w})\right].
\]

Evaluate the Kaiser angular factors and Gaussian at each component's transformed
coordinates. The widths remain fixed numbers, although D changes when its
coordinates are rescaled. Instrument response and sampling noise are applied
separately in observed coordinates.

Supported setups include wiggle scaling with no-wiggle scaling fixed at unity
for a BAO-oriented forecast; independent scaling blocks for split full-shape
forecasts; or tied blocks for common AP scaling. Nuisance and target parameter
selection is explicit in each setup. Do not add an implicit polynomial
peak-extraction derivative or reuse lyaforecast's finite-k BAO derivative.

BAO alpha definitions involving the sound horizon must be documented separately
from geometric AP q definitions. For cosmological forecasts the geometry adapter
derives AP from the varied cosmology unless the external model already does so.
Survey volume, integration weights, selection, and measured coordinates remain
fiducial; changing the model does not change the hypothetical observed dataset.

#### Step 08 reviewed implementation boundary

The [Step 08 handoff](../reviews/step-08.md) documents KaiserModel and
Scaling, a prepared callable and scaling record using the existing P3D provider
contract and original observed
field indices. Prepare one instance per redshift bin with fixed z, G, h_fid, and
widths. Each bias, forest beta, shared galaxy f, and component-scaling coordinate
can be fixed or explicitly bound to an ordered local free slot. Existing global
equality bindings express tracer or component ties; field labels never imply
sharing. Different scaling bases can be used by the two components, but automatic
nonlinear ties between bases are not part of this step.

Reuse Step 06 numerical derivatives with explicit steps for the initial built-in
forecast path. Independent analytic Kaiser/dilation derivatives and Fisher basis
transformations validate it; a complete production analytic Jacobian is deferred.
Both component domains are checked at every parameter state. Step 08 documents
coverage for derivative stencils; generic whole-run domain preflight belongs to
later orchestration. Observed grid cuts and units remain the caller's explicit
fixed convention.

Synthetic forecasts use supplied response/noise and the accepted covariance and
Fisher APIs. Geometry/response/default P1D subsequently arrived in Step 09;
forest-noise preparation and realistic forecast validation remain later work. Individual read-only Vega routine comparisons complement
analytic tests without adopting Vega's full configuration or cache behavior.
The [independent review](../reviews/step-08-review-r1.md) passed with no
required changes; the user's Step 09 request subsequently authorized progression.

### 6.5 Default P1D

Retain the exact default PD2013 P1D formula and redshift-dependent low-k floor.
Apply pixel sinc and resolution kernels outside it. Step 09 implements a zero-local-
parameter callable compatible with evaluate_p1d and empty BoundParameters; custom
P1D remains independently selectable. The floor is the stationary point of the
unclamped fit: power and its first k derivative join continuously, while the
second derivative need not. Test this behavior without changing the formula for
autodiff convenience. The full coefficients and independent checks belong to the
current Step 09 assignment and the reviewed models/p1d.py implementation.

### 6.6 Step 09 reviewed survey primitives

The user confirmed the three choices below on 2026-09-13. Their implementation
satisfies [Step 09 revision 1](../notes/IMPLEMENTATION_STEP.md); the
[handoff](../reviews/step-09.md) and
[independent review](../reviews/step-09-review-r1.md) document APIs,
source identity, numerical validation, and installed-wheel checks.

- Supply independent H(z) and transverse comoving distance D_M(z) callables in
  km/s/Mpc and Mpc, or use a thin optional Astropy FLRW adapter. Explicit h_fid
  stays fixed and need not equal the cosmology's h. The adapter uses unit-aware
  H and transverse distance calls and feeds the same array preparation. No
  named default cosmology, CAMB dependency, or cosmology solver is introduced.
  The cosmology extra includes Astropy/SciPy with lazy imports;
  supplied-background geometry, response, and P1D require only NumPy.
- Integrate V=Omega*h_fid^3*integral[c*D_M(z)^2/H(z) dz] using NumPy
  Gauss–Legendre quadrature with explicit z_order. Geometry retains one explicit
  evaluation redshift, bin, common area/volume, fixed h_fid, and owned immutable
  sampled arrays. Use transverse distance for curved backgrounds. The legacy
  bin-centre/log-velocity-depth approximation is restricted to labeled reference
  checks; no equality to it is required for wide bins. Test analytic flat limits,
  curved Astropy volume checks, quadrature convergence, and unchanged Fourier cuts.
- With a_v=H(z_eval)/[(1+z_eval)*h_fid], response uses observed q=k*mu/a_v and
  W=sinc_unscaled(q*Delta_v/2)*exp[-(q*sigma_v)^2/2]. Delta_v is full top-hat
  pixel width; sigma_v is Gaussian one-sigma, both in km/s. W(0)=1 and signed
  sinc lobes are retained. Per-field settings remain separate from ObservedField
  identity. Required-pair products W_i*W_j multiply intrinsic signal and mean
  derivatives once; noise receives no automatic response factor. P1D consumers
  explicitly apply W^2 and the 1/a_v comoving conversion.
- Wavelength conversions use c*width_lambda/lambda_obs in the local narrow-width
  approximation. Gaussian FWHM converts by 1/[2*sqrt(2*log(2))]. Ordinary
  resolving power defined from FWHM gives sigma_v=c/[R*2*sqrt(2*log(2))]; a
  separately named legacy helper gives c/R as sigma. Use c=299792.458 km/s in
  new code; isolate the reference's 299800 km/s difference in comparisons.
  Ly-alpha examples explicitly use lambda_obs=1215.67*(1+z_eval) Angstrom.

Source/SNR loading, weighting, aliasing/sampling noise, cross-P1D/overlap models,
and survey-to-forecast orchestration remain later work. Step 09 checks primitives
and fixed-covariance composition with supplied noise; no real forecast recapture
is required. Revision 1 passed review with 555 quick tests plus independent
volume/P1D/response/Fisher and installed-wheel checks; the Step 10 request
subsequently authorized progression.

## 7. Noise and Gaussian covariance

Use independent local Fourier blocks, with a common survey volume for all fields
within each redshift bin, as confirmed by the user.

Construct a symmetric field matrix at each integration node:

\[
T_{ij}=W_iW_jP_{ij}+N_{ij}.
\]

W is the observed field response. The default galaxy diagonal noise is 1/n_i.
For a forest field, retain the weighted aliasing and effective pixel-noise
calculation from lyaforecast; schematically in its observed units,

\[
N_{FF}=P_w P_{1D,\mathrm{smoothed}}+P_N^{\mathrm{eff}}.
\]

Convert both terms consistently to the P3D convention. Apply each field's
response once to a cross spectrum and twice to its auto spectrum. Shared-object
shot noise and shared-sightline noise require explicit off-diagonal N_ij.
Independent sampling can justify zero cross-noise; it never justifies dropping
the signal covariance between different spectra.

Expose two noise entry points: full survey/source/instrument preparation, and
directly supplied effective noise arrays for users who already know N_ij. Both
feed the same covariance kernel. Weighting should define each observed field
consistently across its pairs. If separate pair estimators require different
weights, they need a more general estimator covariance than a single T matrix.

For selected pairs A=(i,j), B=(m,n),

\[
C_{AB}=\frac{T_{im}T_{jn}+T_{in}T_{jm}}{N_{\mathrm{modes}}}.
\]

For real, even spectra using positive μ and counting both conjugate hemispheres,
N_modes=V k² Δk Δμ/(2π²), matching the current convention. General quadrature
replaces cell widths with integration weights. Keep normalization in one place
and verify Var(P_ii)=2 T_ii²/N_modes. Complex/odd spectra would require an
explicitly extended estimator contract.

Support arbitrary selected-pair covariance as the baseline. For all pairs with
identical cuts, an optional field-level trace expression can be more efficient:

\[
F_{ab}=\sum_{\mathrm{nodes}}\frac{N_{\mathrm{modes}}}{2}
\mathrm{Tr}[T^{-1}T_{,a}T^{-1}T_{,b}].
\]

Its equivalence requires a complete set of pair observables and consistent mean,
noise, and parameter treatment. Do not use it for an arbitrary subset of spectra.
For the normal bandpower calculation use

\[
F=\sum_{\mathrm{nodes}}J^\mathsf{T}C_{\mathrm{fid}}^{-1}J+F_{\mathrm{prior}}.
\]

The user confirmed fixed fiducial covariance and fixed survey weights on
2026-09-11. Parameter perturbations affect the mean prediction; they do not
redefine measurement weights or append covariance-derivative information.

The mean convention must specify noise subtraction. Default known sampling noise
is excluded from the mean and fixed at its fiducial value. If a stochastic/noise
amplitude is forecast, its residual contribution must appear in the mean
derivative relative to a fixed subtraction; covariance changes alone do not make
it a constrained parameter in this formulation.

Do not automatically append the Gaussian covariance-derivative trace term to the
bandpower Fisher matrix. It can assign spurious additional information to
two-point estimators of Gaussian fields. A different likelihood/information model
needs a separate explicit design. See [Carron (2013)](https://arxiv.org/abs/1204.4724).
The forest noise starting point is the existing implementation of
[McDonald & Eisenstein](https://arxiv.org/abs/astro-ph/0607122).

The initial common-volume assumption excludes general survey-window coupling and
correlated redshift bins. Keep field/sample overlap distinct from volume overlap:
fields in a common volume can still share individual objects or sightlines, which
requires explicit cross-noise. The simple default assumes disjoint sampling and
zero cross-noise; user-supplied cross-noise arrays permit other specified cases.
Correlation functions and their projected covariance are outside FishHighz's
scope, including future scope.

### 7.1 Step 04 covariance boundary

The reviewed implementation accepts supplied observed total powers in
`PairSelection.required_pairs` order and explicit positive mode counts
`V_fid * grid.q_mode`. A small composition helper adds supplied observed signal
and noise with matching shapes. Response and survey/noise preparation remain
separate later modules. Covariance output has shape
`(n_node, n_selected, n_selected)` in the selected-pair order, and reuses the
existing four pair-lookup arrays in a dedicated array-only numerical kernel.

Boundary validation checks the active total field-power matrix for positive
semidefiniteness, using normalized correlations and a documented roundoff
tolerance. Valid singular matrices, including zero-power fields, are allowed
without regularization; the later Fisher step must assess solvability of the
selected covariance. Invalid inputs receive node/field diagnostics. Validation
and eigensystem checks stay outside the kernel, and callers can process bounded
node slices. The [Step 04 review](../reviews/step-04-review-r1.md)
records the passing analytic and synthetic checks; no real reference forecasts
were required for that step.

### 7.2 Step 05 Fisher and result boundary

The reviewed implementation takes supplied selected-mean Jacobians in the global
parameter basis and fixed covariance blocks. It separates reusable Cholesky
factorization from triangular solves and accumulation of `Y.T @ Y`, where
`Y = solve(L, J)`. Covariance already contains mode normalization; accumulation
adds no extra weights. Numerical work lives in dedicated functions with array
inputs. Singular or numerically unresolved selected covariance receives an error
with node diagnostics; no automatic observable projection is introduced.

Named results keep data information separate from prior information and retain
the existing parameter registry. Independent data contributions must use the
same ordered parameter metadata, with zeros for bin-inactive dependencies.
Combine unmarginalized data first and apply a shared prior once. The initial
combination API rejects contributions already carrying nonzero priors. Finite
diagonal or correlated PSD prior information is supported; bounds are not priors.

Singular Fisher information remains inspectable with scale-aware rank and null
directions, while joint marginalized errors require an invertible retained
matrix. Conditional errors, marginalization over other free parameters, and
explicit fixing through a principal information submatrix are distinct
operations. No finite pseudoinverse errors or implicit regularization are used.
The [Step 05 review](../reviews/step-05-review-r1.md) records its scope
and passing analytic checks. Step 06 subsequently implemented external-model
invocation and derivative generation, as described in section 5.3.

### 7.3 Step 10 fixed weighting and noise (revision 2 reviewed)

The user confirmed the preparation boundaries on 2026-09-13. Revision 1 is
implemented. The [revision-2 review](../reviews/step-10-review-r2.md)
closes R1 from independent review: inexact underflow in individual integral
products raises a field-context error before another positive term can conceal
the loss. Exact-zero and exactly representable subnormal products remain valid.
Revision 2 is accepted for progression after the Step 11 request. The scientific
conventions below are unchanged.

Use one prepared source sample per observed forest field/bin: ordered magnitude
nodes and explicit positive quadrature weights, already normalized density rho
in deg^(-2) (km/s)^(-1) mag^(-1), pixel-noise variance s2, one positive velocity
forest length L_v and full pixel width Delta_v, and the accepted field response.
Source redshift describes the density/variance sample; the common bin z_eval
describes the clustering/P1D and geometry. Convert a supplied differential
redshift-density row with (1+z_source)/c explicitly. No target-density or raw
file-normalization policy is inferred. Keep legacy rectangular magnitude sums
available through explicit weights; do not silently change their endpoints.

For a legacy iteration, initialize w=B/(B+Delta_v*s2) and update simultaneously
using w_new=S/[S+s2/(I1_m*L_v/Delta_v)], where I1_m is the **cumulative** weighted
source-density integral. S is the field's positive intrinsic auto-power sampled
at an explicit auxiliary mode, multiplied by its W^2 and a_v/d_deg^2; B is
independently supplied P1D times W^2 at that mode. Use existing callable bindings
and original field indices. Three explicit updates reproduce the legacy recipe;
no convergence or joint optimality is claimed. Supplied nonnegative fixed weights
bypass auxiliary model queries. Zero-density bins have zero iterative weight;
no artificial density or undefined leading-prefix ratios.

With quadrature masses r=rho*q_m, compute cumulative I1=sum(r*w),
I2=sum(r*w^2), and I3=sum(r*w^2*s2). Final coefficients are
A=I2/(I1^2*L_v) in deg^2 and P_pixel=I3*Delta_v/(I1^2*L_v) in deg^2 km/s.
At observed q=k*mu/a_v, forest noise is
[A*P1D(q)*W(q)^2+P_pixel]*d_deg^2/a_v in (Mpc/h_fid)^3. Pixel noise is not
smoothed again. Galaxy Poisson noise is 1/n_bar with supplied comoving density;
a separately labeled local-density helper converts dndzdm at z_eval. Bin-averaged
count/volume inputs are the caller's explicit choice, not silently inferred.

Prepare one weight/noise definition per field, shared across its pairs. Generated
diagonal noise requires an explicit independent-sampling declaration. Alternatively
accept a full replacement in required-pair order; validate noise PSD independently
of signal with the accepted normalized field-matrix conventions. Do not infer
cross-noise from labels or allow positive signal to conceal indefinite noise.
Known noise, weights, response, geometry and covariance remain fiducial; only the
intrinsic mean is differentiated. Step 11 now scopes raw density/SNR readers
and orchestration; overlap-derived noise remains later work. Step 10 completed
a bounded DESI-2 source/variance comparison separately from synthetic tests;
that comparison ran no forecast.

### 7.4 Step 11 implemented orchestration and revision-2 raw-input boundary

The user confirmed on 2026-09-13 a Python API plus standalone raw density/SNR
readers, with INI translation, CLI and serialization deferred. The detailed
[revision-2 assignment](../notes/IMPLEMENTATION_STEP.md) describes the implemented
host-side survey/forecast composition and labeled optional SciPy input adapter.
Revision 2 passed independent review, resolving R1 and R2, and is accepted for
progression after the Step 12 request.

The adapter explicitly interprets legacy three-column density values as cell
counts, applies caller-declared raw magnitude masks and target-normalization
windows, divides row i by Delta_z[i]*dm, and prepares the legacy quadratic
(kx=ky=2) spline. The user confirmed nonuniform raw redshift nodes with explicit
positive per-cell widths in sorted redshift order. Do not infer edges from centres.
An explicitly labeled legacy_first_spacing mode instead uses z[1]-z[0] throughout,
including on irregular reference grids; it is not an automatic physical-width
estimate. Existing uniform-z conversion remains valid on verified uniform grids.
Keep magnitude spacing uniform. Record axes, widths and policy in provenance;
normalization sums raw counts, independently of forecast-bin boundaries.
No field name chooses normalization. SNR files are keyed by validated headers;
explicit wavelength-axis smoothing (legacy sigma=10 samples, reflect boundaries)
precedes linear interpolation. Pixel variance uses SNR per Angstrom and explicit
pixel width/exposure count. Domains and invalid interpolated values fail without
legacy floors, bright clamps, extrapolation, or automatic node removal. Paths,
hashes, normalization measures, smoothing and exposure metadata are retained.
SciPy stays in an optional survey extra; normalized-array orchestration is NumPy-only.

Forecast bins may have unequal widths, non-midpoint evaluation redshifts and
gaps; integrate each actual interval and retain nonoverlap/interior checks. They
are independent of raw source-table cells. The current runner already supports
unequal widths; revision 2 makes this explicit in tests and the example.

Prepared independent bins retain the existing geometry/grid, original pair
selection, one global registry with explicit per-bin bindings, per-field response,
and either generated noise with declared independence or a full replacement.
Caller-specified source redshifts/lengths and local versus supplied galaxy density
remain explicit. Fixed weights/noise, modes and covariance factors are reused
across intrinsic-mean derivative runs. Built-in and external models share the
same execution path; node batches preserve full pair covariance. Sum bin data
Fishers before adding one prior or marginalizing; singular bin information is
allowed. Retained factors have O(n_node*n_selected^2) storage, while transient
Jacobians are batched. No disk cache or performance claim is introduced.

Ownership copies must not coerce invalid complex/Boolean/string/object scientific
arrays into valid float64 inputs before validation. Preserve metadata dtypes;
accepted numerical state remains immutable float64 (R2).

Acceptance includes generated-file normalization/interpolation oracles, a two-bin
shared-parameter analytic Fisher, batch/permutation equivalence, fixed-state call
counts, strict failure contexts, new-reader comparisons at Step 10's bounded
five-population input point, and installed-wheel/blocked-optional checks. No real
full forecast is authorized by the Step 11 assignment.

### 7.5 Step 12 validation, comparison profiles and compatibility boundary

Step 12 revision 1 implemented explicit seven-case recipes, seven synthetic
selections, one real 15x2pt bin, opt-in legacy input policies and an intrinsic
lyaforecast P3D bridge. Review R1 concerns scientific evidence validation, not
an observed failure of the saved one-bin Fisher. Preserve original reports and
artifacts. Revision 2 implemented the expanded comparison. Its independent
review verifies the final numerical contractions and closes the five original
R1 probes, but finds R2 metric/trial binding missing. Weighting stabilization
fails in all 34 forest bins, with magnitude convergence also failing in 23;
46/468 diagnostics are unavailable. The user first narrowed revision 3 to
15x2pt, then archived it during performance-only revision 4. Revision 5 implemented the restored one-case assignment, retaining the reviewed
optimizations. Its six accuracy bins remain unconverged; the fixed-weight
quadrature and cumulative-weight diagnosis are independently verified. The r5 R2 defect allowed false convergence for a weak spectrum because aggregate
norms concealed its changed operands. Revision 6 implements per-spectrum binding
and trial-derived metrics/verdicts and passes independent synthetic/offline review. Weighting/input choices, new real
runs and acceptance remain with the user.

The user confirmed wiggle-only ap/at per bin with fixed smooth scaling, biases,
beta/f/G, widths, covariance, weights and observed response. The accuracy profile
uses the accepted full KaiserModel wiggle mapping, including volume prefactor,
remapped RSD and fixed-width damping. The validation-only compatibility profile
instead reproduces the legacy peak estimator and upstream model/J/mode/pair
inputs, then independently assembles Wick covariance and FishHighz Fisher
information. Never install legacy estimator switches into default scientific
APIs or reuse legacy final F as the independent result. Diagnose any legacy
field-PSD inconsistency; a literal validation-only Wick matrix may be used if
selected C is SPD, without relaxing production validation or adding jitter.

The user chose a converged existing-model accuracy configuration: integrated
volumes, refined k/mu/magnitude quadrature, full wiggle derivatives, physical
FWHM-to-sigma resolution, CAMB power growth, consistent per-field noise and
accepted cross damping. Keep original samples, cuts and bin edges; retain the
explicit single-z_eval approximation and independent sampling assumption.
Assess stabilization of the existing cumulative weights; no new optimizer or
claim of joint optimality. New overlap noise, within-bin evolution and other
unimplemented realism remain outside this step. Accuracy is conditional on the
existing physical model and raw data, not an absolute survey-realism guarantee.

Step 11 strict readers remain default. The separate opt-in adapter retains
out-of-magnitude density 1e-20 and labeled spline redshift-domain extension, SNR
bright clamp, out-of-domain RMS 1e10 (variance 1e20), and post-pixel/exposure SNR
floor 1e-10. Negative interpolated density flooring to 1e-20 is an explicitly
approved extension beyond actual legacy behavior; exact zeros/positive values
stay unchanged. Malformed inputs/nonfinite arithmetic still fail. Record raw
and effective coordinates/values, masks, policy constants and source identities.
The user approved retaining these policies and legacy_first_spacing conversion
where calibrated coverage/physical widths are unavailable, with sensitivity
tests. Diagnostic perturbations are not calibrated systematic uncertainties;
never infer physical cell edges from centres or silently change primary samples.

The intrinsic-P3D adapter remains a thin wrapper around a caller-prepared object
with explicit domains, pair routes and h units. Its real amplitude example is
separate from BAO validation, with independent P1D and no eager reference imports.
Pair-specific legacy behavior belongs only to the compatibility validation path,
not this adapter or the physical per-field survey API.

The user explicitly requested all seven real benchmark cases and two FishHighz
profiles per case. Revision 2 produced a fresh full legacy reference plus all 39
original bins in each FishHighz profile, selected-pair diagnostics, numerical
convergence, retained-input sensitivity, plots and controlled explanations of
material differences. The scientific acceptance gate fails as described above;
revision 5 supplies revision 3's 15x2pt-only inventory: 12 case/bin/profile
records, 180 individually selected-spectrum results and 72 diagnostic requests
(60 completed, 12 unavailable). Six compatibility records pass; no accuracy
record converges. The legacy reference is reused; the other six cases and
reference recapture are excluded. Revision 6 completes and passes the offline/synthetic evidence repair; any new
real run requires an explicit user request. Preserve
historical baseline/one-bin evidence separately from new runs. Require semantic
case/profile/bin/pair/parameter/node validation and independent saved-array
consistency checks, not only hashes and self-reported shapes. The implementation
must not claim full completion for missing/failed records or unexplained material
differences. General INI translation and production CLI/serialization remain
deferred. R3 transfers to the user-requested Step 13 investigation. See the
[current assignment](../notes/IMPLEMENTATION_STEP.md).

### 7.6 Step 13: normalized weight-limit investigation and review

The user selected option 2 on 2026-09-14: investigate the existing cumulative
rule before choosing a different estimator or redefining its convergence
requirement. Step 13 revision 1 is implemented and requires corrections after
independent review; revision 2 is proposed for user approval/dispatch. The
production weighting method remains unchanged.

With r_i=rho_i*q_i, alpha=L/Delta_v and a_i=variance_i/S, the accepted map is
w_i(next)=N_i/(N_i+a_i), N_i=alpha*sum(j<=i,r_j*w_j). Its noise ratios
A=sum(r*w^2)/(L*sum(r*w)^2) and
P_pixel=Delta_v*sum(r*w^2*variance)/(L*sum(r*w)^2) are invariant under common
weight rescaling. The map itself is not. Any stable equivalent representation
must retain amplitude as well as shape; setting the maximum weight to one at
each original update changes the finite recurrence.

The r1 implementation provides an exact amplitude/shape diagnostic. Independent
checks verify finite-step coefficients and all sixty fixed-grid limits. The
measured P_pixel refinement exponent is 0.9615–0.9798 over orders 16/32/64,
but a continuum theorem for all real populations is not established by these
three-order screens. Revision 2 requires a controlled asymptotic argument or an
explicitly unresolved continuum conclusion, source/attempt binding and scalar
range corrections. A bounded normalized shape can concentrate on
vanishing source measure and yield divergent noise coefficients; normalization
alone is not a convergence proof. The constant-coefficient linearized continuum
control in the planning note illustrates this obstruction without asserting it
for the full nonlinear DESI-2 rule. Small toy checks during planning passed;
no actual forecast or new sample preparation ran.

Use saved arrays for both forest populations in all six 15x2pt bins. Prototype
numerics remain validation-only, with preserved zero/range semantics and
provenance. Report a supported finite grid-independent limit, a demonstrated
failure of that condition, or an explicitly unresolved bounded investigation.
A later user decision governs production integration, a replacement prescription
if needed, and any real forecast validation. R2 remains closed; R3 and scientific
acceptance remain open. No density floor, support, P1D, response, covariance,
parameter or existing convergence threshold is changed by this investigation.

## 8. Performance and compilation

The user brought targeted performance work forward to Step 12 revision 4 on
2026-09-14. It is now implemented and independently reviewed for its optional
compiled path; revision 5 retains it as the implementation baseline. Matrix checks/factorization batch
at most 256 cells; failed or near-threshold batches replay scalar diagnostics.
The compiled Fisher contraction retains one-cell solve workspace and finite
checks. NumPy remains unconditional/default; the optional compiled extra adds
Numba/SciPy lazily. Accuracy reuses prepared factors and at most three study
payloads while retaining a direct NumPy oracle and the trial inventory.
At 2048 synthetic cells the composed gains are 5.77x NumPy / 19.21x compiled;
matrix gains are 2.77x / 17.56x, leaving the NumPy-only 5x matrix target unmet.
The handoff includes subsequent profiles and conditional runtime predictions;
1167.37 s accuracy and 12.61 s compatibility historical records cover unequal
workloads. No real forecast ran under the r4 performance assignment; revision 5
subsequently used the compiled path for its bounded scientific assignment.
Revision 6 preserves these optimizations without a new forecast. See the r4 handoff and
independent review for measurements, numerical checks and extrapolation limits.

Optimize data flow before adding decorators:

- Prepare geometry, quadrature, density/S/N tables, field responses, and integer
  pair maps once. Evaluate all needed pairs in batches where supported.
- Read the matter template and prepare signed-component interpolation once.
  Keep template I/O separate from the rescaling, Kaiser, and Gaussian kernels.
  Precompute fixed pair-width combinations, but recompute damping when rescaled
  coordinates change; caching by fixed widths alone would give stale derivatives.
- Cache fiducial covariance factors for repeated derivative and prior studies.
  Separate model, derivative, noise, and factorization timing.
- Process bounded batches of Fourier nodes; solve C X=J using Cholesky when
  positive definite, then accumulate JᵀX. Avoid explicitly forming C inverse.
- Reuse a factorization across all parameter columns. Stream parameter blocks
  if a large Jacobian exceeds the memory budget.
- Use float64 and diagnose non-positive covariance, near-degenerate fields,
  zero-density samples, and unconstrained parameter combinations explicitly.
- Resolve cuts into static index groups before compiled execution. Remove
  inactive observables before factorization; a zero integration weight does not
  make a singular covariance safe to invert.
- Keep aggressive floating-point reassociation/fastmath off initially.

With s selected spectra, dense covariance factorization costs roughly O(s³) per
node. Since s grows as T², the field trace route becomes useful for many tracers,
but fifteen spectra is a modest dense block. External cosmological models may
dominate runtime; compiling the Fisher contraction cannot remove that cost.

The end goal is to compile FishHighz-owned intensive computations. Compilation
of arbitrary external Python models is outside FishHighz's control. Those models
can still be batched, cached, or replaced by a user-supplied compiled model.

## 9. NumPy/Numba versus JAX

| Aspect | NumPy with optional Numba kernels | JAX kernels and optional native models |
| --- | --- | --- |
| Arbitrary external packages | Natural host-array interface | Evaluate externally, then transfer arrays, or supply a JAX-compatible adapter |
| Derivatives | Supplied derivatives or finite differences | Same options, plus autodiff through compatible numerical paths |
| Intensive loops | Numba-compatible standalone functions; compiled linear algebra where appropriate | Pure functions, static shapes, batching and JIT compilation |
| Dynamic setup and file inputs | Ordinary Python | Ordinary Python outside transformations |
| Numerical precision | Explicit float64 arrays | Explicitly enable/verify x64 before constructing numerical state |
| Repeated same-shape calculations | Low orchestration overhead after compilation | Compilation may amortize well; changed shapes can trigger recompilation |
| Accelerators | CPU-first approach | Accelerator execution possible, subject to workload, transfer, and float64 performance |
| Maintenance | Simplest initial interoperability path | Extra precision, compilation, device, and derivative validation |

Numba supports a subset of Python/NumPy, so arbitrary object methods and SciPy
interpolators cannot simply be decorated. Prepare tables outside the kernels and
use supported array operations inside them. Profile before replacing efficient
NumPy/BLAS calls. [Numba performance guidance](https://numba.readthedocs.io/en/stable/user/performance-tips.html).

JAX changes implementation discipline more than the scientific architecture:
numeric parameter arrays or pytrees, functional updates, fixed shapes, static
pair structure, and compatible control flow inside transformed functions.
Initialize and validate its float64 mode explicitly rather than silently changing
global JAX configuration at package import.
[JAX constraints and precision](https://docs.jax.dev/en/latest/notebooks/Common_Gotchas_in_JAX.html).

The useful JAX modes are distinct:

1. **External models plus compiled numerical kernels:** external P and derivatives
   are evaluated on the host and transferred in batches. This preserves the
   confirmed interoperability requirement.
2. **Native differentiable models:** compatible model, AP, response, and derivative
   calculations can be compiled together. For many model outputs and relatively
   few parameters, forward-mode Jacobians or blocked JVPs are reasonable starting
   points; benchmark their memory and runtime.
   [JAX Jacobian guidance](https://docs.jax.dev/en/latest/301/cookbook.html).

Calling external NumPy through `pure_callback` does not compile the external
code or supply its derivative. Custom derivative rules require additional work
and still depend on a trustworthy derivative implementation. Treat callbacks as
an adapter-specific option, not the default architecture.
[JAX external callbacks](https://docs.jax.dev/en/latest/external-callbacks.html).

**How much harder?** This is a qualitative engineering estimate, not a measured
schedule. Porting a small, already pure covariance/Fisher kernel is moderate
additional work. Building all internal noise/weighting and interpolation paths
for JAX from day one is a larger but manageable design investment. Requiring
end-to-end autodiff through arbitrary external packages is substantially harder
and sometimes unavailable without rewriting those packages or supplying custom
derivatives. Supporting two complete production backends also adds a continuing
test and maintenance burden.

The user confirmed NumPy with dedicated Numba-friendly kernels as the starting
backend on 2026-09-11. Implement and validate the reference numerical path with
NumPy, then compile measured bottlenecks in separately reviewed work. Optional
derivative adapters remain compatible with this choice. Preserve the possibility
of JAX with pure array boundaries, but a JAX backend is not an initial delivery
requirement. A thin `xp` alias alone would not address compilation, interpolation,
control flow, and external model boundaries.

Compare cold startup/compilation, steady-state runtime, memory, and complete
forecast time. JAX timings must synchronize completion and separate transfers
from execution. Do not assume a GPU helps small dense covariance blocks or a
single forecast. [JAX benchmarking guidance](https://docs.jax.dev/en/latest/benchmarking.html).

## 10. Validation and staged implementation

The authoritative step sequence is in the separate
[implementation roadmap](FISHHIGHZ_IMPLEMENTATION_PLAN.md). Detailed instructions
and acceptance tests are written only for the current step. The user reviews
plans, dispatches implementation agents, reviews their results, requests an
independent agent review, and decides whether to revise the same step or advance.

Scientific validation must cover the following across those steps:

- Capture current, reproducible outputs for all seven DESI-2 examples before
  numerical implementation. Save effective inputs, versions, hashes, logs, Git
  state, outputs, and timings in the same environment used for comparisons.
- Check analytic covariance/Fisher limits, mode normalization, pair closure,
  parameter ordering/sharing, priors, marginalization, and degeneracies.
- Test template metadata, units, signed wiggles, redshift amplitudes, interpolation
  accuracy, rescaling-domain coverage, and fixed observed cuts.
- Verify equivalent model predictions and transformed Fishers in all three
  scaling bases; check component ties, zero-width damping, cross-width
  normalization, and growth derivatives with fixed damping and template amplitude.
- Compare survey/noise intermediates, density normalizations, responses, and the
  default P1D against the reference before assessing complete forecasts.
- Demonstrate an external model through the same pipeline and establish
  convergence with model-grid resolution and parameter-derivative steps.
- Explain differences from lyaforecast's polynomial-extraction BAO recipe rather
  than forcing agreement between different models or estimators.
- Profile representative workloads and validate numerical agreement after any
  targeted compilation or optimization.

Step 02 has completed reference capture and tool validation as recorded in its
review reports. The numerical checks above remain requirements for the later
implementation steps, not claims of completed FishHighz scientific tests. No
scientific source, example configurations, or reference outputs are changed by
these planning updates.

## 11. Change record

- 2026-09-10: initial source-grounded design; user confirmed arbitrary external
  models first and optional automatic differentiation.
- 2026-09-10: user confirmed independent redshift bins and common volume,
  Fourier-space only with no current or future correlation-function support,
  and cosmological Fisher forecasts without a sampling likelihood. Updated the
  architecture and removed the conditional projection/likelihood branches.
- 2026-09-10: user specified Vega-format decomposed file input, fixed per-tracer
  Gaussian BAO widths, separate component rescaling in three parameterizations,
  forest beta versus galaxy growth-rate Kaiser parameters, and fixed user-chosen
  k cuts. Replaced the generic built-in-model proposal and expanded kernel,
  caching, and validation steps. User confirmed mean-squared widths for cross
  spectra, separately parallel and transverse to the line of sight.
- 2026-09-10: user selected separate implementation agents and a user-controlled
  review cycle. Created the high-level roadmap, package directory, single-step
  assignment for the scaffold, and standing AGENTS.md; moved sequencing out of
  this design to prevent competing plans.
- 2026-09-11: user confirmed NumPy with dedicated Numba-friendly kernels as the
  initial backend. Added a source-grounded comparison of 2D spectra and
  multipoles, recommending the 2D observable pending the user's decision; recorded
  lyaforecast's actual grid and the need for angular convergence checks.
- 2026-09-11: user confirmed the 2D P(k,mu) approach. Synchronized the roadmap
  and agent guidance and checked Step 01 readiness. The scaffold remains the only
  detailed implementation assignment; no package code has been implemented.
- 2026-09-11: reviewed the implemented Step 01 revision-3 scaffold and its report.
  User requested a direct fishhighz/ layout. Step 01 revision 4 requests that
  migration and fresh-environment/artifact validation on every repeat run;
  original artifact installs and source checks passed independent review.
- 2026-09-11: Step 01 revision 4 passed independent review and the user requested
  Step 02. Proposed a fresh seven-case lyaforecast reference capture, with source,
  input, environment, and result provenance plus a separate 15x2pt repeat. This
  changes implementation status only; scientific design decisions are unchanged.
- 2026-09-11: reviewed Step 02 reference capture and proposed revision 2 for the
  user's speed feedback and four evidence-validation gaps. Routine checks stay
  synthetic; a quick full-resolution 15x2pt run compares with saved full evidence,
  and complete captures are reserved for relevant review checkpoints. Scientific
  model and grid choices are unchanged; no parallel reference execution is added.
- 2026-09-11: user directed that agents run quick tests by default and execute
  full mode only on explicit request. Step 02 revision 3 changes the planned CLI
  default to quick and removes automatic full runs from completion/acceptance
  checks. This policy applies across future steps and supersedes the previous
  checkpoint cadence; existing full evidence remains available for comparison.
- 2026-09-11: user requested Step 03 after Step 02 passed review and was committed.
  Prepared the next single-step assignment for fields/pairs, parameter bindings,
  integration grids, and callable array contracts. The user confirmed
  Gauss–Legendre quadrature with custom-grid support and fixed fiducial covariance
  and survey weights. Step 03 uses synthetic quick checks and a fresh installed
  wheel; no real reference forecast is required.
- 2026-09-11: Step 03 revision 1 implemented and independently reviewed with no
  required changes. All 146 tests and installed-wheel checks passed, with
  independent pair/Jacobian/quadrature checks. Core contracts now exist; covariance,
  Fisher, and scientific model implementation remain later steps.
- 2026-09-11: user requested Step 04, accepting progression after the Step 03
  review. Prepared the Gaussian covariance assignment using supplied observed
  field powers, existing pair lookups, and fixed mode counts. It separates
  physical input validation from the numerical kernel and tests valid singular
  cases without introducing inversion or regularization. Implementation remains
  pending user dispatch; quick synthetic and fresh-wheel checks are sufficient.
- 2026-09-11: Step 04 revision 1 implemented and independently reviewed with no
  required changes. All 186 quick tests, 63 independent quadratic-form cases,
  and source/installed-wheel checks passed. Covariance from supplied observed
  powers now exists; Fisher implementation remains a later step. Await the user's
  next-step instruction before preparing Step 05.
- 2026-09-11: user requested Step 05 after the Step 04 review passed. Prepared
  the assignment for Fisher accumulation, reusable factors, named information,
  independent-bin combination, explicit priors, and uncertainty/degeneracy
  reporting. Synthetic quick checks and one fresh installed-wheel probe are
  required; no real reference forecasts. Implementation awaits user dispatch.
- 2026-09-11: Step 05 revision 1 implemented and independently reviewed with no
  required changes. All 261 quick tests, 36 independent assembly/result cases,
  8 singular-diagnostic cases, and source/installed-wheel checks passed. Fisher
  assembly and named results now exist for supplied arrays. Await the user's
  instruction before preparing Step 06 for external models and derivatives.
- 2026-09-11: user requested Step 06 after the Step 05 review passed. Prepared
  the assignment for explicit callable routes, mixed supplied/numerical global
  derivatives, tied bindings, bounds-aware stencils, and opt-in convergence
  studies. Independent P1D invocation and a fixed-covariance synthetic forecast
  are included. Quick synthetic and fresh-wheel checks suffice; implementation
  awaits user dispatch.
- 2026-09-12: reviewed Step 06 revision 1. All 327 quick tests, 9 independent
  mapping/bounds cases, and artifact checks passed, but an independent reproducer
  exposed false convergence when h and h/2 share identical float64 stencils.
  Prepared revision 2 requiring full-study preflight rejection and focused tests.
  No source fixes or real reference forecasts were performed by the reviewer;
  Step 07 remains pending successful repair and user progression.
- 2026-09-12: Step 06 revision 2 passed independent review; the repeated-stencil
  finding is resolved. All 332 quick tests, 16 independent zero-call rejections,
  6 valid one-point refinements, and artifact/installed-wheel checks passed.
  No further repair is required. Await the user's instruction before Step 07.
- 2026-09-12: user requested Step 07 after the Step 06 revision-2 review. Drafted
  template loading/interpolation with explicit units, growth amplitude, signed
  components, derivatives, and domain checks. Reinspected both Vega reference
  files and recorded hashes; negative PKSB samples require signed interpolation.
  Asked the user about optional Astropy/SciPy preparation versus custom NumPy
  spline preparation. The draft awaits that choice; no implementation started.
- 2026-09-12: user confirmed Astropy plus SciPy for template preparation. Resolved
  the dependency gate and synchronized the Step 07 assignment, roadmap, and agent
  guidance. Optional preparation dependencies feed NumPy evaluation kernels;
  the base external-model installation stays NumPy-only. No implementation started.
- 2026-09-12: Step 07 revision 1 passed independent review with no required
  changes. All 408 quick tests, lint/format, 24 independent analytic cases, both
  read-only FITS reference checks, and artifact/installed-wheel checks passed.
  Template loading/interpolation is implemented; the Kaiser/BAO model remains
  future work. Await the user's request before preparing Step 08.
- 2026-09-12: user requested Step 08 after the Step 07 review passed. Prepared
  the intrinsic Kaiser/BAO callable assignment with fixed/free slots, existing
  equality ties, all three component scaling bases, fixed per-field damping,
  strict domains, and numerical derivatives through the existing engine.
  Acceptance includes analytic derivative/basis checks, bounded Vega convention
  comparisons, a synthetic fixed-covariance forecast, and one installed-wheel
  probe. No real forecast recapture is required. Implementation has not started.
- 2026-09-12: Step 08 revision 1 passed independent review with no required
  changes. All 468 quick tests, lint/format, 18 independent five-field model
  cases and 216 derivative columns, and artifact/installed-wheel checks passed.
  The intrinsic Kaiser/BAO model is implemented. No real reference forecasts
  ran. Await the user's request before preparing Step 09.
- 2026-09-12: added a fresh-agent handover with the implemented API map,
  survey-integration contracts, remaining design choices, checkout/environment
  and evidence navigation, and a reusable planning/review prompt. This is a
  documentation handoff only; Step 08 remains review-passed and Step 09 awaits
  the user's progression request.

- 2026-09-13: user requested Step 09 after Step 08 passed review, and confirmed
  optional Astropy background preparation, integrated volume, and explicit
  instrument conversions with labeled legacy compatibility. Prepared revision 1
  with bounded survey primitives/default P1D scope, analytic and convergence
  oracles, optional-dependency/wheel checks, and a stop-for-review handoff.
  Synchronized status and roadmap; no implementation or agent dispatch. Live
  source identity and installed Astropy API availability were checked; historical
  numerical-suite results were not rerun or relabeled as new evidence.

- 2026-09-13: Step 09 revision 1 passed independent review with no required
  changes; the user supplied no comments. All 555 quick tests, lint/format,
  48 independent volume cases, 272 high-precision P1D points, signed-response
  covariance/Fisher checks, source hashes, and the installed-wheel probe passed.
  Geometry/response/default P1D and optional Astropy preparation are implemented.
  No real reference forecast ran. Preserve revision 1 and wait for the user's
  acceptance/progression instruction before planning Step 10.

- 2026-09-13: user requested Step 10 following the successful Step 09 review,
  and confirmed prepared-array inputs, explicit legacy cumulative iterations or
  supplied weights, and declared independence or full supplied noise. Prepared
  revision 1 with fixed per-field auto weighting, explicit density/response/unit
  ownership, independent noise PSD checks, scalar/matrix acceptance oracles,
  bounded DESI-2 input comparisons, and installed-wheel checks. No implementation,
  dispatch, forecast execution, commit, or push; await user plan review.

- 2026-09-13: reviewed Step 10 revision 1 after the user reported no comments.
  The 628-test quick suite, 256 independent Decimal cases, field-matrix Fisher,
  81 bounded reference comparisons, and the existing installed-wheel probe passed.
  A separate mixed-underflow regression failed on source and wheel: supplied
  weights can lose an I3 term and violate normalization invariance. Recorded P2
  finding R1 and prepared revision 2 repair instructions; no production changes
  or Step 11 work. User retains repair dispatch, acceptance, and progression.

- 2026-09-13: Step 10 revision 2 passed independent review with R1 resolved.
  All 638 quick tests, original R1 source/wheel replays, 12 independent boundary
  rejection cases, six zero controls, one exact subnormal control, 96 Decimal
  weighting comparisons, and the existing installed-wheel probe passed. All
  34 reference hashes match; the 81 earlier bounded comparisons remain historical.
  No further repair or Step 11 planning; await user acceptance/progression.

- 2026-09-13: user requested Step 11 after Step 10 revision 2 passed review,
  accepting it for progression. User confirmed Python API plus raw readers with
  INI/CLI/serialization deferred, and optional explicit legacy-compatible SciPy
  density/SNR preparation with strict domain/value errors. Prepared revision 1
  for fixed multi-bin orchestration, normalization/interpolation provenance,
  batching, independent tests and bounded reader comparisons. Planning verified
  all 28 reviewed production modules against the Step 10 wheel; no implementation,
  test-suite/forecast execution, dispatch, commit or push. Await user plan review.

- 2026-09-13: Step 11 revision 1 reviewed; changes required. User requested
  nonuniform redshift bins and confirmed explicit raw redshift-cell widths plus
  labeled legacy first-spacing conversion. Unequal forecast bins already pass
  independent analytic checks; the raw reader blocks three LBG/LAE comparisons
  under the contradictory revision-1 uniform-grid plan. R2 reproduces 20 dtype
  coercion bypasses in ForestInput. New review: 723 quick tests, 32-module
  existing-wheel identity/probe, five examples and unequal-bin checks passed;
  34 reference hashes match, 17 reader comparisons pass and three are blocked.
  Prepared revision-2 instructions only; no production changes, real full run,
  dispatch, acceptance, progression, commit or push.

- 2026-09-13: Step 11 revision 2 passed independent review, closing R1/R2 with
  no further required changes. All 833 quick tests, existing 32-module installed
  wheel probe, five examples, independent nonuniform density/count measures,
  legacy-width oracles and original 20 dtype-bypass rejection cases passed.
  All 34 reference hashes match; 33 new bounded comparisons pass with zero
  blocked, including both known negative-spline rejections. No real full run,
  production fixes or Step 12 planning. Await user acceptance/progression.

- 2026-09-13: user requested Step 12, accepting Step 11 revision 2 for
  progression. Confirmed wiggle-only BAO ap/at, opt-in legacy input floors/clamps,
  explicitly labeled negative-density floor extension and a real lyaforecast
  intrinsic-P3D adapter. Prepared revision-1 instructions with seven synthetic
  selections, bounded one-bin real checks and separate explicit full execution.
  Planning verified 32 modules against the reviewed wheel and offline-checked
  both legacy bundles; preserved schema-v1 limitations. No implementation,
  numerical-suite/forecast execution, agent dispatch, commit or push.

- 2026-09-13: reviewed Step 12 revision 1. All 1017 quick tests, existing
  40-module installed probes, 104 new raw-reference comparisons and independent
  reconstruction of the saved one-bin F pass. R1: checker accepts five invalid
  scientific payloads with matching hashes. Prepared revision 2 to repair R1
  and deliver the user's requested seven-case full reference/two-profile
  comparison, plots, sensitivities and attribution. User confirmed validation-only
  literal legacy algebra, converged existing-model accuracy and retained-input
  sensitivity policies. No production fixes, new full forecast, agent dispatch,
  acceptance, progression, commit or push during this review.

- 2026-09-14: Step 12 revision 2 independently reviewed. All 1059 quick tests,
  78 saved primary C/F reconstructions, 39 primary/trial bindings, existing
  installed probes and 720 plotted series pass. The five original R1 probes
  reject. R2 demonstrates accepted false convergence from unbound metric/trial
  operands. R3 retains the correctly reported 34 unconverged forest accuracy
  bins and 46 unavailable diagnostics; compatibility passes all 39 bins.
  Revision 3 specifies bounded evidence repair and saved-array diagnosis,
  reserving weighting/input choices and any further full run for the user.
  Production source and historical evidence remain unchanged. No new full
  forecast, installation, dispatch, acceptance, progression, commit or push.

- 2026-09-14: user restricted Step 12 revision 3 real validation to
  `lya_qso_lbg_lae_15x2pt`, using all six bins and its 15 individual-spectrum
  forecasts to diagnose and resolve convergence. The other six cases are
  excluded. Reuse the saved actual reference; bounded accuracy/convergence
  work is authorized when the user dispatches implementation. Weighting/input
  choices remain subject to approval. The revision number stays 3. Planning
  inspected saved timing records only and executed no forecast or benchmark.

- 2026-09-14: user archived Step 12 revision 3 for later scientific repair and
  requested performance optimization first. The exact narrowed r3 instructions
  are preserved in `lib/fishhighz/reviews/step-12-instructions-r3.md`. Revision 4
  proposes matrix/loop and repeated-work optimizations, subsequent profiling and
  a detailed runtime prediction report. Short synthetic assessment completed;
  no actual 15x2pt, production changes or acceptance. The earlier real-run
  authorization is suspended with r3; R2/R3 remain open. Step 13 stays deferred.

- 2026-09-14: Step 12 revision 4 independently reviewed against its exact
  performance-only assignment. No required code corrections found; optional
  compilation meets both targets, while the default NumPy matrix target remains
  unmet. New review passes 1163 ordinary tests (25 compiler skips), 358 installed
  regressions, 129 performance regressions outside the checkout, independent
  C/F/range/layout checks, five baseline-matched studies and installed examples.
  Source/wheel identity and short benchmark timings agree. No real forecast,
  scientific acceptance or progression; user acceptance pending. R2/R3 remain
  archived and no repair revision or next step was drafted.

- 2026-09-14: user requested Step 12 revision 5 from the archived revision-3
  instructions. Restored evidence binding repair and 15x2pt-only convergence,
  sensitivity/attribution and handoff requirements, retaining all six bins,
  15 individual spectra and the joint forecast. Preserve r4 optimizations and
  all 1188 current tests. Scientific alternatives still require user approval;
  bounded runs become active on plan approval/dispatch. No new forecast, code
  edit, test-suite execution, dispatch or scientific acceptance during planning.

- 2026-09-14: Step 12 revision 5 independently reviewed: 1183 quick tests
  (25 optional skips), 191 installed regressions, all 12 saved C/F reconstructions,
  180 individual-spectrum checks, 40 coupled-weight trial comparisons and 240
  plotted series pass. A new synthetic weak-spectrum mutation nevertheless
  falsely certifies a 1.9419% error change; R2 remains incomplete. All six actual
  accuracy reports truthfully fail, with 12 diagnostic variants unavailable.
  Fixed-weight quadrature and cumulative-weight decay diagnosis are verified;
  R3 still requires a user scientific decision. Revision 6 proposes bounded
  synthetic/offline evidence repair with per-spectrum binding and trial-derived
  metrics/verdicts. The exact r5 instructions, handoff and historical forecasts
  are preserved. No production fix, real forecast or scientific acceptance
  during review; no Step 13.

- 2026-09-14: Step 12 revision 6 passed independent review of the bounded
  evidence repair. R2 is resolved: the original weak-spectrum false-pass probe
  now rejects through validator, writer and offline reader. New review evidence
  includes 1243 quick tests/25 skips, 251 installed regressions, six examples,
  21 Decimal norm controls, four mixed-scale binding mutations and all 12 joint
  plus 180 individual-spectrum saved-array checks. All source/wheel/installed
  modules match; physical kernels and historical evidence are unchanged. R3
  retains six unconverged accuracy bins and 12 unavailable diagnostics; the user
  retains scientific decisions and acceptance. No new real run, code fix, next
  revision or Step 13 during review.

- 2026-09-14: the user selected option 2 for R3 and requested a new Step 13.
  Revision 1 proposes investigation of the existing rule's normalized limit,
  with exact finite-step equivalence, discrete/continuum asymptotics, diagnostic
  prototypes and bounded saved-input checks. This authorizes progression to the
  investigation, not acceptance of the failed accuracy forecasts. The previous
  r6 instructions are archived exactly. Fifteen scalar and seven continuum
  planning controls pass; no production implementation, real forecast or agent
  dispatch. Later prospective performance/distribution steps shift to 14/15.

- 2026-09-14: Step 13 revision 1 independent review requires S13-R1/R2/R3
  corrections. Finite-step and fixed-grid numerics are independently verified;
  a continuum verdict cannot follow from the current empirical screen alone.
  Swapped populations and omitted checkpoints pass the submitted finalizer;
  a finite-positive scalar range control fails. Revision 2 retains diagnostic-only
  scope. No production change, real run, estimator choice, dispatch or progression.

- 2026-09-16: W08 revision 2 passed independent scientific review. The W08-R1
  positive-product underflow reproducer now rejects with field context, while
  exact-zero and exactly representable subnormal controls remain valid. Focused
  analytic, legacy/supplied, noise, survey and saved-reference checks pass; no
  consequential correction is required. User acceptance remains pending. No
  profile adoption, real forecast, package Step 13 change or progression.

- User-requested progression after W08 review: W09 revision 1 is proposed for
  approval and dispatch. It reconstructs the fixed-reference noise/covariance/
  Fisher calculation on one saved accuracy QSO-forest auto-spectrum/bin and
  compares with W07. Only planning documents changed; no numerical replay,
  production prescription/default change, forecast acceptance or dispatch.

- User-requested progression after W09 review: W10 is proposed for approval and
  dispatch, with coefficient/refinement checks for two additional saved samples
  and conditional historical BAO comparisons. The remaining diagnostic stages
  are now outlined at the user's request; profile adoption and broader runs
  retain separate decision points. No numerical/production work during planning.

- User-requested progression after W10 review: W11 is proposed as a bounded
  reference-mode sensitivity test on the original QSO sample. Two alternatives
  span the inherited P1D floor; early stopping follows a 0.1% error-change
  screen. This does not choose a reference convention or adopt a profile.
  W10's supported two-sample result and weak-LBG limitations remain recorded.

- User decision after W11 review: retain q_star=0.00035 s/km and plan accuracy-
  profile adoption of fixed inverse-variance weights. W12 defines necessary
  host-side method/refinement-contract changes and one saved three-spectrum
  check. W11 sensitivity is recorded, not reclassified as quadrature error.
  Compatibility, physical input policies and historical verdicts are preserved;
  implementation approval/dispatch and full-run execution remain separate.
