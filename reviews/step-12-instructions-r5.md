# Step 12: 15x2pt convergence diagnosis and evidence repair

Revision: 5 (2026-09-14; user restored the archived revision-3 assignment).
State: proposed for user plan approval and dispatch. Revision-2 review requires convergence-evidence repair (R2); the forest
accuracy comparison remains unconverged (R3). User controls approval, scientific
decisions and implementation-agent dispatch. Step 12 is not accepted.
Implementation root: lib/fishhighz.

Read [AGENTS.md](AGENTS.md), [design section 0](../../FISHHIGHZ_DESIGN.md), the
[roadmap handover](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md), the
[revision-2 implementation](reviews/step-12.md), and
[independent review](reviews/step-12-review-r2.md), plus the
[revision-4 handoff](reviews/step-12-performance-r4.md) and
[revision-4 review](reviews/step-12-review-r4.md). Implement this step only when
dispatched by the user. Complete the checks and comparison below, write the
handoff, and stop for review. Do not start Step 13, dispatch agents, commit or push.

## Restoration and revision-4 baseline

The user requested resumption of the exact archived revision-3 scientific
assignment as revision 5. The archive `reviews/step-12-instructions-r3.md`
remains byte-for-byte unchanged (SHA256
`f74db52795b8b27cb3fd9839406e60236a4edde4ea6223df45df03d86f9bc752`).
Its scope, scientific conventions, one-case inventory, convergence thresholds,
R2/R3 findings and user scientific-decision requirements are restored below.
Revision 4's performance-only execution restriction does not apply to this
restored assignment once the user approves and dispatches it. This planning
request itself executes no forecast and dispatches no agent.

Retain the reviewed revision-4 implementation: batched matrix validation and
factorization, optional compiled contraction, prepared-factor reuse, independent
direct NumPy checking and bounded study-payload caching. Do not restore the old
production source together with these instructions. Evidence repair must bind
cached results to the correct actual controls/trials without reintroducing
repeated expensive evaluations or hiding failures. Preserve scalar error/range
semantics and the numerical equivalence of NumPy and compiled paths. No new
performance optimization or compiler redesign is requested; the disclosed
NumPy-only matrix target remains unmet and is not waived by this restoration.

The current baseline contains 1188 ordinary parametrized cases: revision-4
review passed 1163 with 25 optional-compiler skips, plus 358 installed compiled
regressions and 129 performance tests outside the checkout. Those are historical
checks, not new revision-5 passes. Retain all original tests and the 129 added
performance cases; rerun affected tests with both backends where available.
Use the existing optional compiled path for the bounded accuracy study when
available, recording actual backend, versions and cold/warm timing separately;
confirm a short NumPy/compiled numerical comparison before the real study.
Do not change the NumPy default, silently install into a sibling environment,
or rerun a real case solely to benchmark the other backend. Keep the existing
separate direct-solve oracle for independently checking reused factors/Jacobians.

Preserve `reviews/step-12-performance-r4.md`, `reviews/step-12-review-r4.md`,
`.validation/step12-r4-20260914T181308/` and `.validation/step12-review-r4/`.
The previous active revision-4 document is snapshotted in
`.validation/step12-plan-r5/IMPLEMENTATION_STEP.md`. Whenever an archive target
below already exists, verify its identity and retain it; never overwrite it.

## Revision 5 repair scope and remaining scientific decision

This revision runs only `lya_qso_lbg_lae_15x2pt`, retaining all six original
redshift bins and all 15 individually selected spectra as well as their joint
forecast. Do not run any of the other six DESI-2 cases, including as a fallback,
reference recapture, sensitivity study or final acceptance check. Ordinary
self-contained unit regressions remain required; this restriction concerns the
seven real survey benchmark cases, not generic synthetic array tests.

On user dispatch, the 15x2pt accuracy calculation and its bounded convergence,
sensitivity and attribution studies are authorized for diagnosis and validation
of an approved repair. Reuse the saved 15x2pt lyaforecast reference and verify
the compatibility calculation from its upstream arrays; no fresh NewForecast
capture is required or authorized. Necessary model/background setup for this
one case is in scope. This planning update does not dispatch or execute a run.

Sections 1–9 retain the scientific contracts and acceptance thresholds with
the one-case inventory below. The earlier seven-case requirement is historical.
The exact reviewed revision-2 assignment is preserved in
`.validation/step12-review-r2/IMPLEMENTATION_STEP-r2.md`. Its implementation,
full comparison and plots are under `.validation/step12-r2-implementation/`.
The seven-case run was already attempted in revision 2. Preserve all of it;
new real numerical work and new comparison plots are restricted to 15x2pt.

The original five R1 payload probes now reject. Historical revision-2 review passed 1059 quick
tests, all 78 independent primary C/F reconstructions, all 39 primary/trial
bindings, existing installed probes and the 720-series plot check. Nevertheless:

- **R2, convergence evidence:** the writer/checker accepts a copied failed
  accuracy bin whose metric operands have been replaced by repeated copies of
  the higher-resolution result. Recomputed metrics are zero while the saved
  original trial matrices still imply a 20.07% combined-error change. See
  `.validation/step12-review-r2/convergence-probe.json` and `independent.py`.
- **R3, scientific acceptance:** only 5/39 accuracy bins converge. Weighting
  fails in all 34 forest bins; magnitude integration fails in 23. Of 468
  diagnostic requests, 46 are unavailable through guarded underflow. These
  results are correctly disclosed in the handoff, and the full gate fails.

For R2, extend only the validation evidence layer and relevant tests/scripts:

1. Require an explicit inventory of study trials, actual controls, successful or
   failed outcomes, and metric-to-trial references for real accuracy records.
   Bind the primary result to its final trial and the combined lower result to
   its declared lower controls. Bind combined F, selected-pair F and V together.
   Validate shapes, unique identities, ordering/coverage and exact control
   associations before calculating any convergence metric. Predeclare required
   control families and allowed bounded levels; record attempted extensions and
   failures rather than relying solely on an unbound failure list.
2. Derive each metric from its referenced trials, requiring distinct actual
   refinement settings and consistency with final other controls. Preserve the
   existing finite-difference representability checks. Do not reject identical
   numerical matrices at distinct valid controls: constant/analytic limits are
   legitimate convergence controls. Missing operands, failed required trials,
   repeated trial IDs, or references to unrelated trials cannot establish a pass.
3. Reproduce R2 using self-contained tiny real-shaped evidence and a copied
   failed 15x2pt accuracy record with the same duplicate-operand mutation as the
   preserved review fixture. Both writer and offline checker must reject after
   hashes/inventories are consistently updated. Mutate wrong/repeated trial
   selection, missing/reordered trials, control identities, primary/final F,
   pair F, volume, combined lower controls and suppressed failed-trial records.
   Keep finite, null-direction, roundoff and identical-result/distinct-control
   successes. Retain all original five R1 regressions and all 1188 current tests, including revision-4 performance regressions.
4. Preserve schema-2 historical bundles byte-for-byte. If the required contract
   changes, version it and make legacy limitations explicit. Derive any new
   checked evidence into an exclusive directory from existing saved numerical
   inputs, recording producer and checker identities. All 78 historical primary
   C/F arrays remain immutable. Revalidate only the six original 15x2pt accuracy
   trial bindings and corresponding compatibility records; their historical
   failures remain failures. New primary coverage is 12 case/bin/profile records
   and 180 individual-spectrum/bin/profile results. The retained 12 diagnostic
   variants yield 72 requests in this case, including explicit zero-mask proofs
   and failures. Preserve the original seven-case inventory as historical;
   do not require its other six cases for this revision's acceptance.

For R3, do not choose a new scientific prescription without user approval.
Start with the saved 15x2pt density/variance/quadrature/weight samples and small
self-contained examples to diagnose the cumulative update's fixed points,
quadrature dependence and first unrepresentable operations. Separate finite
iteration changes from magnitude integration and demonstrate any explanation
against its recorded trial outputs. The 15x2pt-only refinement study may use
the bounded levels in section 5; do not run the other six cases. Isolate magnitude
quadrature first with fixed weights, then examine the coupled calculation.
Report the
diagnosis and consequential alternatives to the user before changing weighting,
sample support, density/SNR floors, arithmetic guards or acceptance criteria.
There is no authorization for a replacement optimizer or revised physics here.
The original converged-accuracy requirement remains unmet pending that decision.

Preserve the current handoff as `reviews/step-12-r2.md` before an implementation
agent replaces `reviews/step-12.md`. Preserve all prior assignment snapshots,
reports, figures and failed evidence; snapshot current dirty/untracked files
before edits. Keep production scientific modules, strict readers, adapters,
examples, reference assets and environments unchanged for this repair.

After repairs run the existing quick runner and `git diff --check`, focused
semantic mutations, the offline primary/trial binding check, and affected plot
checks. Build/install one exact new wheel in a fresh isolated environment and
run changed evidence regressions and retained examples outside the checkout
with `-I`, recording module/metadata identity. The revised real validation and
its plots cover 15x2pt only. Report R2 closure separately from R3; retain links
to historical figures without rerunning or regenerating the other six cases.
Stop for independent/user review. Scientific changes, dispatch and Step 13
remain with the user; do not re-ask permission for the bounded 15x2pt runs
already included in this assignment once dispatched.

## Scope, authority and preservation

The user originally requested all seven real benchmark cases for revision 2.
Revision 5 is explicitly restricted to 15x2pt, using its 15 individual-spectrum
forecasts to diagnose the problem. Both profile definitions remain unchanged:

- **Lyaforecast maximum compatibility**: a validation-only path reproduces the
  legacy model, BAO derivative, mode counts and pair-specific inputs; FishHighz
  independently assembles covariance and Fisher information.
- **FishHighz maximum accuracy**: a converged configuration within the existing
  model, with integrated volumes, refined k/mu/magnitude quadrature, full wiggle
  derivatives, physical FWHM-to-sigma resolution, CAMB growth and per-field noise.
  Survey samples and per-bin wiggle-only ap/at targets stay fixed.
- Retain approved legacy floors/clamps and first-spacing raw density conversion
  where calibrated coverage or physical cell widths are unavailable, with
  sensitivity tests. Do not invent replacement cell boundaries or input data.

Historical revision 2 explicitly authorized a fresh full seven-case lyaforecast
reference capture, both full FishHighz profiles, and bounded convergence,
sensitivity and attribution studies. That work has been attempted and preserved.
Revision 5 authorizes only the bounded 15x2pt work specified above. Ordinary quick tests remain
self-contained. Follow package login-node/thread guidance; no Slurm.

The revision-1 implementation passed 1017 quick tests and installed probes in
independent review. The saved real one-bin Fisher was independently reconstructed
within 5.02e-15 relative; 104 newly evaluated raw-reference checks passed.
Historical **R1** concerned invalid dimensions, negative Fisher information,
inconsistent errors, incorrect payload pairs and missing scientific pass status.
Revision 2 rejects those five probes. R2 above identifies the remaining
convergence-to-trial binding defect.
The seven-case comparison was revision-2 scope, not a retroactive r1 failure.

Preserve the dirty/untracked checkout at HEAD
0d69786a06d5d676564a51fad14a7156951c2228. Record fresh source/configuration hashes
and Git status, including untracked files, before edits. Revision-1 assignment
and governance snapshots and failing review probes are in
`.validation/step12-review-r1/`. Preserve `.validation/step12-r1-implementation/`,
all prior baselines and reports. Copy the current revision-2 `reviews/step-12.md`
byte-for-byte to `reviews/step-12-r2.md` before replacing it; the r1 handoff
already exists and must not be overwritten. Use exclusive fresh
output directories; never overwrite or silently resume historical evidence.

Extend the existing small validation/recipe layer and its scripts. Keep the two
implemented input/P3D adapters and six examples. Use existing scientific engines
for the accuracy profile. Validation-specific legacy estimator code belongs in
`fishhighz/validation/` or maintained validation scripts, never default kernels,
strict readers or KaiserModel. Reusing read-only legacy upstream P3D/J inputs is
allowed; reusing its final Fisher as a FishHighz result is not.

Excluded: the other six real DESI-2 benchmark cases, new physical models,
overlap-derived noise, cross-bin covariance,
covariance derivatives, nuisance/full-shape forecasts, general INI translation,
production CLI/serialization, new optimization algorithms, further compilation
work and Step 13. Retain the already reviewed optional compiled backend. Preserve the existing recipe library; execute only its 15x2pt recipe and
associated validation JSON/NPZ/plot work. NumPy is the sole unconditional
dependency. Optional SciPy/Astropy/CAMB,
lyaforecast and plotting dependencies remain lazy and outside base imports.
Preserve licenses/provenance; do not distribute sibling raw assets or modify
sibling source, shared environments, authoritative INIs or FITS templates.

## 1. Explicit legacy input compatibility

Keep new compatibility entry points distinct from DensityReader/SNRReader. They
consume the same validated raw tables, metadata and normalization policies, with
owned snapshots. Bad shapes, non-real/nonfinite raw inputs, negative raw counts,
malformed SNR headers, inconsistent grids and invalid exposures remain errors.
Do not turn parser or arithmetic failures into physical fallback values.

Density compatibility must reproduce the actual reference in
`lyaforecast/tracer.py:Tracer.get_dn_dzdm` and the setup methods:

- Preserve explicit raw magnitude masks, target normalization/strict threshold,
  quadratic kx=ky=2 spline and Step 11 width policy. Select legacy_first_spacing
  explicitly for irregular LBG/LAE reference tables; do not infer physical widths.
- A queried magnitude strictly outside the raw table range returns 1e-20 in
  deg^-2 redshift^-1 mag^-1. Endpoints are inside. This is applied after raw
  normalization/interpolation and is not a new normalization of returned rows.
- Reproduce the reference spline's behavior for redshift queries outside its
  tabulated range, labeled as a legacy domain extension. Verify it directly
  against the installed reference spline, including both redshift boundaries;
  do not claim it is a physically validated extrapolation. Strict queries still
  fail. Redshift and query coordinates must be finite and physically valid.
- Apply the separately confirmed negative-density policy confirmed above, retaining raw
  interpolated values and identifying any extension beyond reference behavior.
- Keep exact zeros and positive in-domain densities unchanged unless an explicit
  named policy requires otherwise. No implicit max(density, floor) everywhere.

SNR compatibility follows `Spectrograph.get_pixel_rms_noise`, individually for
each source population. Keep header-based setup and wavelength-only legacy
smoothing from Step 11. For requested pixel width Delta_lambda and exposure
count N_exp, with file N_file:

    if m > m_max or z_source outside [z_min,z_max]
       or lambda_obs outside [lambda_min,lambda_max]: variance = 1e20
    otherwise:
        m_query = max(m, m_min)
        S = SNR_per_Angstrom(m_query,z_source,lambda_obs)
            * sqrt(Delta_lambda) * sqrt(N_exp/N_file)
        variance = 1 / max(S, 1e-10)^2

The 1e10 reference sentinel is RMS; return variance 1e20. Apply the SNR floor
AFTER pixel/exposure scaling; out-of-domain sentinel branches are independent
of that scaling. Preserve closed endpoints and reject incompatible explicitly
supplied exposure time. Validate representability; do not mask NaN/overflow
using fmax or a sentinel. No geometric-mean SNR between the two forests, no
pair-dependent response, and no fallback P1D; those would violate accepted
per-field ownership in this adapter. Literal legacy pair-specific behavior belongs
only in the separate validation profile described below. Generic unrestricted SNR extrapolation is outside scope.

Each sampled result must carry the original and effective query coordinates,
policy names/constants, reason masks/counts for density floor, negative handling,
redshift extension, SNR bright clamp, out-of-range sentinel and SNR floor, plus
source hashes/units/normalization metadata. Persist those diagnostics into the
plain provenance passed to ForestInput and validation outputs. Do not mutate
reader counters during derivative runs. A forecast using compatibility must
identify it prominently, with no claim of unchanged strict-domain physics.

## 2. Complete case and result inventory

Verify `../lyaforecast/examples/desi2/lya_qso_lbg_lae_15x2pt.ini` read-only against
its explicit recipe. Preserve each original field index, exact selected pair order,
covariance closure, area, exposure, density normalization, magnitude domain,
redshift edges, observed k/mu cuts and reconstruction inputs. A number in a case
name is insufficient to identify its selection.

| Case ID | Selected pairs | Original bins |
| --- | ---: | ---: |
| lya_qso_lbg_lae_15x2pt | 15 | 6 |

Required primary coverage is six original bins per profile: six preserved
actual legacy reference results and 12 FishHighz results across both profiles.
Each profile also contains 15 individually selected-spectrum forecasts per bin
(90 per profile, 180 total). The fields are Ly-alpha from QSO spectra, QSO,
LBG, LAE and Ly-alpha from LBG spectra. Preserve all 15 auto/cross selections,
including the forest-forest cross, and the original bins over z=2.0–3.41.
Refinements and sensitivity runs are additional labeled evidence. Reuse the
reference capture; do not invoke the seven-case capture protocol or its repeat.
Never infer completion from a convenient intersection of available results.
Original recipes happen to have uniform edges; preserve Step 11's arbitrary,
nonuniform BinSpec support and its regression tests.

For every case/bin/profile save combined-selected-pair Fisher, parameter
covariance, marginalized sigma_ap/sigma_at and correlation coefficient, and the
corresponding forecasts using each individually selected pair. Compute the
combined forecast with the full selected covariance, not by summing single-pair
information. Unselected zero-filled legacy output placeholders are not forecasts
and must not enter ratios or plots. Use named per-bin ap/at blocks in the same
order; do not invert unrelated zero columns or tie distance parameters across
redshift. No priors or nuisance marginalization in the primary comparison.

## 3. Preserved 15x2pt reference and maximum-compatibility profile

Use `.validation/step12-r2-implementation/reference-01/lya_qso_lbg_lae_15x2pt/`
and its parent manifest as historical actual `NewForecast` evidence. Check this
case's configuration, source/resource identities, all six bin arrays and results
before using them. Retain original/effective INIs, logs, versions, source/import/
resource hashes, output dictionaries and timings. Do not label these results
newly executed or mutate the old environment. If the reference is incompatible,
report the specific mismatch; do not recapture the seven-case suite or silently
substitute a different reference.

The compatibility profile must reproduce the actual legacy observable and
numerical estimator, not merely set the production model's inputs similarly.
Inspect live `forecast_new.py`, `fisher.py`, `covariance.py`, `weights.py`,
`power_spectrum.py`, `tracer.py` and `spectrograph.py`. Record and match:

- Exact original k/mu nodes, endpoint handling and weights; bin-centre volume
  and explicit mode counts. Legacy endpoint sums need not integrate to the
  nominal cut span. Do not force them through a grid validator by silently
  renormalizing weights or changing cuts; the validation path can pass explicit
  node/mode arrays to the existing low-level Fisher routines.
- Actual mean and covariance evaluation redshifts. The current new runner uses
  arithmetic survey centres for cached mean P3D and geometric centres inside
  covariance preparation; reproduce and label that distinction. Preserve actual
  growth, bias, damping, template/power interpolation and all constants.
- Pair-specific weighting/source/noise/response, original magnitude rectangular
  sums including endpoints, three weight updates, legacy resolution convention
  and speed-of-light constants at each reference call site. Do not replace these
  with approved physical per-field defaults while calling the result maximum
  compatibility. Capture exactly the total-power array used by legacy Fisher.
- Legacy degree-8 log-polynomial BAO peak extraction, its weights/regularizing
  constant, backward log-k derivative and zero first node, mu factors, damping
  and pair-specific reconstruction. Reference means already include instrument
  smoothing; do not multiply the Jacobian by response a second time.

Reuse the saved upstream arrays for this case; do not recapture or rerun
NewForecast. A missing required array must be reported before proposing any
additional reference execution. Do not edit the sibling checkout.
Save mean, derivative, required total powers, selected pairs and modes before
Fisher assembly. Independently form each Wick element
`C_(ij)(mn) = (T_im*T_jn + T_in*T_jm)/N_modes` and calculate
`F = sum_nodes J.T @ solve(C, J)` using FishHighz Cholesky/Fisher machinery.
Test this implementation against a separate direct-solve oracle. Use legacy
Fisher/output only as a comparator; it cannot supply the alleged independent F.

Do not weaken default PSD, dtype, range or response-ownership checks. Literal
legacy total powers can be pair-dependent or reflect negative raw spline values.
Record field-matrix eigenvalues and legacy inconsistencies. If such an array is
outside the physical field-PSD contract but its selected Wick covariance is SPD,
a clearly labeled validation-only Wick assembly may feed existing factor/Fisher
routines directly. Explain that condition and why this is algebraic legacy
reproduction, not a physically consistent production model. If selected C is
not SPD or named information is unconstrained, record the failure/null result;
no jitter, absolute values, floors or pseudoinverse finite-error substitution.
The opt-in physical adapter above still rejects negatives unless floor_negative
is explicitly selected. It does not acquire a pass-negatives mode.

For well-conditioned matched arrays require rtol <=5e-12 (explicit absolute
scales only for zeros/cancellation); for final matched finite uncertainties
require relative agreement <=1e-6 with the preserved actual legacy calculation. Report
Fisher norm differences, correlations and condition numbers too. Investigate
any breach, even below the 1% scientific-reporting threshold. Conditioning may
explain a discrepancy but is not permission to silently weaken acceptance.

## 4. Maximum-accuracy profile within the existing model

Use the existing prepare_bin/run_bin/run_forecast and KaiserModel scientific
engines through the optimized revision-4 accuracy recipe; retain factor reuse
and its separate direct NumPy oracle. “Maximum
accuracy” means demonstrated numerical convergence under the following fixed
physical model and supplied inputs, not a claim of all survey realism.

- Use the same bin edges, samples, normalization targets and observed cuts as
  the corresponding reference. Set a single explicit evaluation redshift per
  bin, retaining r1's `sqrt((1+z_lo)*(1+z_hi))-1` convention for all mean/noise
  preparation. Keep arithmetic centres separately as comparison plot labels.
  Integrate `V = Omega * integral D_M(z)^2*c/H(z) dz * h_fid^3` with consistent
  physical Mpc background quantities and explicit h_fid; refine volume order.
  Within-bin signal/density evolution remains the declared z_eval approximation.
- Use Gauss-Legendre k/mu quadrature on the fixed observed domain, refining k
  intervals and mu order. Mode counts use integrated V and the accepted measure.
  Do not alter k_min/k_max to obtain convergence or mask source-domain failures.
- Use the explicit Vega K/PK/PKSB template (the reviewed local Planck18_z_2.406
  template is available). PKSB stays smooth, PK-PKSB stays signed wiggle. Record
  its cosmology, redshift, h conversion and hash. Replace the labeled EdS power
  growth with CAMB linear growth: power multiplier
  `G = [sigma8(z_eval)/sigma8(z_template)]^2`, with f at z_eval. Prepare CAMB
  outputs at the actual redshifts, including z_template and damping reference
  z_ref, or independently demonstrate interpolation convergence. Report any
  template/background mismatch; do not regenerate PKSB or retune it for agreement.
- Wiggle-only ap/at at fiducial (1,1), with independent named parameters per bin;
  identity smooth scaling, fixed biases/beta/f/G and auto widths. Differentiate
  all existing wiggle-mapping terms: inverse coordinates, `Q=1/(ap*at^2)`,
  remapped RSD and damping at fixed widths. No covariance, noise, geometry,
  response or weight derivatives. Keep full fiducial signal in covariance.
- Damping uses explicit accepted auto widths from the reference sigma8 ratio:
  Sigma_perp = 3.26*sigma8(z_eval)/sigma8(z_ref)/sqrt(r),
  Sigma_parallel = (1+f)*Sigma_perp; r=1 for forests and the original configured
  reconstruction factor for galaxies. Cross squared widths are the mean of
  squared auto widths. This differs from the legacy pair reconstruction rule.
- Interpret ordinary resolving power R as FWHM: Gaussian sigma_v is
  `c/(R*2*sqrt(2*ln(2)))`, with c=299792.458 km/s and explicit pixel-Angstrom
  conversion. Preserve signed sinc/Gaussian field transfers and apply them once
  to signal/mean derivatives. Supplied pixel/Poisson noise is not W-multiplied.
- Prepare density, SNR and weights separately per observed field, using each
  forest's own auxiliary auto-power and independent default P1D including its
  low-k treatment. Declare independent sampling when using diagonal noise.
  Do not infer overlap noise from shared tracer labels or invent cross-noise.
- Refine the fixed-domain magnitude integral with explicit nodes and measures,
  including cumulative weighting integrals. Partition at known interpolation,
  raw-support and clamp/sentinel boundaries so discontinuities are resolved.
  Legacy rectangular endpoint sums remain compatibility-only. No hidden extra
  renormalization to restore a desired density after changing quadrature.
- Assess the existing cumulative weight update at explicit counts 3, 6 and 12
  (up to 24 if needed), recalculated once per field then held fixed. Use its
  stabilized result if demonstrated; do not call three updates converged or
  claim the fixed point is jointly optimal. If the rule does not stabilize,
  retain diagnostic results and report an unresolved limitation for review;
  do not introduce a new optimization algorithm or silently call it converged.

Primary raw inputs retain the explicitly approved out-of-domain/clamp policies,
negative-density extension and legacy_first_spacing where physical cell widths
are unknown. Record raw/effective samples, reason masks and contributions for
every field/bin. Numerical refinement does not make these policies physical.
No calibration, source-overlap model, distortion/mode-mixing physics, redshift
errors, new P1D physics or within-bin evolution model is authorized here.

## 5. Convergence, sensitivity and causal comparison

Run only the 15x2pt profiles over their six original bins. For accuracy, demonstrate separate
k, mu, magnitude, volume, finite-difference and weight-iteration sensitivity;
change one numerical control at a time with all others at final settings.
Cache fixed preparation appropriately; batch nodes and do not reread raw files
or rebuild CAMB in derivative loops. Compare results on the same physical domain.

Initial levels: k intervals 32/64/128 with order 4; mu and volume orders 8/16/32;
ap/at steps 1e-3, 5e-4, 2.5e-4. For magnitude, use at least three levels of
composite quadrature on a declared common breakpoint partition, e.g. orders
4/8/16. These are starting checks, not universal convergence guarantees. Confirm
the selected final combination against a combined refinement as well as isolated
controls. Where necessary allow at most two further factor-two quadrature
refinements per control; report any remaining failures without weakening tests.

For each bin and combined forecast require final-two relative Fisher Frobenius
change <=1e-3 and relative finite marginalized-error change <=1e-3; volume
relative change <=1e-6. Apply the same information/error criteria to stabilization
of the existing weighting rule. Store actual controls, arrays, conditions and
metrics, not only booleans. For selected single-pair curves, check their reported
errors converge to the same 1e-3 budget, or explicitly mark unresolved curves.
If no stable finite-difference window exists, report truncation/roundoff and
conditioning rather than selecting whichever step best matches legacy.

Quantify retained input-policy sensitivity at least in every affected case/bin:
negative/out-of-magnitude density floors at 1e-22, 1e-20 and 1e-18 (diagnostic
variants, central primary remains 1e-20); contribution-removal controls for
bright-clamped and out-of-domain/SNR-sentinel samples; and sensitivity to the
unknown raw redshift-cell conversion. For the last, keep legacy_first_spacing
in the primary and use explicitly artificial +/-10% cell-width perturbations
at fixed normalized cell counts as local sensitivity diagnostics. Do not label
these widths measured or infer physical edges from centres. Recompute derived
noise/weights consistently. If a policy is never exercised, save the zero-mask
proof; no empty sensitivity run is needed. State which tests change effective
sample support and why they are diagnostics, not alternate primary surveys.
No perturbation represents a calibrated systematic uncertainty or error bar.

Compare by stages: raw/effective density and variance; geometry/modes;
response/noise; intrinsic and observed required-pair powers; selected C and
its conditioning; Jacobians; F, sigma_ap/sigma_at and correlation. Save aligned
intermediate comparisons at common nodes where grids differ. Use relative
Frobenius metrics and explicit zero-safe scaled differences rather than dividing
by signed/zero wiggle or cross-power values blindly.

For reporting, flag >=1% changes in either finite uncertainty or Fisher norm,
and >=0.01 absolute changes in ap/at correlation. These are material forecast
changes, not observational statistical significance. Every compatibility
numerical-tolerance breach is a finding regardless of these thresholds.
Document the sign and size of BOTH FishHighz/reference differences and the
accuracy/compatibility difference for every bin and selected-pair diagnostic.

For every case with a flagged primary difference, perform controlled swaps at
least at its largest-discrepancy bin and at any bin showing a different dominant
cause. Cover model/template/peak extraction, full versus legacy BAO derivative,
growth/redshift, volume/modes/quadrature, physical resolution, reconstruction,
per-field versus pair-specific noise, magnitude/weight convergence and raw
policies as relevant. Label the starting profile, changed setting and held-fixed
quantities for each diagnostic. A sequential chain must reproduce its endpoints;
show the residual and explain order/interaction dependence, without claiming
independent additive effects. If a swap is not physically representable, use
clearly labeled supplied-array diagnostics. A list of possible causes alone is
not a comprehensive explanation. Unexplained significant differences remain
open handoff findings, not a passed comparison.

## 6. Plots and comparison artifacts

Generate reproducible scientific plots using ordinary plotting tools (e.g.
optional Matplotlib with a noninteractive backend). Save a maintained plotting
script plus exact machine-readable JSON/CSV data and source bundle hashes.
Plotting must work offline without rerunning forecasts. Do not use raster image
generation for scientific curves.

Required outputs:

1. **One 15x2pt case figure**, showing sigma_ap(z), sigma_at(z) and their
   correlation, with three labeled curves: actual lyaforecast, FishHighz maximum
   compatibility, and FishHighz maximum accuracy. Include fractional-difference
   panels to legacy (absolute difference for correlation). Show all original
   bins, with original centres and bin widths; record actual evaluation redshifts
   separately. Use consistent fraction/percent units and no misleading offsets.
2. Supplemental plots or a multipage PDF for all 15 individually selected
   pair forecasts, each across all six bins, with the same three profiles and
   reference ratios. Exclude
   unselected legacy placeholders; show failed/unconstrained entries explicitly.
3. A 15x2pt summary plot/table of both profiles' uncertainty changes and
   convergence status over all six bins, plus attribution/sensitivity figures.
   Do not regenerate the six excluded cases. Do not hide failures by plotting
   only common successes.

Provide PNG plus PDF or SVG for each main figure, descriptive legends/captions,
and a manifest of filenames/hashes. Verify all plotted values match the checked
numerical tables, all six bins/both profiles/15 selected pairs are present, ratios use the correct
reference, and visually inspect representative/all main figures for labeling,
clipping, zero ranges and overlaid-curve visibility. If optional ellipses are
shown, specify joint confidence level and Delta-chi-squared convention explicitly.

## 7. Real external intrinsic-P3D adapter

Provide a thin callable around a caller-prepared object exposing
compute_p3d_hmpc(z,k,mu,corr). Require explicit field-index-to-reference-label or
pair routing, source h units and a supported k/z domain. Validate paired-node
shape, pair bounds/coverage, finite real outputs and domain before invoking a
provider that clamps internally. Preserve arbitrary requested pair order and
signed crosses. Do not call smoothed/observed P3D entry points or add noise.

The reference get_pk_lin table is returned in its native h units. If different
h_fid is supported, explicitly map k_source=k_fid*h_fid/h_source and
P_fid=P_source*(h_fid/h_source)^3; alternatively reject unequal h with a clear
contract. No silent conversion. The wrapper does not infer cosmology or rebuild
CAMB on parameter perturbations. A fixed external object can use zero local
parameters; demonstrate a separate explicit amplitude wrapper A*P_external to
produce a real Fisher calculation with an analytic amplitude oracle. Do not
claim this amplitude demonstration implements the primary wiggle-only BAO model.

The actual installed/read-only lyaforecast PowerSpectrum demonstration already
ran in revision 2. Preserve that evidence and adapter regressions; a new external
amplitude forecast is not required in this 15x2pt-only assignment. Record
import origins, versions, source hashes and caller setup; no neighboring import
or environment modification on normal FishHighz import. Repeated calls must be
deterministic; P1D remains the separate accepted provider. Do not introspect or
serialize arbitrary object internals as reproducible state.

## 8. R1 repair: scientific evidence validation

Version the validation-only evidence contract. Distinguish real BAO (two named
parameters per bin), synthetic amplitude (one parameter), external amplitude,
compatibility/accuracy profile and sensitivity/attribution records. Each must
have an explicit required payload schema, not arbitrary worker-reported shapes.
Preserve historical schema-1 artifacts and their documented weaker checks;
legacy inspection may return an explicitly limited result but must not satisfy
the new 15x2pt acceptance gate. Do not silently rewrite old bundles.

At production of evidence and in an independent offline checker:

- Bind exact case/profile/bin bounds, original field identities, selected and
  required pairs (including payload integer contents/order/closure), parameter
  IDs/order, nodes/measures and aligned array dimensions to the declared request.
  Require the exact primary inventory and separate diagnostic inventory. Catch
  missing, duplicate, extra and swapped records, including plausible same-shaped
  case/bin/profile swaps. Match original and effective configuration hashes.
- Validate finite real arrays and appropriate integer index arrays. Fisher must
  have the declared square parameter shape, be symmetric within a documented
  numerical tolerance, and have valid information sign/rank. Recompute named
  covariance/errors/correlation using rank-aware semantics and check supplied
  results. Do not allow finite error bars for unconstrained directions; encode
  unavailable results explicitly rather than forcing all JSON numbers finite.
- For real primary records, retain enough total power, modes, selected C or
  factors and observed J to reconstruct C and F independently offline. Make
  response ownership explicit for each profile so verification never applies
  it twice. Validate F against the saved node contraction. Stream chunks if
  needed; hashes alone do not establish mathematical consistency.
- Require `passed` to be the boolean True for a successful record and evaluate
  saved convergence/comparison metrics against recorded predeclared thresholds.
  `complete` as accepted evidence must require successful scientific checks;
  distinguish execution finished from scientific validation passed. Record
  failed/incomplete tasks and logs without passing a partial run.
- Record actual imported FishHighz/lyaforecast/CAMB source paths and versions,
  FishHighz source and exact wheel identities, resolved resource paths/hashes,
  plot/data provenance, policy masks and all settings. Confirm imported reference
  code/resources are the ones inventoried, not merely files beneath a requested
  checkout while Python uses a different installed copy. No pickle/model imports
  when reading evidence; validate paths and inventories before loading arrays.

Reproduce all five original review probes from
`.validation/step12-review-r1/independent.py`; each must now reject. Also mutate
otherwise valid new synthetic and tiny-real-shaped bundles, updating hashes and
reported inventories consistently so tests reach semantic validation: wrong
F/error/correlation, pair IDs, parameter order, node alignment, missing/nonboolean
pass, and failed convergence disguised as complete. Include a correct finite
case, a legitimate null direction, and roundoff-level controls to avoid rejecting
valid evidence. Missing old schema guarantees must be visible to callers.

## 9. Acceptance, packaging and handoff

Preserve all 1188 current tests and their numerical contracts. Add focused
self-contained checks for R1/R2, both profile dispatches, exact one-case/six-bin/
12-primary inventory and 180 individual-spectrum results,
selected-pair output, numerical controls, sensitivity labeling, failure handling,
and plot-data alignment. Ordinary pytest must not import machine-local assets,
solve CAMB or execute real forecasts. Test full orchestration with synthetic
workers. Existing synthetic seven-case inventory tests may remain; they execute
no real survey forecast. Keep quick as the default CLI mode. Explicitly select
`--cases lya_qso_lbg_lae_15x2pt` and all six bins for real validation; never invoke
an unfiltered full command. If the CLI uses `--suite full` to select all bins,
the explicit single-case filter is mandatory. Check the task inventory before
expensive setup and add a synthetic dispatch test proving no other case runs.
The existing seven-case `check_desi2_full.py` gate remains a historical full-suite
contract; add an explicitly scoped 15x2pt gate rather than making it accept a
one-case bundle as a completed seven-case suite.

Retain adapter endpoint/negative/zero/dtype/exposure/domain checks and real
Tracer/Spectrograph branch comparisons. Keep arbitrary external P3D/P1D and
optional-dependency blocker tests. Demonstrate exact single-application W/noise
ownership, fixed covariance/weights during derivatives, bin independence and
nonuniform BinSpec support. Preserve the independent five-field signed-cross
Wick oracle and full-model analytic wiggle derivative oracle:
`dP_w/dap = -P_w-k*mu^2*P_w'`,
`dP_w/dat = -2*P_w-k*(1-mu^2)*P_w'`
for isotropic factors and zero damping. Include a nonzero-width anisotropic
finite-difference check and literal legacy derivative/endpoint fixture tests.

Run `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` and `git diff --check` after
changes. Use OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and MKL_NUM_THREADS=1 for
all numerical runs. Keep the six 15x2pt bins serial or conservatively bounded on the
login node under package guidance. This authorization does not cover Slurm.
If convergence work becomes materially more intensive than these bounded
settings, report the unresolved limit rather than launching other resources.

Build one exact wheel into a fresh directory and install that artifact in a
fresh isolated environment with needed extras. Outside the checkout, with `-I`,
verify all source/wheel/installed modules, METADATA/WHEEL bytes, wheel hash and
import origins. Run six examples and new profile/evidence regressions there.
Run real profiles with this exact implementation in an identified environment;
if reference dependencies require a separate interpreter, verify the same wheel
payload and document that split. Do not mutate historical baseline environments.
Optional plotting does not become a mandatory runtime dependency.

Required revision-3 scientific acceptance evidence: the preserved actual 15x2pt
reference, both FishHighz profiles for all six bins, all 15 individual selected-
pair diagnostics per bin, matched compatibility checks,
accuracy convergence, retained-input sensitivities, causal comparisons and plots.
The retained 12 diagnostic variants have exactly 72 case/bin/variant records,
with inactive policies proved by zero masks and unsuccessful variants explicit.
Do not require fresh evidence from the six excluded cases for this revision.
No reduced bin/pair/grid subset can be relabeled complete, and no 15x2pt-only
result can be described as new seven-case acceptance. A failed bin or
unexplained significant difference remains explicit and prevents an unqualified
pass; preserve its evidence and report it for review without silently changing
scientific assumptions or tolerances.

Update README with actual Python/validation entry points, profile definitions,
explicit 15x2pt-only invocations, plot regeneration commands and limitations. Historical
r1 single-bin evidence stays separately labeled. In `reviews/step-12.md`, record:

- Revision-2 preservation and R2 closure, changed-file manifest, exact commands,
  interpreter/source/resource/wheel identities, logs, timings and failure history.
- A concise 15x2pt comparison table and complete six-bin/15-pair
  numerical tables with both FishHighz/legacy ratios and accuracy/compatibility
  differences; link its main figure and supplementary/attribution plots.
  Clearly separate historical seven-case evidence from new 15x2pt-only results.
- Every per-bin profile setting, physical/evaluation-redshift conventions,
  convergence levels/metrics, floor/clamp/extension counts and sensitivities.
- Significant differences ranked by size, controlled evidence explaining each,
  residuals and limitations. Separate numerical reproducibility from physical
  adequacy, and present missing sample-overlap/coverage/evolution information
  honestly. “Maximum accuracy” is bounded by this existing-model contract.
- Actual completed/failed coverage and remaining findings; no invented pass,
  acceptance, or equivalence claim. Stop for independent/user review. No Step 13.
