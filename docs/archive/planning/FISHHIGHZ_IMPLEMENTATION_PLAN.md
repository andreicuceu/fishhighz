# FishHighz implementation roadmap

Updated: 2026-09-21. Status: **S1–S5 complete with passing independent reviews**.
The five-stage compatibility weighting study and additional grid-to-BAO check
are complete. The user selected **early lyaforecast weights as the baseline**,
with McDonald weights retained as an explicit alternative. S1 implements profile
integration; S2 qualifies the revised accuracy calculation within its stated
finite-refinement tests; S3 completes the individual and joint comparison.
S4 completes the controlled difference attribution. On 2026-09-21 the user
selected the tested S2–S4 accuracy recipe as the research baseline, retaining
its mixed-pair damping and full wiggle AP derivatives. S5 completed its three
implementation parts; non-critical sensitivities are deferred. The user
authorized S5 implementation with Sol (medium) agents
and Sol (high) scientific, software and documentation reviews. S5 is now
complete with a passing independent review and coordinator verification.

This document records the requested direction, recommendations and confirmed choices.
The execution record below covers completed S1–S5. Step 13, W12 and the old W01–W13 diagnostic
plan remain closed. Their historical conclusions and numerical evidence stand.
The preceding long roadmap remains in the
[2026-09-17 archive](FISHHIGHZ_IMPLEMENTATION_PLAN_ARCHIVE_2026-09-17.md).

## Scientific basis and current state

The authoritative weighting choice and equations are in
[FOREST_WEIGHTING_DECISION.md](../../research/FOREST_WEIGHTING_DECISION.md).
The [completed weighting study](FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md)
and its linked reports supply the calculation and independent reviews.
Early weights reached finite, nonzero iteration convergence for both forest
populations in all six bins on all three tested magnitude grids. The final
213-to-425-node refinement changed individual BAO errors by less than 0.1%
and joint errors by at most 0.00562% in **bins 2–6 only**. S2 subsequently added
the selected three-spectrum bin-1 check: the final refinement changes individual
and joint errors by at most 0.00413% and 0.00152%, respectively. These are
finite-refinement results, not continuum bounds.

| Work | Disposition |
| --- | --- |
| Steps 01–11 | Accepted for progression: Fourier covariance/Fisher core, external models, templates/Kaiser/BAO, geometry/response, weighting/noise and Python survey orchestration. |
| Step 12 | Retain validation, reviewed performance work and revision-6 evidence repair. Historical unconverged accuracy forecasts remain scientifically unaccepted. |
| Step 13 | Done and closed by user decision. Preserve the unresolved historical continuum result; no further assignment. |
| W12 | Done and closed. Preserve the fixed-reference inverse-variance implementation and bounded passing review as an explicit alternative and historical comparison. |
| Old W01–W13 plan | Closed; deferred repairs and the former W13 proposal remain retired. |
| Compatibility weighting stages 1–5 | Completed study and accepted baseline choice, as recorded in the weighting decision. S1 integrated the accepted full-sample recurrences and stopping into a shared numerical kernel. |
| New three-profile comparison | S1 integration, S2 numerical qualification and S3 comparison are complete with passing independent scientific reviews. Scientific acceptance remains with the user; S4 attribution is complete with a passing independent scientific review. |
| Package completion | S5 completes the research-ready Python API documentation/example milestone with passing scientific, software and documentation review. Numerical assumptions and finite qualification remain explicit. A broader CLI/configuration and serialization release is deferred. |

## 1. Define the three profiles

Use these names in new configuration, reports and figures. Preserve historical
records carrying `compatibility` or `accuracy` under their original meanings;
identify new results by profile name **and recipe revision**. If `compatibility`
is retained as an API alias, it must mean `full-compatibility` only and the
resolved name must be recorded. Never reinterpret an old cache as a new profile.

| Choice | full-compatibility | fixed-compatibility | accuracy, revised |
| --- | --- | --- | --- |
| Purpose | Reproduce the captured modern lyaforecast choices and implementation | Isolate the accepted weighting change | Numerically qualified forecast within the existing physical model |
| Forest weights | Literal prefix/intrinsic recurrence, exactly three updates | Early lyaforecast full-sum recurrence, confirmed convergence | Same early recurrence and stopping, using this profile's fiducial field inputs |
| Magnitude integration | Original rectangular grid, 107 nodes for this case, full endpoint weights | Identical original grid and measure | Composite Gauss–Legendre integration on the existing support/interpolation partition; demonstrate refinement |
| Density/SNR policies | Literal saved compatibility inputs, including signed interpolants | Identical saved inputs; no new floor or clipping | Retain existing explicitly labeled accuracy adapters and sensitivity tests; strict public defaults stay unchanged |
| Mean, derivatives and response | Captured legacy estimator, redshift conventions and resolution | Identical to full-compatibility | Existing KaiserModel wiggle AP mapping, consistent evaluation redshift, physical FWHM-to-sigma response |
| Geometry, growth and integration | Captured volumes, growth, Fourier nodes and mode counts | Identical to full-compatibility | Integrated volume, CAMB power growth and refined k/mu quadrature |
| Covariance | Independent reconstruction of the literal legacy Wick covariance | Rebuild all affected entries after replacing forest-auto noise | Consistent per-field powers/noise and full inter-spectrum covariance |

“Full” denotes the captured lyaforecast revision and supported benchmark, not a
promise to reproduce every upstream option. Bind that revision and reuse the
immutable reference inputs; a fresh upstream run is unnecessary for this task.
The shared sample definition starts from the six-bin DESI-2 case with fields
lya(qso), lya(lbg), QSO, LBG and LAE, subject to the confirmed selection below.
Preserve the survey area, magnitude limits, bin edges, observed k cuts, template
input, parameter definitions and priors.
Existing profile differences in redshift evaluation and model treatment must be
listed explicitly, rather than accidentally equated when holding inputs fixed.

### Confirmed selection for all future forecasts

On 2026-09-18 the user excluded **every bin-1 correlation involving LBG, LAE
or lya(lbg)** from future forecasts because this redshift bin is too weak for
those samples. Bin 1 means the lowest-redshift bin (mean z = 2.1175; internal
index 0). This is an observable-selection rule, not merely a plotting cut.

- Bin 1: retain only lya(qso) auto, QSO auto and lya(qso) × QSO. Compute the
  joint result using the covariance of these three selected spectra.
- Bins 2–6: retain all 15 spectra and the full 15-spectrum joint covariance.
- Apply the identical selection to all three profiles, all refinements and
  all attribution/sensitivity forecasts. This gives 78 individual forecasts
  and six joint results per profile when every selected result is available.
- Mark excluded entries as excluded by the scientific selection, distinct from
  solver failure or nonconvergence. Do not compute excluded individual forecasts
  or include them in the joint information. Prepare only the forest fields and
  required covariance spectra needed by the selected observables.
- Recompute the bin-1 joint Fisher from its selected covariance/Jacobian; neither
  a sub-block of the old joint Fisher nor the sum of three individual Fishers
  gives this result. Existing selected individual results remain useful controls.

Historical files and their all-15 bin-1 results remain unchanged. The selection
change applies to full-compatibility as well: its numerical conventions reproduce
lyaforecast for the selected observables, while its new bin-1 joint result is not
expected to equal the historical 15-spectrum joint forecast.

The primary fixed-compatibility result must retain 107 nodes: silently adopting
425 nodes would mix weighting and quadrature changes. Use a separately labeled
**fixed-compatibility magnitude-refinement diagnostic** at 107/213/425 nodes,
including only the three retained spectra in bin 1. This is not a fourth
primary profile. Likewise,
retain exactly-three-update early weights only as a diagnostic separating the
change of recurrence from convergence of that recurrence.

For fixed-compatibility, follow the reviewed reassembly: replace the required
forest-auto noise contributions (QSO forest in bin 1; both forests in bins 2–6),
retain intrinsic and cross powers, galaxy
noise, mean Jacobians and modes, and rebuild the full Wick covariance. Extra
historical pair-specific preparations remain diagnostic records; do not assign
independent physical weights to each cross-spectrum or introduce a new overlap
noise prescription. Galaxy-only individual forecasts must be unchanged, though
the full joint result can change through its covariance with forest spectra.

## 2. Integrate the baseline into accuracy

### Weight equation and ownership

For each observed forest field and redshift bin, define full-sample scalar
moments using its prepared magnitude measure, r_i = n_i q_i:

    I1 = sum(r_i w_i)
    I2 = sum(r_i w_i^2)
    I3 = sum(r_i v_i w_i^2)
    w_i(0) = B / (B + lp v_i)
    w_i(t+1) = (P L I1(t) + B) / (P L I1(t) + B + lp v_i)
    A = I2 / (L I1^2)
    P_pixel = lp I3 / (L I1^2)

P and B are fixed representative, response-smoothed intrinsic auto-P3D and P1D;
only the weight moments are updated. Both are evaluated from the accuracy
profile's fiducial model/response, rather than copied from compatibility.
Use the full fiducial auto-P3D provider, consistent with existing
`sample_auxiliary`; retain the independently specified default P1D and its
low-k convention. No integration of P3D to replace P1D is proposed.

Use the existing angular/velocity preparation and coordinate conversions, or
an explicitly equivalent comoving preparation. In angular/velocity units,
`P = P3D_comoving * W^2 * a_v / d_deg^2`, `B = P1D_velocity * W^2`, and density,
L and lp must use matching units. Check equivalence of the resulting observed
noise after conversion. Never mix comoving P with angular density moments.

Final noise remains `A * P1D_resp(k_parallel) + P_pixel`, with appropriate
coordinate conversion. The early approximation `1/(L I1)` belongs only to the
weight update; it does not replace final A. Apply squared field response once
to P3D/P1D and do not smooth additive pixel noise again.

Compute weights once per observed forest field at fiducial parameters. QSO and
LBG forests have separate densities, variances and weights even when sharing
intrinsic physics. Reuse each field's weights in all auto/cross covariances,
and hold weights, covariance, response and damping widths fixed during BAO
finite differences. Changing a diagnostic's fiducial model or response requires
recomputing its weights; perturbing AP for a derivative does not.

### Representative mode: confirmed initial choice

User confirmed on 2026-09-18: initially retain `(k_t_deg, k_p_velocity) = (2.4, 0.00035)` for both revised
profiles. This holds the mode convention fixed while testing the accuracy
changes; the corresponding comoving k and mu vary with redshift and must be
recorded. Accuracy's B still depends on its own physical response, and its P
on its own fiducial auto model. W12's B-only preparation is insufficient for
this recurrence: a representative P3D query is also required.

A fixed comoving mode, for example k = 0.07 h/Mpc and mu = 0.5, is a useful
separate sensitivity calculation, with consistent conversions of both transverse
and radial coordinates and P1D evaluation. Do not change only q_star while
claiming to have fixed both k and mu. Choice of this mode is not numerical
convergence or a claim of BAO-optimal weighting. The user also selected a
separate fixed-comoving-mode sensitivity test; its exact mode and bounded
execution scope should be specified in that later assignment.

### Solver, quadrature and implementation boundary

Carry over the tested stopping rule: rtol = 1e-4 for amplitude, normalized
shape, the full vector, A and P_pixel; at least three updates and three stable
transitions; confirmation through the actual doubled candidate count, checking
all intervening states; cap 96. Return the confirmed state. Preserve amplitude
inside the nonlinear recurrence and introduce no absolute convergence floor.
Record counts, residuals and explicit converged/capped/arithmetic/ineligible
status. A failed solve cannot supply a qualified forecast or trigger a hidden
fallback. Retain fixed-count execution as an explicitly labeled option.

Generalize eligibility from the validation module's two named forest autos to
a matching observed forest field with valid auto inputs. Keep strict public
nonnegative density/variance validation and positive finite signal/reference
requirements. Signed compatibility inputs continue through the labeled
validation path; do not relax production validation to accommodate them.
Analytic limiting-case tests may cover P = 0 without extending runtime eligibility.

Use accuracy's composite quadrature and positive quadrature weights directly in
the scalar moments. Recompute the converged solution at every magnitude
refinement. Do not interpolate a saved 107/425-node weight vector onto new
nodes or transfer the 425-node qualification to Gauss–Legendre integration.

Keep one small array-level implementation of the early and McDonald equations,
with separate strict public and literal compatibility preparation where needed.
Integrate through `weights.py`, `survey.py`/`forecast.py` and the accuracy recipe
without making the production package import validation machinery. Provisional
public method names are `early_lyaforecast` and `mcdonald`, with recorded mapping
to validation identifiers `sum_historical` and `sum_aliasing`. Preserve existing
`legacy`, `inverse_variance` and supplied-weight behavior. Do not duplicate the
nonlinear solver inside an accuracy-only recipe or silently change old API calls.

Extend profile selection, convergence studies, evidence records, cache identity
and plotting for the new methods. Record equation/method, reference mode,
response/input policy, actual quadrature, stopping controls and result. Use a
new trial-contract version for adaptive weighting; retain historical versions
1 (legacy) and 2 (W12 fixed reference) unchanged. The current cache check on
weight method alone is insufficient for distinguishing these new recipes.

## 3. Proposed bounded implementation and testing sequence

**New naming sequence:** the steps below start at **S1** and continue through
**S5**. Use these identifiers and names in subsequent assignments and reports.
This sequence is separate from the historical numbered implementation steps,
W-series diagnostics and compatibility-weighting stages; it does not renumber
or reopen that work.

### Execution authorization, 2026-09-18

The user explicitly authorized the coordinator to dispatch **Astra (light)**
implementation agents and **Sol (high)** scientific review agents, check their
results and progress through **S1–S3 inclusive**. Reviews should propose changes
only when likely to affect the scientific conclusions. This includes the scoped
real forecasts in S2/S3 and supersedes the earlier requirement for a separate
user dispatch between these steps. The user subsequently authorized **S4** with the same
agent roles and review criterion. On 2026-09-21 the user began S5 jointly with
the planning agent: summarize unresolved numerical/input sensitivity, interpret
the BAO results, and propose the next work. That initial request authorized
assessment and plan updates only; the later S5 execution authorization below
superseded that scope.

The user also explicitly directed that these forecasts run on **login nodes**,
because they are short enough not to require compute nodes. Use one numerical
thread, bounded calculations and applicable saved results. Submit no Slurm jobs.
Keep implementation/review evidence separate from scientific acceptance, and
stop after the currently authorized step. Retired current-step files remain historical.

| Step | Execution status |
| --- | --- |
| S1 | Implemented by Astra; [handoff](../reviews/s1.md), [comparison outputs](../reviews/s1-comparison-outputs.md) and [passing Sol scientific review](../reviews/s1-review.md) complete. No consequential correction requested. Coordinator reran all 10 new integration tests and verified the installed wheel's 61 modules outside the checkout. |
| S2 | Complete; [handoff](../reviews/s2.md) and [passing Sol scientific review](../reviews/s2-review.md). All six accuracy studies pass with no unresolved controls; 78 saved trials independently checked. Maximum individual/joint BAO-error refinement changes are 0.0011411%/0.0012685%. The shared three-profile run saved all 18 records in 3m51s; selected bin-1 grid refinement also passes. Saved run (local-only path: `lib/fishhighz/.validation/s2/profiles-r1/`); bin-1 diagnostic (local-only path: `lib/fishhighz/.validation/s2/bin1-grid-bao.json`). |
| S3 | Complete: [handoff](../reviews/s3.md), [passing Sol scientific review](../reviews/s3-review.md), summary (local-only path: `lib/fishhighz/.validation/s3/summary.json`) and 12 figures. Reuses the frozen S2 bundle without another forecast. All 252 selected individual/joint matrices and three comparison ratios independently reconstructed. Each profile has 78 individual and six joint results; bin-1 exclusions are retained. |
| S4 | Complete; [passing Sol scientific review](../reviews/s4-review.md), with no consequential correction outstanding. Astra produced 364 direct calculations across six bins, with forward/reverse effects for 15 settings, targeted interactions and refinement controls. [Handoff](../reviews/s4.md), [joint impact table](../reviews/s4-impact-table.md), [source inventory](../reviews/s4-profile-inventory.md) and numerical summary (local-only path: `lib/fishhighz/.validation/s4/attribution-summary.json`). Maximum componentwise endpoint mismatch is 0.000129%; maximum tested error refinement is 0.010783%. |
| S5 | Complete: Sol (medium) implementation, [handoff](../reviews/s5.md) and [passing Sol (high) scientific/software/documentation review](../reviews/s5-review.md). [Research baseline](../../research/RESEARCH_BASELINE.md), updated API documentation and standalone example (local-only path: `lib/fishhighz/examples/research_bao_forecast.py`) are available. 101 focused tests and 2 installed checks pass; independent direct Wick contraction agrees to 4.50e-17 relative. Numerical code and retained S1–S4 evidence are unchanged. D1/D2 remain deferred. |

Measured S3 endpoint differences: accuracy/fixed joint radial errors decrease
by 2.6693–3.8610%, and transverse errors by 0.6665–1.6089%. The largest
individual reduction is 8.8562%, in bin-1 lya(qso) × QSO radial error. Changing
full-compatibility to fixed-compatibility alone changes individual errors from
−44.7950% to +4.6131%. These are measured differences, with no causal
attribution or scientific acceptance implied. Numerical accuracy qualification
is the finite-refinement result recorded under S2.

S4 identifies mixed forest–galaxy damping as the dominant isolated reduction:
introducing accuracy's mean squared auto-width convention lowers joint radial
errors by 2.376–4.105% and transverse errors by 0.954–2.383%. The grouped AP
derivative change opposes it, raising radial errors by 1.020–1.174% and transverse
errors by 1.123–1.301%. These are forward controls relative to fixed-compatibility;
the linked table also gives reverse controls. Effects depend on the other
settings and are not additive. This quantifies the reconstruction/damping
assumption; it does not establish that its physical treatment is preferred.
The weighting recurrence and representative mode are shared by the endpoints.

S2 execution bounds retain the existing refinement inventory: k intervals
32/64/128/256/512 with k order 4; mu orders 8/16/32/64/128;
magnitude orders 4/8/16/32/64 on the existing composite partition; volume
orders 8/16/32/64/128; and derivative steps 1e-3 down to 6.25e-5 by successive
halving. Compare adaptive rtol 1e-4 and 1e-5 with cap 96, including actual
per-field stopping states, and check a combined refinement. Stop refinement
within these bounds and report any unresolved numerical criterion. The
fixed-compatibility bin-1 diagnostic uses 107/213/425 magnitude nodes and only
the three retained spectra. Reuse applicable bins 2–6 refinement evidence.

### Steps

1. **S1 — Profile definitions and baseline-weight integration.** Use the confirmed choices,
   implement explicit profile/method selection and per-field preparation, and
   connect records and plots. Check scalar one/two-cell updates, homogeneous
   fixed points, zero-noise and small-P limits, full-sample order invariance to
   roundoff, and agreement with saved early/McDonald trajectories on unchanged
   inputs. Test that amplitude rescaling changes an update while leaving final
   A/P_pixel invariant. Exercise cap/failure handling, response/units, field
   reuse and frozen BAO derivatives. Reproduce the legacy control, W12 explicit
   option and accepted fixed-compatibility evidence for retained selections with
   focused tests. Verify that the bin-1 selection removes all 12 excluded spectra
   before the Fisher contraction and that its three-spectrum joint agrees with
   an independent covariance solve.
2. **S2 — Weight convergence and numerical refinement.** For the QSO forest
   in bins 1–6 and LBG forest in bins 2–6, confirm iteration stability on accuracy's
   magnitude grids and check tighter stopping, e.g. rtol = 1e-5. Separately refine
   magnitude, k, mu, volume and derivative step, then check a combined
   refinement. Test effects on each individual and joint Fisher/error result;
   a dominant spectrum must not hide an unconverged weak one. Complete the
   fixed-compatibility grid diagnostic for the retained bin-1 spectra, omitting
   LBG-forest preparation there, and reuse existing bins 2–6
   evidence where the inputs are unchanged. Bound all refinements in the
   execution assignment and report unresolved controls without extending them
   automatically. Saved-input coefficient/forecast checks can precede the
   separately authorized real accuracy runs.
3. **S3 — Three-profile BAO forecast comparison.** Once implementations
   pass focused review and the real run is explicitly assigned, recompute the
   FishHighz comparison from applicable immutable reference inputs, including
   the three selected bin-1 spectra, all 15 spectra in bins 2–6, and the
   corresponding joint contractions. Numerical qualification and primary calculations may share actual trial outputs to
   avoid repeated runs. Deliver convergence status alongside every result;
   unresolved accuracy bins remain diagnostic, not accepted primaries.
4. **S4 — Attribution of fixed-compatibility/accuracy differences.** Use the
   predeclared threshold and the tests below. Include interactions where
   isolated switches do not explain the endpoint difference. Stop with the
   scientific interpretation for user review before selecting further changes.
5. **S5 — Package completion reassessment.** After reviewing these results, identify
   remaining scientific checks and document the chosen public methods and
   limitations. Review measured performance with the existing optimizations,
   and packaging/examples only where changed APIs require it. Do not resume
   old numbered assignments or expand the reference survey suite by default.

Retain the existing numerical criteria: fractional joint Fisher, joint-error
and per-spectrum error refinement changes <= 1e-3, and volume <= 1e-6. Confirm
both radial and transverse errors separately, including points omitted from
plots, and perform a combined refinement. These are finite-resolution targets,
not total error bounds. full-compatibility and primary fixed-compatibility are
reproducibility controls and need not meet accuracy-profile refinement targets.

For implementation checks, direct auto variance `2 T_FF^2 / N_modes` and cross
variance `(T_FF T_GG + T_FG^2) / N_modes`, plus independent selected-cell joint
solves, should verify the new noise propagation. Keep existing conditioning and
positive-definiteness checks; do not add jitter or regularize a discrepancy.

## 4. Comparison products and attribution

For every profile, bin and spectrum, retain Fisher matrices, marginalized
sigma(alpha_parallel), sigma(alpha_perp), their correlation coefficient,
rank/availability, and numerical status. A 2D ellipse-area ratio is a useful
additional diagnostic of changes concealed by one-dimensional errors.

The primary joint result uses the full covariance of the selected spectra:
three in bin 1 and 15 in bins 2–6. It is not a sum of individual Fisher matrices.
Keep independent redshift bins and
the existing separate AP parameters per bin; a single cross-redshift two-
parameter forecast would require a further shared-parameter assumption.
Report fractional changes fixed/full, accuracy/full and accuracy/fixed, alongside
weights versus magnitude, A, P_pixel and signal/noise comparisons needed for
interpretation. Preserve the existing 0.2 plot cut on uncertainties and paired
component masking; retain all valid values in tables and in scientific checks.
A plot cut never removes a spectrum from the joint calculation.

**Attribution trigger, confirmed by the user on 2026-09-18:** an absolute
fractional change >= 1% in either BAO uncertainty for any individual spectrum or joint bin, including
valid hidden points. Numerical refinement effects must be demonstrably smaller
(target 0.1%); otherwise resolve or report that ambiguity first. This threshold
is a practical diagnostic choice, not statistical significance. Always inspect
rank/availability changes and conspicuous ellipse-orientation changes, even if
marginal errors happen to cancel.

Start from an explicit inventory of differences, then test forward switches
from fixed-compatibility and reverse switches from accuracy:

| Change | Required controlled comparison |
| --- | --- |
| Magnitude integration | Original rectangular versus refined rectangular versus composite Gauss–Legendre at fixed support/policy; reconverge weights at each grid. |
| Density/SNR treatment | Isolate the existing negative-density extension and relevant floor/clamp/width policies from quadrature; preserve literal signed inputs only in diagnostics. |
| Fourier integration | Change k/mu nodes and measures with fixed observed cuts, evaluating the same model at the new nodes. |
| Volume | Centre approximation versus integrated volume with common powers/Jacobians; check the expected sigma proportional to V^(-1/2) limit. |
| Evaluation redshift | Separate legacy mean/covariance redshifts versus consistent z_eval, tracking all model, response, source-query and conversion dependencies. |
| Instrumental response | Legacy resolving-power convention versus physical FWHM conversion; account for its effects on mean, covariance and representative P/B. |
| Growth and damping | Isolate power-growth law, reconstruction/damping and cross-damping conventions after confirming which actually differ at the endpoints. |
| BAO derivatives | Legacy peak estimator versus full wiggle AP mapping, including volume prefactor, remapped RSD and fixed-width damping. |
| Field/pair construction | Identify any remaining legacy pair-power/noise inconsistency versus accuracy's per-field construction without inventing overlap noise. |
| Representative mode | Only if selected conventions differ, or as a separately labeled sensitivity after the primary comparison. |

For a physical switch that affects P/B or the magnitude inputs, the principal
comparison recomputes weights consistently. Where it helps interpret the
result, also hold numerical weights fixed on common magnitude nodes to separate
the direct change from its induced weighting change. Do not call such hybrids
new physical profiles. Dependencies that cannot be varied consistently alone
must be grouped and their limitation stated.

Check that all compatibility settings recover the fixed-compatibility endpoint
and all accuracy settings recover the accuracy endpoint. Compare forward and
reverse effects; if interactions exceed numerical uncertainty and matter for
the 1% trigger, test the implicated two-change combinations. Do not assume that
fractional error changes add or run a full factorial study by default.
The existing array-swap attribution chain can localize differences, but its
interpolation and grouped total-power/Jacobian replacements do not establish
physical attribution. Prefer direct model evaluation on common nodes; quantify
interpolation effects if captured legacy arrays force their use.

## 5. Deferred studies and decision record

Non-critical tests are tracked in
[DEFERRED_SCIENTIFIC_TESTS.md](../../research/DEFERRED_SCIENTIFIC_TESTS.md).
D1 is the user-agreed reference-mode study, explicitly deferred on 2026-09-21.
D2 is the optional matched early-versus-W12 inverse-variance comparison; it is
not an implementation requirement or a reopening of W12. Exact reference-mode
choices and study cases belong to those future assignments, not S5.

McDonald remains an explicit optional method with its known failures preserved.
No new McDonald survey or additional modeling study is required for S5.
The confirmed bin-1 exclusions, 1% attribution trigger and 0.1% numerical target
remain in force for any future forecast.

## S5 — Research-ready Python API plan, 2026-09-21

This plan reuses S2–S4 evidence and current source. The user confirmed a
research-ready Python API first, adopted the tested accuracy recipe as the
research baseline, and deferred the reference-mode study. These decisions do
not require changing the numerical recipe or its identity. S5 implementation, bounded validation and independent review are complete. Step 13, W12 and the old W-series remain closed.

### Numerical convergence and remaining sensitivity

There is no unresolved convergence failure in the tested revised-accuracy
baseline. All 11 required forest preparations converge, and all six bins pass
isolated and combined k, mu, magnitude, volume, derivative-step and weight-tolerance
refinements. S2 maxima are 0.0011411% individually and 0.0012685% jointly.
S4 extends the tested Fourier resolution to 256 intervals and mu order 64:
its largest individual/joint error refinements are 0.0102175%/0.0107825%, still
below the 0.1% target. These are finite tests on the selected DESI-2 case.

| Issue | What is established | What remains limited or open |
| --- | --- | --- |
| Early-weight iteration | All required populations converge with doubled-count confirmation; tightening tolerance has negligible BAO effect. | No theorem for arbitrary survey inputs or claim of globally optimal multi-tracer BAO weighting. |
| Fixed-compatibility magnitude grid | Final 213→425 refinement changes individual/joint errors by at most 0.0843233%/0.00562013% in bins 2–6, and 0.0041283%/0.0015122% in selected bin 1. | Bin-1 increments are nonmonotonic; no continuum bound. Primary 107-node compatibility grid stays fixed for reproduction. |
| Accuracy quadrature | All prescribed refinements pass; composite magnitude order 8→16 is at roundoff for this interpolation partition. | Qualification does not transfer automatically to other surveys, cuts or interpolation policies. |
| Legacy Fourier estimator | Exact legacy choices remain reproducible; S4 resolves effects of its specified extension to new nodes. | Native rectangular integration and grid-dependent peak extraction are not certified as a converged accuracy method. |
| Attribution interactions | 58 two-switch controls; largest interaction is 0.414 in 100 times the log-error ratio, individually, and 0.281 jointly. | Not every interaction was separately refined; no unique additive decomposition or exhaustive factorial study. |
| Weighting reference mode | Both endpoints retain (2.4 deg^-1, 0.00035 s/km). | The separate fixed-comoving-mode sensitivity is deferred to D1; it is not an implementation prerequisite. |
| Input policies | S4 isolates negative-density flooring and magnitude integration. | Shared SNR/clamp/sentinel and source-cell-width assumptions have not been newly qualified by these comparisons; sharing a policy does not establish its physical accuracy. |
| McDonald alternative | Explicit method and failure behavior retained. | Historical low-z LBG-forest nonconvergence includes bin 2, which remains in the forecast. This is not an unresolved failure of the accepted early-weight baseline. |

Historical prefix-weight continuum questions stay documented and closed; resolving
them is not a prerequisite for the accepted replacement's use.

### BAO constraint changes

These ranges are signed percentage changes in uncertainties, spanning all retained
bins, spectra and both AP components. Negative means smaller uncertainties. Bin 1
contains the three QSO-related spectra; bins 2–6 contain all 15, with separate AP
parameters per bin and full inter-spectrum covariance.

| Comparison | Individual error range (%) | Joint error range (%) |
| --- | ---: | ---: |
| Fixed/full compatibility: weighting replacement alone | -44.795 to +4.613 | -6.688 to +1.699 |
| Accuracy/fixed compatibility | -8.856 to +0.790 | -3.861 to -0.666 |
| Accuracy/full compatibility | -45.790 to +5.097 | -10.255 to +0.347 |

The largest weighting-only reduction is the bin-2 LBG-forest transverse auto
error; the largest accuracy/fixed reduction is the bin-1 QSO-forest × QSO radial
error. Accuracy/fixed joint radial reductions span 2.669–3.861%, transverse
reductions 0.666–1.609%, and ellipse areas shrink 2.55–4.55%. The 1% trigger is
met by 56 of 78 individual results and all six joint results. No selected result
loses rank. The [S3 report](../reviews/s3.md) preserves exact values.

The principal lessons are: weighting can materially change a forest auto forecast
while having a smaller effect on the covariance-coupled joint result; revised
baseline numerics are stable at the chosen target; and the remaining profile
differences predominantly reflect modeling/estimator assumptions. S4 identifies
mixed forest–galaxy reconstruction/damping as the largest isolated joint-error
reduction, opposed by the full AP derivative prescription. Smaller errors do not
establish a physically preferred reconstruction model or more optimal weights.

### Three implementation parts

On 2026-09-21 the user removed the sensitivity study from S5. The former S5.2
is now D1/D2 in the deferred-study file; former S5.3 and S5.4 are renumbered
S5.2 and S5.3 below. The overall S1–S5 sequence is unchanged.

1. **S5.1 — Document the chosen scientific forecast prescription.** Record the
   tested S2–S4 accuracy recipe as the recommended research baseline, as selected
   by the user. Retain early weights, the present angular/velocity reference,
   mixed forest–galaxy mean squared damping widths, and full wiggle AP
   derivatives. Explain fixed fiducial weights/covariance, nuisance parameters
   and priors, input/response conventions and finite numerical qualification.
   Preserve full-compatibility and fixed-compatibility as reproduction and
   weighting-isolation controls, with their profile and recipe identities.
   Distinguish adoption of the stated assumptions from independent physical
   validation of reconstruction. Reuse S4 evidence; no new physical prescription
   or broadband/nuisance model is introduced by this step.
2. **S5.2 — Align Python API documentation and examples with the prescription.**
   Audit public method names/options, convergence failure behavior, strict input
   handling versus validation-only legacy extensions, and available result
   metadata. The weighting methods already exist; keep method selection explicit
   and avoid a second weighting API. Refresh the README and affected docstrings
   and examples: the current README incorrectly says survey preparation and P1D
   are unimplemented and compiled support is only planned. Provide a compact
   reproducible example using early weights and a joint covariance, separating
   standalone package inputs from neighboring-checkout validation resources.
   Document the three profile recipes without promising that captured legacy
   validation is a general-purpose public profile factory. A broader CLI/config
   and general serialization release remains outside this milestone.
3. **S5.3 — Validate the changes and record completion limits.** Run focused tests
   and installed-package/example checks appropriate to the actual edits. Reuse
   frozen S2–S4 evidence for the unchanged numerical prescription; identify any
   necessary numerical correction before selecting a bounded scientific rerun.
   Do not expand by default to seven reference cases, fresh upstream capture,
   deferred sensitivities, or new performance infrastructure. The measured S2
   three-profile/refinement run took 3m51s on one thread, including about 125s of
   shared CAMB preparation; broader optimization needs a measured workload.
   Deliver an implementation/review record identifying the documented baseline,
   working examples, validation results and retained scientific limitations.

### Decisions recorded before implementation

- **Milestone:** research-ready Python API first, confirmed by the user.
- **Research baseline:** retain the tested S2–S4 accuracy recipe, including its
  mixed-pair damping convention and full wiggle AP derivatives, confirmed by
  the user on 2026-09-21. Preserve both compatibility controls.
- **Deferred tests:** D1 is a separate later study; D2 remains an optional
  proposal. Neither blocks S5. Their future mode/case choices remain in the
  [deferred-study file](../../research/DEFERRED_SCIENTIFIC_TESTS.md).
- **Remaining scientific decisions before S5:** none identified within this
  scope. Routine documentation, example and validation choices can follow the
  tested implementation. A consequential numerical/model change discovered
  during implementation would require a specific discussion.

Execution authorization, 2026-09-21: the user requested S5 implementation using
Sol (medium) implementation agents and Sol (high) review agents, with coordinator
dispatch and verification. Reviews cover scientific correctness, software and
documentation. The authorized work includes standalone synthetic examples,
focused tests and appropriate installed-package checks. D1/D2 and new real-survey
studies remain outside S5; no commits, pushes or Slurm actions are authorized.
Preserve the numerical recipe and existing evidence; discuss any consequential
model or numerical change discovered during implementation.

### S5 completion record, 2026-09-21

All three S5 parts are complete. The [research guide](../../research/RESEARCH_BASELINE.md)
records the adopted assumptions and exact evidence controls; the refreshed README,
public docstrings and standalone early-weight joint-BAO example expose the current
Python API. The guide and examples are included in the source distribution.
[Independent review](../reviews/s5-review.md) passes scientific,
software and documentation checks, including direct covariance/Fisher reconstruction.

The coordinator verified the installed example outside the checkout, exact
agreement with its source output, all 61 installed modules against the wheel and
source, and preservation of the numerical implementation and 18 retained evidence
files. Changes to package code are docstrings and three error-message occurrences.
The [handoff](../reviews/s5.md) records commands, artifacts, 101 focused
and two installed tests, and limitations. No new real-survey forecast was needed.
D1/D2, a broader release, and further scientific prescriptions remain outside
this completed assignment. No subsequent implementation step is authorized.

## Working conventions

Use AGENTS.md (local-only path: `lib/fishhighz/AGENTS.md`) for standing rules and
[FISHHIGHZ_DESIGN.md](FISHHIGHZ_DESIGN.md) for scientific/API context. This dated
revision supersedes their older scheduling statements, not historical evidence.
The user retains scientific acceptance. The explicit S1–S5 authorizations above
permitted coordinator dispatch and review within their stated scopes; they do
not authorize commits, pushes or Slurm actions. The requested comparison uses the DESI-2
15×2pt case with the confirmed bin-1 exclusions; other reference cases and fresh
upstream captures are outside scope.

Preserve the dirty checkout, saved inputs, original reports and literal legacy
arithmetic. Use focused analytic and saved-input checks before real runs. The
user directed login-node execution for these short forecasts; keep numerical
threads limited to one. Keep implementation and numerical status distinct from
scientific acceptance. Update this roadmap as the user resolves choices and
reviews results. Root planning documents are outside the independent
`lib/fishhighz` Git checkout and must be preserved separately when transferring it.
