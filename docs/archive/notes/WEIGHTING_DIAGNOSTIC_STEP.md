# W12 — done and closed

W12 revision 1 is done and closed at the user's direction on 2026-09-17.
Its bounded independent scientific review passed. The entire W-series plan is
closed; the proposed W13 continuation is withdrawn.

The current direction is the [compatibility weighting test plan](../planning/FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md)
and [implementation roadmap](../planning/FISHHIGHZ_IMPLEMENTATION_PLAN.md).
The assignment below is retained as historical text only; its dispatch and
pending-work instructions are inactive.

---

# W12 — adopt the declared fixed reference in the accuracy profile

Revision: 1. Status: proposed for user approval and implementation dispatch.
After W11 detected reference-mode sensitivity, the user explicitly selected:
**"Retain baseline and plan adoption (recommended)"**. The requested convention
is fixed inverse-variance weighting with q_star=0.00035 s/km, explicitly recorded.
This authorizes planning adoption, not implementation dispatch, acceptance of old
forecasts or a full suite. Package Step 13 remains separate.

## Scientific question and scope

Can the accuracy profile use the declared fixed reference consistently through
forest preparation, applicable refinement controls and a small Lyalpha–QSO
three-spectrum calculation, while preserving historical cumulative evidence and
compatibility behavior?

W11 gives a reason to specify the convention, not to optimize it. At q_star=0.001
s/km radial/transverse errors changed by -0.00147679%/-0.0000622378%; at 0.003
by +0.116214%/+0.220774%, relative to 0.00035 in one QSO sample. These are model/
convention sensitivities, not quadrature errors. Their existence does not change
any accuracy tolerance. No additional mode scan is assigned.

W12 has three inseparable parts: connect the existing weight method to the
accuracy recipe; remove the inapplicable iteration test from its new refinement
schedule without weakening historical validation; and validate one saved bin
including a cross-spectrum. Do not broaden this into a controller redesign.

Read workspace/package AGENTS.md, design section 0, main-roadmap handover and
progress register, and IMPLEMENTATION_STEP.md for separate package context.
Then read the [diagnostic roadmap](../planning/FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W09 review](../reviews/weighting-diagnostics-w09-review-r1.md),
[W11 review](../reviews/weighting-diagnostics-w11-review-r1.md), the live
`validation/accuracy.py`, `study.py`, `trials.py`, relevant `schema.py` checks,
and the existing `ForestInput`/weight/noise interfaces. Inspect cached-evidence
reuse in `validation/profiles.py`; preserve its scientific-source checks.

## 1. Accuracy-profile preparation and declared convention

On implementation, make the new accuracy recipe select
`method="inverse_variance"` by default. Retain an explicitly selected legacy
accuracy route for historical replay/tests; compatibility continues its existing
finite cumulative calculation. The core public weighting API/default contracts,
legacy kernels and compatibility recipe are not changed.

For each forest field at its fixed fiducial bin redshift:

```text
q_star = 0.00035 s/km
B_star = P1D(z_eval,q_star)*W_field(q_star)^2
nu_i = B_star/(B_star+l_p*v_i)
```

Use the retained independent P1D model and field-specific response once. Do not
query P3D to prepare these weights. Pass direct `alias=B_star` and omit legacy
signal, iteration count and auxiliary coordinates in the `ForestInput`. Do not
use legacy zero iterations as an implementation shortcut. Final forest noise
still uses intrinsic mode-dependent P1D(q), with response on aliasing alone:

```text
N_F(q) = [A*P1D(q)*W_F(q)^2 + P_pixel]*d_deg^2/a_v.
```

Keep q_star fixed within/across the profile's numerical-refinement trials;
B_star is field/bin-specific and may change with redshift or response. Keep
weights/noise/covariance fiducial during mean-parameter differentiation. The
legacy route retains its original P3D/P1D sampling and explicit iteration count.

Record the chosen method, q_star and units, intrinsic reference P1D, response
factor and B_star per forest in new settings/provenance. Do not imply that
q_star was fitted, optimized or rendered irrelevant by numerical convergence.
Preparation/cache identity must distinguish method and reference convention;
never reuse legacy prepared weights or an old cached primary as the new method.
New metadata must state that iterations are inapplicable, not converged.

Use the smallest host-side implementation that both `AccuracyRecipe.prepare`
and the saved-input check can exercise. A small forest-input preparation helper
inside the accuracy module is acceptable; no new public survey API or generic
configuration framework is needed. Keep the existing model, density/SNR,
response, pixel/forest scales and independence assumptions unchanged.

## 2. Applicable refinement controls and historical evidence

The live `study.py` unconditionally tries 3/6/12(/24) updates; `trials.py` and
`schema.py` require that family. Leaving those trials as repeated fixed-reference
calculations would create false iteration-convergence evidence. Make only the
method-dependent changes necessary to remove that problem.

- New fixed-reference studies retain **k, mu, magnitude, volume, derivative-step
  and combined** checks, with the existing levels, thresholds, per-spectrum
  assessment and trial/operand binding. There is no `weights` metric, no
  iteration sweep and no fake duplicate baseline standing in for one.
- The active fixed-reference numerical controls exclude `iterations`. Reject
  an explicit iteration request for that method rather than silently ignoring
  it. Legacy accuracy retains the existing six-control schedule and verdicts.
- Give the new fixed-reference trial contract a distinct explicit version
  (version 2); bind its method and fixed reference convention to the report.
  Preserve version-1 historical semantics, including legacy reports without a
  new method field. Absence of metadata in old evidence means historical legacy,
  never automatic reinterpretation as inverse variance. New fixed-reference
  evidence must explicitly identify its method/convention.
- Validation and offline controller replay must select the correct active
  families from that contract and verify matching settings. Do not make
  `weights` optional for arbitrary old reports or allow other families to vanish.
  The combined lower-control trial varies only applicable numerical controls.
- A changed reference mode is a scientific sensitivity experiment, not another
  convergence family. Keep q_star out of the numerical refinement schedule.
- The existing `three_weights` diagnostic must remain explicitly a legacy-method
  comparison (with truthful metadata), or reject that invocation as inapplicable
  to the new method. It cannot claim to vary iterations of fixed weights. Prefer
  a clear rejection here over adding a new attribution experiment.
- Preserve scientific-source/input checks on cached full-profile reuse. A cached
  historical cumulative primary cannot be relabelled as the new accuracy result.
  Do not regenerate old artifacts, alter their verdicts or migrate their files.

Version 2 is this narrow semantic distinction, not a new general evidence system.
No relaxation of tolerances, failed-trial handling, weak-spectrum checks or
scientific acceptance criteria is authorized. Galaxy-only numerical results must
remain unchanged; label inapplicable forest weighting honestly. Broader controller
performance and real convergence results belong to later work.

## 3. Small saved-input three-spectrum validation

Use only accuracy bin 0 ([2,2.235]), fields F=`lya(qso)` and G=`qso`, spectra
FF, FG and GG, parameters (`ap_0`,`at_0`). Select/remap all pairs by metadata from
the saved 15-spectrum inputs before constructing any covariance. This is a new
calculation of a two-field subset, not a subblock of a joint inverse. No LBG,
LAE, other redshift or real-survey driver is needed.

All inputs are relative to `.validation/step12-r5-20260914T191855Z/`:

| Input | SHA-256 |
| --- | --- |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |

First verify identity, selection, parameter order and the native order-32 sample
anchor. Then perform exactly **two** saved-input evaluations: reference forest
magnitude orders 32 and 64. Use saved masses/quadrature/variance for the forest,
and hold saved observed signal, mean Jacobian, Fourier nodes, mode counts and
GG noise fixed between the two. Keep FG noise zero under the existing independent
sampling assumption; do not introduce a physical overlap-noise model. This
isolates forest quadrature, not a re-integration of the quasar number density.

Exercise the actual accuracy-profile forest-input construction from part 1,
then public weight/noise/covariance/Fisher functions. A script that constructs
an independent copy of the profile weighting decision is insufficient to test
its adoption. W09's constant-background adapter may construct the public
context; its artificial volume must be discarded in favor of saved mode counts.
Use intrinsic mode-dependent P1D and saved response conventions; no P3D, CAMB,
reader, derivative, new Fourier grid or physical volume calculation.

For T_FF=S_FF+N_F, T_GG=S_GG+N_G and T_FG=S_FG, independently reconstruct
in selected order (FF,FG,GG):

```text
M*C = [[2*T_FF^2,       2*T_FF*T_FG,             2*T_FG^2],
       [2*T_FF*T_FG,    T_FF*T_GG + T_FG^2,     2*T_GG*T_FG],
       [2*T_FG^2,       2*T_GG*T_FG,             2*T_GG^2]]
F_ab = sum_n J_na^T C_n^{-1} J_nb.
```

Do not divide a cross-spectrum variance by another factor of two. Use an
independent direct NumPy solve for the 3x3 blocks, not production factorization
or contraction helpers. Report the joint 2x2 Fisher/errors and each individual
spectrum's 2x2 Fisher/errors. These individual results use their own marginal
spectrum variance, not a subblock of the joint parameter covariance. Use analytic
2x2 inversion where positive definite; no priors, regularization or pseudoinverse.
Stop and report rank/conditioning limitations if a requested error is unavailable.

## Minimal acceptance checks

1. **Profile construction.** A tiny self-contained fixture must call the actual
   recipe preparation path. Check direct inverse-variance routing, no legacy
   auxiliary P3D query, per-field response/B_star metadata, fixed derivative
   noise, and cache separation from the explicit legacy route. Use different
   synthetic forest responses if needed to catch accidental shared B_star.
   Reuse established W08 tests; do not repeat its range/keyword campaign.
2. **Controller semantics.** A tiny deterministic evaluator exercises the new
   five-control schedule, combined comparison and version-2 offline replay.
   Confirm no iteration trial or weight metric, and check an explicit iteration
   request is rejected. Two narrowly targeted negative controls must reject
   relabelling a legacy result as fixed-reference and dropping a still-required
   refinement family. No broad mutation campaign or unrelated validator repairs.
3. **Historical compatibility.** Run focused existing controller/trial tests
   under the explicit legacy route. Validate the one saved bin-0 historical
   record with `require_pass=False`; its unresolved status must remain unchanged.
   Do not require acceptance of failed historical forecasts. Keep the real
   compatibility recipe/kernel unchanged and verify its relevant source identity.
4. **Numerical subset.** At each magnitude order independently check coefficients,
   forest noise, every selected covariance block and joint/individual Fisher
   matrices. Verify the FF-only order-32 result reproduces W09's reference
   (errors approximately 0.029619522112556678, 0.027975323168394890). Check the
   cross-only variance uses both auto totals plus the squared cross signal.
5. **Bounded refinement.** Compare order64/order32 minus one for A, P_pixel,
   joint and individual BAO errors, holding the other inputs fixed. Require
   absolute changes <=0.001 (0.1%); report the actual operands/ratios. Numerical
   equivalence uses rtol=5e-12, atol=0 for nonzero quantities and exact expected
   zeros; signed contrasts use absolute 5e-12 in fractional units near zero.
   Neither numerical nor scientific tolerances may be relaxed.

Stop at a consequential implementation failure or an unsupported numerical
constraint. An independently supported refinement failure is a scientific
finding requiring review, not permission to extend the grid. A pass supports
local adoption and this one forest-quadrature comparison. It does not establish
real Fourier/redshift/derivative convergence or accept the full accuracy suite.
Synthetic controller outcomes are not survey convergence evidence.

## Edit boundary, execution and handoff

After user approval/dispatch, allow only the necessary host-side changes in
`validation/accuracy.py`, `study.py`, `trials.py`, and method-dependent checks in
`schema.py`; a minimal routing/cache change in `profiles.py` is allowed only if
required for the stated adoption. Add a short saved-input script (suggested
`scripts/check_accuracy_fixed_reference.py`), focused tests and concise README
method/convention documentation. No core numerical-kernel, compatibility-recipe,
physical input-policy, INI, CLI or general schema redesign. If the integration
needs more, report the specific obstruction instead of expanding the assignment.
Do not change plans, historical reports/artifacts, the scientific note or siblings.

Use one process and OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1.
Cap the saved-input script at 60 seconds and the focused tests at 60 seconds
separately, using the existing `.venv/bin/python -B`. No new environment, wheel,
full pytest suite, full controller on real models, survey driver, Slurm, agent
dispatch, commit or push. No broader execution is authorized by this plan.

Write `reviews/weighting-diagnostics-w12-r1.md` and a fresh
`.validation/forest-weight-diagnostics/w12-r1-<UTC>/` directory. Record the user
selected convention, exact source/input identity, changed files, actual commands,
controller/legacy checks, small matrices/errors and refinement ratios, runtime,
limitations and any early stop. Save small outputs only. Distinguish new local
numerics, synthetic controller evidence and historical results. Stop for review.

## Independent review

Inspect only the adopted convention, profile routing/cache, applicable-family
selection, old/new trial interpretation and the cross-spectrum calculation.
Rerun focused checks and independently reconstruct selected covariance blocks,
one joint 2x2 Fisher sum and all quoted error/refinement ratios. Verify no full
survey-convergence claim follows from the synthetic controller or fixed saved
Jacobian/grid. Check compatibility and historical failed verdicts are preserved.

Begin `reviews/weighting-diagnostics-w12-review-r1.md` with the scientific result
and how any proposed correction could affect it. Require changes only when
scientifically consequential, with the smallest closing test; explicitly state
when none are needed. Revise this same assignment if needed, never silently fix
production code. Passing review does not authorize W13 or a full run. The user
retains dispatch, acceptance and further execution decisions.
