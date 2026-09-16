# W09 — fixed-reference replay of one saved BAO forecast

Revision: 1. Status: proposed for user approval and implementation dispatch.
The user requested progression after W08 revision 2 passed scientific review.
This is the parallel weighting investigation, not package Step 13.

## Scientific question and reason for this step

Does the implemented `method="inverse_variance"`, propagated through the public
forest-noise, covariance and Fisher functions, reproduce the fixed-reference
BAO errors previously inferred from saved matrices in W07?

W07 independently inverted historical Fisher matrices; W08 checked the new
preparation against analytic/saved coefficients and synthetic noise/survey
examples. Neither reconstructed W07's mode sum using the new implementation.
Close this one remaining connection on the same sample before considering any
broader forecast or profile adoption. A pass establishes consistency of this
implemented calculation, not a new optimum or an independent validation of the
saved signal, derivatives, density/SNR policies or mode integration.

Read the workspace/package AGENTS.md, design section 0, the main-roadmap
handover/progress register and package IMPLEMENTATION_STEP.md for context. Then
read the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W07 review](reviews/weighting-diagnostics-w07-review-r1.md),
[W08 r2 review](reviews/weighting-diagnostics-w08-review-r2.md), and only relevant
parts of `scripts/diagnose_desi2_weights.py`,
`scripts/compare_forest_weight_bao_errors.py`, `fishhighz/weights.py`,
`noise.py`, `response.py`, `covariance.py`, and `fisher.py`.

## Bounded comparison and immutable inputs

Use only **accuracy bin 0, `lya(qso)` auto-spectrum, magnitude order 32**:
5184 magnitude nodes and 16384 saved Fourier nodes. Identify the spectrum and
parameter order from metadata; do not assume a column solely from its index.
The bin bounds are [2, 2.235]; parameters are `ap_0`, `at_0`.

Inputs below are relative to
`.validation/step12-r5-20260914T191855Z/`:

| Input | SHA-256 |
| --- | --- |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |

Check these identities and the actual sample/order/field selection. If missing
or inconsistent, report the obstruction; do not regenerate evidence. Use saved
`k`, `mu`, `modes`, the auto column of `observed_signal`, and its `observed_j`
(two parameters). Use saved magnitude, density, quadrature, variance and forest
length, with the already declared nonnegative accuracy density policy. Cross-check
`rho*quadrature` against the saved `order32_field0_masses`. No re-interpolation,
new density flooring, SNR policy change or target renormalization is allowed.

Hold the observed signal, derivatives, modes, redshift, geometry conversion,
pixel width, resolution, source distribution and B_star fixed. The only new
weight preparation is the explicit inverse-variance branch. Historical t=3/t=6
matrices are comparison denominators/numerators, not instructions to re-run
cumulative iterations. No order-64 rerun is assigned: W06/W07 already answered
that bounded refinement question.

## Calculation and scientific conventions

1. Prepare the public weights once, with `alias=B_star` and no signal,
   iterations or auxiliary sampling. Recover B_star from the matching diagnosis
   row; expected value is approximately 16.356748967852234 km/s. It is the
   **already response-smoothed** P1D at the reference parallel mode
   0.00035 s/km. It sets the fixed weights only:
   `nu_i = B_star / (B_star + l_p*v_i)`.
2. Preserve `rho=dN/(dv dm deg^2)` and the actual magnitude quadrature. Expected
   scales are L=44147.09373976799 km/s, l_p=63.328211161339375 km/s,
   a_v=102.61326549228777 (km/s)/(Mpc/h_fid), and
   d_deg=64.41180749124948 (Mpc/h_fid)/degree. Read full values from evidence;
   these printed values are identification checks, not replacement inputs.
3. Construct a public geometry context with the saved z_eval, h_fid, a_v and
   d_deg. If the saved background is insufficient for a full geometry replay,
   a diagnostic adapter through `prepare_geometry` may use constant
   H=a_v*h_fid*(1+z_eval) and D_M=d_deg/(h_fid*pi/180). Label it explicitly:
   its integrated volume is artificial and must never enter this calculation.
   All mode counts/volume factors come exclusively from saved `modes`.
   Do not bypass frozen dataclasses or add a production geometry API.
4. Use the retained analytic `default_p1d` at q=k*mu/a_v and the recorded
   instrumental response. This small analytic evaluation is permitted; no P3D,
   CAMB, new template, reader or derivative calculation is needed. Establish
   that this P1D/response reconstructs the historical calculation by the
   control below; stop if it does not. Match the response to the saved forest
   response column. Pass intrinsic P1D to `forest_noise`:
   `N(k,mu) = (A*P1D(q)*W(q)^2 + P_pixel)*d_deg^2/a_v`.
   Apply response exactly once to aliasing, never to the pixel term. Do not
   substitute constant B_star for the mode-dependent P1D.
5. Use a one-field/one-auto `PairSelection`, `combine_observed_power`,
   `gaussian_covariance` and `fisher_matrix`. Slice and remap the saved arrays
   explicitly. Keep diagonal Gaussian covariance fixed during differentiation;
   the saved Jacobian is the derivative of the mean alone. Weights and noise
   are fiducial and fixed. No joint 15-spectrum inverse or subblock thereof.

## Minimal acceptance checks and early stopping

Run in this order; stop at the first scientifically consequential discrepancy
and report its location. Do not run further samples to mask a failed check.

1. **Historical control and conventions.** With the saved t=6 coefficients,
   independently reconstruct the selected historical noise and its one-spectrum
   Fisher matrix. Compare noise to the appropriate original saved noise column
   after verifying its final controls are order 32/t=6; compare Fisher to
   `order32_iterations6_pair_fisher` in `bin-0.npz`. Check `W` against the saved
   response. This isolates signal/response/unit/pair/mode-count mistakes before
   testing the new prescription. Use scalar formulas, not a new cumulative run.
2. **Reference coefficients and noise.** Compare public A/P_pixel to the saved
   fixed coefficients and the independently reviewed W08 values:
   A=0.029419845663300446 deg^2 and
   P_pixel=0.31435492218062416 deg^2 km/s. Independently sum the three moments
   from saved masses and nu without importing the production integral kernel;
   compute N from the explicit formula above and compare every selected mode
   to `forest_noise`.
3. **Independent Fisher reconstruction.** For this single auto-spectrum,
   compute `C_n=2*(P_observed,n+N_n)^2/M_n` and
   `F_ab=sum_n J_na*J_nb/C_n` independently of the production covariance,
   factorization and contraction helpers. A simple scalar sum with `math.fsum`
   suffices. Compare C and F to the public functions and F to the selected
   `order32_iterations0_pair_fisher`. Check positive definiteness and invert
   the 2x2 F independently, e.g. using its analytic determinant.
4. **Scientific contrast.** Recover radial/transverse marginalized errors
   approximately 0.029619522112556678 / 0.027975323168394890. Report
   `100*(sigma_saved_cumulative/sigma_new_reference - 1)` for each parameter.
   Expected t=3 contrasts are +7.6305714% / +10.3256347%, and t=6 contrasts
   +8.6195623% / +11.7086249%. Reconstruct contrasts from full saved values;
   printed percentages are rounded checks. Report these as historical
   cumulative results relative to a newly reconstructed reference.

Use rtol=5e-12, atol=0 for nonzero coefficients, response, noise, covariance,
Fisher elements and errors. Expected exact zeros must remain exact. For a
ratio formed from nonzero full-precision inputs, apply the same relative
criterion; compare printed contrasts only at their quoted rounding precision.
Do not relax tolerance silently or treat 0.1% scientific refinement tolerance
as a replacement for numerical equivalence. No new convergence/optimality
claim is an acceptance criterion.

## Minimal implementation, execution and handoff

After user approval/dispatch, add at most one short diagnostic script,
`scripts/replay_inverse_variance_bao.py`, and its report. Keep independent
arithmetic in that script or one small separate check. Do not edit production
modules, existing tests, legacy code, profiles, plans or the scientific note.
No new framework, evidence schema, mutation tests, wheel, installation or broad
suite. If a production defect prevents the calculation, report it for a bounded
repair decision; do not silently fix it.

Run one process with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1,
using `.venv/bin/python -B` and a 60-second cap. Select only the one spectrum
before covariance/Fisher construction. No full forecast command, other bins,
fields, integration controls, mode-dependent weights, Slurm or agent dispatch.
If the cap is insufficient, report the measured obstruction instead of expanding
execution. No calculation is authorized merely by this planning document.

Write `reviews/weighting-diagnostics-w09-r1.md` and small outputs in a fresh
`.validation/forest-weight-diagnostics/w09-r1-<UTC>/` directory. Include the
question, held/varied quantities, exact input/source hashes and command, sample
selection, geometry-adapter limitation, independent formulas, worst residuals,
2x2 matrices, two errors and four contrasts. Save only the small outputs needed
to reproduce the conclusion; retain all original inputs. Clearly separate new
mode sums from historical comparisons. No separate general provenance audit.
Stop for user review; do not adopt a profile or plan W10.

## Independent review instructions

Check the exact assignment, decisive source and handoff; rerun the small replay
and independently verify the scalar covariance/Fisher normalization, reference
errors and contrast denominator. Reuse the established W08 coefficient result;
do not repeat its full range/keyword/provider tests. Propose changes only when
likely to alter this scientific conclusion or its supported scope, with the
smallest check closing each finding. Begin the review with a short explanation
of how each proposed change could affect the conclusion, or explicitly state
that no consequential changes are needed. Record the review separately in
`reviews/weighting-diagnostics-w09-review-r1.md`. Revise this same assignment if
needed, never silently fix production code. A pass awaits the user's acceptance
and next-step request; it does not change defaults or accept Step 12 forecasts.
