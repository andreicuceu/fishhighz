# Forest weighting diagnostic W06, revision 1: saved magnitude-refinement comparison

Status: revision 1 implemented and independently reviewed; no scientifically
consequential correction is needed. Awaiting user review and acceptance; no
progression or production adoption is authorized. The scientific assignment
below is unchanged. See [review](reviews/weighting-diagnostics-w06-review-r1.md)
and [original instructions](reviews/weighting-diagnostics-w06-instructions-r1.md).
Package Step 13 and forecast acceptance remain unchanged. Created: 2026-09-15.

## 1. Scientific question and existing evidence

For one saved QSO-forest bin with a nonnegative accuracy-profile measure, how
stable are A, P_pixel and Q=A B+P_pixel under the existing magnitude refinements
when using the explicit inverse-variance reference? How do the original
cumulative weights at fixed counts t=3 and t=6 compare on those same inputs?
Separate refinement at fixed iteration count from changes in iteration count.

W05 establishes that nu=B/(B+l_p v) is the legacy seed and minimizes Q for a
nonnegative measure under the stated diagonal-covariance approximation. Its
signed compatibility example has Q excesses of 14.3618%, 12.8358% and 19.1309%
for prefix t=3 and full-sample fixed/refreshed t=24, respectively. That ranking
alone proves no physical optimum. W04's bounded iteration stability establishes
neither a minimum nor magnitude convergence.

There is already relevant historical evidence: the
[Step 12 r5 review](reviews/step-12-review-r5.md) verifies the fixed-initial-weight
magnitude control, with order-32/64 individual-spectrum error differences at
roundoff. W05 identifies that fixed function with the inverse-variance
reference. Do not present W06 as the first fixed-weight convergence test or
rerun the historical forecasts. W06 connects that evidence to explicit A,
P_pixel, Q and n_eff comparisons, with the cumulative rule as a matched-input
control. Only saved arrays and scalar metadata are needed.

Scope: `lya(qso)`, accuracy bin 0, z=[2.0,2.235], orders 16/32/64, and three
vectors at each order: nu and saved prefix t=3/t=6. Nine coefficient evaluations;
no new recurrence, survey/model evaluation, Fisher sum or forecast. The
reference is a diagnostic candidate, not an adopted prescription.

## 2. Read and identify the saved inputs

Read workspace/package AGENTS.md, design section 0, main roadmap handover/
progress register, the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W05 report](reviews/weighting-diagnostics-w05-r1.md),
[W05 review](reviews/weighting-diagnostics-w05-review-r1.md), and the fixed-weight
paragraphs in the [Step 12 r5 handoff](reviews/step-12-r5.md) and review above.
The [reviewed W05 instructions](reviews/weighting-diagnostics-w05-instructions-reviewed-r1.md)
are archived. `IMPLEMENTATION_STEP.md` remains package Step 13.

Inspect only the relevant source: W05's coefficient/independent-sum helpers;
`scripts/diagnose_desi2_weights.py` for saved masses, seed, smoothing and input
policies; `scripts/diagnose_weight_limit.py::_field_inputs` for scalar mapping;
and `fishhighz/kernels/weights.py::_integrals` for product order. Do not execute
those controllers or import their model/forecast preparation paths.

Authoritative historical inputs, relative to the package root:

```text
.validation/step12-r5-20260914T191855Z/weight-diagnosis/diagnosis.json
SHA-256: 069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008
.validation/step12-r5-20260914T191855Z/weight-diagnosis/bin-0.npz
SHA-256: d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73
.validation/step12-r5-20260914T191855Z/profiles-checked/records-001.report.json
SHA-256: 9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666
```

Read `diagnosis["bins"][0]`, whose NPZ reference/hash must agree above. For each
n in (16,32,64), select the `lya(qso)` field metadata from a row of that order.
Read NPZ keys `order{n}_field0_` followed by `magnitudes`, `measure`, `masses`,
`variance`, `weights_3` and `weights_6`. Verify field identity rather than relying
only on index 0. The node counts are 2592, 5184 and 10368: n is the Gaussian
order per partition interval, not the total node count. Retain the saved common
magnitude partition and support [16.1,26.75]; these Gauss nodes are not nested.

Read the accuracy report's `settings.samples["lya(qso)"]` and model `z_eval`.
Confirm profile, bin, bounds and field identity. Derive the common scalars by
the existing mapping:

```text
L = sample["length_velocity"]                 # 44147.09373976799 km/s
z_source = sample["z_source"]                 # 2.3830099582078716
l_p = c * pixel_width_angstrom / [1215.67 * (1+z_eval)]
c = 299792.458 km/s; pixel_width_angstrom = 0.8
z_eval = 2.1152848986890427
P0 = field_metadata["signal"]                 # 1.7182791088399005 deg^2 km/s
B = field_metadata["alias"]                   # 16.356748967852234 km/s
```

Check P0 and B agree across the three orders. At order 32, independently match
magnitudes, quadrature and variance to the accuracy report; match masses to
`sample.density * (1+z_source)/c * sample.quadrature`. Require rtol=5e-12,
atol=0 with exact reference-zero handling. Check positive quadrature, ordered
in-support nodes, nonnegative masses of positive total, finite nonnegative
variance and the common partition; sum of quadrature weights must match its
interval length within this tolerance.

Record the inherited density `floor_negative` and SNR `legacy_floor_clamp`
policies and relevant scalar settings from source/report. They were applied in
Step 12; W06 must neither apply them again nor modify them. A nonnegative
result here is conditional on that existing input model, not evidence that its
flooring or extrapolation is physically correct. W05's signed 107-node inputs
remain historical evidence; do not compare profile-to-profile differences as
if quadrature alone changed. Missing or inconsistent inputs require a bounded
handoff, not source regeneration or repairs to Step 13 evidence machinery.

## 3. Coefficients and controlled comparisons

Use the saved masses directly as r_i: they already include quadrature and the
redshift-to-velocity density conversion. Do not multiply them by `measure`, dm,
or (1+z_source)/c again. For each order form

```text
D_i = B + l_p v_i
nu_i = B / D_i
I1 = sum r_i w_i
I2 = sum r_i w_i^2
I3 = sum r_i w_i^2 v_i
A = I2 / (L I1^2)
P_pixel = l_p I3 / (L I1^2)
Q = A B + P_pixel.
```

For historical replay preserve `rw=r*w`, `rww=rw*w`, `rwwv=rww*v`, final
cumulative sums and the sequential divisions in `_integrals`. Keep arithmetic
guards; no underflow suppression or regularization. At zero-mass nodes the
choice of nu is immaterial to these sums; do not alter masses or mask positive
support. No displayed normalization feeds any update; no update is performed.

For nu also compute K=sum(r/D), n_eff=L B K and Q=1/(L K)=B/n_eff. P0 is held
fixed but absent from D and Q. B is already smoothed at the historical weighting
mode; l_p v is the individual 1D instrumental power. Apply no extra response or
noise term. A has units deg^2, n_eff deg^-2, and B, l_p and L km/s; Q and
P_pixel have units deg^2 km/s. r has units deg^-2 (km/s)^-1.

1. Reconstruct A and P_pixel for nu and the two saved cumulative vectors at
   each order. Compare nu with metadata `fixed_coefficients`, and prefix values
   with the matching field metadata `rows` entries at counts 3/6. Require their
   recorded availability. These are recomputations of historical coefficients.
2. Report A, A B, P_pixel, Q and (for nu) n_eff by order and vector. The new
   interpretation is the combined-noise comparison tied to W05's reference.
3. For each vector definition, report signed fractional changes 32/16, 64/32
   and 64/16 for A, P_pixel and Q. Label the reference empirically stable over
   the tested orders if all nine absolute changes are <=1e-3. Apply the same
   descriptive test separately to prefix t=3 and prefix t=6. Show individual
   comparisons even if the combined criterion fails; do not infer nonconvergence
   from a failed coarse comparison or a continuum limit from a pass.
4. At each order report prefix Q/Q_nu-1 at counts 3 and 6, and t=6/t=3 changes
   of A, P_pixel and Q. Keep these iteration comparisons separate from the
   order comparisons. Check Q_prefix>=Q_nu within numerical tolerance under the
   verified nonnegative measure. A significant violation signals a mapping or
   arithmetic error requiring investigation, not a superior estimator.

This test cannot attribute compatibility/accuracy differences, establish a
joint limit in order and iteration count, or prove convergence of the cumulative
rule as t grows. It does not compare full-sample W04 branches on accuracy inputs.
A stable reference and sensitive cumulative control support the narrower claim
that resolving a fixed weight function does not resolve the discrete recurrence.
If both are stable over this bounded selection, state that result without
contradicting historical failures at other counts/orders or claiming their cure.

## 4. Minimal acceptance checks and stopping conditions

- Reuse W05's exact two-bin control with nonuniform quadrature: q=(1/2,2),
  rho_v=(2,1/2), hence r=(1,1); v=(1,4), L=l_p=B=1. Obtain nu=(1/2,1/5),
  A=29/49, P_pixel=41/49, Q=10/7 and n_eff=7/10. Check with rational arithmetic;
  this specifically checks that the saved integration measure is used once.
- Independently compute all nine real-input moment/coefficient sets using
  scalar 80-digit Decimal sums of exact saved binary values. For nu, independently
  form weights from B, l_p and v. Do not call the float moment helper. Check
  direct Q=sum(r D w^2)/(L I1^2) against A B+P_pixel, and reference Q against
  B/n_eff. Check the three refinement contrasts and within-order contrasts.
- Use rtol=5e-12, atol=0 for nonzero coefficient/moment comparisons; check
  reference zeros exactly. Fractional contrasts require absolute agreement
  <=5e-12, avoiding relative comparisons of nearly zero differences. Classify
  stability from independently checked operands, not a stored verdict.

Stop on an identity mismatch, changed policy/scalar across orders, missing
checkpoint, nonfinite arithmetic, D<=0, I1=0, or violation of the stated input
domain. Confirm a numerical obstruction independently and retain it; do not
change guards, floor inputs or extend the calculation. A stable or unstable
result is a scientific answer, not a reason to add orders. Stop after 16/32/64;
no adaptive additions, other populations, new trajectories or raw interpolation.

One short standalone script with assertions and a concise report is sufficient;
a separate focused test file is optional. Run only these checks and Ruff on new
Python. Use the existing interpreter, one process and OMP_NUM_THREADS,
OPENBLAS_NUM_THREADS and MKL_NUM_THREADS set to 1. Cap each numerical invocation
at 30 seconds; preserve partial results on a cap and do not auto-extend. No
full test suite, wheel, environment, generalized validator or mutation framework.

## 5. Handoff and independent review

Allowed additions: `scripts/compare_forest_weight_refinement.py`, optional focused
tests, `reviews/weighting-diagnostics-w06-r1.md`, and small outputs under a new
`.validation/forest-weight-diagnostics/w06-r1-<UTC>/`. Preserve existing scripts,
production code, saved artifacts, reports and planning documents during
implementation. Record source/input hashes, exact commands, actual checks,
selected metadata, coefficient tables, independent discrepancies and stop reason.
Reference immutable input vectors rather than copying entire historical bundles.

Lead with the scientific answer and distinguish:

- Historical fixed-seed magnitude stability already reviewed in Step 12.
- W05's identification of that seed and its nonnegative-measure minimum.
- W06's independently reconstructed coefficient/refinement/control comparison.
- Limits imposed by fixed mode, finite orders/counts and the inherited accuracy
  input policy; no physical validation of that policy or forecast improvement.

The reviewer checks the measure/scalar mapping, decisive source and tiny control,
and independently reconstructs the nine coefficient sets and contrasts. Write
`reviews/weighting-diagnostics-w06-review-r1.md`, beginning with **Scientific
conclusion and impact of proposed changes**. Propose corrections only if likely
to affect the scientific conclusion or its justified scope, giving the possible
impact and smallest resolving check. State explicitly when none are needed.
Do not add optional cleanup, hypothetical hardening or unassigned calculations.

Stop for user review. No production adoption, new density/response convention,
new survey/model evaluation, Fisher sum or forecast, full-sample/cumulative
trajectory, further step, package Step 13 repair, Slurm, agent dispatch, commit
or push is authorized. Any consequential scientific choice remains with the user.
