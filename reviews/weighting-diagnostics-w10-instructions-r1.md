# W10 — fixed-reference checks for a second forest population and redshift

Revision: 1. Status: proposed for user approval and implementation dispatch.
The user requested this next assignment after W09 passed independent scientific
review. This belongs to the parallel weighting investigation, not package Step 13.

## Question and bounded scope

Does the fixed inverse-variance reference retain magnitude-refinement stability
when we change the background population or redshift, and what do the saved
individual-spectrum BAO comparisons show in those two cases?

Use exactly two saved accuracy samples, in this order:

| Sample | Bounds | Field | Orders / magnitude-node counts | Parameter order |
| --- | --- | --- | --- | --- |
| A: population extension | [2, 2.235] | `lya(lbg)` | 16/32/64: 2592/5184/10368 | `ap_0`, `at_0` |
| B: redshift extension | [3.175, 3.410] | `lya(qso)` | 16/32/64: 2640/5280/10560 | `ap_5`, `at_5` |

Sample A changes population at the W06–W09 redshift; sample B retains the QSO
population at the highest saved redshift. This is two controlled extensions,
not a population-by-redshift survey. Neither requires new inputs or a mode sum.
For each sample, check coefficients first, then the small saved Fisher matrices
only if the coefficient checks permit proceeding. Stop on the first substantive
failure; do not expand the sample or repair a production module.

Read workspace/package AGENTS.md, design section 0 and main-roadmap handover/
progress register; IMPLEMENTATION_STEP.md is separate package context. Then read
the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W06 review](reviews/weighting-diagnostics-w06-review-r1.md),
[W07 review](reviews/weighting-diagnostics-w07-review-r1.md), and
[W09 review](reviews/weighting-diagnostics-w09-review-r1.md). Inspect only the
relevant parts of `compare_forest_weight_refinement.py`,
`compare_forest_weight_bao_errors.py`, `replay_inverse_variance_bao.py`, and the
public weight preparation/integral code. Preserve all existing work.

## Frozen evidence and sample identity

Paths below are relative to `.validation/step12-r5-20260914T191855Z/`:

| Input | SHA-256 |
| --- | --- |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |
| `weight-diagnosis/bin-5.npz` | `0246bc1093e98e7ac48e859e1e86a575f10a0bf6afbbd29908c791c35558fe13` |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-011.report.json` | `f2585acd91a965d90dd52a97887b3365fe0dbb25c64eacb187f181600506f09f` |

Identify profile, bin, field, physical/background labels, unique auto pair and
parameter order from each report's context; match the diagnosis bin and field
by identity. Field indices are expected to be 4 and 0 respectively, but verify
that mapping before selecting `order{order}_field{index}_...` arrays. Read
`magnitudes`, `measure`, `masses` and `variance` for orders 16, 32 and 64.
Verify the diagnosis NPZ hash binding, positive quadrature, finite nonnegative
masses/variance and nonempty density support. No raw readers or new flooring.

The report provides an independent density/quadrature anchor at its native
order: **bin 0 uses order 32, bin 5 uses order 64**. At that order compare
magnitudes, quadrature and variance, and reconstruct saved masses from
`dN/(dz dm deg^2)*(1+z_source)/c * quadrature`. At other orders use saved
`rho=masses/measure`, checking that the reconstruction does not lose support.
Do not substitute a report's native arrays for a different order.

Planning verified the listed arrays and reference/t=3/t=6 order-32 Fisher
matrices exist for both samples. Reference order-64 matrices also exist. The
bin-0 **order-64/t=6 joint diagnostic is unavailable**; it is not required here.
Do not regenerate it or infer an individual-spectrum failure from that omission.
Input absence or an identity discrepancy is an obstruction, not permission to
recapture or silently select another bin.

## Held quantities and conventions

For each sample separately, hold z_eval, z_source, L, l_p, response, density/SNR
policies and B_star fixed across magnitude orders. Read B_star from the matched
field's diagnosis `alias`; verify it is the same across the three orders and
positive. It is the already response-smoothed P1D at the historical reference
parallel mode 0.00035 s/km. Do not change B_star or recompute S/weights iteratively.

Prepare public `method="inverse_variance", alias=B_star` weights once per order,
with no signal, iteration count or auxiliary call. Reconstruct the geometry/
response context as in W09 from saved a_v, d_deg, h_fid and z_eval. If needed,
use W09's explicitly artificial constant-background adapter; its volume has no
role. Take pixel width from each sample's saved SNR policy and wavelength, and
resolution from the historically identified recipe. No new cosmology or model
calculation is needed, and response is not applied again to B_star.

With r_i=rho_i*dm_i, v_i the pixel variance, and nu_i=B_star/(B_star+l_p*v_i),

```text
I1 = sum r_i nu_i
I2 = sum r_i nu_i^2
I3 = sum r_i nu_i^2 v_i
A = I2/(L I1^2)
P_pixel = l_p I3/(L I1^2)
Q_star = A B_star + P_pixel = B_star/(L I1).
```

A has units deg^2; P_pixel and Q_star have units deg^2 km/s. The last identity
is conditional on the fixed positive B_star, common forest/pixel scales and
nonnegative measure. It verifies this reference's combined noise at the
reference mode; it does not assert a multi-mode Fisher optimum.

All physical input policies remain inherited accuracy policies, including any
negative-density flooring, extrapolation/SNR clamps and legacy source-cell
width convention recorded in these samples. List them; do not validate or
change them by implication. Across orders only magnitude quadrature varies;
across samples the physical populations/redshifts differ and must be labelled.

## Minimal checks, in order for each sample

### 1. Coefficients and magnitude refinement

- Compute public A and P_pixel for the three orders (six preparations maximum
  for the whole assignment). Independently compute moments and coefficients
  from the saved masses using scalar arithmetic and `math.fsum`, without the
  production integral helper. Verify the Q_star identity above.
- Compare A and P_pixel with the matched diagnosis `fixed_coefficients`; compare
  public and scalar values, including Q_star. These bind the implemented
  reference to the historical fixed-seed calculation used below.
- Report signed fractional changes `X_high/X_low - 1` for A, P_pixel and Q_star,
  for 32/16, 64/32 and 64/16. Require each absolute change <=0.001 (0.1%) for
  bounded refinement stability. Finite stability is not a continuum theorem.
- If numerical equivalence fails, or any refinement contrast exceeds 0.1%,
  record the quantity, operands and discrepancy and stop. This is a valid
  diagnostic outcome requiring review, not a reason to increase orders or
  relax tolerances. Do not proceed to another sample or Fisher analysis.

### 2. Saved individual-spectrum BAO comparison

Only after that sample's coefficient checks pass, select its own 2x2
`pair_fisher` matrix for orders 32 and 64, iterations=0 (the fixed initial
reference). Use the metadata `pair_fisher_key` and verified unique selected auto
pair; never take a block of the joint inverse. Invert independently with
`F=[[a,b],[b,d]]`, determinant `a*d-b*b`, and
`sigma=(sqrt(d/det),sqrt(a/det))`. Verify symmetry/positive definiteness and
compare the two errors with the saved `pair_errors` and an ordinary 2x2 solve.
No priors, regularization, pseudoinverse or new mode sum.

Report the two reference errors at both orders and signed 64/32 changes.
Apply the same 0.1% bounded-stability criterion to each error. Stop on a failure,
including a singular/ill-conditioned matrix that prevents the required numerical
agreement. Report a limitation instead of inventing a marginalized constraint.

At order 32 only, invert the saved t=3 and t=6 individual matrices and report
`100*(sigma_cumulative/sigma_reference - 1)` for each BAO parameter. Use the
order-32 reference as denominator and full saved values; check the division
independently. **The sign or size of an improvement is a result, not a pass
condition.** A reversal is scientifically useful evidence to report; it is
not grounds to change the sample or prescription. If a cumulative matrix is
unexpectedly absent/unusable, record that comparison as unavailable and retain
the supported reference findings; do not regenerate cumulative weights.

### Numerical versus scientific criteria

Numerical equivalence uses rtol=5e-12, atol=0 for nonzero moments, coefficients,
errors and independent reconstructions; expected exact zeros must be exact.
For a contrast near zero, compare its independently reconstructed value with
an absolute tolerance of 5e-12 in fractional units (5e-10 percentage points),
while retaining full-precision operands. This handles subtraction of unity
without imposing a relative precision on zero. Do not use the 0.1% stability
criterion to excuse a numerical mismatch. State distinct numerical and
scientific outcomes, and distinguish new coefficient sums from historical
Fisher calculations and their new inversions.

## Implementation boundary, execution and handoff

After user approval/dispatch, add one short standalone script,
`scripts/check_inverse_variance_samples.py`, a concise handoff and small outputs.
Reuse established conventions without refactoring existing diagnostics. No
production, existing test, profile, plan or manuscript edits. No generic evidence
framework, mutation campaign, repeated W08/W09 test suite, installation or wheel.
The identity checks and scalar reconstructions above are sufficient; add another
check only to resolve a concrete scientific ambiguity.

Use one process with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1,
`.venv/bin/python -B` and a 60-second cap for the complete assignment. No new
P3D/P1D, derivative, covariance or Fourier integration; no full forecast, other
bin/population, Slurm or agent dispatch. If the cap or a range guard is reached,
report where and stop; no automatic extension. Planning authorizes no execution.

Write `reviews/weighting-diagnostics-w10-r1.md` and a fresh
`.validation/forest-weight-diagnostics/w10-r1-<UTC>/` directory. Record the input
and decisive source hashes, command/runtime, identities/policies/held scalars,
small coefficient and refinement tables, 2x2 matrices/errors/contrasts if reached,
independent residuals, early stopping and unavailable comparisons. Preserve
original arrays; avoid large copies and mode-level output. Summarize what changed
relative to W06–W09 and what these two samples cannot establish. Stop for review.

## Independent review

Review only the exact W10 assignment, decisive source and evidence. Rerun the
bounded diagnostic and independently check one scalar coefficient reconstruction
per sample reached, all refinement ratios/criteria and the 2x2 error/contrast
arithmetic. Use a separate calculation, not an import of the submitted script.
Check native-order anchoring (32 versus 64), sample/pair identity, nonnegative
measure and the distinction between historical matrices and new calculations.
No need to repeat W09's mode sum or W08's arithmetic-range campaign.

Begin `reviews/weighting-diagnostics-w10-review-r1.md` with the scientific finding
and a short explanation of how any proposed changes could affect it. Propose
changes only if likely to affect the conclusion or its supported scope; give the
smallest closing check, or state that no consequential correction is needed.
A supported instability or reversal is a scientific result, not automatically
an implementation defect. Revise this same assignment only for consequential
repairs; never silently fix production code. A review pass does not adopt a
profile or authorize W11. The user controls acceptance and further progression.
