# Forest weighting diagnostic W07, revision 1: one-spectrum BAO error comparison

Status: proposed for user approval and implementation-agent dispatch. W06 passes
independent scientific review; the user requested the next diagnostic. This is
an analysis of four existing Fisher matrices, not a new forecast run or adoption
of a weighting prescription. Package Step 13 and forecast acceptance remain
unchanged. Created: 2026-09-15.

## 1. Scientific question and minimum scope

How do the explicit inverse-variance reference and cumulative t=3/t=6 weights
affect marginalized radial and transverse BAO errors for one QSO-forest
auto-spectrum in the same accuracy-profile redshift bin?

W06 established empirical magnitude stability of the reference over orders
16/32/64. The cumulative t=3 coefficients also pass its 0.1% criterion; t=6
fails only P_pixel at 64/16 (+0.129453%), while adjacent-order changes pass.
At fixed order, t=3 to t=6 raises Q at the weighting mode by 1.7039–1.7533%.
Cumulative Q exceeds the reference by about 14.7% and 16.7%, respectively.
These are not BAO-error percentages. The reference minimum is conditional on
the nonnegative measure and the single-mode covariance approximation.

Step 12 already saved and independently checked the individual-spectrum Fisher
matrices needed here. Reuse them: three weight choices at magnitude order 32,
plus the reference at order 64. W07 reports their previously unisolated
comparison in the diagnostic sequence. No weighting trajectory, coefficient
refinement, model evaluation or Fourier-mode Fisher sum needs to be repeated.
The new arithmetic is four 2x2 inversions, checks and fractional contrasts.

Select only `lya(qso)_lya(qso)`, accuracy bin 0, z=[2.0,2.235], with parameters
`ap_0`, `at_0`. Do not include cross-spectra, other tracers or redshift bins.
A fixed-mode noise minimum does not guarantee a minimum of both integrated BAO
errors when the same weights are used across modes. Report the actual result,
including a trade-off or no improvement if found; no preferred ranking is an
acceptance condition.

## 2. Relevant history, source and immutable inputs

Read workspace/package AGENTS.md, design section 0, main roadmap handover/
progress register, the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W06 report](reviews/weighting-diagnostics-w06-r1.md),
[W06 review](reviews/weighting-diagnostics-w06-review-r1.md), and the fixed-weight
and individual-spectrum paragraphs in the
[Step 12 r5 report](reviews/step-12-r5.md) and
[review](reviews/step-12-review-r5.md). W05 supplies the reference derivation.
The [reviewed W06 instructions](reviews/weighting-diagnostics-w06-instructions-reviewed-r1.md)
are archived. `IMPLEMENTATION_STEP.md` remains the separate package assignment.

Inspect `scripts/diagnose_desi2_weights.py` for the saved fixed-seed/count
variants and their common Fourier grid, J, model and geometry. Inspect
`fishhighz/validation/numerics.py::wick`, `contract` and `summaries` to establish
that `pair_fisher` is the information from each spectrum selected independently.
It is not a block extracted from the joint inverse covariance. Use the recorded
rank convention in `fishhighz/_information.py`; no production edits or imports
of forecast/controller paths are needed.

Inputs relative to the package root, under
`.validation/step12-r5-20260914T191855Z/`:

| File | SHA-256 |
| --- | --- |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |

Check `diagnosis.bins[0]` references the matching NPZ/hash. Verify the report's
accuracy/bin/bounds, field identities, selected-pair ordering and parameter
ordering; compare the primary NPZ `selected_pairs` with that report. Resolve the
index for field pair [0,0] semantically. Its current index is 0; do not select a
matrix solely because it is first. The two parameter names must be ap_0, at_0,
with their order explicitly recorded.

Select exactly these available diagnosis rows and their named `pair_fisher_key`:

| Interpretation | Magnitude order | Stored iterations | Key |
| --- | ---: | ---: | --- |
| Fixed explicit reference nu | 32 | 0 | `order32_iterations0_pair_fisher` |
| Cumulative rule t=3 | 32 | 3 | `order32_iterations3_pair_fisher` |
| Cumulative rule t=6 | 32 | 6 | `order32_iterations6_pair_fisher` |
| Reference refinement check | 64 | 0 | `order64_iterations0_pair_fisher` |

`iterations=0` means the fixed initialization, not zero weights/no noise.
Verify `fixed_initial_weights=True` only for those reference rows, and
`available=True` for all four. These arrays have shape (15,2,2); extract only
the selected auto-spectrum. Associated `pair_covariance`, `pair_errors`,
`pair_correlation`, `pair_constrained` and `pair_rank` share each row's prefix.
Match the order-32/t=6 selected Fisher matrix to `records-001.npz::pair_fisher`
at the same pair index; it is the saved primary (controls order 32, count 6).

The order-64/t=6 joint diagnostic row is unavailable because the LBG forest
coefficient failed historically, although its QSO coefficient is available.
It is deliberately outside W07. Do not substitute another row, reconstruct
missing matrices, repair the LBG calculation or interpret that absence as a
failure of the QSO-only observable.

## 3. What is held fixed and what these matrices mean

The historical comparison holds the observed signal, BAO Jacobian, Fourier
nodes, mode counts, geometry/volume, model parameters, response, input policies
and redshift bin fixed. It changes forest weights and the corresponding noise.
The historical controller also changed the LBG-forest noise; demonstrate from
the independently selected auto-spectrum covariance that this does not enter
the QSO-forest pair Fisher matrix:

```text
T00(k,mu) = P00_observed(k,mu) + N0(k,mu)
C_auto(k,mu) = 2 T00(k,mu)^2 / modes(k,mu)
F_ab = sum_modes J_auto,a J_auto,b / C_auto.
```

This equation establishes the selected information's meaning; no new mode sum
is assigned. Each two-parameter covariance is F^-1. Its diagonal yields errors
marginalized over the other BAO parameter, not 1/sqrt(F_aa). No extra priors,
nuisance parameters, parameter ties, derivative regeneration or covariance
derivatives are introduced. Do not invert the joint 15-spectrum covariance.

The reference uses nu=B_star/(B_star+l_p v) prepared at the single historical
weighting mode and held fixed as a function of magnitude across forecast modes.
The actual mode-dependent noise uses the existing P1D and response. Do not
reinterpret these saved matrices as weights reoptimized separately at each
k_parallel, or rescale Fisher matrices/errors by the W06 Q ratios. The
covariance contains the signal too, and radial/transverse derivatives weight
modes differently.

Retain the existing accuracy input policies from W06, including density
`floor_negative`, SNR `legacy_floor_clamp` and legacy source-cell widths.
This comparison neither changes nor physically validates those assumptions.
The reported errors refer to dimensionless alpha_parallel and alpha_perp;
if quoting percentages of fiducial values, state alpha_parallel=alpha_perp=1.

## 4. Four-matrix calculation and concrete checks

1. Verify hashes, row/pair/parameter identity and primary t=6 matrix agreement.
   Reported constraints must include both parameters with rank 2. Independently
   check finite positive diagonals and symmetry/positive definiteness in the
   diagonally normalized basis, using the existing 64*eps*n rank convention.
   If a selected matrix is singular, indefinite, missing or materially
   inconsistent, hand off that finding. Do not add priors, a pseudoinverse or
   regularization to manufacture finite errors.
2. Compute covariance and errors with a direct float64 solve. Independently use
   80-digit Decimal scalar arithmetic on the exact saved binary matrix entries:

   ```text
   F = [[a,b],[b,d]], determinant = a*d-b*b > 0
   Cov = [[d,-b],[-b,a]] / determinant
   sigma_parallel = sqrt(d/determinant)
   sigma_perp = sqrt(a/determinant)
   correlation = -b/sqrt(a*d).
   ```

   Check symmetry before using this expression; only average sub-tolerance
   off-diagonal roundoff in a local copy, consistent with the existing rule.
   Compare independently computed covariance, errors and correlation with the
   stored pair summaries. Require rtol=5e-12, atol=0 on nonzero quantities and
   exact reference-zero handling. No general evidence validator is needed.
3. Tiny analytic control: F=[[4,1],[1,9]] has determinant 35,
   Cov=[[9,-1],[-1,4]]/35, errors (3/sqrt(35),2/sqrt(35)), and correlation -1/6.
   Verify the solve/scalar calculations and the use of marginalized errors.
   This one control suffices; no broad singular-matrix test framework.
4. Report all four F matrices, sigma_parallel, sigma_perp and correlation.
   For order 32, report signed fractional error changes t=3/nu, t=6/nu and
   t=6/t=3, with the denominator named. A positive change means larger error.
   Compare the reference at order 64 with order 32 using both error changes
   and ||F64-F32||_F/||F32||_F. Classify only this reference pair as empirically
   stable if both error changes and the Fisher discrepancy are <=1e-3.
   Check dimensionless contrasts with independent scalar arithmetic to absolute
   tolerance 5e-12; use absolute changes if a reference is exactly zero.

There is no requested outcome for the error ranking. A finite difference above
0.1% is a sensitivity result, not an implementation failure. A reference
refinement pass reproduces a bounded historical check, not convergence in k,
mu, redshift, all magnitude orders or the full forecast. No prescription is
accepted by these calculations.

Stop after these four matrices and the control. Do not calculate another
spectrum or bin, recover unavailable states, create modes, evolve weights or
rerun a profile/forecast. Existing Step 12 review supplies the independent
validation of the stored Fourier-mode contractions; W07's new independent
checks concern selection, matrix inversion and the specified contrasts. State
this distinction explicitly rather than claiming a fresh end-to-end forecast.

## 5. Minimal implementation, handoff and review

Allowed additions: `scripts/compare_forest_weight_bao_errors.py`, optional
focused tests, `reviews/weighting-diagnostics-w07-r1.md`, and a small output
under a new `.validation/forest-weight-diagnostics/w07-r1-<UTC>/` directory.
One short script with assertions is sufficient. Save the four selected matrices,
row/pair/parameter mapping, derived tables, independent discrepancies, exact
commands and relevant input/source hashes. Reference historical arrays rather
than copying full bundles. Preserve all prior code, reports and artifacts;
do not edit planning documents during implementation.

Run only the bounded calculation/control and Ruff on new Python, with the
existing interpreter, one process and OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and
MKL_NUM_THREADS set to 1. Cap each numerical invocation at 30 seconds and retain
partial results if reached; no automatic extension. No full pytest, wheel,
new environment, generic replay framework, validator hardening or optimization.

The report leads with the effect on each marginalized BAO error and its limits.
Explain whether the one-mode Q ranking carries over to this one-spectrum
forecast, whether t=3/6 dependence matters for its errors, and the reference's
32/64 stability. Distinguish W06 coefficient results, historical Step 12 Fisher
matrices and new W07 derived comparisons. The actual forecast calculations are
historical; a new interpretation is not a fresh survey run or full-suite pass.

The independent reviewer checks the decisive source/mapping, primary match,
tiny control and four-matrix arithmetic. Write
`reviews/weighting-diagnostics-w07-review-r1.md`, beginning with **Scientific
conclusion and impact of proposed changes**. Propose corrections only if likely
to affect the scientific conclusion or its justified scope; state the possible
impact and smallest resolving check. Explicitly state when none are needed.
No optional cleanup or hypothetical hardening list. Stop for user review.

No production weighting or input-policy adoption, mode-dependent reoptimization,
new physical evaluation, Fisher mode sum or forecast execution, additional
population/bin/step, package Step 13 repair, Slurm, agent dispatch, commit or
push is authorized. Consequential scientific decisions remain with the user.
