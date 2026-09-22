# W07 revision 1: one-spectrum marginalized BAO errors

**The fixed inverse-variance reference gives the smallest radial and transverse
BAO errors among the three saved order-32 choices.** For the accuracy-profile
QSO-forest auto-spectrum in bin 0, cumulative t=3 increases these errors by
**7.63057% and 10.32563%**, respectively, relative to the reference; t=6 increases
them by **8.61956% and 11.70862%**. Changing t=3 to t=6 raises them by
**0.91888% and 1.25355%**. Finite-count sensitivity therefore exceeds 0.1% for
both errors. The order-32/64 reference comparison is empirically stable to
roundoff under the assigned 0.1% criterion.

W06's single-mode Q ordering carries over to both marginalized errors in this
particular comparison. Its Q percentages are not BAO-error percentages and were
not used to rescale information or errors. These results establish neither a
multi-mode optimum nor a production prescription. **Implementation is complete
within W07 scope; awaiting independent/user review.**

The forecast mode sums are historical Step 12 calculations. The new W07 work is
matrix selection, inversion and comparison, with independent scalar arithmetic
and one analytic control. No new forecast or Fourier-mode contraction was run.

## Selected information and results

Only `lya(qso)_lya(qso)`, accuracy bin 0, z=[2.0,2.235], is selected. Matrix
rows/columns and errors are ordered `ap_0`, `at_0`, corresponding to dimensionless
alpha_parallel and alpha_perp. Errors below are absolute dimensionless errors;
at fiducial alpha_parallel=alpha_perp=1 their percentage values are 100 sigma.
Each error is marginalized over the other BAO parameter, with no added prior or
nuisance parameter.

| Magnitude order | Weights | F_parallel,parallel | F_parallel,perp = F_perp,parallel | F_perp,perp |
| ---: | --- | ---: | ---: | ---: |
| 32 | reference nu | 1419.7381662460089 | 667.4320627022987 | 1591.5274282829370 |
| 32 | cumulative t=3 | 1227.5500924474698 | 564.8346409624346 | 1309.6749369362570 |
| 32 | cumulative t=6 | 1205.5387504892117 | 553.0960157089298 | 1277.7024563498703 |
| 64 | reference nu | 1419.7381662460130 | 667.4320627023006 | 1591.5274282829423 |

| Order | Weights | sigma_parallel | sigma_perp | correlation |
| ---: | --- | ---: | ---: | ---: |
| 32 | nu | 0.029619522112556678 | 0.027975323168394890 | -0.44401319174889775 |
| 32 | t=3 | 0.031879660901309670 | 0.030863952834406856 | -0.44547153034773057 |
| 32 | t=6 | 0.032172595269768860 | 0.031250848833533480 | -0.44565119739870906 |
| 64 | nu | 0.029619522112556633 | 0.027975323168394843 | -0.44401319174889770 |

Signed fractional changes are sigma(numerator)/sigma(denominator)-1, not
percentages; positive means larger uncertainty.

| Numerator / denominator | Radial change | Transverse change |
| --- | ---: | ---: |
| t=3 / nu, order 32 | +0.07630571418958998 | +0.10325634662463501 |
| t=6 / nu, order 32 | +0.08619562285678639 | +0.11708624938564127 |
| t=6 / t=3, order 32 | +0.009188754214357164 | +0.01253552975545369 |
| nu order 64 / nu order 32 | -1.5543122344752192e-15 | -1.7763568394002505e-15 |

The reference Fisher discrepancy ||F64-F32||_F/||F32||_F is
**3.0790395131128374e-15**. It and both absolute reference error changes are
below 1e-3. This is only a bounded empirical reference refinement check, not
convergence in k, mu, redshift, arbitrary magnitude orders or the full forecast.

## Meaning and historical provenance

Source inspection of `fishhighz/validation/numerics.py::wick`, `contract` and
`summaries` establishes that `pair_fisher` contains independently selected
spectrum information: `contract` divides each spectrum's Jacobian product by
its own covariance diagonal. It does not extract a block of the joint inverse
covariance. For this auto-spectrum,

```text
T00 = P00_observed + N0
C_auto = 2 T00^2 / modes
F_ab = sum_modes J_auto,a J_auto,b / C_auto
Cov_parameters = F^-1
```

Consequently, the historical LBG-forest noise changes do not enter this QSO-only
information. `scripts/diagnose_desi2_weights.py` uses the same saved observed
signal, observed J, Fourier nodes and mode counts while changing forest noise.
Geometry, volume, model parameters, response and input policies are common.
This source argument is not a new execution of the mode sums. Their independent
validation remains the historical [Step 12 r5 review](step-12-review-r5.md),
including its saved individual-spectrum and fixed-weight checks.

The reference nu=B_star/(B_star+l_p v) is the fixed initialization (`iterations=0`),
not zero weights. It is prepared at the historical single weighting mode and
held fixed as a function of magnitude across forecast modes. The saved
mode-dependent noise retains the existing P1D and response; these matrices do
not represent weights reoptimized at each k_parallel.

## Identity and actual checks

Inputs are under `.validation/step12-r5-20260914T191855Z/`. All four prescribed
SHA-256 identities matched before and after execution:

| Input | SHA-256 |
| --- | --- |
| weight-diagnosis/diagnosis.json | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| weight-diagnosis/bin-0.npz | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |
| profiles-checked/records-001.report.json | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| profiles-checked/records-001.npz | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |

The diagnosis bin names the matching NPZ/hash. Report context/settings agree on
accuracy, bounds, fields and parameter ordering. Field 0 is identified by its
`lya(qso)` name, forest kind, physical `lya` and `qso` background. Its [0,0]
auto-pair resolves uniquely to selected-pair index 0; the entire primary NPZ
selected-pair ordering matches the report. No matrix is selected solely by its
position.

Exactly four unique rows were used: `order32_iterations0_pair_fisher`,
`order32_iterations3_pair_fisher`, `order32_iterations6_pair_fisher`, and
`order64_iterations0_pair_fisher`. Each has shape (15,2,2), `available=True`,
and `fixed_initial_weights=True` only for the two reference rows. Only the
selected 2x2 matrix and its pair summaries were evaluated. The order-32/t=6
matrix matches primary `records-001.npz::pair_fisher` with maximum elementwise
relative discrepancy **3.5591e-16**; primary controls specify order 32/count 6.

All four saved ranks are 2, and both parameters are constrained. Independent
checks find finite positive diagonals, exactly symmetric matrices and positive
normalized eigenvalues. The smallest normalized eigenvalue is 0.5543488,
well above the corresponding 64*eps*2 spectral tolerance (about 4.11e-14).
No off-diagonal averaging was needed. No prior, pseudoinverse or regularization
was used.

Direct float64 solves were compared with independent 80-digit Decimal scalar
inversions of exact saved binary entries. For F=[[a,b],[b,d]], this uses
Delta=a*d-b*b>0, covariance=[[d,-b],[-b,a]]/Delta,
errors=(sqrt(d/Delta),sqrt(a/Delta)), and correlation=-b/sqrt(a*d).
Maximum solve/Decimal relative discrepancy is **2.474e-16**. Comparison with
stored pair covariance/errors/correlations gives at most **7.368e-16**.
Both pass rtol=5e-12, atol=0, with explicit exact-reference-zero handling.
All error contrasts agree with independent Decimal ratios to absolute
**2.013e-16**; the Frobenius contrast agrees at float64 rounding precision,
within the assigned absolute 5e-12 tolerance.

The sole analytic control F=[[4,1],[1,9]] reproduces determinant 35,
covariance=[[9,-1],[-1,4]]/35, errors=(3/sqrt(35),2/sqrt(35)) and correlation
-1/6. Solve/scalar/analytic comparisons pass (maximum relative discrepancy
2.221e-16). The marginalized errors differ from (1/2,1/3), explicitly checking
that conditional 1/sqrt(F_aa) errors were not substituted.

## Execution and limits

New standalone script (local-only path: `../scripts/compare_forest_weight_bao_errors.py`) and
summary.json (local-only path: `../.validation/forest-weight-diagnostics/w07-r1-20260916T000045Z/summary.json`)
contain the four matrices, mapping, covariances, errors, correlations, scalar
results, discrepancies, source/input hashes and timing. Full historical arrays
were not copied. Script SHA-256:
`bcb7004e9e3b9be2fd2cbe80df2bbe9eafc09669e2e156b5dfcd2b62ff8e875a`.

Commands from the package root:

```bash
.venv/bin/ruff format scripts/compare_forest_weight_bao_errors.py
.venv/bin/ruff check scripts/compare_forest_weight_bao_errors.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B scripts/compare_forest_weight_bao_errors.py
.venv/bin/ruff format --check scripts/compare_forest_weight_bao_errors.py
.venv/bin/ruff check scripts/compare_forest_weight_bao_errors.py
```

The single calculation/control invocation passed in **0.3724 s** internally on
login32, existing Python 3.13.15 / NumPy 2.5.3, one process and one numerical
thread. A 30-second alarm retains partial results on a domain/check/timeout
exception; the outer timeout also caps the invocation. No stop occurred.
Two preceding read-only JSON metadata inspections used the same interpreter,
thread settings and `timeout 30`; they performed no matrix calculation.
Ruff lint and format checks pass. Exact execution commands are also saved in
`commands.txt` beside the summary.

Retained accuracy policies include density `floor_negative` (floor 1e-20),
`legacy_floor` magnitude support, `legacy_spline_extension` redshift support,
SNR `legacy_floor_clamp` (floor 1e-10, sentinel variance 1e20, four exposures),
and `legacy_first_spacing; unknown physical cells`. W07 neither changes nor
physically validates these assumptions. W06 coefficient evidence, historical
Step 12 Fisher sums and new W07 arithmetic remain distinct.

Order-64/t=6 remains excluded and unavailable because of the historical LBG
coefficient failure; its missing matrix was not reconstructed. That absence
is not a failure of the QSO-only observable. No other spectrum/bin, weighting
trajectory, model evaluation, mode creation, mode sum, forecast, full pytest,
wheel, environment, generalized validator, Slurm action, delegation, commit or
push was performed. No production or planning file was edited; all historical
dirty/untracked work was preserved. `IMPLEMENTATION_STEP.md` retains SHA-256
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.
The home workspace resolves to the requested CFS checkout. The plan's proposed
status is superseded for this execution by the user's explicit W07 r1 approval;
its text was preserved. W01 repair, package Step 13 and scientific acceptance
remain unchanged. Stop for independent/user review.
