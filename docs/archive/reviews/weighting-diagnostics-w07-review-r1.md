# Scientific conclusion and impact of proposed changes

**W07 revision 1 passes independent scientific review. No scientifically
consequential correction is needed.** For the saved accuracy-profile QSO-forest
auto-spectrum in bin 0, the fixed inverse-variance reference has smaller
marginalized radial and transverse BAO errors than either cumulative choice.
Relative to that reference, t=3 increases the errors by **7.63057% and
10.32563%**, and t=6 by **8.61956% and 11.70862%**. Increasing t=3 to t=6 raises
them by **0.91888% and 1.25355%**. Both finite-count differences exceed 0.1%.
The reference order-32/64 refinement passes the assigned empirical criterion
to roundoff. This supports the reported finite comparison, not a general
integrated-BAO optimum or production adoption.

No concrete feedback accompanied the request: its feedback field remained a
placeholder. This review assesses the exact revision-1 assignment and handoff.
No repair revision or implementation change is proposed. A pass is not user
acceptance or permission to advance; stop for user review.

## Selection and scientific meaning

All four prescribed input SHA-256 identities match before and after the
independent check. The diagnosis bin identifies the matching `bin-0.npz` and
hash. Report context/settings identify accuracy bin 0, z=[2.0,2.235], and
parameter order `ap_0`, `at_0`. Field `lya(qso)` is a forest with physical tracer
`lya` and background `qso`; its unique auto-pair resolves to index 0. The
primary NPZ selected-pair ordering agrees with the report. Selection was by
identity, not an assumption about array position.

The four available rows are order32/iterations0, order32/iterations3,
order32/iterations6 and order64/iterations0. Their named `pair_fisher_key`
arrays have shape (15,2,2); only the selected matrix was calculated. The
fixed-initial flag is true only for the references. The selected order-32/t=6
matrix matches the saved primary to **3.5591e-16** maximum elementwise relative
difference; primary controls specify order 32/count 6.

Source inspection of `scripts/diagnose_desi2_weights.py::run` establishes that
`iterations=0` selects `w0=1/(1+pixel*variance/alias)`, equivalent to
nu=B_star/(B_star+l_p v). It does not mean zero weights or zero noise. The
reference is prepared at the historical weighting mode
k_parallel=0.00035 s/km, then held fixed as a function of magnitude. The
mode-dependent noise uses the historical P1D and response with these fixed
coefficients; it is not reoptimized at each Fourier mode.

`fishhighz/validation/numerics.py::wick`, `contract` and `summaries` establish
that `pair_fisher` is computed for each spectrum selected independently.
The `single` contraction uses the reciprocal covariance diagonal, separately
from the joint solve. For the QSO-forest auto-spectrum this gives

```text
T00 = P00_observed + N0
C_auto = 2 T00^2 / modes
F_ab = sum_modes J_auto,a J_auto,b / C_auto.
```

No joint inverse-covariance block or information from another tracer enters
this matrix. In particular, the historical controller's LBG-noise changes do
not enter this selected auto-spectrum. The saved observed signal, J, Fourier
grid, modes, geometry, response and model are held fixed while noise changes.
This is a source-level identification, not a new execution of the contraction.

The [Step 12 r5 review](step-12-review-r5.md), particularly its scientific
assessment and saved-array checks, supplies the existing independent mode-sum
evidence. Its weak-spectrum validator finding did not identify an inconsistency
in these submitted arrays. W07 does not reopen or certify the full forecast.
The [W06 review](weighting-diagnostics-w06-review-r1.md) supplies the distinct
coefficient and conditional single-mode minimum evidence. The observed error
ranking agrees with W06's Q ranking here, but no Q ratio was used to rescale
Fisher matrices or BAO errors. A minimum at one weighting mode alone cannot
establish that ranking for integrated radial and transverse information.

## Independent arithmetic

The independent script imports NumPy and the standard library only; it neither
imports nor executes the implementation script or any forecast/controller.
It checks the four matrices and F=[[4,1],[1,9]] analytic control. All four
matrices are exactly symmetric, finite, have positive diagonals, saved rank 2,
and both parameters constrained. Their normalized minimum eigenvalue is at
least **0.5543488026**, far above the `64*eps*2` spectral threshold (about
4.11e-14), following `fishhighz/_information.py`. No averaging, regularization,
prior or pseudoinverse is needed.

Direct float64 solves agree with independently formed 80-digit Decimal
inverses of the exact saved binary entries. For a=F00, b=F01, d=F11, the check
uses Delta=a*d-b*b>0, covariance=[[d,-b],[-b,a]]/Delta,
sigma=(sqrt(d/Delta),sqrt(a/Delta)), and rho=-b/sqrt(a*d).
Thus errors are marginalized over the other BAO parameter. The analytic
control reproduces determinant 35, covariance=[[9,-1],[-1,4]]/35,
sigma=(3/sqrt(35),2/sqrt(35)), and rho=-1/6, distinct from conditional errors.

| Order / weights | sigma_parallel | sigma_perp | correlation |
| --- | ---: | ---: | ---: |
| 32 / reference | 0.02961952211255668 | 0.02797532316839489 | -0.4440131917488977 |
| 32 / t=3 | 0.03187966090130967 | 0.03086395283440686 | -0.4454715303477305 |
| 32 / t=6 | 0.03217259526976886 | 0.03125084883353349 | -0.4456511973987090 |
| 64 / reference | 0.02961952211255663 | 0.02797532316839485 | -0.4440131917488976 |

These dimensionless errors refer to alpha_parallel and alpha_perp; at fiducial
values of one their percentages are 100 sigma. All four reported matrices
agree exactly with the saved entries. Covariances, errors and correlations
agree with both stored pair summaries and the implementation output; the
maximum relative discrepancy across the checks is **7.368e-16**, below
rtol=5e-12 with atol=0 and exact handling of reference zeros.

| Numerator / denominator | Radial fractional change | Transverse fractional change |
| --- | ---: | ---: |
| t=3 / reference, order 32 | +0.07630571418958998 | +0.10325634662463493 |
| t=6 / reference, order 32 | +0.08619562285678649 | +0.11708624938564127 |
| t=6 / t=3, order 32 | +0.00918875421435737 | +0.01253552975545378 |
| reference 64 / reference 32 | -1.4876851e-15 | -1.6892669e-15 |

The table uses independent Decimal ratios; float64 rounding explains the tiny
last-row differences from the handoff. All contrasts agree within the assigned
absolute tolerance 5e-12. The independently evaluated Fisher Frobenius
relative discrepancy is **3.0790395131128373e-15**. It and both absolute error
changes pass 1e-3. This establishes only the specified reference order-32/64
stability, not convergence in Fourier quadrature, redshift or arbitrary
magnitude orders, or a continuum limit of the cumulative rule.

## Evidence, execution and disposition

Independent check.py (local-only path: `../.validation/forest-weight-diagnostics/w07-review-r1-independent/check.py`)
and result.json (local-only path: `../.validation/forest-weight-diagnostics/w07-review-r1-independent/result.json`)
retain the exact matrices, mapping, input hashes, Decimal results and checks.
Commands from the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B .validation/forest-weight-diagnostics/w07-review-r1-independent/check.py
.venv/bin/ruff check scripts/compare_forest_weight_bao_errors.py .validation/forest-weight-diagnostics/w07-review-r1-independent/check.py
.venv/bin/ruff format --check scripts/compare_forest_weight_bao_errors.py .validation/forest-weight-diagnostics/w07-review-r1-independent/check.py
```

Two bounded independent arithmetic invocations passed in under one second
each, with one process and all three numerical thread limits set to one. The
second followed simplification of the review control assertion. Ruff initially
flagged import ordering in the new review script; it was corrected, and final
lint/format checks pass for both scripts. No implementation repair was made
and no implementation forecast or mode contraction was rerun.

Inherited policies remain density `floor_negative` (1e-20), `legacy_floor`
magnitude support, `legacy_spline_extension` redshift support, SNR
`legacy_floor_clamp` (1e-10 floor, 1e20 sentinel variance, four exposures), and
`legacy_first_spacing; unknown physical cells`. Neither this comparison nor
its numerical pass physically validates these input assumptions. The excluded
order-64/t=6 matrix remains unavailable because of the historical LBG coefficient
failure; it was not recovered and its absence is not a failure of this QSO-only
observable.

Archived the original W07 instructions before status-only synchronization of
the current diagnostic assignment, roadmap, prompts, package guidance and
parallel-diagnostic handovers. The scientific assignment remains revision 1.
Package Step 13, `IMPLEMENTATION_STEP.md`, all historical evidence and forecast
acceptance status remain unchanged. No new inputs, trajectories, modes, joint
calculations, mode sums, forecasts, broad suites, Slurm actions, agent dispatch,
commits or pushes were performed. User approval, acceptance and progression
remain with the user.
