# Scientific conclusion and impact of proposed changes

**W06 revision 1 passes independent scientific review within the exact assigned
scope. No scientifically consequential correction is needed.** The reference
nu and saved cumulative t=3 coefficients satisfy the prescribed 0.1% magnitude
refinement criterion. Saved t=6 fails only the order-64/16 P_pixel comparison,
which changes by +0.129453%. Its adjacent-order changes pass. This finite
comparison establishes neither nonconvergence nor a continuum limit.

At fixed order, increasing the saved count from t=3 to t=6 raises Q by
1.703884–1.753327%. All six cumulative Q values exceed the reference minimum
on the verified nonnegative measure. The conclusion is conditional on the
inherited accuracy input policy and the stated single-mode covariance
approximation; numerical nonnegativity does not establish physical validity
of the density floor or extrapolation. No production prescription follows.

No repair revision or consequential scientific choice is required. **Awaiting
user review; this pass is not acceptance or permission to advance.** Package
Step 13, Step 12 scientific acceptance and deferred W01 findings are unchanged.

## Assignment and evidence inspected

Reviewed [W06 implementation report](weighting-diagnostics-w06-r1.md) against the
[exact archived revision-1 assignment](weighting-diagnostics-w06-instructions-r1.md).
The request supplied no substantive feedback beyond the unfilled placeholder.
Read workspace/package AGENTS.md, design section 0, main roadmap handover and
progress register, diagnostic roadmap, W05 report/review, and the specified
Step 12 r5 fixed-weight evidence. The latter already records fixed-seed
order-32/64 individual-spectrum error changes at most 1.09e-14 in its independent
review. W06 correctly labels those results historical. No new forecast is
needed to establish that previously reviewed result.

Inspected the complete W06 standalone script, W05 moment/scalar helpers,
`diagnose_desi2_weights.py:35–214`,
`diagnose_weight_limit.py:243–305`, the density-to-velocity conversion in
`fishhighz/weights.py:33–49`, and `_integrals` in
`fishhighz/kernels/weights.py:26–40`. The standalone assertions supply the
focused checks; no separate W06 test file or broad suite is necessary.
No implementation helper or production module was imported by the review check.

## Input identity, measure and scalar mapping

All three authoritative input SHA-256 values in the assignment match live
files. All twelve source/document hashes recorded by W06 match before status
synchronization, including the W06 script and package `IMPLEMENTATION_STEP.md`.
The latter retains SHA-256
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.

The diagnosis and accuracy report identify bin 0, z=[2.0,2.235], with field 0
`lya(qso)` (forest, physical Ly-alpha, QSO background). For each selected order,
all repeated field metadata agree, and the t=3/t=6 checkpoints are available.
The common 162-interval partition is identical to the report. Orders 16/32/64
have 2592/5184/10368 ordered interior nodes, positive quadrature and strictly
positive masses. Each interval's quadrature sums to its width, and the total is
10.65 mag on [16.1,26.75]. Variances are finite and nonnegative.

Order-32 magnitudes, quadrature and variance match the accuracy report.
Independently reconstructing its masses as
`density * ((1+z_source)/299792.458) * quadrature` also agrees. The historical
controller explicitly constructs every order's masses by this same operation.
**The saved masses already include both magnitude quadrature and the velocity
density conversion; neither factor is applied twice.** The order-16/64 policy
mapping rests on this shared controller and recorded input identity, not on
new interpolation or independently recorded per-order policy objects.

The scalar mapping gives L=44147.09373976799 km/s,
z_source=2.3830099582078716, z_eval=2.1152848986890427 and
l_p=63.328211161339375 km/s from the 0.8 Angstrom pixel width.
Every order has B=16.356748967852234 km/s and
P0=1.7182791088399005 deg² km/s. Source inspection confirms B includes the
historical response squared at k_parallel=0.00035 s/km. No additional response
belongs in D=B+l_p v; P0 is absent from Q and the reference weight.

The report accurately retains `floor_negative` (1e-20), `legacy_floor` magnitude
support, `legacy_spline_extension`, normalization measure 3789.61349941401,
`legacy_floor_clamp` SNR (floor 1e-10, sentinel variance 1e20, four exposures),
and `legacy_first_spacing; unknown physical cells`. These numerical conventions
are unchanged and have not been physically validated by W06.

## Independent derivation and arithmetic

For r_i>=0, L,B,l_p>0, v_i>=0, finite weights and I1!=0, set
D_i=B+l_p v_i and K=sum(r_i/D_i)>0. Cauchy–Schwarz gives

```text
I1² <= K sum(r_i D_i w_i²),
Q[w] >= 1/(L K).
```

Equality holds for w_i=c/D_i on positive support, c!=0. With c=B,
nu=B/D, n_eff=L B K and Q[nu]=B/n_eff. The report applies this statement only
with common scalars, independent instrumental sightline noise and the stated
diagonal-covariance approximation. Its units and separation of individual 1D
instrumental power l_p v from effective 3D P_pixel are consistent. The signed
W05 grid does not meet the measure assumption and is not used for this test.

Independent exact Fraction arithmetic using q=(1/2,2), rho_v=(2,1/2) gives
r=(1,1), nu=(1/2,1/5), A=29/49, P_pixel=41/49, Q=10/7 and n_eff=7/10.
The nonuniform quadrature enters once.

The review script forms scalar 80-digit Decimal sums from exact saved binary
values, deriving nu independently from B, l_p and v. It reconstructs all nine
I1/I2/I3, A, A B, P_pixel and Q sets, reference n_eff, direct-numerator Q,
Q=A B+P_pixel, Q[nu]=1/(L K)=B/n_eff, and historical fixed/prefix coefficients.
All agree at rtol=5e-12, atol=0 with exact zero handling. The maximum relative
discrepancy is **1.058e-14**. All 42 requested fractional contrasts agree to
absolute **1.525e-14**, below 5e-12. Stability is classified from the independent
operands, not copied from the stored verdict.

| Vector | Largest absolute fixed-count refinement | Criterion |
| --- | ---: | --- |
| nu | roundoff, below 4e-15 in the reported float contrasts | Pass |
| t=3 | 0.0539419% in P_pixel, 64/16 | Pass |
| t=6 | 0.129453% in P_pixel, 64/16 | Fail |

Reference n_eff=20.5598389416803 deg^-2 and Q=0.795567952368185 deg² km/s
are stable to roundoff. All individual refinement contrasts in the implementation
report are supported. Its within-order percentages are independently reproduced:

| Order | Q3/Qnu−1 | Q6/Qnu−1 | A6/A3−1 | P6/P3−1 | Q6/Q3−1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 | 14.669166% | 16.622995% | 1.004695% | 2.163397% | 1.703884% |
| 32 | 14.691813% | 16.683345% | 1.009266% | 2.214121% | 1.736420% |
| 64 | 14.703492% | 16.714618% | 1.011601% | 2.240500% | 1.753327% |

Thus resolving the fixed reference does not remove the measured finite-count
sensitivity of the cumulative rule. Neither the t=3 pass nor the t=6 failure
settles an asymptotic, continuum or joint limit. No cross-profile attribution,
accuracy-profile full-sample-branch comparison or forecast improvement is claimed.

## Execution and disposition

Independent evidence is in
scalar_check.py (local-only path: `../.validation/forest-weight-diagnostics/w06-review-r1-independent/scalar_check.py`)
and result.json (local-only path: `../.validation/forest-weight-diagnostics/w06-review-r1-independent/result.json`).
Commands from the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B .validation/forest-weight-diagnostics/w06-review-r1-independent/scalar_check.py
.venv/bin/ruff check scripts/compare_forest_weight_refinement.py .validation/forest-weight-diagnostics/w06-review-r1-independent/scalar_check.py
.venv/bin/ruff format --check scripts/compare_forest_weight_refinement.py .validation/forest-weight-diagnostics/w06-review-r1-independent/scalar_check.py
```

The single numerical invocation passed in 0.827 s internally, under the
30-second cap and one process with numerical thread limits of one. An earlier
incorrect relative script path failed before script creation or numerical work.
Ruff initially identified two style issues in the new review script; these were
corrected and formatting applied. Final lint and format checks pass for both
scripts. No implementation code was changed or rerun.

For subsequent replay, the review script checks the exact assignment archive.
The original run verified the roadmap hash before synchronization; that
historical hash remains in result.json and is not required of the updated roadmap.

Archived the original W06 instructions before updating their status only.
Synchronized W06's diagnostic roadmap, package guidance and parallel-diagnostic
handover paragraphs. The scientific assignment remains revision 1 with no
repair tasks. Historical evidence, production files and package Step 13 remain
unchanged. No new inputs, trajectories, orders, counts, bins, populations, modes,
Fisher calculations, forecasts, broad tests, Slurm, agent dispatch, commits or
pushes were performed. Stop for user review.
