# Forest-weighting diagnostic W05 revision 1

## Scientific answer

**The explicit inverse-variance reference is already the legacy seed. The
original prefix t=3 weights and both historically stable W04 t=24 vectors have
different shapes and higher combined noise Q at the saved weighting mode.**
Relative to the reference, Q is higher by **14.3618%** (prefix), **12.8358%**
(fixed P0), and **19.1309%** (refreshed aliasing). Refreshed Q exceeds fixed Q
by **5.57897%**, despite its smaller aliasing coefficient A.

For a nonnegative measure, the derivation below establishes the reference's
minimum of Q under the specified covariance approximation. The saved grid has
nine negative densities: its ordering is **algebraic only**, and establishes
neither physical optimality nor a forecast improvement. A recurrence may
stabilize without attaining the positive-measure noise minimum. The observed
signed-grid ordering alone cannot test that physical minimum.

Status: W05 r1 implemented, awaiting independent/user review. The user's explicit
approval supersedes the planning documents' proposed status; those files remain
unchanged. No production prescription is selected.

## Reference, assumptions, and derivation

[McQuinn & White (2011), section 2.2, equations (11)–(13)](https://arxiv.org/pdf/1102.1752)
gives sightline weights proportional to the inverse sum of line-of-sight forest
and instrumental powers. Equation (11)'s B(k_parallel) absorbs the arbitrary
normalization, fixed by their mean-weight convention; it is **not our B**.
Their variance argument neglects subdominant off-diagonal covariance and uses
factorized sightline weights (also stated with equation 2). Equations (12)–(13)
express the additive power as P_los/n_eff with nu=P_los/(P_los+P_N,n).
Source inspected: arXiv:1102.1752v3, 31 December 2014, printed page 4; the PDF
has a later printed title-page date. No equivalence of the ME07/FR14 recurrence
to these weights is inferred.

The following derivation is in the assigned FishHighz units. At this one mode,
all forests share L, l_p and B; instrumental noise is uncorrelated between
sightlines. The supplied P0 and B already contain the instrumental response.
The dimensionless saved pixel variance is v_i, so **l_p v_i is individual 1D
instrumental power**, whereas **P_pixel is effective 3D pixel-noise power**.
No extra smoothing, deconvolution, or coordinate conversion is applied.

Put r_i=rho_i dm and D_i=B+l_p v_i. Then

```text
nu_i = B/D_i = (B/l_p)/(B/l_p+v_i),
I1 = sum r_i w_i,
I2 = sum r_i w_i^2,
I3 = sum r_i w_i^2 v_i,
A = I2/(L I1^2),
P_pixel = l_p I3/(L I1^2),
Q[w] = A B + P_pixel = sum r_i D_i w_i^2 / (L I1^2).
```

The seed equivalence includes amplitude, not just relative shape. This reference
has no iteration count and no dependence on P0. It is not inserted into an
update. P0 is common to all four vectors and is absent from D; aliasing is
included once through A B.

For r_i>=0 with at least one positive mass, B,L,l_p>0, v_i>=0, finite real
weights and I1!=0, let K=sum r_i/D_i>0. Cauchy–Schwarz gives

```text
I1^2 = [sum (sqrt(r_i D_i) w_i) sqrt(r_i/D_i)]^2
     <= (sum r_i D_i w_i^2) K,
Q[w] >= 1/(L K).
```

Equality requires w_i=c/D_i on positive support, with c!=0; weights on zero
measure are immaterial. In particular, nu chooses c=B. Since I1[nu]=B K,

```text
n_eff = L sum r_i nu_i = L B K,
Q[nu] = 1/(L K) = B/n_eff.
```

Equivalently, writing u_i=w_i/I1,

```text
Q[w] - 1/(L K) = (1/L) sum r_i D_i [u_i-1/(K D_i)]^2.
```

The last expression is nonnegative only for a nonnegative measure. The
identities remain algebraically valid for signed r whenever their denominators
are defined, but the inequality and physical effective-density interpretation
do not follow. We therefore call the real-input L B K a *formal n_eff*.

Under the diagonal Gaussian auto-power approximation, variance is proportional
to (P0+Q)^2 at fixed mode count. For physical nonnegative power and measure,
minimizing Q minimizes this variance. This is a single-mode statement, not a
multi-mode BAO or joint-forecast optimum. Neglected covariance terms need not
share this optimum.

For any real c!=0, w→cw sends I1→c I1 and I2,I3→c²(I2,I3), leaving A,
P_pixel and Q unchanged. This is an output normalization identity, not permission
to normalize weights between nonlinear recurrence steps.

| Quantity | Units |
| --- | --- |
| rho | deg^-2 (km/s)^-1 mag^-1 |
| r, I1, I2, I3 | deg^-2 (km/s)^-1 |
| w, nu, v | dimensionless |
| L, l_p, B, D | km/s |
| K | deg^-2 (km/s)^-2 |
| A | deg^2 |
| n_eff | deg^-2 |
| A B, P_pixel, Q, P0 | deg^2 km/s |

## New four-vector measurements

Only saved-vector moments were computed; no trajectory was evolved. The
reference is new W05 output, prefix weights come from the authoritative row,
and fixed/refreshed t=24 weights are W04 checkpoints. The following values
use historical float64 product order and final `cumsum` elements.

| Vector | I1 | I2 | I3 |
| --- | ---: | ---: | ---: |
| nu | 4.6530242201031624e-4 | 2.8149089096121373e-4 | 4.740612787841867e-5 |
| Prefix t=3 | 5.213095963873360e-4 | 2.653738328759833e-4 | 1.0382411040095769e-4 |
| Fixed P0 t=24 | 7.116830099799945e-4 | 4.877411676323469e-4 | 1.909804326682728e-4 |
| Refreshed t=24 | 7.701606407013814e-4 | 5.421287033890218e-4 | 2.518455402224129e-4 |

| Vector | A [deg²] | A B [deg² km/s] | P_pixel [deg² km/s] | Q [deg² km/s] |
| --- | ---: | ---: | ---: | ---: |
| nu | 0.029449686440879 | 0.481006065094750 | 0.314093507491696 | 0.795099572586446 |
| Prefix t=3 | 0.022118395206768 | 0.361263005837993 | 0.548026906249744 | 0.909289912087738 |
| Fixed P0 t=24 | 0.021812385938153 | 0.356264911393947 | 0.540892274852048 | 0.897157186245994 |
| Refreshed t=24 | 0.020702687061427 | 0.338140036237620 | 0.609069323591706 | 0.947209359829326 |

Signed fractional changes use nu as denominator; shape differences are
max_i |w_i/max(abs(w))-nu_i/max(abs(nu))|.

| Vector | ΔA/A_nu | ΔP_pixel/P_pixel,nu | ΔQ/Q_nu | Maximum shape difference |
| --- | ---: | ---: | ---: | ---: |
| nu | 0 | 0 | 0 | 0 |
| Prefix t=3 | -24.894293% | +74.478903% | +14.361766% | 1.0000000000170253 |
| Fixed P0 t=24 | -25.933385% | +72.207404% | +12.835828% | 0.28584207033001047 |
| Refreshed t=24 | -29.701503% | +93.913376% | +19.130910% | 0.3436826871785209 |

Aliasing A B dominates for nu; P_pixel dominates for the other three vectors.
Refreshed-minus-fixed fractional Q is +0.05578974827451. Its reduction in A B
is outweighed by increased P_pixel. The prefix vector contains three negative
weights; the other three have none. Shapes retain these signs.

All coefficients and moments are finite and positive. D ranges from
16.797447129420597 to 6.332980433473594e21 km/s; K=2.8488227951316808e-5
and formal n_eff=20.542266418926324 deg^-2. For nu, the moment Q, 1/(L K),
and B/n_eff agree. No domain stop was triggered.

Signed cancellation does not limit this ordering numerically: the ratios
sum(abs(terms))/abs(sum(terms)) lie between 1.00003136 and 1.00496916 over
the twelve moments, and equal 1.00066092 for K. A B and P_pixel are positive,
so their final addition has no cancellation. These numerical facts do not
make the measure nonnegative.

## Actual independent checks

The exact control has r=(1,1), v=(1,4), L=l_p=B=1. Rational arithmetic gives
nu=(1/2,1/5), K=n_eff=7/10, I1=7/10, I2=29/100, I3=41/100,
A=29/49, P_pixel=41/49, Q=10/7. Uniform weights give Q=7/4>10/7.
Direct-numerator and effective-density expressions agree exactly. Multiplying
nu by 3 preserves A, P_pixel and Q exactly. The float moment helper also
matches all three rational cases.

For the saved inputs, independent scalar sums use 80-digit Decimal arithmetic
from exact binary-float values of rho, dm, v, L, l_p, B and each saved weight.
They do not call the float moment helper. A second reference calculation
computes nu directly in Decimal. Direct sum(r D w²)/(L I1²) is checked against
the decomposition A B+P_pixel for every vector; the reference also passes both
effective-density identities. Historical prefix A/P_pixel match the row's
last `_aliasing_weights`/`_effective_noise_power` elements. Both W04 pairs and
all their moments match its NPZ and summary.

Nonzero comparisons require rtol=5e-12, atol=0, with exact reference-zero
handling and finite operands. All dimensionless fractional contrasts agree
with Decimal to absolute tolerance 5e-12. Maximum relative coefficient/moment
discrepancy across the four vectors is 8.38e-16. The two float seed expressions
differ by at most 2.58e-16 relatively; direct Decimal nu agrees within 1.64e-16.
No tolerance was relaxed and no failed numerical invocation occurred.

## Input/source identity and reproducibility

Authoritative input:
`.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json`,
row `settings.pair_inputs["lya(qso)_lya(qso)"]`, compatibility bin 0,
z=[2.0,2.235]. All 107 magnitudes, densities and variances match W04 exactly;
endpoints are 16.1 and 26.75 mag. Negative-density indices are
0,2,3,76,77,80,81,84,85. The rectangular rule includes both endpoints.

| Held scalar | Exact saved float |
| --- | ---: |
| dm | 0.100471698113207 mag |
| L | 44148.20436604321 km/s |
| l_p | 63.32980433473595 km/s |
| P0 | 1.7189455896517238 deg² km/s |
| B | 16.333147249645222 km/s |

These scalars equal W04's saved values exactly. Response ownership follows
[W02's source audit](weighting-diagnostics-w02-r1.md) and the inspected live
legacy definitions: `weights.py:104–133` samples P0/B and initializes the seed;
`:184–253` defines moment products; `covariance.py:287–294` extracts coefficients
and `:352–363` adds signal, aliasing and pixel noise once. The live FishHighz
kernel was inspected but never imported: it has strict positive-support and
range behavior, so it is not used to reinterpret this signed compatibility row.

New outputs:
summary.json (local-only path: `../.validation/forest-weight-diagnostics/w05-r1-20260915T231931Z/summary.json`)
(9.3 KiB) records the full table, 80-digit values, checks, scalars, timings and
source hashes;
arrays.npz (local-only path: `../.validation/forest-weight-diagnostics/w05-r1-20260915T231931Z/arrays.npz`)
(12 KiB) retains the three input arrays, four weights and display shapes.
The only new Python is
compare_inverse_variance_forest_weights.py (local-only path: `../scripts/compare_inverse_variance_forest_weights.py`).

Decisive SHA-256 identities (all source/input entries checked unchanged across
the numerical invocation):

```text
3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6  authoritative report
7c250829e919388cb6d32eff4b5d22a599ee61b62b5719e5faf624a358f2c895  W04 arrays.npz
cc7bfc622a336f3c47d71085a93c9194d1293a696c1019853a358783435e0bac  W04 summary.json
46752cca5f6dea7a3605bcf7510131164d5e867c467d95a867977dbf2e8055f6  new W05 script
1adfa3354fc887d5008c3a37596b494ba8585c774a9866aaa366bdf8a462b788  W04 script
10f7da1c68eba0cb25a7c77fa1c92b48ae6bc4323eddea680997f29303328431  legacy weights.py
f251921dd700a3cacd1f618426093237a46983937167bc5570c9202f9e596e1d  legacy covariance.py
3867c9eee6d43f22fca878b16c9f38bafc65a55061355d79828fcb135cc5e66f  FishHighz kernels/weights.py
aa8fe72ce592dae51e57dad2a6a07830fd412ac35964ed71048f4c3c88ac080c  WEIGHTING_DIAGNOSTIC_STEP.md
3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310  IMPLEMENTATION_STEP.md
```

The checkout resolves to the requested CFS package through the supplied home
workspace path. HEAD is `0d69786a06d5d676564a51fad14a7156951c2228`; extensive
historical dirty/untracked work was already present and is preserved. Source
hashes, rather than HEAD alone, identify the inspected implementation.

Commands run from the package root:

```bash
.venv/bin/ruff format scripts/compare_inverse_variance_forest_weights.py
.venv/bin/ruff check scripts/compare_inverse_variance_forest_weights.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -B scripts/compare_inverse_variance_forest_weights.py
.venv/bin/ruff format --check scripts/compare_inverse_variance_forest_weights.py
.venv/bin/ruff check scripts/compare_inverse_variance_forest_weights.py
git diff --check
```

The single numerical invocation passed in 0.2480 s internally (2.53 s command
wall time), using the existing Python 3.13.15 / NumPy 2.5.3 interpreter on
login32, one process and one numerical thread. Its 30-second alarm retains
partial arrays/summary on an exception. Ruff lint/format and tracked
`git diff --check` passed. The latter does not check untracked additions.

## Historical evidence and review boundary

[W04's report](weighting-diagnostics-w04-r1.md) and
[independent review](weighting-diagnostics-w04-review-r1.md) establish bounded
0.1% coefficient stability over t=6/12/24 for both full-sample branches. W05
neither repeats those trajectories nor strengthens that stability claim.
Its new evidence is the reference derivation/control and the four-vector
reference/noise comparison. Historical refresh-convention uncertainty remains.

Only this report, the standalone script and the new output directory were added.
No planning or production files were edited. No flooring, smoothing, mode change,
magnitude refinement, new physical inputs, additional bins/counts, optimizer,
validator framework, pytest suite, wheel/environment, Fisher sum, forecast,
Slurm action, delegation, commit or push was performed. W01 evidence repair and
package Step 13 remain unassigned. No continuum or physical-forecast conclusion
follows. Stop for independent/user review.
