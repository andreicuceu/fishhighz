# W06 revision 1: saved forest-weight coefficient refinement

## Scientific answer

**The inverse-variance reference is empirically stable to roundoff over orders
16/32/64 for accuracy-profile QSO-forest bin 0. Saved cumulative t=3 also meets
the prescribed 0.1% coefficient criterion. Saved t=6 does not: its 64/16
P_pixel change is +0.129453%.** Its adjacent-order changes individually remain
below 0.1%; this failure does not establish nonconvergence. Reference n_eff is
20.5598389416803 deg^-2. All cumulative Q values exceed the nonnegative-measure
minimum, by 14.6692–14.7035% at t=3 and 16.6230–16.7146% at t=6.
Increasing t from 3 to 6 raises Q by 1.7039–1.7533% at fixed order.

Resolving the fixed reference function therefore does not remove the observed
finite-count dependence of the cumulative rule. This is an empirical statement
about these three finite orders, two saved counts and one weighting mode.
It establishes no continuum or joint limit, physical validation of the inherited
positive input model, multi-mode optimum, or forecast improvement. No production
estimator or density prescription is selected. **Implementation complete within
W06 scope; stop for independent/user review.**

## Relation to previous evidence

- [Step 12 r5 handoff](step-12-r5.md#r3-fixed-weight-quadrature-versus-the-cumulative-update)
  already reported fixed-seed magnitude stability. Its
  [independent review](step-12-review-r5.md) found order-32/64 individual-spectrum
  error changes at most 1.09e-14. Those historical forecasts were not rerun here.
- [W05](weighting-diagnostics-w05-r1.md) and its
  [review](weighting-diagnostics-w05-review-r1.md) identify the legacy seed with
  nu=B/(B+l_p v), including amplitude, and establish its nonnegative-measure
  minimum. W05's signed compatibility inputs are not substituted here; their
  comparison cannot isolate quadrature effects relative to this accuracy input.
- W06 reconstructs nine coefficient sets from immutable accuracy arrays,
  connects the fixed seed to n_eff and Q, and separates finite-order contrasts
  from changes between the two saved cumulative counts. No trajectory is evolved.

## Definition, assumptions and units

Use the saved masses r directly: they already include quadrature and
(1+z_source)/c. They receive no second quadrature or velocity factor.

```text
D_i = B + l_p v_i,      nu_i = B/D_i
I1 = sum r_i w_i,       I2 = sum r_i w_i^2,       I3 = sum r_i w_i^2 v_i
A = I2/(L I1^2),        P_pixel = l_p I3/(L I1^2), Q = A B + P_pixel
K = sum r_i/D_i,        n_eff = L B K
Q[nu] = 1/(L K) = B/n_eff.
```

For r>=0, Cauchy–Schwarz implies I1² <= K sum(r D w²), hence
Q[w]>=1/(L K), with equality for weights proportional to 1/D on positive
support. This is conditional on common L, l_p, B, independent sightline
instrumental noise and the stated diagonal-covariance approximation. P0 is
held fixed but absent from D and Q. B already contains the historical response;
no extra smoothing or noise is applied. l_p v is individual 1D instrumental
power, whereas P_pixel is effective 3D pixel-noise power.

A is in deg²; A B, P_pixel, Q and P0 are in deg² km/s; n_eff is in deg^-2;
L, l_p and B are in km/s; r is in deg^-2 (km/s)^-1. Weights and v are
dimensionless. The float reconstruction preserves rw=r*w, rww=rw*w,
rwwv=rww*v, final cumulative sums and sequential divisions from `_integrals`.
All NumPy arithmetic uses raised floating-point exceptions, including underflow.

## Coefficients

| Order | Vector | A [deg²] | A B [deg² km/s] | P_pixel [deg² km/s] | Q [deg² km/s] | n_eff [deg^-2] |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 16 | nu | 0.0294198456633005 | 0.481213030187562 | 0.314354922180623 | 0.795567952368185 | 20.5598389416803 |
| 16 | t=3 | 0.0221183694718741 | 0.361784617029651 | 0.550486516042509 | 0.912271133072160 | — |
| 16 | t=6 | 0.0223405915716006 | 0.365419448129986 | 0.562395725305585 | 0.927815173435570 | — |
| 32 | nu | 0.0294198456633004 | 0.481213030187562 | 0.314354922180624 | 0.795567952368186 | 20.5598389416803 |
| 32 | t=3 | 0.0221174079156969 | 0.361768889096642 | 0.550682420406181 | 0.912451309502824 | — |
| 32 | t=6 | 0.0223406314618806 | 0.365420100605283 | 0.562875195600552 | 0.928295296205835 | — |
| 64 | nu | 0.0294198456633004 | 0.481213030187561 | 0.314354922180623 | 0.795567952368184 | 20.5598389416803 |
| 64 | t=3 | 0.0221169109895518 | 0.361760761000431 | 0.550783458732569 | 0.912544219733001 | — |
| 64 | t=6 | 0.0223406459783711 | 0.365420338047874 | 0.563123761917365 | 0.928544099965239 | — |

### Refinement at fixed vector definition

Entries are signed fractional changes (higher order / lower order minus one),
not percentages. Each vector is classified from all nine changes in this table:
all absolute changes must be <=1e-3. The operands are independently checked.

| Vector | Orders | ΔA/A | ΔP_pixel/P_pixel | ΔQ/Q |
| --- | --- | ---: | ---: | ---: |
| nu | 32/16 | -7.772e-16 | +3.553e-15 | +8.882e-16 |
| nu | 64/32 | -2.109e-15 | -3.331e-15 | -2.665e-15 |
| nu | 64/16 | -2.998e-15 | +2.220e-16 | -1.665e-15 |
| t=3 | 32/16 | -4.34731945e-5 | +3.55874954e-4 | +1.97503159e-4 |
| t=3 | 64/32 | -2.24676484e-5 | +1.83478395e-4 | +1.01824864e-4 |
| t=3 | 64/16 | -6.59398661e-5 | +5.39418644e-4 | +2.99348133e-4 |
| t=6 | 32/16 | +1.78555165e-6 | +8.52549679e-4 | +5.17476739e-4 |
| t=6 | 64/32 | +6.49779774e-7 | +4.41601120e-4 | +2.68022213e-4 |
| t=6 | 64/16 | +2.43533258e-6 | +1.29452729e-3 | +7.85637647e-4 |

nu and t=3 pass; t=6 fails only the last P_pixel comparison. The t=3 pass is
restricted to this selection and does not resolve historical failures at other
orders/counts. No convergence theorem follows from either pass.

### Within-order comparisons

All entries below are percentages. The last three columns compare t=6 to t=3;
the first two compare each cumulative Q with the reference at the same order.

| Order | 100(Q3/Qnu−1) | 100(Q6/Qnu−1) | 100(A6/A3−1) | 100(P6/P3−1) | 100(Q6/Q3−1) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 16 | 14.669166 | 16.622995 | 1.004695 | 2.163397 | 1.703884 |
| 32 | 14.691813 | 16.683345 | 1.009266 | 2.214121 | 1.736420 |
| 64 | 14.703492 | 16.714618 | 1.011601 | 2.240500 | 1.753327 |

Both A and P_pixel increase between these saved counts. All six cumulative
vectors satisfy the nonnegative-measure lower bound in float and Decimal.

## Input identity and retained policies

Selected only `lya(qso)`, accuracy bin 0, z=[2.0,2.235], from
`.validation/step12-r5-20260914T191855Z/`. Verified the following exact SHA-256
identities before and after the comparison:

| Relative input | SHA-256 |
| --- | --- |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |

The diagnosis bin's NPZ name/hash agree. Report context and settings identify the
accuracy profile, bin bounds and QSO forest as field 0. Selected metadata rows
at each order agree on that identity, P0/B and available t=3/6 checkpoints.
The common partition exactly matches the report: 162 intervals on [16.1,26.75].
Orders 16/32/64 mean nodes per interval (2592/5184/10368 total), not total
node counts. These Gaussian nodes are not nested. Each interval has the expected
number of ordered interior nodes and positive quadrature weights summing to
its width; total quadrature sums to 10.65 mag. All masses are strictly positive,
with total about 0.001439259272800875; variances are finite and nonnegative.

Order-32 magnitudes, quadrature and variance independently match the accuracy
report. Its masses match `density*((1+z_source)/299792.458)*quadrature` with
rtol=5e-12, atol=0, including exact reference-zero handling.

| Common scalar | Value |
| --- | ---: |
| L | 44147.09373976799 km/s |
| z_source | 2.3830099582078716 |
| z_eval | 2.1152848986890427 |
| pixel width | 0.8 Angstrom |
| l_p=c*0.8/[1215.67*(1+z_eval)] | 63.328211161339375 km/s |
| P0 | 1.7182791088399005 deg² km/s |
| B | 16.356748967852234 km/s |

P0 and B agree exactly across orders. L, redshifts and pixel width use the
single shared report mapping inspected in `_field_inputs`; the historical
controller uses the same settings and policies for every order. Policies are
shared source/report metadata, not separately recorded per-order policy objects.
No order-dependent replacement is introduced here. D is positive, from
16.821037167341775 to 6.332821116133938e21 km/s; K is about
2.8472172163562443e-5 deg^-2 (km/s)^-2.

Retained density policies: `floor_negative`, density floor 1e-20,
`magnitude_domain=legacy_floor`, `redshift=legacy_spline_extension`, with the
negative-density treatment explicitly marked a reference extension. The saved
normalization measure is 3789.61349941401, using the recorded legacy reduction.
Retained SNR policy: `legacy_floor_clamp`, SNR floor 1e-10, sentinel variance
1e20, exposure count 4, pixel width 0.8 Angstrom. Source widths retain
`legacy_first_spacing; unknown physical cells`. These choices were already
applied in Step 12. W06 neither reapplies nor physically validates them.
No response change, raw interpolation or provider call occurs.

## Actual checks and reproducibility

New script: [compare_forest_weight_refinement.py](../scripts/compare_forest_weight_refinement.py).
New output: [summary.json](../.validation/forest-weight-diagnostics/w06-r1-20260915T234041Z/summary.json).
It contains all float moments and coefficients, independent 80-digit values,
direct-numerator Q, contrasts, domain summaries, retained scalar policy settings,
input/source hashes and elapsed time. Historical arrays are referenced, not copied.

The exact rational control uses q=(1/2,2), rho_v=(2,1/2), hence r=(1,1),
v=(1,4), L=l_p=B=1. It gives nu=(1/2,1/5), A=29/49, P_pixel=41/49,
Q=10/7 and n_eff=7/10. This explicitly forms the nonuniform measure once;
the float helper also matches these rational moments and coefficients.

All nine real-input sets were independently computed with scalar 80-digit
Decimal sums of exact binary input values. The reference weights were formed
independently in Decimal from B, l_p and v; cumulative weights were read directly
from the saved checkpoints. The scalar calculation does not call the float
moment helper. Direct Q=sum(r D w²)/(L I1²) agrees with A B+P_pixel for every
vector, and reference Q agrees with both 1/(L K) and B/n_eff. Historical
fixed_coefficients and available t=3/6 A/P_pixel entries all match.

Maximum relative discrepancy across moment, coefficient, mapping and n_eff
comparisons: **1.059e-14**, below 5e-12 (zero absolute tolerance, exact reference
zeros). All refinement, t=6/t=3 and cumulative/reference contrasts agree with
independent Decimal ratios to absolute **1.525e-14**, below 5e-12. No domain,
nonfinite, underflow, missing-checkpoint or timeout stop occurred. No tolerance
was relaxed, and the single coefficient-check invocation passed.

Commands from the package root:

```bash
.venv/bin/ruff format scripts/compare_forest_weight_refinement.py
.venv/bin/ruff check scripts/compare_forest_weight_refinement.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -B scripts/compare_forest_weight_refinement.py
.venv/bin/ruff format --check scripts/compare_forest_weight_refinement.py
.venv/bin/ruff check scripts/compare_forest_weight_refinement.py
```

The script enforces a 30-second alarm and writes partial results on an exception.
The numerical comparison passed in **0.764 s** internally on login32, existing
Python 3.13.15 / NumPy 2.5.3, one process with all three numerical thread limits
set to one. Two earlier read-only metadata inspections used the same interpreter
and thread settings under `timeout 30`; they computed no coefficient sets.
Ruff lint and format checks pass.

The script SHA-256 is
`d86e110951f18042ddbf6a2cd004f65dbda4d4d2c5636a78817da0947d2f1057`.
The summary also records hashes of inspected coefficient/mapping sources, W05
and Step 12 reports, the diagnostic plan and assignment. All recorded sources
were unchanged across execution. `IMPLEMENTATION_STEP.md` retains SHA-256
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.
The supplied home workspace resolves to the requested CFS checkout.

Only the new standalone script, this report and the new numerical-output
directory were added. Existing dirty/untracked work and planning/production files
were preserved. No full pytest, wheel/environment work, optimizer, generalized
validator, new trajectory, other order/bin/population, Fisher sum, forecast,
Slurm action, delegation, commit or push was performed. W01 repair and package
Step 13 remain unassigned. Stop for independent/user review.
