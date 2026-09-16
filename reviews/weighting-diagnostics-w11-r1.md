# W11 revision 1: fixed-reference-mode sensitivity

**Reference-mode sensitivity is detected at the assigned second alternative.**
For accuracy bin 0 and the `lya(qso)` auto-spectrum at magnitude order 32,
changing `q_star` from 0.00035 to 0.001 s/km changes the marginalized radial
and transverse BAO errors by -0.00147679% and -0.0000622378%, respectively.
Both are below the 0.1% screen, so the conditional 0.003 s/km calculation was
required. At 0.003 s/km the changes are **+0.116214% and +0.220774%**, both
above the screen. The calculation therefore stops at 0.003 s/km.

This is the specified sensitivity finding, not a production defect, an
optimization result, or evidence for selecting another reference convention.
No numerical check failed. The result awaits independent and user review; W12
was not started.

## Identity, intervention and held calculation

Saved metadata identifies profile `accuracy`, bin 0 with bounds [2, 2.235],
field `lya(qso)` (physical tracer `lya`, background `qso`), its unique selected
and required auto-pair at index 0, and parameter order (`ap_0`, `at_0`). The
calculation uses 5,184 magnitude nodes and 16,384 saved Fourier nodes. The four
input identities below were verified before and after execution; paths are
relative to `.validation/step12-r5-20260914T191855Z/`.

| Input | SHA-256 |
| --- | --- |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |

The native order-32 magnitudes, quadrature and variance agree exactly between
the report and diagnosis arrays. With
`rho=dN/(dz dm deg^2)*(1+z_source)/c`, `rho*quadrature` agrees exactly with the
saved masses. Fixed scalars are:

| Quantity | Value |
| --- | ---: |
| `z_eval` | 2.1152848986890427 |
| `z_source` | 2.3830099582078716 |
| `h_fid` | 0.6736 |
| `a_v` [(km/s)/(Mpc/h_fid)] | 102.61326549228777 |
| `d_deg` [(Mpc/h_fid)/degree] | 64.41180749124948 |
| `L` [km/s] | 44147.09373976799 |
| `l_p` [km/s] | 63.328211161339375 |
| Gaussian `sigma_v` [km/s] | 50.92405402826606 |

The W09 constant-H and constant-D_M adapter constructs only the public geometry
context and reproduces `a_v` and `d_deg` exactly. Its artificial integrated
volume is discarded; the calculation uses the saved mode counts.

Across trials, the saved observed signal, saved mean Jacobian, k, mu, mode
counts, geometry conversions, response parameters, magnitude nodes, masses,
variance, quadrature, L and pixel width are unchanged. The mode-dependent
physical P1D and response arrays are evaluated on the saved q nodes and held
fixed. The only intervention is

```text
q_star -> B_star=P1D(q_star)*W(q_star)^2 -> fixed inverse-variance weights.
```

Noise and therefore fiducial covariance are recomputed between trials. Within
each trial, those weights and that covariance remain fixed while contracting
the unchanged saved mean Jacobian. No covariance derivative is used. Constant
`B_star` never replaces physical `P1D(q_n)`, and no mode-dependent weights are
introduced.

The exact source manifest was also checked after calculation. The decisive
entries are:

| Source | SHA-256 |
| --- | --- |
| `WEIGHTING_DIAGNOSTIC_STEP.md` | `84dbdd5d07d48b2cf780600661575aa5b0c3b450c0cdfc6256e20b27ac21c162` |
| `scripts/check_forest_reference_mode.py` | `211c482cec8c0e6f4e1d538603e1017e27761686860fb93c7a800613370850ed` |
| `fishhighz/weights.py` | `5e9110b63d13eded2cb84dbd0e82b958418ed59602d910bfeb1428a27d67f663` |
| `fishhighz/noise.py` | `9cb580b1cb52860d3b65572bb18dfb707795b30a9646ff92778b930b429ef7b6` |
| `fishhighz/response.py` | `ea6c61048ea8a9eae809acd304b66d2b7e97aea5f4817b704fd37255024137b4` |
| `fishhighz/covariance.py` | `0c9cfbe1508aecfd5074342f799599e569dab91b2cf6c1e48db8277c85a9e4f2` |
| `fishhighz/fisher.py` | `459ba376ba1382071939d95c98e67b82ba82c1f9092274c825dc31c2b26186a8` |
| `fishhighz/models/p1d.py` | `907023a259cb0bbe47994b6de5c52f9010106f0c77d31cd40c23beaa2e0ce32b` |

The result artifact contains the full 17-file manifest, including the numerical
kernels and context dependencies.

## Reference points, coefficients and scalar checks

The independently reconstructed P1D floor is
0.0008164315784041575 s/km. The saved physical-mode range is
[1.3686646589112106e-7, 0.004863411498638425] s/km. Thus the baseline lies below
the inherited stationary floor, both alternatives lie above it, and all three
points lie inside the saved mode domain.

The public and independent scalar formulas agree at every reference point. The
table gives public values; `W` is the field response and `B_star` is already
response-smoothed exactly once.

| `q_star` [s/km] | intrinsic P1D [km/s] | W | W^2 | `B_star` [km/s] |
| ---: | ---: | ---: | ---: | ---: |
| 0.00035 | 16.362615788177436 | 0.9998207087918288 | 0.9996414497289949 | 16.356748967852234 |
| 0.001 | 16.295450037915643 | 0.9985373329203274 | 0.9970768052356408 | 16.247815263681932 |
| 0.003 | 13.813255102182207 | 0.9869123551380103 | 0.9739959967240542 | 13.454055171253588 |

For every trial, the public preparation is called once with
`method="inverse_variance", alias=B_star`, without signal, iterations or
auxiliary sampling. Independently, Python scalar arithmetic and `math.fsum`
evaluate

```text
nu_i = B_star/(B_star+l_p*v_i)
I1 = sum r_i nu_i
I2 = sum r_i nu_i^2
I3 = sum r_i nu_i^2 v_i
A = I2/(L I1^2)
P_pixel = l_p I3/(L I1^2).
```

The mode noise is then independently reconstructed with the fixed physical
P1D and response arrays. Pixel noise is checked separately to be unsmoothed.
Scalar `C_n=2*(P_observed,n+N_n)^2/M_n` and all four Fisher entries are summed
without importing the production integral, covariance or Fisher kernels.

| `q_star` [s/km] | A [deg^2] | `P_pixel` [deg^2 km/s] | `Q_0=A*B_0+P_pixel` [deg^2 km/s] | `A*B_star+P_pixel` [deg^2 km/s] |
| ---: | ---: | ---: | ---: | ---: |
| 0.00035 | 0.029419845663300446 | 0.31435492218062416 | 0.7955679523681858 | 0.7955679523681858 |
| 0.001 | 0.029476112426414035 | 0.3134376504937113 | 0.7955710220007555 | 0.7923600798896060 |
| 0.003 | 0.031136479966332096 | 0.2888850209773887 | 0.7981766075292429 | 0.7077969402830528 |

Here `B_0=16.356748967852234 km/s` is fixed in every Q_0 comparison. Both
alternative Q_0 values satisfy the assigned lower bound relative to the
baseline. The final column is not used as a between-trial physical-mode
comparison; within each row it agrees with `B_star/(L*I1)` within 3.45e-15
relative.

The baseline public coefficients reproduce W09, and its Fisher matrix and
errors reproduce the saved `order32_iterations0_pair_fisher` result. Across all
public/scalar checks, the worst relative residual is 1.0295e-14, from the
0.003 covariance array. This is well below rtol=5e-12. The largest absolute
discrepancy in a signed sensitivity contrast is 3.77e-15 in fractional units,
also below the assigned 5e-12 tolerance. Every Fisher matrix is exactly
symmetric and positive definite under the analytic 2x2 checks; no prior,
regularization or pseudoinverse is used.

## Fisher matrices, errors and stopping decision

Matrices use parameter order (`ap_0`, `at_0`). Each public matrix below is a
fresh one-spectrum mode sum.

```text
q_star = 0.00035 s/km:
[[1419.738166245999,   667.4320627023011],
 [ 667.4320627023011, 1591.527428282942 ]]

q_star = 0.001 s/km:
[[1419.7803744056564,  667.4427262366473],
 [ 667.4427262366473, 1591.5297167799017]]

q_star = 0.003 s/km:
[[1416.5066441868166,  665.2780056341063],
 [ 665.2780056341063, 1584.5933127900716]]
```

All sensitivity ratios use the full-precision baseline errors in the first row
as their denominator:

| `q_star` [s/km] | sigma_parallel | sigma_perp | Delta_parallel | Delta_perp | Decision |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0.00035 | 0.02961952211255682 | 0.027975323168394888 | 0 | 0 | Baseline reproduced |
| 0.001 | 0.02961908469573798 | 0.027975305757164514 | -1.47678553751307e-5 | -6.223781676695239e-7 | Both below 0.001; continue |
| 0.003 | 0.029653944084125597 | 0.028037085296255797 | +0.0011621379790656139 | +0.0022077359925081197 | Both above 0.001; sensitivity detected and stop |

The signs are retained: the first alternative produces a very small decrease,
whereas the second produces increases above 0.1%. The criterion is based on
BAO errors, not Q_0 or the within-trial Q identity. Since 0.003 s/km is the last
assigned point, no prescribed reference trial is omitted; execution stops there
without interpolation or extension.

## Execution, artifacts and limitations

The sole bounded calculation was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/check_forest_reference_mode.py
```

It completed in 3.2873 s with one Python process, all numerical thread limits
set to one, the NumPy backend, Python 3.13.15 and NumPy 2.5.3. Exactly three
public weight preparations and three one-spectrum mode sums were performed.
The checkout remained at HEAD `0d69786a06d5d676564a51fad14a7156951c2228`
with extensive pre-existing modified/untracked work preserved.

The [script](../scripts/check_forest_reference_mode.py) has SHA-256
`211c482cec8c0e6f4e1d538603e1017e27761686860fb93c7a800613370850ed`.
The only numerical output is
[result.json](../.validation/forest-weight-diagnostics/w11-r1-20260916T044318380696Z/result.json),
SHA-256 `68c97083fd08cde9732dc77a4d1557a90d12a49a69d46b9bd089b446d32f4091`.
It contains full-precision public and independent matrices, covariance inverses,
errors, baseline-denominator sensitivity ratios, Q values, residuals, fixed
scalars and complete input/source manifests. The original saved evidence was
not copied or changed.

Final focused checks were Ruff lint and formatting for the new script plus
`git diff --check`; they pass. No pytest or broad suite was needed or run.

No cumulative comparison, magnitude refinement, other bin, population or
spectrum, joint covariance, P3D call, new derivative, raw-reader resampling,
physical geometry-volume calculation, full forecast, installation, Slurm
action, agent dispatch, commit or push was performed. Production modules,
profiles, tests, plans, the scientific note, sibling packages and historical
evidence are unchanged.

The saved calculation retains the inherited `floor_negative`, `legacy_floor`,
`legacy_spline_extension`, `legacy_floor_clamp` and
`legacy_first_spacing; unknown physical cells` policies. The result does not
validate those physical assumptions, the empirical P1D outside its calibration,
the saved signal/Jacobian/mode construction, a general reference-scale optimum,
or behavior in another spectrum. It only detects sensitivity above 0.1% at the
tested 0.003 s/km point under the stated held inputs. Stop here for independent
and user review; do not proceed to W12.
