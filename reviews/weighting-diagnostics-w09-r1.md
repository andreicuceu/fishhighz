# W09 revision 1: fixed-reference BAO replay

**The public fixed-reference calculation reproduces the saved W07 result for
accuracy bin 0, `lya(qso)` auto-spectrum, magnitude order 32.** All assigned
checks pass at rtol=5e-12, atol=0, with exact handling of expected zeros. The
largest relative residual is 8.489e-15 (public versus scalar Fisher). This is
implementation evidence awaiting independent/user review and acceptance.

Only one saved spectrum/bin was calculated. This establishes numerical
consistency with W07, not a general multi-mode optimum, physical validation of
the retained inputs, adoption of a profile default, or acceptance of the
unconverged Step 12 forecasts. Package Step 13 remains separate.

## Selection, identities and held quantities

Metadata identifies accuracy bin 0, bounds [2, 2.235], the forest field
`lya(qso)` with physical tracer `lya` and background `qso`, and parameter order
`ap_0`, `at_0`. Its unique auto-spectrum resolves to field index 0 and
selected/required pair index 0. The primary NPZ pair arrays agree with that
metadata. The script explicitly slices the signal, noise and mean Jacobian,
then remaps to a one-field `PairSelection` with pair (0,0), before any public
covariance or Fisher construction. There are 5,184 magnitude nodes and 16,384
saved Fourier nodes. Primary final controls explicitly specify order 32/t=6.

The four input identities below were checked before and after calculation.
Paths are relative to `.validation/step12-r5-20260914T191855Z/`:

| Input | SHA-256 |
| --- | --- |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |

The matching diagnosis bin binds its NPZ hash; the selected rows identify
order32/iterations0,3,6 and their `pair_fisher_key`. Only iterations0 is marked
as fixed initial weights. No cumulative iteration was executed.

Saved observed signal, mean Jacobian, k, mu, mode counts, magnitude distribution,
quadrature, variance, source redshift, forest length and conversion factors
remain fixed. No derivative, raw-reader interpolation, renormalization or new
flooring was performed. The report density is converted explicitly as
rho=dN/(dz dm deg²)*(1+z_source)/c. Its product with the saved quadrature agrees
exactly with `order32_field0_masses`; magnitudes, quadrature and variance also
agree exactly between the two saved inputs.

| Fixed quantity | Value |
| --- | ---: |
| z_eval | 2.1152848986890427 |
| z_source | 2.3830099582078716 |
| h_fid | 0.6736 |
| a_v [(km/s)/(Mpc/h_fid)] | 102.61326549228777 |
| d_deg [(Mpc/h_fid)/degree] | 64.41180749124948 |
| L [km/s] | 44147.09373976799 |
| l_p [km/s] | 63.328211161339375 |
| Gaussian sigma_v [km/s] | 50.92405402826606 |
| B_star [km/s], already response-smoothed | 16.356748967852234 |

B_star comes from the matching diagnosis field, at the historical reference
parallel mode 0.00035 s/km. The pixel width comes from the recorded 0.8 Angstrom
width and observed wavelength. Resolving power R=2500 comes from the retained
case recipe; its full source hash matches the report's historical producer
hash. The report specifies the FWHM convention, giving
sigma_v=c/(2 R sqrt(2 ln 2)).

The saved background does not supply H(z) and D_M(z) across the bin. The allowed
diagnostic adapter calls `prepare_geometry` with constant
H=a_v*h_fid*(1+z_eval) and D_M=d_deg/(h_fid*pi/180), saved bin bounds and area,
and a one-node preparation rule. It reproduces both conversion factors exactly.
**Its integrated volume is artificial and never enters this calculation.**
Only saved `modes` supplies volume/mode normalization. No frozen dataclass is
bypassed and no production geometry API is added.

## Independent formulas and ordered checks

The historical control ran first, before public reference weight preparation.
For q=k*mu/a_v, the retained analytic `default_p1d` and response give

```text
W(q) = sinc(q*l_p/(2*pi)) exp[-(q*sigma_v)^2/2]
N(q) = [A P1D(q) W(q)^2 + P_pixel] d_deg^2/a_v
C_n = 2 [P_observed,n + N_n]^2 / modes_n
F_ab = sum_n J_na J_nb / C_n.
```

The scalar C/F calculation uses Python scalar arithmetic and `math.fsum`, with
no production covariance, factorization or contraction helper. Saved t=6
coefficients A=0.022340631461880645 deg² and
P_pixel=0.5628751956005517 deg² km/s reproduce the original noise and the selected
historical Fisher matrix. W is checked against the saved field response. This
control fixes the response, unit, pair and mode conventions before proceeding.

Public weights are prepared once with `method="inverse_variance", alias=B_star`;
no signal, iteration count or auxiliary sampling is passed. Independently,
with r_i=rho_i*quadrature_i and nu_i=B_star/(B_star+l_p*v_i), the script sums

```text
I1 = sum_i r_i nu_i
I2 = sum_i r_i nu_i^2
I3 = sum_i r_i nu_i^2 v_i
A = I2/(L I1^2),  P_pixel = l_p I3/(L I1^2).
```

These sums do not import the production integral kernel. They yield
I1=0.0004657121726488611, I2=0.0002816940591039239,
I3=0.00004752918223323785. Public coefficients are
A=0.029419845663300446 deg² and P_pixel=0.31435492218062416 deg² km/s,
identical to the reviewed W08 printed values. Independent scalar coefficients
are 0.029419845663300394 and 0.3143549221806231.

Intrinsic, mode-dependent P1D is passed to `forest_noise`; response is applied
once to aliasing and never to the pixel term. Constant B_star sets the weights
only. All selected modes are checked against the independent noise formula,
including a separate unsmoothed-pixel check. Public `combine_observed_power`,
`gaussian_covariance` and `fisher_matrix` then operate on the one-spectrum
arrays. Covariance, weights and noise remain fiducial; J differentiates the
saved mean alone. No joint inverse or joint inverse subblock is used.

For each 2x2 F=[[a,b],[b,d]], the independent inverse is
[[d,-b],[-b,a]]/(a*d-b²), with marginalized errors
(sqrt(d/(a*d-b²)), sqrt(a/(a*d-b²))). Exact symmetry, positive diagonals and
positive determinant are checked; no prior, regularization or pseudoinverse
is introduced. The new reference determinant is 1814086.6742375866 and its
normalized minimum eigenvalue is 0.5559868082510999.

| Check | Worst elementwise relative residual |
| --- | ---: |
| Historical W / saved response | 2.274e-16 |
| Historical t=6 noise / saved noise | 4.233e-16 |
| Historical scalar F / diagnosis t=6 F | 5.659e-16 |
| Historical scalar F / primary pair F | 7.119e-16 |
| Reference moments / independent sums | 1.514e-15 |
| Reference coefficients / independent sums | 3.356e-15 |
| Reference coefficients / saved fixed coefficients | 1.766e-16 |
| Reference noise / independent formula | 2.872e-15 |
| Unsmoothed pixel term / independent formula | 3.355e-15 |
| Reference C / scalar C | 5.573e-15 |
| Reference F / scalar F | 8.489e-15 |
| Reference F / saved W07 F | 6.887e-15 |
| Reference errors / scalar errors | 5.740e-15 |
| Reference errors / W07 errors | 4.803e-15 |
| t=3 and t=6 error ratios / historical-reference ratios | 4.907e-15 |

All use rtol=5e-12, atol=0 for nonzero operands; expected zeros must be exact.
The printed contrast checks use their quoted rounding precision, 5e-8 percentage
points. No tolerance was changed and no consequential failure occurred.

## New mode sums and historical comparisons

Matrices below use parameter order (ap_0, at_0). The first three are newly
summed; the remaining three are saved historical matrices.

```text
Historical t=6 scalar control, newly summed:
[[1205.538750489211,  553.0960157089297],
 [ 553.0960157089297, 1277.7024563498699]]

New public fixed reference:
[[1419.738166245999,  667.4320627023011],
 [ 667.4320627023011, 1591.527428282942]]

New independent scalar fixed reference:
[[1419.7381662460111, 667.4320627023002],
 [ 667.4320627023002, 1591.5274282829423]]

Saved W07 fixed reference:
[[1419.7381662460089, 667.4320627022987],
 [ 667.4320627022987, 1591.527428282937]]

Saved cumulative t=3:
[[1227.5500924474698, 564.8346409624346],
 [ 564.8346409624346, 1309.674936936257]]

Saved cumulative t=6:
[[1205.5387504892117, 553.0960157089298],
 [ 553.0960157089298, 1277.7024563498703]]
```

The newly reconstructed public reference has marginalized errors
**sigma_parallel=0.02961952211255682** and
**sigma_perp=0.027975323168394888**. Its parameter covariance is

```text
[[ 0.0008773160901762424, -0.000367916303107626 ],
 [-0.000367916303107626,   0.0007826187063761317]].
```

Each contrast is **100*(sigma_saved_cumulative/sigma_new_reference - 1)**.
The denominator is the new public result above, not the old saved reference,
not a Q ratio, and not the cumulative result.

| Historical numerator | sigma_parallel | sigma_perp | Radial contrast [%] | Transverse contrast [%] |
| --- | ---: | ---: | ---: | ---: |
| t=3 | 0.03187966090130967 | 0.030863952834406856 | +7.630571418958487 | +10.3256346624635 |
| t=6 | 0.03217259526976886 | 0.03125084883353348 | +8.619562285678107 | +11.708624938564128 |

These historical cumulative comparators have larger errors in this sample.
They were inverted from full saved matrix entries; their cumulative weights
were not recomputed. All errors are dimensionless and marginalized over the
other BAO parameter.

## Execution, artifacts and limitations

The single numerical replay completed in **1.2747 s**, under both the external
60-second timeout and internal 60-second alarm, on login32 with one process,
all three numerical thread limits equal to one, and the NumPy backend.
Interpreter: Python 3.13.15; NumPy 2.5.3. The live checkout is branch `main`,
HEAD `0d69786a06d5d676564a51fad14a7156951c2228`, with extensive pre-existing
modified/untracked work. It was inspected before editing and preserved.
The machine UTC timestamp placed this run on 2026-09-16.

Commands from the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 .venv/bin/python -B scripts/replay_inverse_variance_bao.py
.venv/bin/ruff check scripts/replay_inverse_variance_bao.py
.venv/bin/ruff format --check scripts/replay_inverse_variance_bao.py
git diff --check
```

All checks passed. Only the new script was formatted before execution. No
pytest suite was needed or run for this diagnostic-only assignment.

The [script](../scripts/replay_inverse_variance_bao.py) has SHA-256
`135cd6bbe06736a7b042b9aff80c643f0a0fbabd1ab6c81828829b1e7309a811`.
The small [result.json](../.validation/forest-weight-diagnostics/w09-r1-20260916T015605976471Z/result.json)
contains all full-precision matrices, covariance inverses, errors, contrast
numerators/denominators, residuals, fixed scalars, exact input hashes, and the
full SHA-256 manifest of the script and decisive source files. Those source
hashes were rechecked after execution. The adjacent
[git-status.txt](../.validation/forest-weight-diagnostics/w09-r1-20260916T015605976471Z/git-status.txt)
records the live dirty/untracked state at execution. Original evidence is
unchanged; no large arrays or copies of input bundles were written.

The only additions are this report, the replay script, and that fresh artifact
directory. No production code, profiles, existing tests, planning documents,
scientific note or other package was edited. No full forecast, other bin or
spectrum, magnitude refinement, new physical integration, new derivative,
installation, Slurm action, agent dispatch, commit or push was performed.
The diagnostic geometry preparation described above is not a physical volume
calculation and its volume was discarded.

The inherited nonnegative accuracy inputs retain `floor_negative`,
`legacy_floor` magnitude support, `legacy_spline_extension` redshift support,
`legacy_floor_clamp` SNR, and `legacy_first_spacing; unknown physical cells`.
This numerical agreement does not validate those policies or the saved
signal, derivatives and Fourier quadrature independently. It adds no
convergence, continuum-limit or general optimality claim. Stop here for user
review and acceptance; no subsequent step is planned or executed.
