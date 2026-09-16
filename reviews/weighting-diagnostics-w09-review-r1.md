# Scientific finding and impact of proposed changes

**W09 revision 1 passes independent scientific review. No scientifically
consequential change is needed.** The public fixed-reference weighting, noise,
one-spectrum Gaussian covariance and mean-only Fisher calculation reproduces the
saved W07 fixed-reference result for accuracy bin 0 and the `lya(qso)`
auto-spectrum. The historical t=3/t=6 comparisons therefore retain their stated
sign and magnitude. A correction to sample/parameter identity, response ownership,
the Gaussian factor of two, mode normalization, fixed-covariance treatment or the
percentage denominator could have changed that conclusion; the bounded rerun and
an independent scalar reconstruction found none.

This is a pass of the assigned numerical connection only. It is not user
acceptance, profile-default adoption, physical validation of the inherited
density/SNR policies, a general optimality result, or acceptance of the Step 12
forecasts. No revision of `WEIGHTING_DIAGNOSTIC_STEP.md` is warranted.

## Identity and held calculation

All four assigned immutable-input SHA-256 identities match. The public rerun also
rechecked its decisive source manifest after calculation. Metadata, rather than
column position alone, identifies:

- profile `accuracy`, bin 0, bounds [2.0, 2.235];
- field `lya(qso)`, physical tracer `lya`, background `qso`;
- its unique auto-spectrum at selected and required pair index 0;
- magnitude order 32 with 5,184 nodes and 16,384 saved Fourier nodes; and
- parameter order (`ap_0`, `at_0`).

The primary pair arrays agree with the report metadata, and final primary controls
identify order 32/t=6. The saved signal, mean Jacobian, `k`, `mu`, mode counts,
source distribution, response settings and geometry conversion remain fixed. The
replay slices and remaps the one auto-spectrum before public covariance/Fisher
assembly at `scripts/replay_inverse_variance_bao.py:210-213` and
`scripts/replay_inverse_variance_bao.py:319-324`. No joint inverse or covariance
subblock is used.

The artificial constant-background geometry adapter is used only to construct a
valid public context and recover `a_v` and `d_deg`; its volume is discarded. The
calculation uses only the saved `modes` for normalization, as recorded at
`scripts/replay_inverse_variance_bao.py:151-202`.

## Response, noise and fixed-information ownership

Source inspection and the historical control establish the intended ownership:

- `fishhighz/weights.py:268-299` prepares
  `nu=B_star/(B_star+l_p*v)` once, with no signal, iterations or auxiliary
  sampling, then forms the fixed A and P_pixel coefficients;
- `fishhighz/noise.py:52-84` evaluates the intrinsic mode-dependent P1D and
  applies W squared only to `A*P1D`; the P_pixel term is unsmoothed and both terms
  receive `d_deg^2/a_v` once;
- `fishhighz/covariance.py:151-174` uses the saved counts in
  `C_AB=(T_im*T_jn+T_in*T_jm)/modes`, which reduces to
  `C_n=2*T_n^2/modes_n` for this auto-spectrum; and
- `fishhighz/fisher.py:1-6` and `fishhighz/fisher.py:206-220` contract the saved
  mean Jacobian against that fixed fiducial covariance without covariance
  derivatives or another mode factor.

Before preparing the reference weights, the replay reconstructs the historical
t=6 response, noise and Fisher matrix from its saved coefficients. The independent
review rebuilt W and the PD2013 P1D formula directly rather than importing the W09
script or production response/P1D/covariance/Fisher functions. Residuals were
2.274e-16 for W, 2.166e-16 for historical noise, 5.658e-16 for the diagnosis
Fisher matrix and 7.118e-16 for the primary Fisher matrix. This control rules out
a response, unit, spectrum or mode-count convention change at the precision
relevant here.

## Independent scalar reconstruction

From the saved nonnegative density, quadrature and variance, the independent
fixed-reference sums give

```text
A       = 0.029419845663300394 deg^2
P_pixel = 0.3143549221806231   deg^2 km/s
```

They agree with the reviewed W08/public values within 3.36e-15 relative. Using
Python scalar arithmetic and `math.fsum` for

```text
C_n  = 2*(P_observed,n + N_n)^2/modes_n
F_ab = sum_n J_na*J_nb/C_n
```

gives

```text
historical t=6 control:
[[1205.538750489211,   553.0960157089297],
 [ 553.0960157089297, 1277.7024563498699]]

independent fixed reference:
[[1419.7381662460114, 667.4320627023003],
 [ 667.4320627023003, 1591.5274282829425]]
```

The fixed-reference Fisher matrix agrees with the saved W07 matrix within
3.43e-15 and with the new public rerun within 8.65e-15 relative. The independent
determinant is 1814086.6742376077, so the matrix is positive definite. Analytic
2x2 inversion gives

```text
sigma_parallel = 0.02961952211255665
sigma_perp     = 0.027975323168394843
```

These agree with the public rerun within 5.74e-15 relative. No prior,
regularization or pseudoinverse enters the calculation.

For both historical comparisons, the public artifact stores its denominator as
the new public fixed-reference error vector exactly. The independent check then
recomputed

```text
100*(sigma_saved_cumulative/sigma_new_reference - 1)
```

and found:

| Historical numerator | Radial contrast | Transverse contrast |
| --- | ---: | ---: |
| t=3 | +7.630571418959109% | +10.325634662463678% |
| t=6 | +8.619562285678750% | +11.708624938564327% |

Thus the denominator is neither the cumulative result nor a Q ratio. The largest
independent relative residual is 8.15e-14, arising in a percentage comparison
after subtracting unity; all primary response/noise/Fisher/error residuals are at
most 8.65e-15. Both are well below the assigned 5e-12 numerical criterion.

## Execution and disposition

The public replay was rerun from the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/replay_inverse_variance_bao.py
```

It passed under the 60-second cap and wrote
`.validation/forest-weight-diagnostics/w09-r1-20260916T021526420747Z/`;
its internal elapsed time was 2.030 s. The independent calculation is in
`.validation/forest-weight-diagnostics/w09-review-r1-independent/check.py` with
full-precision output in the adjacent `result.json`. It imports neither the W09
replay nor the production numerical helpers named above.

Final bounded checks passed:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B \
  .validation/forest-weight-diagnostics/w09-review-r1-independent/check.py
.venv/bin/ruff check scripts/replay_inverse_variance_bao.py \
  .validation/forest-weight-diagnostics/w09-review-r1-independent/check.py
.venv/bin/ruff format --check scripts/replay_inverse_variance_bao.py \
  .validation/forest-weight-diagnostics/w09-review-r1-independent/check.py
git diff --check
```

No full forecast, W08 test campaign, other bin/spectrum, magnitude refinement,
production repair, profile change, plan update, Slurm action, agent dispatch,
commit or push was performed. W09 revision 1 is ready for user review; approval,
acceptance and any progression remain with the user.
