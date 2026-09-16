# Scientific finding and effect of possible corrections

**W11 revision 1 passes independent scientific review in its exact bounded
scope. Reference-mode sensitivity is detected, and no scientifically
consequential correction is needed.** For the accuracy bin-0 `lya(qso)`
auto-spectrum at magnitude order 32, changing `q_star` from 0.00035 to
0.001 s/km changes the radial/transverse marginalized BAO errors by
-0.00147679%/-0.0000622378%, both below the 0.1% screen. The required
conditional point at 0.003 s/km changes them by **+0.116214%/+0.220774%**,
both above the screen. The stated classification and stop at 0.003 s/km are
therefore supported.

A correction to the sample or parameter identity, P1D floor/domain, response
ownership, common-mode noise comparison, covariance treatment or sensitivity
denominator could have changed that finding. The live rerun and an independent
raw-array reconstruction found no such error. No assignment revision or
production correction is warranted. The detected sensitivity is evidence that
the reference convention matters at the stated threshold in this calculation;
it does not select a new convention or imply that the sign must be positive.

This is a one-spectrum, one-bin diagnostic under inherited physical input
policies. It is not a reference-scale optimization, profile adoption, physical
validation of those policies, a broader population statement, acceptance of the
Step 12 forecasts, or authorization to proceed to W12.

## Identity, floor and reference domain

All four assigned saved-input SHA-256 identities match before and after the
fresh rerun and independent check. Metadata, not column position alone,
identifies:

- profile `accuracy`, bin 0 and bounds [2.0, 2.235];
- field `lya(qso)` with physical tracer `lya` and background `qso`;
- its unique selected and required auto-pair at index 0;
- parameter order (`ap_0`, `at_0`); and
- magnitude order 32 with 5,184 magnitude nodes and 16,384 Fourier nodes.

The report-native magnitudes, quadrature and variance match the diagnosis arrays
exactly. Independently reconstructed masses,
`density*(1+z_source)/299792.458*quadrature`, match within `2.22e-16` relative.
The W11 script binds these identities at
`scripts/check_forest_reference_mode.py:154-218`; the independent review repeats
the selection at
`.validation/forest-weight-diagnostics/w11-review-r1-20260916T045000Z/check.py:153-214`.

The independent PD2013 formula gives a stationary floor of
`0.0008164315784041575 s/km`; the saved parallel-mode domain is
`[1.3686646589112106e-7, 0.004863411498638425] s/km`. Thus 0.00035 lies below
the inherited floor, 0.001 and 0.003 lie above it, and all three lie inside the
saved domain. This verifies the assigned bracket without asserting empirical
P1D validity outside its calibration range.

## Intervention, response ownership and held information

Live source inspection supports the handoff's ownership statement:

- `fishhighz/models/p1d.py:33-56` returns intrinsic node-wise P1D without
  response or noise;
- `fishhighz/response.py:82-101` constructs the field transfer once;
- `fishhighz/weights.py:268-299` uses scalar `B_star` only to prepare the fixed
  inverse-variance weights;
- `fishhighz/noise.py:52-84` applies node-wise `W(q)^2` only to
  `A*P1D(q)`, leaves `P_pixel` unsmoothed, and converts both terms once;
- `fishhighz/covariance.py:151-218` rebuilds the one-spectrum Gaussian
  covariance from each trial's changed total power and the unchanged saved mode
  counts; and
- `fishhighz/fisher.py:206-220` contracts the unchanged saved mean Jacobian
  against that fixed fiducial covariance, with no covariance derivative.

The submitted script evaluates the physical node-wise P1D and response arrays
before the trial loop at `scripts/check_forest_reference_mode.py:311-350`. Within
the loop it changes `q_star -> B_star`, makes exactly one public weight
preparation, reconstructs noise and then rebuilds covariance and Fisher at lines
366-503. It never substitutes constant `B_star` for physical `P1D(q)` and does
not introduce mode-dependent weights.

The independent calculation constructs the physical P1D/response arrays once at
`.validation/forest-weight-diagnostics/w11-review-r1-20260916T045000Z/check.py:238-256`
and reuses them for all three scalar mode sums. The alternative covariance arrays
differ from the baseline by maximum pointwise relative amounts
`8.02364e-4` and `1.75064e-2`, respectively, confirming that covariance was not
incorrectly held identical between trials. Within each run one fixed weight
vector supplies A and P_pixel for all modes and the saved Jacobian is unchanged.

## Independent coefficients, Fisher sums and common-mode comparison

The independent calculation imports neither the W11 script nor FishHighz
weight, P1D, response, noise, covariance or Fisher helpers. It evaluates the
P1D/response formulas with Python scalar arithmetic, forms the weights and
moments with `math.fsum`, constructs
`C_n=2*(P_observed,n+N_n)^2/M_n`, sums all four Fisher entries, and performs the
positive-definite 2x2 inversion analytically. Its principal results are:

| `q_star` [s/km] | `B_star` [km/s] | A [deg^2] | `P_pixel` [deg^2 km/s] | common-mode `Q_0` [deg^2 km/s] | `sigma_parallel` | `sigma_perp` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00035 | 16.356748967852234 | 0.029419845663300394 | 0.3143549221806231 | 0.7955679523681838 | 0.02961952211255665 | 0.027975323168394843 |
| 0.001 | 16.247815263681932 | 0.029476112426413896 | 0.3134376504937124 | 0.7955710220007544 | 0.02961908469573796 | 0.027975305757164424 |
| 0.003 | 13.454055171253588 | 0.031136479966332276 | 0.2888850209773896 | 0.7981766075292469 | 0.02965394408412571 | 0.028037085296255898 |

Here every `Q_0=A*B_0+P_pixel` uses the fixed baseline
`B_0=16.356748967852234 km/s`. Both alternatives satisfy the assigned lower
bound relative to the baseline. The within-run identity
`A*B_star+P_pixel=B_star/(L*I1)` also passes for every point. The very small
error decreases at 0.001 are therefore allowed: the common-mode Q bound does
not determine the full mode-summed Fisher response.

The independent baseline Fisher matrix is

```text
[[1419.7381662460114, 667.4320627023003],
 [ 667.4320627023003, 1591.5274282829425]]
```

and reproduces the W09/saved order-32 fixed-reference matrix within
`3.43e-15` relative. All three independent matrices reproduce both the original
W11 artifact and the fresh public rerun within the assigned `5e-12` tolerance.
The largest discrepancy in the complete review is `9.55e-15`, from the 0.003
sensitivity fraction. Every matrix is symmetric and positive definite; no
prior, regularization or pseudoinverse is used.

All signed sensitivities use the independent full-precision baseline error
vector as denominator:

| `q_star` [s/km] | `Delta_parallel` | `Delta_perp` | Screen |
| ---: | ---: | ---: | --- |
| 0.001 | -1.47678553701347e-5 | -6.22378169334859e-7 | both below 0.001; continue |
| 0.003 | +0.00116213797907516 | +0.00220773599251323 | both above 0.001; sensitivity detected |

The first decision correctly permits the conditional second point. The second
decision crosses the screen in both parameters and ends the sequence. Because
0.003 is the final assigned point, there is no omitted prescribed trial; no
interpolation or additional scan was performed.

## Execution, evidence and disposition

The final live-source bounded rerun was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/check_forest_reference_mode.py
```

It completed successfully in 3.34 s wall time and wrote
`.validation/forest-weight-diagnostics/w11-r1-20260916T044755730574Z/result.json`
(SHA-256 `bd7a1d06829fe0bc9c3200fa61d30fe885c96f3872aadd410d925eb25fbac71f`).
The three error vectors and contrasts exactly reproduce the handoff run.

The separate independent calculation is
[check.py](../.validation/forest-weight-diagnostics/w11-review-r1-20260916T045000Z/check.py),
with full-precision output in
[result.json](../.validation/forest-weight-diagnostics/w11-review-r1-20260916T045000Z/result.json).
Their SHA-256 values are respectively
`85d5ad54f9152e2819f96cc1cae6637ba24d51c1f9121b706440830f49a61304`
and `351755e925fcfbb4e7382d8efb4f2064f458a0e6ccb6f1cb6dfd1c8c9e3589e5`.
The independent check completed successfully in 0.84 s wall time.

Focused lint/format checks for the submitted and independent scripts and
`git diff --check` pass. No broad pytest suite was needed or run. No production
source, test, profile, plan, assignment, historical evidence or sibling package
was changed. No extra mode, population, magnitude refinement, cumulative
comparison, P3D/model call, derivative, joint covariance, full forecast,
installation, Slurm action, dispatch, commit or push was performed.

W11 revision 1 is ready for user review. Scientific acceptance, selection of a
reference convention, profile adoption and any progression remain with the
user.
