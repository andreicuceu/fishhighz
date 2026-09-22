# W11 — sensitivity to the fixed reference parallel mode

Revision: 1. Status: proposed for user approval and implementation dispatch.
The user requested progression after W10 passed independent scientific review.
This is the parallel weighting sequence; package Step 13 remains separate.

## Scientific question and minimal comparison

How sensitive are the BAO errors to the parallel reference mode used to set
B_star in the fixed inverse-variance weights? Does the inherited reference
convention matter at the 0.1% level in one controlled saved-input calculation?

Return to **accuracy bin 0, `lya(qso)` auto-spectrum, magnitude order 32**,
parameters (`ap_0`, `at_0`), 5184 magnitude nodes and 16384 Fourier nodes.
W09 already verified its public noise/covariance/Fisher chain, so no additional
population, redshift, magnitude-refinement or cumulative comparison is needed.
W10 extends the empirical reference stability to two other samples; it does
not choose B_star or establish a universal optimum.

Use the following three reference modes, in this order, with early stopping:

| Role | q_star [s/km] | Purpose |
| --- | ---: | --- |
| Baseline | 0.00035 | Reproduce the W09 inherited reference convention. |
| First alternative | 0.001 | Probe just beyond the model's low-q plateau. |
| Conditional second alternative | 0.003 | Probe farther along the declining P1D, only if the first alternative is insensitive by the criterion below. |

These are diagnostic probes, not candidate defaults to rank or optimize. The
retained analytic model imposes a stationary low-q floor at approximately
0.00081643 s/km at this redshift. Halving/doubling 0.00035 would keep both
intrinsic P1D samples below that floor and mostly test response smoothing.
The chosen points bracket that feature and lie inside the saved q range
(approximately 1.37e-7 to 0.0048634 s/km). Check the exact floor and range from
inputs, and retain the floor as an inherited modelling assumption. No claim
about empirical P1D validity follows from this choice.

Read workspace/package AGENTS.md, design section 0, main-roadmap handover and
progress register, and IMPLEMENTATION_STEP.md as separate package context.
Then read the diagnostic roadmap (local-only path: `../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md`),
W09 review (local-only path: `reviews/weighting-diagnostics-w09-review-r1.md`),
W10 review (local-only path: `reviews/weighting-diagnostics-w10-review-r1.md`), and only relevant
parts of `replay_inverse_variance_bao.py`, public weight/noise/covariance/Fisher
functions, `models/p1d.py` and response functions. Do not dispatch other agents.

## Immutable inputs and held calculation

Paths below are relative to `.validation/step12-r5-20260914T191855Z/`:

| Input | SHA-256 |
| --- | --- |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |

Verify identities and metadata for profile, bin bounds [2,2.235], forest/background
labels, unique selected auto pair and parameter order. Anchor the order-32
magnitude, variance and quadrature arrays in the report, and check
`rho*quadrature` against saved masses using `rho=dN/(dz dm deg^2)*(1+z_source)/c`.
Use only the selected column of saved observed signal and mean Jacobian, and
saved k, mu and mode counts. Do not construct a joint 15-spectrum covariance.

Hold geometry conversions, redshift, density/SNR policies, L, l_p, response,
signal, mean derivatives, integration nodes and mode counts fixed across runs.
Use the W09 geometry context adapter if needed; its artificial volume must never
enter the calculation. Inherited flooring, support, SNR and source-cell-width
policies are unchanged and remain unvalidated physical assumptions.

The only intervention is q_star -> B_star -> fixed weights. Within each run,
weights/noise/covariance remain fixed when interpreting the saved mean Jacobian.
Across the runs, covariance is recomputed from the changed noise. Do not confuse
fixed covariance during differentiation with holding it identical across tests.

## Equations and response ownership

At z_eval=2.1152848986890427, use the same retained intrinsic analytic P1D and
instrument response as W09:

```text
q_n = k_n*mu_n/a_v
W(q) = sinc(q*l_p/(2*pi))*exp[-(q*sigma_v)^2/2]  # NumPy sinc convention
B_star(q_star) = P1D(z_eval,q_star)*W(q_star)^2
nu_i(q_star) = B_star/(B_star+l_p*v_i)
I1 = sum r_i nu_i; I2 = sum r_i nu_i^2; I3 = sum r_i nu_i^2 v_i
A = I2/(L I1^2); P_pixel = l_p I3/(L I1^2)
N_n = [A*P1D(z_eval,q_n)*W(q_n)^2 + P_pixel]*d_deg^2/a_v
C_n = 2*(P_observed,n + N_n)^2/M_n
F_ab = sum_n J_na*J_nb/C_n.
```

Here r_i=rho_i*dm_i is nonnegative, B_star has units km/s, A deg^2,
P_pixel deg^2 km/s, and N comoving P3D units. Record intrinsic P1D, W squared
and B_star separately at each reference mode. Do not smooth B_star twice.
The physical P1D(q_n) and W(q_n) arrays at all forecast nodes are evaluated once
and held fixed; constant B_star never replaces them in the noise formula.
Pixel noise is not response-smoothed. There are no mode-dependent weights.

To compare combined noise at a common physical mode, report
`Q_0(q_star)=A(q_star)*B_0+P_pixel(q_star)`, where B_0 is the baseline B_star.
Do not compare `A*B_star+P_pixel` between different q_star as though the physical
mode were held fixed. The latter quantity may be used solely for the independent
identity `A*B_star+P_pixel=B_star/(L*I1)` within each run. For nonnegative measure,
Q_0 is minimized by the baseline reference; that does not imply its BAO errors
are minimal when information from all q_n is combined.

## Ordered checks, scientific criterion and early stopping

1. **Baseline.** Reconstruct the baseline B_0 from the analytic model/response
   and compare to saved alias 16.356748967852234 km/s. Prepare public weights,
   noise, one-spectrum covariance and Fisher, reproducing W09's coefficients
   and `order32_iterations0_pair_fisher`. Expected marginalized errors are
   approximately (0.029619522112556678, 0.027975323168394890). Use full saved
   values for numerical tests. No need to repeat W09's historical t=6 control.
   Stop if identity, units, response or baseline equivalence fails.
2. **First alternative.** Compute B_star at q_star=0.001, prepare weights once
   with public `method="inverse_variance", alias=B_star`, then noise/covariance/
   Fisher using the unchanged mode inputs. Do not pass signal, iterations or
   auxiliary sampling. Record A, P_pixel, Q_0, 2x2 Fisher and the two errors.
3. **Independent checks for every computed run.** Reconstruct scalar P1D and W
   at the reference points from the retained formula, without calling their
   production helpers. Sum moments with `math.fsum`, reconstruct mode noise
   and `C_n=2*T_n^2/M_n`, and sum the four Fisher entries with scalar arithmetic.
   Compare to public outputs and invert the positive-definite 2x2 F analytically.
   Do not import production integral/covariance/contraction helpers as the oracle.
   Verify the within-run Q identity, and Q_0 >= baseline Q_0 within the numerical
   tolerance. If that bound fails beyond rounding, diagnose before interpreting
   Fisher changes; do not require either BAO error to increase.
4. **Sensitivity decision.** Report signed fractional changes
   `Delta_a(q_star)=sigma_a(q_star)/sigma_a(q_baseline)-1`, for radial and
   transverse errors. If either absolute Delta exceeds 0.001 (0.1%), record
   **reference-mode sensitivity detected** and stop; omit the second alternative.
   This is a valid scientific result, not an implementation failure. Otherwise
   run q_star=0.003 once with the same checks and decision. If both alternatives
   stay within 0.1%, report **no sensitivity above 0.1% at the tested points**,
   not stability throughout the interval or a general reference-scale optimum.

The 0.1% threshold is an operational sensitivity screen at the existing forecast
accuracy scale. It is not a new accuracy-profile acceptance rule, convergence
test or authority to choose a preferred q_star. Significant changes of either
sign require a user decision before profile adoption. Do not extend the scan,
interpolate an optimum, compare other samples or proceed to W12 automatically.

For nonzero numerical quantities use rtol=5e-12, atol=0, with exact expected zeros.
For signed contrasts near zero use an absolute 5e-12 tolerance in fractional
units (5e-10 percentage points), retaining full-precision operands. Apply the
Q_0 bound with a relative 5e-12 rounding allowance. Do not relax equivalence
criteria to 0.1%; distinguish numerical failures from detected physical-model
sensitivity. Positive definiteness and ordinary 2x2 inversion are required;
no regularization, prior or pseudoinverse.

## Minimal implementation, execution and handoff

After approval/dispatch, add at most one small standalone script,
`scripts/check_forest_reference_mode.py`, a concise handoff and small numerical
outputs. Reuse W09 conventions without editing/refactoring the historical script.
Do not change production modules, tests, profile recipes, plans, the scientific
note or either scientific package. No new framework, evidence schema, broad
suite, repeated W08 range campaign, installation or wheel is needed.

At most three public preparations and three single-spectrum mode sums; one
process, OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1,
`.venv/bin/python -B`, 60-second cap for the entire assignment. Only small
analytic P1D/response evaluations and saved-input sums are authorized. No
P3D, CAMB, raw reader, derivative or new grid/geometry-volume evaluation; the
explicitly labelled W09 context adapter is allowed. No full forecast, other
bin/spectrum, magnitude refinement, Slurm, dispatch, commit or push. Stop and
report an obstruction if an identity, arithmetic or runtime bound fails.

Write `reviews/weighting-diagnostics-w11-r1.md` with outputs in a new
`.validation/forest-weight-diagnostics/w11-r1-<UTC>/` directory. Include the
question, intervention/held quantities, exact source/input identities and command,
plateau/domain check, compact reference/P1D/response/coefficient/error table,
2x2 matrices, Q_0 and error contrasts with baseline denominators, worst residuals,
actual stopping point and omitted trials. Distinguish fresh mode sums from
historical baseline evidence. Save only what is needed to reproduce the result,
not large copied inputs. State limitations and stop for user review.

## Independent review

Review only the decisive source/evidence, rerun the bounded assignment and
independently reconstruct reference B_star and coefficient/error contrasts, plus
one alternative's scalar noise/Fisher calculation if reached. A separate short
calculation should verify the intervention affects weights alone, the baseline
matches W09, final P1D/response remain mode-dependent and unchanged, and the
sensitivity/early-stop decision follows the declared denominator and threshold.
No additional modes, populations, test campaigns or optimization.

Begin `reviews/weighting-diagnostics-w11-review-r1.md` with the scientific finding
and how any proposed corrections could affect it. Propose changes only when
likely to affect the conclusion or supported scope, with the smallest closing
check; explicitly state when none are needed. Detected sensitivity is not a
reason to repair production or change the reference automatically. Revise this
same assignment only for consequential corrections. The user controls acceptance,
reference-convention choices, profile adoption and any W12 request.
