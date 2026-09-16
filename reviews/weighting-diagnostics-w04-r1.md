# Forest-weighting diagnostic W04 revision 1

## Scientific answer and limits

**Both full-sample branches satisfy the prescribed 0.1% coefficient-stability
criterion over t=6/12/24. Additional iterations alone reduce the fixed-intrinsic
branch's original drift; aliasing is not necessary for this bounded stability.**
Explicitly refreshing S_t=P0+A[w_t]B reduces the measured iteration drift and
changes the coefficient levels: at t=24, A is 5.08747% lower and P_pixel is
12.60455% higher than with fixed P0. This comparison does not determine which
weight preparation is physically appropriate or optimal.

Neither branch passes the three-checkpoint test at (3,6,12), because comparisons
involving t=3 exceed 0.1%. This requires the permitted extension to 24; both
then pass all three comparisons (12/6,24/12,24/6), and execution stops.
These are fixed-grid finite-count observations, not a proof of an asymptotic or
magnitude-continuum limit. The nine signed negative-density nodes remain literal
historical inputs, not a physical nonnegative sampling measure.

Status: implemented; awaiting independent/user review. The user's explicit
approval supersedes the assignment's awaiting-approval text. Planning documents
and package Step 13 are unchanged.

## Checkpoints and matched-count signal effect

A has units deg² and P_pixel has units deg² km/s. All entries below are newly
computed; fixed t=3/6 additionally replay the historical W03 checkpoints.
W03 did not execute full-sample t=12/24. No prefix trajectory was run here.

| t | A, fixed P0 | A, refreshed | P_pixel, fixed P0 | P_pixel, refreshed |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 0.022021379679306282 | 0.020832476681501746 | 0.5298546380249203 | 0.6001368886956997 |
| 6 | 0.021818865253617570 | 0.020705076728801980 | 0.5405425825641871 | 0.6089022322878723 |
| 12 | 0.021812392226710760 | 0.020702687879214025 | 0.5408919352169831 | 0.6090692663929413 |
| 24 | 0.021812385938153173 | 0.020702687061427468 | 0.5408922748520478 | 0.6090693235917064 |

Matched-count differences are refreshed minus fixed, with fixed as the
fractional reference:

| t | ΔA | ΔA/A_fixed | ΔP_pixel | ΔP_pixel/P_fixed |
| ---: | ---: | ---: | ---: | ---: |
| 3 | −0.001188902997805 | −5.398858% | +0.07028225067078 | +13.264440% |
| 6 | −0.001113788524816 | −5.104704% | +0.06835964972369 | +12.646487% |
| 12 | −0.001109704347497 | −5.087495% | +0.06817733117596 | +12.604612% |
| 24 | −0.001109698876726 | −5.087471% | +0.06817704873966 | +12.604552% |

Each branch's signed fractional changes use the earlier count as reference:

| Counts | Fixed ΔA/A | Fixed ΔP/P | Refreshed ΔA/A | Refreshed ΔP/P |
| --- | ---: | ---: | ---: | ---: |
| 6 vs 3 | −0.919626% | +2.017147% | −0.611545% | +1.460557% |
| 12 vs 6 | −0.0296671% | +0.0646300% | −0.0115375% | +0.0274320% |
| 12 vs 3 | −0.949021% | +2.083080% | −0.623012% | +1.488390% |
| 24 vs 12 | −0.0000288302% | +0.0000627917% | −0.00000395015% | +0.00000939118% |
| 24 vs 6 | −0.0296959% | +0.0646928% | −0.0115415% | +0.0274414% |

The largest absolute fractional change in the final triple is 0.0006469283
(fixed) and 0.0002744140 (refreshed), both below 0.001. Signed absolute changes
for every comparison are retained in `summary.json`. There are no zero-reference
coefficients; the script separately handles exact zero-to-zero and
zero-to-nonzero comparisons without dividing by zero.

## Explicit-refresh assumption and held quantities

Each independent branch starts at w0=(B/lp)/(B/lp+v). Full-sample moments use
W03's float64 product order and the last `cumsum` element:

```text
I1 = cumsum((rho*w)*dm)[-1]
I2 = cumsum((rho*w**2)*dm)[-1]
I3 = cumsum(((rho*w**2)*v)*dm)[-1]
A = I2/(L*I1**2); P_pixel = lp*I3/(L*I1**2)
N = I1*(L/lp)
S = P0                  [fixed]
S = P0 + A*B            [refreshed, current pre-update weights]
w_new = S/(S + v/N)     [simultaneous over all nodes]
```

This is the approved diagnostic refresh convention. W02's historical source
assessment leaves initialization and refresh of an aliasing-inclusive signal
unresolved; this experiment does not resolve that historical uncertainty.
No pixel-noise term enters S, and no second aliasing term is added to final
noise. Neither branch borrows A from the other or rescales its evolving weights.

| Held quantity | Value |
| --- | ---: |
| Magnitude nodes | all 107, endpoints 16.1 and 26.75 mag |
| Rectangular dm, including endpoints | 0.100471698113207 mag |
| L | 44148.20436604321 km/s |
| lp | 63.32980433473595 km/s |
| Intrinsic P0 | 1.7189455896517238 deg² km/s |
| B | 16.333147249645222 km/s |

The signed rho array, dimensionless variance v, grid, response, population,
redshift bin and seed are fixed. Negative-density indices are
0,2,3,76,77,80,81,84,85; no floor or new physical input is used.
The integrated moments have units deg⁻² (km/s)⁻¹, so A has units deg² and
A B has units deg² km/s, matching P0 and v/N. The saved P0 and B already include
the appropriate instrumental smoothing, as established in the W02 source
assessment; B is held at its saved weighting mode and is not smoothed again.
This is a dimensional/response-ownership check, not a new model evaluation.

The fixed signal is P0 for every executed update. Refreshed S0 is
2.199951654746474 and S23 is 2.057085625889348 deg² km/s. All S_t for t=0..23
are saved, alongside independent 80-digit values. The refreshed A_t B/P0 is
0.197946875842, 0.196736341780, 0.196713643312 and 0.196713635541 at
checkpoints 3,6,12,24 respectively. The last ratio characterizes the final
state; no update from t=24 was performed.

## Weight amplitudes, signs and arithmetic

| Branch | t | max(abs(w)) | I1 |
| --- | ---: | ---: | ---: |
| Fixed | 3 | 0.991089953400174 | 0.000701684611540801 |
| Fixed | 6 | 0.991464603678096 | 0.000711368652103311 |
| Fixed | 12 | 0.991476463994446 | 0.000711682704738839 |
| Fixed | 24 | 0.991476475513314 | 0.000711683009979995 |
| Refreshed | 3 | 0.993190554981609 | 0.000762801495437410 |
| Refreshed | 6 | 0.993401619144651 | 0.000770023770181012 |
| Refreshed | 12 | 0.993405551073653 | 0.000770160593852928 |
| Refreshed | 24 | 0.993405552419528 | 0.000770160640701381 |

All checkpoint weights, moments and coefficients are positive and finite.
Every update passes checks for positive N/S, nonnegative coefficients, defined
normalization, finite arithmetic and nonzero denominators. No domain/arithmetic
obstruction occurred in either float64 or Decimal. Output-only signed shapes
w/max(abs(w)) are saved; they never enter an update. Both trajectories retain
48 explicit update-attempt records in total and their final finite states.

## Independent checks and input/source identity

One-bin analytic control: with rho*dm=L=lp=P0=B=v=1, A=1 and N=w.
The intrinsic map is w/(1+w), giving 1/3 from w0=1/2. The refreshed map is
2w/(1+2w), leaving 1/2 fixed with S=2. Both float updates match.

Two-bin analytic control: from (1/2,1/5), A0=29/49 and S0=78/49 give
w1=(39/74,39/179). Writing the components over a common denominator gives
A1=(179²+74²)/253²=37517/64009 and S1=101526/64009. Substituting
N1=39*253/(74*179) yields
w2=(1979757/3655376,1979757/8682233). Exact Fraction arithmetic verifies these
expressions and explicitly distinguishes w2 from freezing S at 78/49.
Both float updates and their signals agree with the rational values.

An independent scalar loop uses 80-digit Decimal arithmetic, converting exact
saved binary floats with `Decimal.from_float`, forming direct scalar sums and
recomputing its own signal before each update. It calls neither the float
moment nor update function. It checks the seed, all eight checkpoint weight
vectors and [I1,I2,I3,A,P_pixel], and all 48 executed signals. Maximum relative
discrepancy is 1.2721e-15 for fixed and 1.0940e-15 for refreshed. No precision
increase was needed. Finite operands and rtol=5e-12, atol=0 are required,
including exact handling of reference zeros.

The authoritative source is
`.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json`,
`settings.pair_inputs["lya(qso)_lya(qso)"]`. Its SHA-256 is verified:
`3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6`.
Compatibility profile, bin 0, z=[2.0,2.235], field selection, finite arrays,
107 nodes and nine negative densities are explicitly checked. Input arrays and
seed match W03 exactly. Fixed t=3/6 weights and all moments/coefficients pass
the W03 NPZ replay and the assignment's printed coefficient references.

New evidence directory:
`.validation/forest-weight-diagnostics/w04-r1-20260915T225726Z/`.
`arrays.npz` contains inputs, seed, checkpoint weights, output-only shapes,
moments/coefficients, float-converted Decimal weights and final finite states.
`summary.json` contains held scalars, identity, every signal/update attempt,
Decimal coefficient/signal values, differences, both stability assessments,
stop reason, command, versions and source hashes. No generalized validation
framework or additional test file was introduced.

Relevant SHA-256 identities:

```text
1adfa3354fc887d5008c3a37596b494ba8585c774a9866aaa366bdf8a462b788  scripts/compare_forest_weight_signals.py
c5b16414ed881e9a46f1f423c6aa5ed9b290e684912a588907492e6463918bf5  W03 arrays.npz
41a9a4097bdfb51364e1d0ac24f036a4b71d888a25250f18fa1dd6057737e979  W03 summary.json
829450f5f01d1faddae417f9544e4b6cf6c31141ab7306952864fc91dac1479e  scripts/compare_forest_weight_integrals.py
2dd75313f84c404a159bb16872b567d3fece57207583698dd61ae0516c55eb02  WEIGHTING_DIAGNOSTIC_STEP.md
3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310  IMPLEMENTATION_STEP.md
```

## Actual execution and disposition

From `lib/fishhighz`, using the existing Python 3.13.15 / NumPy 2.5.3 interpreter:

```bash
.venv/bin/ruff format scripts/compare_forest_weight_signals.py
.venv/bin/ruff check scripts/compare_forest_weight_signals.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -B scripts/compare_forest_weight_signals.py
.venv/bin/ruff format --check scripts/compare_forest_weight_signals.py
.venv/bin/ruff check scripts/compare_forest_weight_signals.py
git diff --check
```

The single numerical invocation passed all assertions in 0.4003 seconds, using
one process and one numerical thread with a 30-second alarm. The alarm raises
an exception so partial arrays/summary are written if the cap is reached.
Ruff lint/format and `git diff --check` pass. A subsequent one-thread inline
Python command only read/printed the saved JSON fields to prepare this report;
it did not rerun either recurrence. No failed numerical invocation occurred.

Only the standalone W04 script, this report and the new small output directory
were added. Historical dirty/untracked work, W01–W03 scripts, production code,
planning files and IMPLEMENTATION_STEP.md were preserved. No broad pytest,
wheel/environment, physical-input preparation, forecast, Fisher sum, other
population/bin, prefix variant, magnitude refinement, Slurm action, agent
dispatch, commit or push occurred.

The remaining scientific question is whether these finite-count plateaus
correspond to well-defined physically admissible limits; this signed-grid
experiment does not establish that. No production prescription or subsequent
assignment is selected. Stop for independent/user review.
