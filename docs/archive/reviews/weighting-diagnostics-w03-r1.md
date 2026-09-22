# Forest-weighting diagnostic W03 revision 1

## Scientific answer and scope

**Replacing magnitude-prefix feedback by the common full-sample integral does
not remove the measured t=3/6 coefficient sensitivity.** On the exact saved
107-node signed QSO-forest example, the absolute fractional change in A increases
from 0.8616% to 0.9196% (and reverses sign); the change in P_pixel increases from
1.1781% to 2.0171%. Both full-sample changes exceed 1e-3 and are independently
confirmed with 80-digit Decimal arithmetic. The prescribed early stop therefore
occurs at t=6 in both branches. **Full t=12 was not attempted.**

This is a finite-count comparison for compatibility bin 0, z=[2.0,2.235],
population lya(qso), at fixed S, B and magnitude sampling. No asymptotic behavior,
continuum limit, optimality, physical admissibility or production prescription
is established. Even a finite fixed-grid plateau would not prove a continuum
limit. Status: ready for independent/user review, not accepted.

The user's explicit W03 r1 approval supersedes the awaiting-approval text in the
assignment. Planning documents were preserved. Package Step 13 and W01 evidence
repair were not assigned or performed.

## Coefficient levels and iteration changes

New float64 results (A in deg²; P_pixel in deg² km/s):

| Feedback | t | A | P_pixel |
| --- | ---: | ---: | ---: |
| Prefix | 3 | 0.022118395206768268 | 0.5480269062497443 |
| Prefix | 6 | 0.022308973469159235 | 0.5544833486755244 |
| Full sample | 3 | 0.022021379679306282 | 0.5298546380249203 |
| Full sample | 6 | 0.021818865253617570 | 0.5405425825641871 |

Signed changes use the earlier checkpoint or matched-count prefix as the
reference. Absolute differences retain the units above.

| Comparison | ΔA | ΔA/A_ref | ΔP_pixel | ΔP_pixel/P_ref |
| --- | ---: | ---: | ---: | ---: |
| Prefix 6 minus 3 | +0.000190578262391 | +0.8616278921% | +0.006456442425780 | +1.1781250797% |
| Full 6 minus 3 | −0.000202514425689 | −0.9196264205% | +0.010687944539267 | +2.0171465478% |
| Full minus prefix, t=3 | −0.000097015527462 | −0.4386191971% | −0.018172268224824 | −3.3159445307% |
| Full minus prefix, t=6 | −0.000490108215542 | −2.1969106567% | −0.013940766111337 | −2.5141902177% |

The full-sample coefficients are lower at each matched count, but their
iteration sensitivity is larger. These are separate observations.

## Held and changed quantities

The only branch-dependent operation is selecting `cumsum((rho*w)*dm)[-1]`
instead of the full prefix array for feedback. The same function then forms
N=I1*(L/lp) and w_new=S/(S+v/N). Each branch has a separate evolving weight
array copied from the common seed `(B/lp)/(B/lp+v)`. The input arrays are shared
read-only; an assertion checks that the evolving arrays do not share storage.
There is no inter-update rescaling, density floor, support change or aliasing
addition to S. The existing response is already represented in S, B and v.

| Fixed scalar | Value | Units |
| --- | ---: | --- |
| dm | 0.100471698113207 | mag |
| L | 44148.20436604321 | km/s |
| lp | 63.32980433473595 | km/s |
| S | 1.7189455896517238 | deg² km/s |
| B | 16.333147249645222 | km/s |

All 107 magnitude nodes from 16.1 to 26.75 retain the endpoint-inclusive
rectangular dm. Density is per deg² per km/s per mag; v is dimensionless.
The nine negative-density indices (zero-based) are
0, 2, 3, 76, 77, 80, 81, 84, 85. All are retained literally.

For both branches, coefficient extraction uses the last cumulative element of
`(rho*w)*dm`, `(rho*w**2)*dm`, and `((rho*w**2)*v)*dm`, preserving W01's product
and summation order. With these full-sample moments, A=I2/(L I1²) and
P_pixel=lp I3/(L I1²). I1/I2/I3 have units deg⁻² (km/s)⁻¹. No sum replacing the
cumulative reduction, new physical scale, signal or model evaluation was used.

## Weight amplitude and shape

| Feedback | t | max(abs(w)) | I1 | Negative weights |
| --- | ---: | ---: | ---: | ---: |
| Prefix | 3 | 0.7640425843083741 | 0.000521309596387336 | 3 |
| Prefix | 6 | 0.7309209913729792 | 0.000495570419532441 | 0 |
| Full sample | 3 | 0.9910899534001736 | 0.0007016846115408007 | 0 |
| Full sample | 6 | 0.9914646036780957 | 0.0007113686521033106 | 0 |

The signed relative shapes w/max(abs(w)) are saved in the checkpoint NPZ; they
are output-only normalizations. Their maximum absolute change between t=3/6 is
0.34769487588670567 for prefix and 0.010735585908529188 for full sample. Thus a
smaller change in this shape norm does not imply smaller coefficient sensitivity.
Full-sample amplitudes rise slightly over these checkpoints; neither branch's
amplitude trend is extrapolated. All evaluated coefficients and moments are
finite and positive. Signed input density nevertheless prevents interpreting
this result as convergence of a physical nonnegative sampling measure.

## New checks and minimal identity

Source: `.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json`,
`settings.pair_inputs["lya(qso)_lya(qso)"]`.
Its hash, compatibility profile, bin, bounds, population selection, 107-node
finite arrays and nine negative densities were checked directly.

- Prefix t=3 weights reproduce the saved `_w_lya`; coefficients reproduce the
  final saved `_aliasing_weights` and `_effective_noise_power` elements.
- Prefix t=3/6 coefficients reproduce the historical reference numbers at
  rtol=5e-12, atol=0, with explicit finite checks. The t=6 differences from the
  printed historical values are at rounding precision.
- The exact two-bin Fraction calculation gives prefix (1/3,7/47) and full
  (7/17,7/47); separate float64 updates match, including the common last bin.
- Independent scalar Decimal recurrence and direct sums at 80 digits reproduce
  all four checkpoint weight vectors, moments and coefficients. Inputs are
  `Decimal.from_float` conversions of the exact saved float64 values (including
  the legacy float64 dm). Largest relative weight discrepancy: 6.25e-16;
  largest moment/coefficient discrepancy: 1.28e-15. No precision increase was
  needed. No singular or nonfinite arithmetic occurred.

New evidence:
`.validation/forest-weight-diagnostics/w03-r1-20260915T224027Z/`.
`arrays.npz` saves the exact inputs, captured references, seed, four checkpoint
weights, signed shapes, and [I1,I2,I3,A,P_pixel] arrays; `summary.json` gives the
numbers, Decimal values, differences, source hashes and stop. Only one standalone
script and this report were added outside that directory.

SHA-256:

```text
3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6  saved report
829450f5f01d1faddae417f9544e4b6cf6c31141ab7306952864fc91dac1479e  scripts/compare_forest_weight_integrals.py
efc3218677af21705a9170e715c4df0e91101fc8a2cadd0fe350c72a6bf1d223  WEIGHTING_DIAGNOSTIC_STEP.md
239ee57d07ff485bc2ac77907ba104b5ed694cf7eb56e1b7b550c3f0847e4ed9  scripts/diagnose_legacy_forest_iterations.py
10f7da1c68eba0cb25a7c77fa1c92b48ae6bc4323eddea680997f29303328431  ../lyaforecast/lyaforecast/weights.py
```

The last two files were inspected read-only. The unmodified package
IMPLEMENTATION_STEP.md retains SHA-256
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.
Extensive initial dirty/untracked package work and historical reports remain.

## Actual execution

From `lib/fishhighz`, using the existing Python 3.13.15 / NumPy 2.5.3 environment:

```bash
.venv/bin/ruff format scripts/compare_forest_weight_integrals.py
.venv/bin/ruff check scripts/compare_forest_weight_integrals.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python scripts/compare_forest_weight_integrals.py
.venv/bin/ruff format --check scripts/compare_forest_weight_integrals.py
.venv/bin/ruff check scripts/compare_forest_weight_integrals.py
sha256sum IMPLEMENTATION_STEP.md WEIGHTING_DIAGNOSTIC_STEP.md scripts/diagnose_legacy_forest_iterations.py ../lyaforecast/lyaforecast/weights.py
git diff --check
```

The numerical script passed all assertions in 0.293 s, below the 30 s cap;
Ruff lint and final formatting checks passed. `git diff --check` passed for the
pre-existing tracked changes; it does not inspect the new untracked files.
An additional one-thread inline NumPy read of the saved NPZ computed
`max(abs(shape_6-shape_3))` and each shape's min/max for the descriptive checks
above; no recurrence was rerun. Its wrapper emitted a MUNGE socket message
after successful output, with exit status 0. No scheduler command was issued.

No full pytest, wheel, environment creation, production-module edit, W01/W02
script/test change, planning-file edit, survey/model preparation, NewForecast,
Fisher sum, another population/bin, magnitude refinement, Slurm, agent dispatch,
commit or push was performed.

## Historical evidence and remaining uncertainty

W01's report/review and saved summary/array inventory were read as historical
context. W01's coefficient replay is newly checked here; its forecast projections,
evidence validator and outstanding R4/R5 findings were neither rerun nor repaired.
W02's passing equation/source audit is reused as history; its literature and C++
search were not reopened. The exact two-bin calculation above is a new check.

The smallest unresolved scientific question is whether the full-sample
coefficients approach a finite limit at this fixed signed grid beyond the
permitted counts. The present early stop does not answer it or authorize further
iterations. Magnitude-continuum behavior and an aliasing-inclusive signal are
also outside this result. No next step is prepared. Stop for independent/user
review.
