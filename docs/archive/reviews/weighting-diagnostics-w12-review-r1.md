# Scientific finding and effect of possible corrections

**W12 revision 1 passes independent scientific review in its assigned bounded
scope. No scientifically consequential correction is needed.** The accuracy
profile now routes forest preparation through fixed inverse-variance weights at
the user-selected baseline `q_star=0.00035 s/km`. Its new refinement contract
correctly treats the weight-iteration family as inapplicable rather than as a
repeated calculation that appears converged. The historical cumulative contract,
including its weight family and failed saved bin-0 verdict, remains intact.

A correction to the reference P1D/response ownership, method/cache identity,
trial-family interpretation, cross-spectrum Gaussian variance, or saved pair and
parameter selection could have changed this conclusion. Focused reruns, source
inspection and an independent raw-array reconstruction found none. In particular,
the cross variance contains `T_FF*T_GG + T_FG**2` with no extra factor of two,
and the independently summed joint Fisher matrix agrees with the fresh W12
calculation to `9.27e-15` relative or better. No revision of
`WEIGHTING_DIAGNOSTIC_STEP.md` is warranted.

This result supports local adoption of the declared convention and roundoff-level
magnitude refinement in the assigned saved bin. It is not full-survey convergence,
an optimization of `q_star`, validation of the inherited physical input policies,
acceptance of historical failed forecasts, or authorization for W13.

## Profile, metadata and cache routing

The live host-side route has the required ownership:

- `fishhighz/validation/accuracy.py:38-113` evaluates intrinsic P1D at the fixed
  reference mode, applies the field response exactly once in
  `B_star=P1D(q_star)*W_field(q_star)**2`, and passes direct `alias=B_star` with
  neither iterations nor auxiliary P3D coordinates. The default recipe method is
  inverse variance at lines 206-233.
- `AccuracyRecipe.prepare` rejects an iteration control for the fixed method,
  includes the method in its prepared-object cache key, and records the common
  convention plus each field's P1D, response and B_star at lines 489-605. The
  focused two-forest fixture obtains distinct B_star values for distinct responses
  and makes no auxiliary P3D call.
- `fishhighz/validation/profiles.py:95-147` refuses cached accuracy results whose
  recorded method differs from the requested route, then retains the pre-existing
  scientific-source/reference/resource identity checks. The full-profile driver
  defaults to inverse variance at lines 186-242; an explicit legacy route remains.
- `three_weights` is excluded from fixed-reference diagnostic requests at
  `profiles.py:215-223` and rejects direct fixed-reference use at
  `accuracy.py:679-705`. It remains a truthful three-iteration cumulative
  diagnostic for explicit legacy accuracy.

The fresh saved-input calculation rechecked the historical r5 hashes for the
compatibility recipe and noise/covariance/Fisher kernels. The live weight source
instead matches the reviewed W08-r2/W09 SHA-256
`5e9110b63d13eded2cb84dbd0e82b958418ed59602d910bfeb1428a27d67f663`,
so the later range repair is not misrepresented as an r5 byte match. Fixed-noise
mean differentiation is covered by the focused frozen-provider regression.

## Refinement contracts and historical semantics

`fishhighz/validation/study.py:22-27` removes `iterations` before any fixed-reference
trial. Its fixed route executes only k, mu, magnitude, volume and derivative-step
families plus the combined lower-control comparison; the legacy-only iteration
schedule remains at lines 92-129 and 185-212. There is no duplicate fixed-weight
baseline standing in for convergence.

`fishhighz/validation/trials.py:38-72` binds these routes to distinct contracts:
version 1 retains the six cumulative controls, while version 2 contains the five
applicable controls and the exact fixed-reference convention. Validation at lines
75-106 rejects a version-1 report relabelled as inverse variance and requires
version-2 method/reference/inapplicable-iteration metadata. Lines 107-302 retain
the saved trial identities, operands, failed attempts, replayed metrics and final
verdict checks. `schema.py:396-413` requires all six version-1 families or all
five version-2 families plus `combined`; it does not make arbitrary old families
optional.

The immutable saved bin-0 report independently remains:

```text
trial contract: version 1
metrics: k, mu, magnitude, volume, step, weights, combined
successful weight levels: 3, 6
final iterations: 6
passed: false
```

Thus fixed-reference iteration convergence is not fabricated, and the historical
failed verdict is not promoted or rewritten. The synthetic version-2 controller
pass establishes contract behavior only; it supplies no survey-convergence claim.

## Independent saved-input reconstruction

All four assigned input SHA-256 identities match. Metadata selects accuracy bin 0
`[2,2.235]`, fields `lya(qso)` and `qso`, full-array indices `(0,1,5)` for FF/FG/GG,
and parameter order `(ap_0,at_0)`. Saved signal, Jacobian, Fourier nodes, mode
counts and GG noise are held fixed; FG noise is exactly zero.

Without importing FishHighz weight, P1D, response, noise, covariance or Fisher
helpers, the independent reconstruction gives

```text
P1D(q_star) = 16.362615788177436 km/s
W_F(q_star) = 0.9998207087918288
B_star      = 16.356748967852234 km/s

order 32:
A           = 0.029419845663300394 deg^2
P_pixel     = 0.3143549221806231 deg^2 km/s
F_joint     = [[2298.665078858900, 1279.499826988574],
               [1279.499826988574, 4120.953027390534]]
sigma_joint = [0.022933139181031716, 0.017127837096461495]
```

For example, the independently reconstructed selected covariance at saved node
8192 is

```text
[[     68.07833824486617,    -139.88253822770648,       287.42071274780756],
 [   -139.88253822770648,   24072.76556337791,       -98335.44901484299   ],
 [    287.42071274780756,  -98335.44901484299,     33643575.790014274     ]]
```

The complete 16,384-node selected covariance was reconstructed, not sampled for
the Fisher sum. Its order-32 SHA-256 is
`1fc7c1b26660fddc196c40d9403557026d7664a34945c909e46bc9eb322b0264`.
The maximum cross-variance identity residual is `2.03e-16`. Direct block solves
followed by independently accumulated node sums reproduce the implementation's
joint Fisher matrix within `7.50e-15` at order 32 and `9.27e-15` at order 64.

The independent order-32 individual error vectors are

| Spectrum | sigma_parallel | sigma_transverse |
| --- | ---: | ---: |
| FF | `0.029619522112556654` | `0.027975323168394847` |
| FG | `0.026119928157594930` | `0.020368103767303090` |
| GG | `0.039557672476964570` | `0.027016471119289592` |

The FF result reproduces the W09 reference. Each row uses its own marginal
spectrum variance, not a subblock of the joint parameter covariance.

Operation-order differences between scalar `math.fsum` and the public NumPy path
are confined to roundoff. The order-64/order-32 minus-one contrasts are:

| Quantity | Fresh public W12 | Independent raw-array reconstruction |
| --- | ---: | ---: |
| A | `-2.109e-15` | `+2.220e-16` |
| P_pixel | `-3.331e-15` | `0` |
| joint errors | `(-1.110e-15,-7.772e-16)` | `(0,+2.220e-16)` |
| FF errors | `(-9.992e-16,-1.221e-15)` | `(-1.110e-16,-1.110e-16)` |
| FG errors | `(-7.772e-16,-1.332e-15)` | `(0,0)` |
| GG errors | `(0,0)` | `(0,0)` |

All are far below the assigned absolute `0.001` refinement threshold. This is a
forest-magnitude result with fixed saved information, not evidence for the other
full-profile convergence families.

## Execution, evidence and disposition

The exact focused regression selection from the handoff was rerun with one
thread/process and its 60 s cap:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B -m pytest -q \
  tests/test_weighting_w12.py \
  tests/test_inverse_variance_weights.py::test_rational_coefficients_and_legacy_seed \
  tests/test_inverse_variance_weights.py::test_mode_dependent_noise_scalar_response \
  tests/test_inverse_variance_weights.py::test_survey_equivalence_and_frozen_provider_counts \
  tests/test_step12_performance.py::test_independent_accuracy_check \
  tests/test_step12_performance.py::test_study_reuse_controls_failures_and_order \
  tests/test_step12_revision5.py::test_distinct_controls_identical_information \
  tests/test_step12_revision5.py::test_failed_attempt_cannot_be_suppressed \
  tests/test_step12_revision5.py::test_rank_and_roundoff_trial_controls \
  tests/test_step12_revision6.py::test_valid_weak_converged_roundoff_and_unconverged \
  tests/test_step12_revision6.py::test_null_and_partially_constrained_trials \
  tests/test_step12_revision2.py::test_reuse_rejects_scientific_changes_and_preserves_producer
```

Result: `19 passed in 14.30 s`.

The assigned saved calculation was rerun unchanged:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/check_accuracy_fixed_reference.py
```

It passed and wrote
`.validation/forest-weight-diagnostics/w12-r1-20260916T053232979880Z/result.json`
(SHA-256 `a2e43830bad420f9317730c34ebb969b8b4b6649b1328bc7cac7b2f208b38b4c`).

The independent check and full-precision output are
`.validation/forest-weight-diagnostics/w12-review-r1-20260916T053314Z/check.py`
and `result.json`, with SHA-256 values
`69356040f3f6338ff9a1610d1478f01b7515314ccb05a9bd9ef1e8f49c735535`
and `e471a73cd2fc9f39000ef2312200c6e42ebace4889085dbe1f68ab89d9059676`.
An initial independent invocation stopped before producing evidence because the
review script had transcribed the resolving power as 3200; correcting it to the
saved 15x2pt recipe value 2500 yielded the final passing calculation above. No
W12 implementation source was changed.

Focused Ruff lint and format checks pass for the seven W12 source/test/script
files and the independent check; `git diff --check` also passes. The only review
changes are this report and the ignored independent evidence directory. Existing
dirty/untracked work and historical artifacts were preserved.

No real-model controller, reference scan, full survey run, broad pytest suite,
production repair, plan revision, sibling edit, Slurm action, dispatch, commit,
push or W13 planning occurred. W12 revision 1 is ready for user review; the user
retains implementation acceptance and all further execution decisions.
