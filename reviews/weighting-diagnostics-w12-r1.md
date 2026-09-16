# W12 revision 1 implementation handoff

Status: implemented and ready for independent/user review. This is not scientific
acceptance, a full-profile convergence result, or authorization to proceed to W13.

## Selected convention and implemented identities

The accuracy recipe now defaults to fixed inverse-variance forest weighting at
the user-selected

```text
q_star = 0.00035 s/km
B_star(field,z) = P1D_intrinsic(z,q_star) W_field(q_star)^2
nu_i = B_star/(B_star + l_p v_i).
```

`AccuracyRecipe.prepare` constructs each `ForestInput` through one shared
host-side helper. It passes `method="inverse_variance", alias=B_star`, with no
signal, iteration count, auxiliary coordinate or P3D query. Settings and forest
provenance record the method, q_star/units, intrinsic reference P1D, field response,
B_star/units, and `iterations={applicable: false, status: inapplicable}`. The
profile cache key includes the method. During mean differentiation, the prepared
weights/noise/covariance remain fixed under the existing `prepare_bin` contract.

The final mode-dependent forest noise remains

```text
N_F(q) = [A P1D_intrinsic(q) W_F(q)^2 + P_pixel] d_deg^2/a_v.
```

An explicit `weight_method="legacy"` accuracy route retains the historical P3D/
P1D auxiliary sampling, `(2.4,0.00035)` coordinates and explicit cumulative
iteration count. The compatibility recipe, weighting/noise/covariance/Fisher
kernels and public core weighting defaults were not changed.

## Controller and historical semantics

The new fixed-reference trial contract is version 2. Its ordered metrics are
`k`, `mu`, `magnitude`, `volume`, `step`, and `combined`; neither an iteration
trial nor a `weights` metric is present. Combined lower controls vary only those
five applicable numerical families. The contract binds
`method="inverse_variance"` and the fixed q_star convention to the report.

Version 1 retains exactly the historical six-control cumulative schedule plus
combined. A version-1 report without `forest_weighting` metadata is interpreted
as legacy. Focused negative controls reject (1) relabelling a version-1 legacy
report as fixed-reference and (2) removing magnitude, or any other still-required,
version-2 family. Cached accuracy primaries must match the requested method;
a cumulative cached primary cannot be reused as fixed-reference. `three_weights`
is scheduled only for explicit legacy accuracy and direct fixed-reference
invocation rejects it as inapplicable.

The saved bin-0 record validates with `require_pass=False` under version-1
semantics and remains `passed=false`. Its unresolved historical result was not
altered. The live compatibility source, noise/covariance/Fisher sources and all
relevant numerical kernels retain their saved r5 hashes. `weights.py` instead
matches the later reviewed W09/W08-r2 source hash
`5e9110b63d13eded2cb84dbd0e82b958418ed59602d910bfeb1428a27d67f663`;
the W08 range repair is not falsely described as an r5 byte match.

## Saved-input FF/FG/GG calculation

The calculation used Git HEAD `0d69786a06d5d676564a51fad14a7156951c2228`
and verified these immutable inputs before and after evaluation:

| Input | SHA-256 |
| --- | --- |
| `profiles-checked/records-001.report.json` | `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666` |
| `profiles-checked/records-001.npz` | `e48c71b7a34abeadda3f8c9fc1dd6ee5bb3a7b13476f68626b0b252eacc2c4c0` |
| `weight-diagnosis/diagnosis.json` | `069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008` |
| `weight-diagnosis/bin-0.npz` | `d1dc7a4d65f73d23bf3cf1c1aba47d2861401fd39215717d1ae9260bff446f73` |

Metadata selected accuracy bin 0 `[2,2.235]`, fields F=`lya(qso)` and G=`qso`,
full-array pair indices `(0,1,5)` for FF/FG/GG, and parameters
`(ap_0,at_0)`. The calculation remapped these to `(0,0),(0,1),(1,1)` before
constructing a new covariance. Saved signal, Jacobian, k, mu, modes and GG noise
were held fixed; FG noise was exactly zero. No P3D, CAMB, reader, derivative,
new grid or physical-volume evaluation was performed.

The actual profile helper gave, at both forest magnitude orders,

```text
P1D_intrinsic(q_star) = 16.362615788177436 km/s
W_F(q_star)           = 0.9998207087918288
B_star                = 16.356748967852234 km/s
```

The public covariance agreed exactly with an independently assembled 3x3 Wick
array. In particular, the FG variance obeyed
`M C_FG,FG = T_FF T_GG + T_FG^2`, with no extra factor of two. Public Fisher
contractions agreed with direct batched `numpy.linalg.solve` results. Maximum
relative residuals were `8.61e-15` for joint Fisher and `8.52e-15` for individual
Fisher; the cross-variance identity agreed within `2.12e-16`.

### Order 32

```text
A       = 0.029419845663300446 deg^2
P_pixel = 0.31435492218062416 deg^2 km/s

F_joint = [[2298.6650788588945, 1279.4998269885684],
           [1279.4998269885684, 4120.953027390565 ]]
sigma_joint = [0.022933139181031712, 0.017127837096461408]

F_FF = [[1419.738166245999 , 667.4320627023011],
        [ 667.4320627023011, 1591.527428282942 ]]
sigma_FF = [0.02961952211255682, 0.027975323168394888]

F_FG = [[1785.185998454070 , 968.4203279473514],
        [ 968.4203279473514, 2935.798406105992 ]]
sigma_FG = [0.026119928157594924, 0.020368103767303074]

F_GG = [[763.5144269090425 , 451.36081791688116],
        [451.36081791688116, 1636.8974279164906 ]]
sigma_GG = [0.03955767247696439, 0.027016471119289617]
```

The FF errors reproduce the W09 reference values
`(0.029619522112556678,0.027975323168394890)` within `4.81e-15` relative.
Each individual result uses its own marginal spectrum variance, not a parameter-
covariance subblock from the joint calculation.

### Order 64

```text
A       = 0.029419845663300383 deg^2
P_pixel = 0.3143549221806231 deg^2 km/s

F_joint = [[2298.6650788588995, 1279.499826988571 ],
           [1279.499826988571 , 4120.953027390571 ]]
sigma_joint = [0.022933139181031688, 0.017127837096461394]

sigma_FF = [0.02961952211255679,  0.027975323168394853]
sigma_FG = [0.026119928157594903, 0.020368103767303046]
sigma_GG = [0.03955767247696439,  0.027016471119289617]
```

Order-64/order-32 minus one is:

| Quantity | Fractional change |
| --- | ---: |
| A | `-2.1094237467877974e-15` |
| P_pixel | `-3.3306690738754696e-15` |
| joint `(parallel,transverse)` errors | `(-1.1102230246251565e-15,-7.771561172376096e-16)` |
| FF errors | `(-9.992007221626409e-16,-1.2212453270876722e-15)` |
| FG errors | `(-7.771561172376096e-16,-1.3322676295501878e-15)` |
| GG errors | `(0,0)` |

Every assigned refinement contrast is below `0.001` in absolute value.

## Changed files

- `fishhighz/validation/accuracy.py`: default fixed-reference forest preparation,
  per-forest convention metadata, explicit legacy route and legacy-only diagnostic.
- `fishhighz/validation/study.py`: method-dependent applicable control schedule.
- `fishhighz/validation/trials.py`: narrow version-2 contract and version-aware replay.
- `fishhighz/validation/schema.py`: version-aware complete-family validation.
- `fishhighz/validation/profiles.py`: method-bound cache reuse and diagnostic routing.
- `scripts/check_accuracy_fixed_reference.py`: bounded saved-input FF/FG/GG check.
- `tests/test_weighting_w12.py`: self-contained profile/controller/cache fixtures and
  the two assigned negative controls.
- `README.md`: concise fixed-reference convention and v1/v2 semantics. This file
  was already modified before W12; only the accuracy-method paragraphs were added.
- `reviews/weighting-diagnostics-w12-r1.md`: this handoff.

All other dirty/untracked files and historical evidence were preserved. No core
kernel, compatibility recipe, physical policy, tolerance, INI/CLI, plan,
manuscript or sibling package was changed.

## Commands and evidence

Saved calculation, one process/thread and 60 s cap:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/check_accuracy_fixed_reference.py
```

Result: passed in `7.599148422013968 s`. Evidence is in
`.validation/forest-weight-diagnostics/w12-r1-20260916T051738345689Z/`:

- `result.json`, SHA-256
  `7adda00c51046cb72ee7358af8a77fb57ad22a1188011d61fd67c365f06053eb`;
- `focused-tests.xml`, SHA-256
  `e3a8ee356606b6139d50c5f45ca415745226eeb40ed529b5e0c0110770609b64`;
- `git-status.txt` records the preserved dirty tree.

Focused tests, one process/thread and 60 s cap:

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

Result: `19 passed in 13.95 s`. These cover W08 fixed-weight arithmetic/noise,
frozen differentiation inputs, W12 construction and cache separation, version-2
controller/replay, exact version-1 legacy behavior and prior per-spectrum binding.

Static checks:

```bash
.venv/bin/ruff check fishhighz/validation/accuracy.py \
  fishhighz/validation/study.py fishhighz/validation/trials.py \
  fishhighz/validation/schema.py fishhighz/validation/profiles.py \
  scripts/check_accuracy_fixed_reference.py tests/test_weighting_w12.py
.venv/bin/ruff format --check fishhighz/validation/accuracy.py \
  fishhighz/validation/study.py fishhighz/validation/trials.py \
  fishhighz/validation/schema.py fishhighz/validation/profiles.py \
  scripts/check_accuracy_fixed_reference.py tests/test_weighting_w12.py
git diff --check
```

Result: Ruff passed; seven files were already formatted; `git diff --check`
passed. An earlier four-file pytest selection reached the 60 s timeout before a
summary, so no pass/fail claim is made from it. The narrower assigned regression
set above completed. The first saved-script attempt stopped before numerics on an
overly strict assertion that live `weights.py` must equal the older r5 producer;
the stopped artifact
`.validation/forest-weight-diagnostics/w12-r1-20260916T051657887523Z/` is
preserved. The final check correctly distinguishes the reviewed W08/W09 weight
source from unchanged historical compatibility/kernel sources.

## Scientific limits and stopping point

The synthetic controller establishes contract behavior, not survey convergence.
The saved calculation establishes local consistency and magnitude refinement for
one bin, two fields, three spectra, a fixed saved signal/Jacobian/grid/mode count,
and fixed quasar noise. It does not establish Fourier, redshift, volume,
derivative-step or full-survey convergence; validate density/SNR or overlap-noise
policies; optimize q_star; accept historical failed forecasts; or resolve package
Step 13. No real-model controller, full forecast, broad suite, installation,
Slurm action, dispatch, commit or push was performed.

W12 revision 1 stops here, ready for independent and user review. W13 has not
been started.
