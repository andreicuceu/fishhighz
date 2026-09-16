# W08 revision 1: explicit fixed inverse-variance weights

**`method="inverse_variance", alias=B_star` now prepares
nu=B_star/(B_star+l_p v) once on positive density support, with zero weights
elsewhere.** The focused checks establish the analytic coefficients, agreement
with independently supplied weights and the legacy zero-iteration seed,
retained arithmetic rejection, and fixed weight/noise ownership in synthetic
survey differentiation. All **109 focused tests pass**. One saved order-32
QSO-reference coefficient check agrees with independent scalar 80-digit Decimal
arithmetic to **3.18e-15** relative error. Implementation is ready for
independent/user review; no profile default is changed.

## Changes and scientific contract

Only existing `fishhighz/weights.py` and the relevant README section were edited.
Added `tests/test_inverse_variance_weights.py`, this report and the small evidence
directory linked below. No additional numeric helper was needed.

The new branch requires positive finite scalar `alias=B_star` in km/s, already
including the response squared at the caller's fixed weighting mode. `weights`,
`iterations` (including zero), `signal` and `auxiliary` must be `None`. It makes
no provider calls and never calls `_iterate`. Zero variance on positive support
gives exactly one. The existing `ForestWeights` stores the method and alias,
`None` for signal/iterations/auxiliary, an empty `weight_changes` array, owned
immutable arrays and the existing context tuple.

The common input validation, quadrature-mass checks, `_integrals`, coefficient
checks and immutable construction remain unchanged. The branch rejects zero or
nonfinite weights on positive support; raised arithmetic exceptions reject
overflow and the existing integral-product underflow cases. No rescaling,
flooring or relaxed tolerance was introduced. `_iterate` and `_integrals` are
byte-for-byte unchanged, as are the legacy/supplied branches.

`ForestInput.weight_options` already carries the direct scalar option through
`prepare_bin`; its existing validation rejects auxiliary coordinates for this
method. No survey, forecast, noise, profile, reader or provider implementation
needed modification. There was no out-of-scope obstruction.

At each forecast mode q=k mu/a_v, noise still uses

```text
N_F = [A P1D(q) W(q)^2 + P_pixel] d_deg^2/a_v.
```

Here P1D is intrinsic and mode dependent; the fixed B_star is used only to
prepare weights. The response acts once on aliasing, and pixel noise remains
unsmoothed. Weights, coefficients and fiducial noise stay fixed during mean-model
differentiation. The README states the conditional minimum of
Q=A B_star+P_pixel for nonnegative measure, common L/l_p and independent sightline
noise under the stated diagonal-covariance approximation. No multi-mode BAO
optimality follows.

## Focused tests

The rational control uses q=(1/2,2), rho=(2,1/2), v=(1,4), L=l_p=B_star=1.
Independent Fraction arithmetic gives r=(1,1), nu=(1/2,1/5), I1=7/10,
I2=29/100, I3=41/100, A=29/49, P_pixel=41/49 and Q=10/7. The new method agrees
with these values, independently supplied nu and unchanged legacy count zero.
An invocation trap verifies that the new branch calls neither `_iterate` nor
the provider/auxiliary helpers.

Additional tests cover both QSO and LBG identities, zero-density/zero-variance
support, immutable input snapshots/results/metadata, field/geometry/response
mismatch, conflicting settings, invalid scalar B_star, computed-weight underflow,
denominator/product overflow, quadrature-mass loss and mixed integral-product
underflow concealed by another positive contribution. Existing cumulative
counts 0/1/2/3/6 and weight/range/model/noise regressions pass.

The direct noise check uses two distinct velocity wavenumbers, 0.003 and
0.02 s/km, and intrinsic P1D values 4 and 9 km/s, distinct from B_star=2.
A scalar sin(x)/x and exponential response calculation independently checks
aliasing and unsmoothed pixel noise; supplied-nu noise agrees.

The existing tiny synthetic forest fixture is used in a new equivalence test.
Each method performs exactly one P3D call and one P1D call during preparation,
with no auxiliary query or extra weight-generation call. The new and supplied
methods agree in weights, all prefixes, A, P_pixel, noise and the synthetic
information matrix. During `run_bin`, traps reject repeated weight, response,
noise or P1D preparation. Each bin makes five mean-model calls for the existing
two-parameter central differences; the two preparations plus two runs make
12 P3D calls in total, while P1D remains at two. Byte snapshots and coefficient
comparisons verify fixed weights/noise. The existing legacy/supplied
`test_generated_spies` cases also pass. No complete forecast test suite was run.

## Saved order-32 coefficient check

Evidence directory:
[W08 evidence](../.validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/).
The [check script](../.validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/check_saved.py)
and [result](../.validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/saved-result.json)
record the calculation and exact input/source identities.

Only `.validation/step12-r5-20260914T191855Z/profiles-checked/records-001.report.json`
was loaded. SHA-256 matched before and after the successful calculation:
`9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666`.
The report identifies accuracy bin 0, bounds [2.0,2.235], forest `lya(qso)`,
physical model `lya` and background `qso`. The selected array mapping is
`settings.samples["lya(qso)"]`. Magnitude order is 32 per interval: 5,184 nodes
in 162 intervals, each checked against its partition bounds.

The calculation uses the saved density, quadrature, variance and
z_source=2.3830099582078716, with rho_v=density*(1+z_source)/299792.458,
L=44147.09373976799 km/s, l_p=63.328211161339375 km/s and
B_star=16.356748967852234 km/s. Density/SNR policies were already applied to
these saved arrays; they were neither reapplied nor physically validated.

The coefficient-only context is explicitly synthetic: z_eval=2.1152848986890427,
h_fid=0.7, constant H=200 and transverse distance=1000, area=1 deg², z_order=2,
and response `(pixel_width_velocity=l_p, gaussian_sigma_velocity=0)`.
These backgrounds and sigma are not a reconstruction of the historical
geometry/instrument. Only L/l_p and the arrays enter A/P_pixel, while B_star is
supplied already smoothed. No raw reader, interpolation, cosmology adapter,
P3D/P1D provider or model was called by this check.

The oracle converts each saved binary float exactly to Decimal, independently
forms the density conversion and nu, and sums all three moments in scalar
80-digit arithmetic. It calls no implementation weight or integral helper.
Every prefix and weight is compared, as well as the coefficients. Finite values
are required with rtol=5e-12, atol=0 and exact reference-zero handling.

| Quantity | New method | Independent Decimal, rounded |
| --- | ---: | ---: |
| A [deg²] | 0.029419845663300446 | 0.029419845663300400 |
| P_pixel [deg² km/s] | 0.31435492218062416 | 0.31435492218062315 |
| Q [deg² km/s] | 0.7955679523681858 | 0.7955679523681840 |

Maximum relative discrepancies are 2.85e-16 for weights, 2.16e-15 for prefixes,
3.18e-15 for coefficients and 1.54e-15 against the W06 printed reference
coefficients. No tolerance was relaxed.

## Commands, failures and preservation

Commands from the package root, using the existing `.venv` interpreter:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 .venv/bin/python -B -m pytest -q tests/test_inverse_variance_weights.py tests/test_weights.py tests/test_weight_range.py tests/test_weight_models.py tests/test_noise.py tests/test_forecast.py::test_generated_spies
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B .validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/check_saved.py
.venv/bin/ruff check --no-respect-gitignore fishhighz/weights.py tests/test_inverse_variance_weights.py .validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/check_saved.py
.venv/bin/ruff format --check --no-respect-gitignore fishhighz/weights.py tests/test_inverse_variance_weights.py .validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/check_saved.py
git diff --check
```

The sole focused pytest invocation passed **109 tests in 10.62 s**;
[pytest.log](../.validation/forest-weight-diagnostics/w08-r1-20260916T002747Z/pytest.log)
is retained. The successful saved check took **0.347 s** internally. Both
numerical checks used one process and numerical thread limits of one on login32,
Python 3.13.15 and NumPy 2.5.3. Neither cap was reached.

The first saved-check invocation failed after 0.249 s: the new check script
passed the serialized key `physical` directly to `ObservedField`, whose API uses
`physical_model`. The script's explicit field mapping was corrected; its failed
result/log remain as `saved-result-attempt1.json` and `saved-check-attempt1.log`.
This occurred before coefficient preparation and required no production change.
Initial Ruff findings were test import ordering and a lambda assignment in the
check script; both were corrected. Final Ruff lint/format and `git diff --check`
pass. There are no remaining numerical failures or obstructions.

The supplied home path resolves to the requested CFS checkout. HEAD remains
`0d69786a06d5d676564a51fad14a7156951c2228`; the implementation resides in
pre-existing uncommitted/untracked work, so HEAD alone is insufficient identity.
`before.json` and `preservation.json` record a 180-file initial snapshot across
source, tests, reports and planning documents. Only `fishhighz/weights.py` and
`README.md` changed among those files; initial copies of both are retained to
separate this change from earlier work. Final decisive source hashes include:

```text
36d4d489d5e268617c4cac782c76aaa7da85c614a4b456b9d5e2954021d51a16  fishhighz/weights.py
3867c9eee6d43f22fca878b16c9f38bafc65a55061355d79828fcb135cc5e66f  fishhighz/kernels/weights.py (unchanged)
648145349a639c4362983015c8645b40353527e46ae2af421f12e65069a88c38  tests/test_inverse_variance_weights.py
3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310  IMPLEMENTATION_STEP.md (unchanged)
```

The result JSON records further noise/survey/forecast/README/assignment and check
script hashes. Tracked `git diff --check` does not cover untracked additions;
Ruff explicitly checked both new Python files, including the ignored evidence
script. Planning/governance files and unrelated dirty/untracked work are preserved.

## Limits and review boundary

[W07's BAO results](weighting-diagnostics-w07-r1.md) remain historical evidence:
for its one accuracy QSO auto-spectrum/bin, cumulative t=3 and t=6 errors were
larger than the reference. W08 has measured no new survey improvement and has
changed no profile result. It repeats neither historical trajectories nor Fisher
comparisons. The synthetic information equivalence tests only the unchanged
orchestration of the explicit option.

There is no per-mode reoptimization, new physical input policy, automatic P1D
sampling or profile-default adoption. No real forecast, full suite, new
wheel/environment, performance study, W01 repair, package Step 13 work, subsequent
step, Slurm action, agent dispatch, commit or push was performed. Independent/user
review and scientific acceptance remain pending.
