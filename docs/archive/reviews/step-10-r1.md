# Step 10 handoff — revision 1

Implemented `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-13), following the
user's explicit dispatch. Ready for user and independent review; acceptance and
progression are not claimed. No Step 11 work was started. Base commit remains
`0d69786a06d5d676564a51fad14a7156951c2228`.

## Changes and preservation

Added:

- `fishhighz/weights.py`: normalized density conversion, auto-only auxiliary
  callable sampling, immutable forest inputs/results and explicit preparation.
- `fishhighz/kernels/weights.py`: array-only cumulative iteration and integrals.
- `fishhighz/noise.py`: forest aliasing/pixel terms, local galaxy density,
  Poisson noise and complete required-pair noise packing/PSD validation.
- `tests/test_weights.py`, `tests/test_noise.py`, `tests/test_weight_models.py`:
  73 new synthetic acceptance tests.
- `examples/weighted_noise.py`: two distinct forests and a galaxy, using fixed
  generated noise with an independent covariance/amplitude Fisher oracle.
- This report. README.md documents signatures, units, shapes and limitations.

Evidence is in `.validation/step10-r1-implementation/` (`$RUN` below).
`before.json` and `git-before.txt` retain pre-edit hashes and Git state for all
80 maintained preexisting files. `preservation.json` confirms 79 are unchanged;
only README.md changed. Existing scientific modules/tests, all three earlier
examples, pyproject.toml, governance/planning documents and prior reports are
unchanged. Existing Steps 03–09 remain uncommitted/untracked. No sibling file or
previous reference/planning bundle was modified. In particular,
`.validation/step10-plan-r1-20260913T021406Z/` remains intact. Final maintained
hashes/Git state and evidence hashes are recorded separately in `$RUN`.

## Public behavior

`prepare_forest_weights` requires explicit field, BinGeometry,
InstrumentResponse, source redshift greater than evaluation redshift, ordered
magnitudes, positive quadrature weights, normalized nonnegative source density,
pixel variance, positive forest length, and weighting method. There is no raw
survey reader, renormalization, inferred magnitude spacing or area factor.
`density_per_velocity` uses c=299792.458 km/s at the explicit source redshift.

`method="supplied"` accepts nonnegative weights with positive weighted support,
without model queries or legacy settings. `method="legacy"` requires an explicit
nonnegative integer iteration count and positive S/B samples. Initialization is
followed by exactly the requested number of simultaneous cumulative updates.
Zero-density samples get zero iterative weight; positive-density zero-variance
samples get one. Results expose weights, I1/I2/I3 prefixes, A in deg²,
P_pixel in deg² km/s and successive maximum absolute weight changes. Zero
prefixes and pixel noise are valid. Invalid/unrepresentable arithmetic fails
without artificial density, clipping, jitter, or a stopping-rule change.

`sample_auxiliary` reuses the original P3D callable/binding with an auto-only
selection retaining original field order. It never requests the negative
forest/galaxy cross for weighting. Independently bound P1D is queried once;
no model evaluation occurs in magnitude/iteration loops. S is the explicitly
chosen full intrinsic auto, transformed with W² and a_v/d_deg², rather than a
requirement for a smooth-only provider. Auxiliary coordinates and fiducial states
are retained, with field/provider context on domain/nonpositive-power failures.
The examples `(2.4, .00035)` and `(7, .001)` are documented, not defaults.
Auxiliary queries outside forecast cuts add no modes or alter those cuts.

Owned arrays are bytes-backed immutable float64. Preparation retains field
identity, source/evaluation redshifts, h, conversion factors, speed convention,
response, input arrays, method, iterations, and optional auxiliary samples/states.
It retains no mutable provider or parameter-dependent cache. Noise evaluation
checks matching field/evaluation geometry/h/response. Area alone can change
without changing local noise.

`forest_noise` consumes independent intrinsic P1D samples at observed
q=k*mu/a_v and returns aliasing, unsmoothed pixel and total arrays in
(Mpc/h_fid)^3. W² applies only to aliasing; both terms get d_deg²/a_v once.
`galaxy_noise` is 1/n_bar. `local_galaxy_density` explicitly evaluates the local
normalized dndzdm row at z_eval; it does not use source redshift, forest length,
survey area, integrated counts or bin volume.

`prepare_noise` generates diagonal noise only with explicit
`independent_sampling=True` and exact active-field coverage. Unused fields need
no evaluations. Alternatively, full packed noise replaces all generated terms;
ambiguous combinations/shapes are rejected. Noise PSD is checked independently
of signal by reusing `_validate_field_power`, with its existing normalized
64*eps64 tolerances and exact zero-row rule. Signed crosses and singular/zero
noise are valid. Covariance code itself was not changed. Noise is added once
through `combine_observed_power`, excluded from the default mean, and fixed
along with weights, responses, modes and factors across intrinsic derivatives.

## Acceptance results

- Focused suite: **73 passed in 0.82 s**, `focused-final.log`.
- Complete quick runner: **628 passed in 21.88 s**, Ruff lint passed and
  **82 files** formatted, `checks-final.log`. All previous 555 tests remain.
- **81 bounded legacy comparisons passed**, `reference-final-limits.log` and
  `reference.json`, covering synthetic arrays and all five real populations.
- Installed exact wheel: scalar weight/noise oracle, full-noise packing and PSD,
  generated-noise amplitude Fisher, all four examples, optional Astropy path,
  and actively blocked Astropy/SciPy subprocess passed. All **28 package Python
  modules** match source/wheel/installed bytes; METADATA and WHEEL match too.
  See `wheel-probe.log`, `wheel-probe.json`, `blocked-optional.log`.
- Git whitespace, preexisting-file preservation and final hash checks passed.

Tests use independent scalar prefix sums, one-population closed forms and
matrix covariance/Fisher oracles. They distinguish cumulative updates from
total-density and in-place alternatives. Covered limits include zero density,
zero variance, deliberate supplied zero weights, nonuniform quadrature,
normalization and length scaling, h³ scaling, zero/sinc-null/negative-lobe
response, source/evaluation separation, independent P1D substitution, built-in
and external auto routes, template domain failures, zero auxiliary response,
active-field subsets, permutations/slices, signed/singular/full noise and invalid
noise matrices. Parameter perturbations reuse the prepared factors and leave
fixed arrays unchanged. Doubling area doubles information with identical noise.

Direct scalar comparisons use rtol at most 5e-14 in ordinary tests, with explicit
absolute allowances at zeros. Successive weight changes subtract nearly equal
weights; their scalar-oracle check uses atol=2e-16 to allow one-ulp cancellation.
The synthetic sinc null uses atol=1e-28 in power units. Numerical amplitude
Fisher comparisons use rtol=3e-12 for central finite differences at step scales
1, 1/2 and 1/4; the corresponding linear intrinsic amplitude makes truncation
zero, leaving subtraction roundoff. No old tolerance was changed.

For the prescribed r=[1,2,4], s²=[1,4,9], S=3, B=2, L_v=10, Delta_v=.5 fixture,
maximum absolute changes for successive updates are:

| Update | Maximum absolute weight change |
| --- | --- |
| 1 | 0.6451492884382849 |
| 2 | 0.025322847776604984 |
| 3 | 0.0004103760471942941 |
| 4 | 0.000006592134646044023 |
| 5 | 0.00000010695465846310981 |
| 6 | 0.0000000017490269232922628 |

Initialization and 0/1/2/3/6-update results/prefixes are retained in reference
JSON. These diagnostics do not assert convergence or general monotonicity.

The new example gives amplitude Fisher values
**7324.260999453555, 7324.260999453395, 7324.260999451951**, and amplitude error
**0.011684714192096686**. The one-parameter information is positive, with no
null direction and no added prior. Its two forest update-change vectors are
`[0.2844736029140229, 0.005398595780757499, 0.00025247037862134904]` and
`[0.10553572481392659, 0.0003813913857891364, 1.5223734949661476e-6]`.
`example-clean.json` retains results and fixed arrays. These are synthetic
normalization checks, not physical survey predictions.

## Bounded reference interpretation

`$RUN/reference.py` extracts named legacy Weights/Covariance methods through AST
and supplies minimal objects with matched rho, constant endpoint dm, variance,
L_v/Delta_v and positive auto S=3/B=2. It compares initialization, each of three
updates, additional iteration diagnostics, prefixes, the actual legacy
coefficient method, forest total-noise conversion and discrete Poisson noise.
The maximum relative discrepancy is **1.6515801740230874e-15**; maximum absolute
discrepancy is **7.275957614183426e-12** for a larger noise value. Every comparison
passes rtol=5e-13, atol=2e-16; no physical cross-based weighting is compared to a
forest auto as if they were identical models.

The separate DESI-2 input/intermediate check uses the unchanged authoritative
15x2pt INI, named Survey/Tracer/Spectrograph input routines and local resources.
No NewForecast, CAMB, legacy Fisher, capture or full-mode execution occurs.
Reference source/resource hashes (34), Git state, normalization measures,
normalized density rows, pixel variances, source/evaluation redshifts,
magnitudes/quadrature, conversion factors, S/B and compared outputs are retained.
GPLv3/source provenance accompanies the adapted formulas; no raw assets were
copied into tests or package data.

The common bin is [2.47, 2.705], z_eval=2.5875. Both forests use
z_source=2.895806842635155, mapping the chosen observed wavelength to the geometric
mean of their rest-frame limits. Geometry is explicitly synthetic supplied
H=250 km/s/Mpc and D_M=5500 Mpc, h_fid=.7; it is not a Planck18 reconstruction.
Each forest uses its own SNR files and exposure/pixel settings, without pair-wise
averaging. Four original Survey magnitude nodes are sampled per population:
QSO/forest-QSO at 21.02311321–21.32452830; LBG/LAE/forest-LBG at
23.73584906–24.03726415. The exact arrays, not rounded values, are in JSON.
These bounded subgrid integrals do not represent complete population counts.

The legacy normalization measures remain unchanged: forest magnitude cuts precede
normalization over raw cell counts at z>2.15; discrete QSO also uses z>2.15;
other discrete populations use the full table. QSO raw normalization measure is
3789.6134994140098, LBG is 380, and LAE is 430.00000000000006, before the original
target factors. Actual raw dz/dm and table shape are recorded; the reference
methods divide cell counts by dz*dm before interpolation. None of these policies
was moved implicitly into the new array API.

Initial candidate magnitudes 24.33867925–24.64009434 yielded negative spline
outputs in both forest density rows (QSO down to -0.00161094; LBG down to
-8.56997373). Both failed checks are preserved in `reference.log` and
`reference-second.log`, and exact rejected rows are retained in final JSON.
The final documented samples lie in the density/SNR domains with positive
interpolated density and SNR; no clipping, floor activation, normalization change
or input-file edit was used. This exposes a limitation of the reference
interpolation rather than hiding it in production weighting.

Matched-array calculations use legacy c=299800 only in reference preparation.
The production density conversion ratio is quantified independently as
299800/299792.458, with c_new/c_legacy also recorded. Nonzero response comparisons
match; at q=0 the old sin(x)/x routine yields NaN while the new response gives
one. Final response/aliasing unit limits are independently covered by synthetic
tests. There is no claim that these intermediates form a physical DESI-2 forecast.

## Commands, environment and artifact

Commands run from the package root unless stated otherwise. Numerical scripts
and tests use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` on the
login node. `$ROOT` is the component root and `$RUN` its evidence directory.

```bash
.venv/bin/python -m pytest tests/test_weights.py tests/test_noise.py \
  tests/test_weight_models.py -q
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
.venv/bin/python "$RUN/reference.py"
.venv/bin/python examples/weighted_noise.py
.venv/bin/python -m build --wheel --outdir "$RUN/dist"
.venv/bin/python -m venv "$RUN/wheel-env"
"$RUN/wheel-env/bin/python" -I -c \
  "import importlib.util; assert importlib.util.find_spec('fishhighz') is None"
"$RUN/wheel-env/bin/python" -m pip install \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl[templates,cosmology]"
# From /tmp, all paths absolute:
"$RUN/wheel-env/bin/python" -I "$RUN/probe.py" \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" "$RUN"
git diff --check
```

One exact wheel was produced, installed into one fresh environment. Initial build
and dependency-install attempts failed at sandbox package-index access; approved
retries succeeded (`build-approved.log`, `wheel-install-approved.log`), with no
prior wheel artifact to replace. No second environment or wheel/sdist matrix was
used. `fresh-environment.log` verifies initial package absence. Byte comparison
and installed origins prevent a stale same-version package from passing.

Wheel SHA-256:
`e379ddee37694daabb1acc9220d32d6ba1f9f7ec6df1a593f08cdd61e47a2cd2`.
Python 3.13.15, FishHighz 0.1.0.dev0, NumPy 2.5.3, Astropy 8.0.1, SciPy 1.18.1,
pytest 9.1.1, Ruff 0.16.7, build 1.6.1, pip 26.2.1. `versions.json` and the
installed probe record environments independently. NumPy remains the only
unconditional requirement; no new extra was added. Installed origins are under
the new `wheel-env/lib/python3.13/site-packages/fishhighz/`.

## Failures, repairs and limits

The initial preservation scan unnecessarily descended into ignored directories;
it was interrupted before source edits and replaced with a bounded maintained
file inventory. A wildcard file-discovery command was blocked by the workspace
traversal hook; a bounded directory inventory succeeded without bypassing the
hook. The first scalar diagnostic check needed the 2e-16 absolute subtraction
allowance described above. Ruff corrections were restricted to new files. A
later test import cleanup was followed by the final complete runner. No accepted
scientific module or numerical tolerance changed.

Intermittent NERSC MUNGE socket messages appeared after successful commands; no
Slurm action occurred. One redirected example stdout file contained trailing
MUNGE text and could not be parsed as JSON. It is preserved as `example.json`;
`example-clean.json` was written directly by the example's `run()` result and
parses successfully. Packaging network failures and reference interpolation
failures were repaired/disclosed as above. All required bounded checks now pass;
no acceptance evidence is knowingly missing.

No real full forecast suite was requested or run under the user's quick-only
policy; it is not an outstanding acceptance gate. No raw production readers,
overlap-derived noise, new optimization, P1D integration from P3D, model physics
change, covariance derivatives, serialization, JIT or performance work is added.
No commit, push, delegation, governance revision or Step 11 work occurred. Stop
for user and independent review.
