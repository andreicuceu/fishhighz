# Step 11 handoff — revision 2

Implemented the dispatched revision-2 repairs for R1 and R2. All required quick,
bounded-reference and fresh installed-wheel checks pass. Ready for user and
independent review; acceptance and progression are not claimed.

The revision-1 handoff is preserved byte-for-byte as
[step-11-r1.md](step-11-r1.md). Its three blocked density comparisons are historical;
all three are closed by newly evaluated revision-2 readers. Previous implementation
and independent-review evidence were not overwritten.

## Changes and preservation

Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228`, with existing
uncommitted work retained. New evidence is in
`.validation/step11-r2-implementation/` (called `V` below). `V/before.json` records
102 pre-edit maintained-file hashes, including untracked work;
`status-before.txt` records Git status. `preservation.json` checks those hashes,
the preserved r1 handoff, and the absence of missing files. `after.json` identifies
the final maintained state. `production.diff` compares the two changed production
modules to the exact r1 wheel.

Changed production files are only `fishhighz/adapters/legacy_inputs.py` and
`fishhighz/survey.py`. `forecast.py`, all 28 pre-Step-11 production modules,
pre-Step-11 tests/tolerances, all four older examples, package metadata/extras,
and governance documents are unchanged. Other edits are README documentation,
`examples/survey_forecast.py`, and the obsolete redshift-spacing rejection in
`tests/test_legacy_inputs.py`, now a nonuniform-magnitude rejection. Added
`tests/test_step11_repairs.py`, the preserved r1 handoff and this replacement
handoff. No revision-1 coverage was removed except the obsolete expectation.

## R1: explicit cell measures and unequal forecast bins

The density reader now accepts:

```python
DensityReader(path, *, semantics, target_density, z_norm_min,
              magnitude_bounds=None, redshift_widths=None,
              width_policy=None, label="density")
```

- `redshift_widths`: explicit real, finite, positive `(n_z,)` physical cell
  measures aligned with the reconstructed **sorted redshift axis**. This selects
  the recorded `explicit` policy. Widths are owned immutable float64 arrays.
- No widths and no policy: retain the existing constant-spacing conversion only
  for uniform z. Explicit `width_policy="uniform"` selects the same requirement.
- Explicit `width_policy="legacy_first_spacing"`: use z[1]-z[0] for every row,
  including irregular z. This reproduces the legacy divisor; it does not infer
  physical cell measures. It is never selected automatically.

Explicit widths and a compatibility policy conflict. Missing widths/policy for
irregular z produce an actionable error. Magnitude spacing must remain uniform
within the unchanged 64-epsilon coordinate tolerance. There is no edge inference,
rebinning, extrapolation, clipping or change to spline order.

Raw magnitude masking and the strict z>threshold raw count sum still precede
normalization. Each row is then divided by `redshift_widths[i]*dm`. Thus selected
`sum(density[i,j]*Delta_z[i]*dm)` recovers the target; the later magnitude
quadrature remains separate. Provenance records axes, effective widths, ordering
and width-policy label, alongside the existing normalization metadata and units.
The scalar `dz` provenance remains available for compatibility modes and is None
for physical explicit widths. SNR smoothing/exposure behavior is unchanged.

Forecast-bin geometry was already correct and required no runner edit. New
oracles test unequal widths [2,2.7] and [2.7,2.9], optional gaps, distinct areas,
non-midpoint evaluation redshifts, reversed caller order and four batch sizes.
Volumes are integrated analytically for constant H and D_M proportional to 1+z;
the k shell is integrated analytically without using production volume/q_mode in
the expected information. Both scalar bin Fishers remain singular; their sum is
full rank and one b prior is added once. Existing touching/overlap rejection
coverage remains. A redshift-dependent model/P1D spy observes each explicit
z_eval; the response spy verifies the corresponding a_v, with no repeated fixed
preparation during derivative execution.

The generated-reader example now uses [2.2,2.5] and [2.5,2.7], with explicit
evaluation redshifts 2.3 and 2.62. It reports its bounds and evaluation redshifts.
Raw source cells and forecast-bin widths remain separate.

## R2: dtype validation is no longer bypassed

`survey.freeze` now copies non-object ndarray data into immutable byte-backed
arrays **without casting their dtype**. Scientific validation in the unchanged
Step 10 APIs sees the original complex/Boolean/string values and rejects them.
Accepted prepared scientific arrays remain immutable float64. Python complex
values in lists are retained until the same validator rejects them; no imaginary
part is discarded and tests turn warnings into errors.

Object arrays are rejected at snapshot setup because object references cannot be
made immutable through a byte-backed snapshot. Numeric-object scientific arrays
therefore fail consistently with the direct API. Typed metadata arrays retain
values and dtypes, including integer indices and Boolean flags. Plain nested
metadata remains supported; the object-array restriction is documented.

New tests compare direct `prepare_forest_weights` with `ForestInput` through
`prepare_bin` for all five fields (magnitudes, quadrature, rho, variance, supplied
weights). They cover nonzero-imaginary complex, Boolean, numeric-string and
numeric-object arrays, corresponding lists where their dtype information exists,
and NaN/infinity. A numeric list made from an object array no longer has object
dtype; ordinary integer/float list controls test that valid case. Accepted
integer/float, exact-zero, caller-mutation and write-enable rejection controls
verify the immutable float64 preparation boundary. Metadata controls verify
faithful Boolean, int32, float32, string and complex dtypes.

## New checks and numerical evidence

All numerical commands used OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and
MKL_NUM_THREADS=1. Lightweight work ran on the login node, with no Slurm actions.

| Command from package root unless noted | Result / evidence |
| --- | --- |
| `.venv/bin/python -m pytest tests/test_step11_repairs.py tests/test_forecast.py tests/test_legacy_inputs.py tests/test_forecast_imports.py -q` | **195 passed** in 2.96 s; `V/focused.log` |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | **833 passed** in 28.21 s, Ruff lint and format passed; `V/quick.log` |
| `git diff --check` | Passed; `V/diff-check.log` |
| `.venv/bin/python V/reference.py` | **34 matching hashes, 33 passing new comparisons, zero blocked**; `V/reference.json` / `.log` |
| `.venv/bin/python examples/survey_forecast.py` | Unequal-bin raw-reader example passed; `V/example.json` |
| `.venv/bin/python -I V/base_oracle.py` | Analytic unequal-bin checks, normalized forest, 35 direct/orchestrated dtype rejection controls and missing-extra errors passed; `V/repair-source.log` |
| `.venv/bin/python -I V/density_probe.py` | Physical-width and legacy-first-spacing scalar oracles passed; `V/density-source.log` |
| `.venv/bin/python V/range_regression.py source` | All four retained Step 10 range cases passed; `V/range-source.json` |
| Fresh wheel probe from `/tmp` with `-I` | Passed; `V/wheel-probe.json` / `.log` |

There are 110 added repair tests; 55 parameterized cases compare direct and
orchestrated rejection of invalid scientific inputs. The independent physical
cell-count and unequal-bin algebra use rtol=5e-13, atol=0; off-grid quadratic
values use rtol=5e-12, atol=0. No existing tolerances changed. Initial development
checks also passed; final logs above reflect the complete non-vacuous matrix of
list/ndarray cases.

## Bounded reference closure

The new readers were evaluated at the exact four original magnitude nodes and
source/evaluation redshifts from the saved Step 10 reference JSON for all five
populations. The supplied background remains H=250 km/s/Mpc, D_M=5500 Mpc,
bin [2.47,2.705], z_eval=2.5875. All 34 original source/resource/INI hashes match.

The QSO reader keeps uniform compatibility. LBG forest, LBG and LAE explicitly
select `legacy_first_spacing`, recorded in each reader's provenance. All former
blocked density normalization/dndzdm/local-density or forest-weight/noise checks
now pass. Both QSO and LBG SNR/variance paths are newly evaluated. Maximum relative
discrepancy over the 33 comparisons is **3.513e-15**, with zero absolute allowance.
Both known QSO/LBG negative spline rows fail at their unchanged original query
coordinates. There is no remaining blocked comparison.

`V/reference.json` retains exact queries, policies/effective widths, raw file and
source hashes, normalization/exposure settings, intermediate arrays, discrepancies
and negative-query errors. Physical-width correctness uses the independent
synthetic oracle; no physical widths were invented for real irregular tables.
Downstream legacy equality retains the recorded c=299800 pixel/length/rho inputs
explicitly, with separate new-c density-ratio checks. The public reader composition
continues using c=299792.458. No floors, bright clamps or alternate query nodes
were used, and no NewForecast/CAMB or real reference forecast was run.

## Fresh installed wheel

One successful wheel was built in fresh `V/dist/` and installed into fresh
`V/wheel-env` with templates, cosmology and survey extras:

- `fishhighz-0.1.0.dev0-py3-none-any.whl`
- SHA-256: `940eae9b77b9b4f9a1b07ce728a5d9addf25ee46a2e03228f4d60d004a36dce6`
- Python 3.13.15; NumPy 2.5.3; SciPy 1.18.1; Astropy 8.0.1.

The exact artifact, all **32 source/wheel/installed production modules**, and
wheel/installed METADATA/WHEEL bytes match. Origins are inside the fresh
environment; NumPy remains the only unconditional dependency. No system site
packages or previous FishHighz installation were used.

From `/tmp` with `-I`, the installed probe runs all four retained examples and the
updated raw-reader example, scalar/analytic unequal-bin Fisher checks, both
cell-width policies, 35 R2 rejection controls, valid integer/float/zero/metadata
controls and the retained Step 10 range regression. Generated density fixtures
live outside the checkout. The active SciPy/Astropy blocker completes normalized
external multi-bin/forest orchestration and confirms actionable survey-extra
errors for both readers. Probe scripts and `blocked-multibin.log`,
`density-wheel.log`, `range-wheel.json` preserve the evidence.

Initial build/install attempts failed because sandbox DNS could not reach PyPI;
approved retries succeeded. Both failed and successful logs are retained as
`build.log`, `build-approved.log`, `wheel-install.log` and
`wheel-install-approved.log`. NERSC MUNGE messages appear in test logs but the
quick runner exited zero. No scheduler operation was attempted.

The updated synthetic example reports F_AA=116.95133768634004 and prior-assisted
marginalized errors [0.0924692646, 0.4999395062, 0.4999544824] for [A,b_low,b_high].
Per bin, preparation still makes three P3D calls (two auxiliary autos plus one
fiducial), and each forest makes two P1D calls. Twelve nodes in batches of five
use 15 numerical P3D calls per bin/run and zero Jacobian callbacks. These are
synthetic checks, not DESI-2 forecasts or performance claims.

## Review boundary

R1 and R2 are implemented and their required evidence passes. Fixed-state reuse,
shared registry, pair order, one combined prior and O(n_node*n_selected²)
in-memory factor storage remain unchanged. INI/CLI/serialization and Step 12
remain deferred. No full forecast was requested or run, and it is not an
outstanding gate. No commit, push, governance edit or agent dispatch occurred.
Stopping for user and independent review.
