# Step 12 handoff — revision 1

Implemented the dispatched DESI-2 validation, opt-in input compatibility and
external intrinsic-P3D bridge. Required quick tests, bounded real validation,
reference comparisons and fresh installed-wheel checks pass. Ready for independent
and user review; this is not user acceptance or authorization for Step 13.

## Changes and preservation

Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228`. Evidence is in
`.validation/step12-r1-implementation/` (`V` below). `before.json` and
`status-before.txt` captured maintained tracked/untracked files before edits.
`after.json`, `preservation.json` and `status-after.txt` identify the final state.
All earlier production modules, tests, five examples, reports, package metadata,
governance documents, raw assets and baseline bundles are retained. README is
the only previously existing maintained file changed. No sibling source edits,
commits, pushes, agent dispatch or Slurm actions occurred.

Added production modules:

- `adapters/legacy_compat.py`: `LegacyDensity(reader, negative_policy)`,
  `LegacySNR(reader)`, `sample_legacy_forest(...)`, plain provenance conversion.
- `adapters/lyaforecast.py`: `IntrinsicP3D(provider, *, routes, n_fields,
  h_source, h_fid, k_domain, z_domain)`.
- `validation/{__init__,_recipes,cases,evidence,synthetic,real}.py`: seven explicit
  recipes, exact INI verification, synthetic covariance oracle, exclusive evidence
  orchestration/checker and `RealRecipe` preparation/refinement/amplitude example.

Added `scripts/desi2_validation.py`, `scripts/legacy_policy_checks.py`,
`examples/desi2_synthetic.py`, `tests/test_step12.py`, README documentation and
this report. No accepted numerical engine or scientific tolerance changed.
GPLv3/reference provenance is recorded for adapted input logic and recipes;
no raw reference assets are distributed. NumPy remains the only unconditional
runtime dependency, with no eager reference/CAMB/Astropy/SciPy imports.

## Input compatibility and ownership

Strict Step 11 readers are unchanged. Compatibility consumes their validated
raw snapshots and prepares owned interpolators. It does not reread inputs in
numerical batches or mutate counters during derivatives. Tests actively mutate
caller interpolators and confirm the compatibility snapshots remain unchanged.

Density magnitude queries strictly outside the raw axis become 1e-20. Closed
endpoints remain inside. Redshift uses the actual reference spline extension;
original/effective coordinates and extension masks are retained. Exact zero and
positive values below the floor remain unchanged. Explicit `floor_negative`
replaces only negative interpolants; `reject` still fails. Negative flooring is
an extension beyond legacy behavior, never reported as reference equality.
All real populations explicitly select `legacy_first_spacing`; no physical
redshift-cell widths are inferred for irregular tables.

An initial endpoint comparison exposed a reduction-order difference at nearly
zero spline values: the strict forest normalizer sums selected rows, while the
reference sums the complete zero-masked flat table. Compatibility now preserves
the reference arithmetic order in its own preparation. The strict reader and
its tolerances are unchanged. The failed initial log is `policy-01.log`; newly
evaluated original-table comparisons pass in `policy-final.json`/`.log`.

SNR independently clamps bright magnitudes, returns exactly 1e20 variance for
faint/z/wavelength sentinel branches, and applies the 1e-10 SNR floor after
pixel/exposure scaling. Exposure time must match the file. Invalid/nonfinite
inputs and arithmetic failures remain errors. Unevaluated raw SNR placeholders
are zero with an explicit evaluated mask. Provenance contains source hashes,
units/normalization, original/effective coordinates, raw interpolants, policies,
constants, masks and counts. `sample_legacy_forest` carries plain diagnostics
into `ForestInput`; each forest owns its source/SNR/response/P1D separately.

## Actual scientific coverage and settings

Only the original **15x2pt case, bin [2.47,2.705]**, ran as a real FishHighz
forecast. `real-02/` is the corrected final numerical study; `real-01/` retains
the earlier successful development run. Other real cases/bins were not run.
All seven cases and all 39 original case/bin tasks are tested synthetically.
The full real path is available with explicit `--suite full --allow-real-full`
and optional case/bin selection; it requires a separate user request.

The seven original INIs are verified read-only against explicit recipes,
including every option. Selected counts are 3,3,6,2,4,8,15 in the assignment's
order. Pair order follows the reference's original-field upper triangle, rather
than the textual INI correlation-list order. Both forests retain their original
indices and all covariance-required spectra remain present.

For the real bin:

| Setting | Value/convention |
| --- | --- |
| Evaluation redshift | 2.5855752676523194, geometric mean of 1+edge redshifts |
| Historical dictionary label | 2.5875000000000004, arithmetic centre |
| Area and integrated volume | 5000 deg²; 4.598178233605516e9 (Mpc/h_fid)³ |
| Observed cuts | [0.01,0.5] h_fid/Mpc, unchanged in every refinement |
| Free parameters | ap_2, at_2; fiducials (1,1); no nuisances or priors |
| Fixed biases, in original order | forest −0.16753635432407413; QSO 3.937751453604524; LBG 3.3; LAE 1.8729593533010616; second forest same bias |
| Forest beta / galaxy f | 1.45 / 0.9762242438440839 |
| Power growth G | 0.9023428889225946, [(1+z_template)/(1+z_eval)]² |
| Damping sigma8 ratio | 0.9227188220114164 relative to reference z=2.3 |
| Forest widths (parallel, transverse) | (5.9446077385713005, 3.0080633597572173) Mpc/h_fid |
| Galaxy widths | (4.203472443437794, 2.1270219999231177), reconstruction factor 2 |
| Magnitudes | All 107 original nodes, 16.1–26.75; endpoint-inclusive rectangular weights |
| Forest source/length | lambda_obs=1215.67(1+z_eval); source z from rest limits 1040/1205; L_v=c log(1205/1040) |
| Response | 0.8 Angstrom pixels, explicitly labeled legacy c/2500 velocity sigma; c=299792.458 |
| Forest weighting | Three cumulative legacy updates, each field's auxiliary auto at (2.4,0.00035); independent default P1D |
| Galaxy noise | Each population's local z_eval density; explicit independent sampling |

The external Vega-format `Planck18_z_2.406.fits` hash is
`b4a73103e1105b7f9cdb59bbe8133ff0f752b05bc27580dd526f80b9c2fdc26a`.
Its z_ref, h_template, h_fid, complete header, raw source/resource hashes,
versions, settings, parameter IDs and overrides are recorded in the manifest.
PKSB remains smooth and PK−PKSB signed wiggle. Identity smooth scaling is fixed;
wiggle derivatives include inverse coordinates, volume prefactor Q, remapped
RSD and damping with fixed widths. Geometry, noise, observed response, cuts and
covariance factors remain fixed during derivative refinement.

Actual fallback counts among each population's 107 magnitude samples:

| Population | Magnitude density floor | Negative-density extension | Bright SNR clamp | SNR sentinel |
| --- | ---: | ---: | ---: | ---: |
| QSO-source forest | 19 | 9 | 32 | 20 |
| QSO galaxies | 19 | 7 | — | — |
| LBG galaxies | 57 | 12 | — | — |
| LAE galaxies | 57 | 12 | — | — |
| LBG-source forest | 57 | 13 | 32 | 20 |

No real-bin redshift extension or post-scaling SNR floor occurred. These counts
make the compatibility choice prominent; this is not unchanged strict-domain
physics. Three updates and fixed magnitude nodes are not convergence claims.

## Numerical checks and results

Every numerical command used OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and
MKL_NUM_THREADS=1 on the login node. No full real suite or new legacy forecast
capture was requested or executed.

| New check | Result and evidence |
| --- | --- |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | **1017 passed in 26.55 s**, Ruff lint/format pass; `quick-final.log` |
| New focused tests | 184 new cases in `tests/test_step12.py`; original 833 retained |
| Seven tiny synthetic forecasts | Independent five-field covariance slices and correlated Fisher matches; exact field/pair inventories |
| Full-mode synthetic workers | All 39 tasks once; failures remain visible; exclusive output directories and corruption tests pass |
| Wiggle BAO oracles | Analytic isotropic oscillatory derivative/Fisher and independent fourth-order anisotropic differences pass; per-bin global null columns and independence verified |
| Original-table actual reference methods | **104 comparisons**, max relative 3.1516201241062103e-16; all 34 previous reference identities match; `policy-final.json` |
| Generated-table actual initialized Tracer/Spectrograph methods | **112 comparisons**, including endpoints, extension, negatives, clamp/sentinel/scaling/floor; `generated-reference.json` and generated fixtures |
| Offline historical bundle checks | Full and quick pass; `baseline-full.log`, `baseline-quick.log`; schema-v1 limitation retained |
| Corrected real bin | All convergence thresholds pass; `real-02/manifest.json`, `case-000.npz`, `real-02.log`; 160.25 s including setup/amplitude |
| Independent selected covariance / Fisher | Relative Frobenius discrepancies 1.2264e-16 / 5.0167e-15 |
| Actual external amplitude bridge | Direct intrinsic P3D equality, deterministic repeated calls, independent full-field trace Fisher relative discrepancy 2.7756e-15 |
| Final-wheel real preparation / independent W and noise | Fixed arrays match `real-02` bitwise; W scalar oracle passes; noise relative 3.54175e-19 using actual reference P1D; `final-preparation.json`, `final-weights.npz`, `noise-oracle.json` |
| Offline new bundle checker | Pass; `check-real.log` |
| Fresh installed wheel | All 40 production modules and METADATA/WHEEL match; five retained examples, synthetic seven-case example, compatibility/bridge oracles and active optional blockers pass |
| `git diff --check` | Pass; `diff-check.log` |

Final-two-level refinement discrepancies, holding other settings at their final
levels (k has order 4 per subinterval):

| Axis and final comparison | Fisher Frobenius relative | Largest error relative | Volume relative |
| --- | ---: | ---: | ---: |
| k intervals 64→128 (also ran 32) | 9.92391e-6 | 2.76228e-6 | 0 |
| mu order 16→32 (also ran 8) | 5.45898e-8 | 1.38151e-8 | 0 |
| geometry z order 16→32 (also ran 8) | 5.10971e-12 | 2.55496e-12 | 5.10947e-12 |
| ap/at step 5e-4→2.5e-4 (also ran 1e-3) | 1.24166e-5 | 8.11295e-6 | 0 |

Limits remain 1e-3 Fisher, 5e-3 errors and 1e-6 volume. No further refinements
or weakened thresholds were necessary. Final two-parameter condition number is
3.10209, with errors **sigma_ap=0.0174151762792** and
**sigma_at=0.0122614424805**. NPZ preserves actual grids/modes, intrinsic and
observed powers, response, noise, factors, Jacobians, selected/required pairs,
independent covariance/Fisher and every refinement matrix/error vector.

## Reference comparison and interpretation

Historical saved errors for the same bin are ap=0.0179967890622 and
at=0.0123719968890. New errors differ by **−3.23176%** and **−0.893586%**.
`historical-comparison.json` retains both values and original redshift label.
These are unmatched forecasts: there is no equality tolerance and these final
numbers cannot identify each convention's causal contribution.

Known non-equivalences are explicitly retained:

- PKSB decomposition versus degree-8 log-polynomial peak extraction.
- Full template/remapping/Q/RSD/damping derivatives versus backward log-k
  differences, zero first node and the legacy radial-derivative-only recipe.
- Fixed per-tracer widths with mean squared cross widths versus legacy pair
  damping; forest–galaxy reference pairs skip galaxy reconstruction.
- Integrated volume versus centre-bin approximation; Gauss–Legendre k/mu
  integration versus reference sums.
- Fixed per-field weights/source/SNR/response versus pair-specific reference
  choices, including mixed-forest source handling.
- Explicit input floors/clamps and the separately confirmed negative extension.
- c=299792.458 versus 299800; explicitly labeled response conversion conventions.
- External template/header versus CAMB Planck18.ini cosmology/decomposition,
  template-z EdS power growth versus independent sigma8 damping growth, and
  recorded differences between numerical environments.

The actual external object is `lyaforecast.power_spectrum.PowerSpectrum`, prepared
once with `CosmoCamb` and explicit configured bias functions. No NewForecast or
smoothed/observed power entry point is called. Source h units and domains are
explicit; unit conversion/order/sign/domain/error contracts also have mock
oracles. The real amplitude wrapper uses the same per-field response exactly
once and fixed noise, leaving default P1D independent. It does not implement
wiggle-only BAO dilation.

The corrected bounded command, from the package root, was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH="$PWD" ../lyaforecast/.validation/dev-env/bin/python -u \
  scripts/desi2_validation.py run --reference ../lyaforecast \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --output .validation/step12-r1-implementation/real-02 \
  --negative-policy floor_negative
```

The same thread/environment prefix ran `scripts/legacy_policy_checks.py` with
`--reference ../lyaforecast`,
`--saved .validation/step10-r1-implementation/reference.json` and
`--output .validation/step12-r1-implementation/policy-final.json`.
`V/generated_reference.py` and `V/noise_oracle.py` are retained bounded oracle
scripts. `V/final_preparation.py` uses the exact wheel path as PYTHONPATH instead
of the checkout and preserves origin/bitwise-comparison evidence.

## Packaging and reproducibility

Final exact wheel:

- `V/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`
- SHA-256 `b240003ad6497d77fe850107152765750065b027c18fb79f92ccfe157f00847a`
- Fresh `V/wheel-env`, no system site-packages: Python 3.13.15, NumPy 2.5.3,
  SciPy 1.18.1, Astropy 8.0.1; templates/cosmology/survey extras installed.

From `/tmp` with `-I`, `probe.py` verifies exact source/wheel/installed module
inventories and bytes, metadata and artifact hash, then runs the retained
examples and prior normalized-array oracles. `step12_probe.py` adds compatibility
scalar algebra, external routing/units and all seven synthetic cases. Separate
`desi2_synthetic.py` execution is saved in `synthetic-wheel.json`. Subprocesses
actively block optional dependencies while running normalized external forecasts;
missing-extra reader/cosmology errors remain actionable.

The real reference demonstrations instead use the existing read-only
`../lyaforecast/.validation/dev-env/bin/python`: NumPy 2.3.5, SciPy 1.15.3,
Astropy 8.0.1, CAMB 2.0.1, lyaforecast 0.1.0. Actual import origins and source
hashes are recorded. No historical environment was installed into or modified.
The final preparation replay uses the exact wheel on that environment's import
path and retains its actual origin, independently checks W, and compares fixed
arrays against `real-02` before saving weight snapshots in `final-weights.npz`.

Build/install failures are retained: the no-isolation attempt lacked setuptools;
the sandboxed isolated build and dependency installation could not reach the
package index. Approved retries succeeded (`build-approved.log`,
`wheel-install-approved.log`). Exactly one successful wheel artifact was built.
MUNGE diagnostics occurred in ordinary logs without failing checks or implying
any scheduler action.

The bundle checker is validation-only JSON/NPZ tooling, not stable production
serialization. The historical full schema-v1 bundle proves recorded content
integrity but does not gain source/input/configuration inventory preservation
from an offline pass. Full real execution is unrequested, not a failed or
outstanding completion gate. General INI translation, production CLI and result
serialization remain deferred. Stop here for user and independent review.
