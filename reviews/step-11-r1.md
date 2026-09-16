# Step 11 handoff — revision 1

Implemented the dispatched Python survey orchestration and explicit raw-input
adapter assignment, revision 1 (2026-09-13). Ready for user and independent
review, **with a blocked portion of the bounded legacy density comparison**
described below. This is not a claim of acceptance or permission to progress.

## Scope and preservation

Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228` with the pre-existing
dirty/untracked work preserved. Only two pre-existing files changed during this
assignment: `README.md` (API/convention documentation) and `pyproject.toml`
(the lazy `survey = ["scipy"]` extra). All 28 pre-existing production modules,
638 accepted tests and their tolerances, four examples, prior handoffs and
existing evidence remain unchanged. Governance documents were not edited.

Added:

- `fishhighz/survey.py`: `ForestInput`, `BinSpec`, `PreparedBin`, ownership helpers.
- `fishhighz/forecast.py`: preparation, bin execution and independent combination.
- `fishhighz/adapters/__init__.py`, `adapters/legacy_inputs.py`: standalone strict
  readers and reader-to-normalized-array composition.
- `examples/survey_forecast.py`: generated fixtures and a two-bin demonstration.
- `tests/test_forecast.py`, `test_legacy_inputs.py`, `test_forecast_imports.py`:
  85 new checks.
- This report.

New evidence is under `.validation/step11-r1-implementation/` (abbreviated `V`
below). `before.json` records 89 pre-edit source/test/example/review/root-file
hashes and `status-before.txt` records Git status. `preservation.json` verifies
those inputs. A separate `planning-manifest-check.json` checks all 92 entries in
the saved planning manifest: only its prior authorized AGENTS/assignment changes
and this assignment's README/pyproject changes differ; no entries are missing.
`after.json` identifies the final maintained files, and `step11.diff` captures
tracked differences (which also contain older uncommitted steps).

## APIs and fixed-state conventions

Full signatures and runnable usage are in README's “Python survey orchestration
and raw inputs” section and public docstrings. Principal entry points:

```python
prepare_bin(spec: BinSpec) -> PreparedBin
run_bin(prepared, *, batch_size=None, steps=None, step_scale=1.0,
        numerical=False) -> BinRun
run_forecast(bins, *, prior_fisher=None, batch_size=None, steps=None,
             step_scale=1.0, numerical=False) -> ForecastRun
```

`BinSpec` contains the existing geometry/grid, PreparedP3D and explicit response
mapping. PreparedP3D owns the only pair selection and original field order.
Generated inputs are active forest `ForestInput` records plus explicit galaxy
n_bar and `independent_sampling=True`. Alternatively, full packed noise replaces
all generated inputs and independence settings; it is PSD-validated before
signal addition. Unused forests require no source/P1D/weight evaluations.

`ForestInput` accepts the normalized keyword inputs of
`prepare_forest_weights`, an independent P1D callable/binding/state, optional
auxiliary `(k_t_deg, k_p_velocity)`, and plain provenance. Legacy callable
weights query the matching auto and P1D once; direct S/B is also supported.
Supplied weights bypass auxiliaries. Final-grid P1D is evaluated in one vector
call per active forest. Source z and evaluation z remain distinct.

Preparation calls the accepted response, noise, mode-count, P3D, covariance and
factor APIs. `M=V*q_mode`, `R=W_i*W_j`, and `T=R*P_fid+N`. Immutable owned arrays
are `k/mu/modes(node)`, `response(node,field)`,
`products/noise/power/total(node,required_pair)` and
`factors(node,selected,selected)`. Structural provider/pair/binding arrays are
also copied into immutable byte-backed arrays. Weight results retain their
existing immutable arrays/context. Plain reader provenance crosses the boundary;
live interpolators and open files do not. Models remain callable references and
must be deterministic and unmutated; arbitrary callable internals are not hashed.

Execution delegates scheduling and derivatives to `evaluate_derivatives`, applies
R once, gathers `selected_to_required` and calls `fisher_from_factors`. No
covariance/noise derivatives, extra mode factor, inversion, cache refresh or
regularization is introduced. Consecutive node batches bound transient Jacobians;
full factors remain in memory at O(n_node*n_selected²) storage. Provider errors
retain their cause and report bin/global node slice plus provider/stencil details.

One exact registry is required across bins, with matching frozen fiducials,
h_fid and reused field identities. IDs are unique and redshift interiors cannot
overlap. Caller order is retained. Each `BinRun.result` has zero prior;
`combine_results` adds the supplied prior once. Singular individual Fishers are
valid. Prior changes reuse preparation. `BinRun.columns/calls/node_slices` and
`PreparedBin.diagnostics/weights` expose execution and preparation evidence.

## Raw-reader conventions

`DensityReader(path, *, semantics, target_density, z_norm_min,
magnitude_bounds=None, label="density")` requires explicit cell-count semantics.
Inclusive raw magnitude masks precede the raw-cell normalization sum, using
strict z>threshold or the whole grid. Explicit target=None skips renormalization.
Only then is dz*dm divided out. Reconstructed coordinate ordering tolerates row
permutations. Uniform spacing permits an absolute coordinate roundoff allowance
of `64*eps64*max(1,max(abs(axis)))`. Quadratic RectBivariateSpline uses kx=ky=2,
s=0 and the full bounding box. Closed-domain queries reject negative overshoot.

`SNRReader(paths, *, smoothing, label="SNR")` parses BAND/MAG/EXPTIME/NEXP and
Wave source-redshift headers, sorts by header magnitude and validates identical
metadata/grids. Explicit `legacy` smoothing is sigma=10 sample indices on the
wavelength axis, reflect boundaries, truncate=4; `none` disables it. Interpolation
is linear. Variance is `1/(SNR²*Delta_lambda*N_exp/N_exp_file)`; zero SNR,
incompatible explicit exposure time and unrepresentable arithmetic fail.

Both readers preserve paths/hashes, owned input tables, axes, settings and units.
`sample_forest_readers` samples source density and SNR at z_source, while using
`1215.67*(1+z_eval)` for observed wavelength. It checks consistent velocity and
Angstrom pixel widths using c=299792.458 km/s. The caller supplies magnitude
quadrature and forest length separately. The local galaxy helper is explicitly
selected and contains no area multiplier. GPLv3/lyaforecast source provenance
for adapted interpolation conventions is retained in the module.

## New acceptance evidence

All numerical commands used OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and
MKL_NUM_THREADS=1 on the login node. No Slurm action or real full forecast ran.

| Check | New result / log |
| --- | --- |
| `.venv/bin/python -m pytest tests/test_forecast.py tests/test_legacy_inputs.py tests/test_forecast_imports.py -q` | 85 passed; `V/focused.log` |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | 723 passed in 41.20 s, Ruff and format passed (95 files); `V/quick.log` |
| `git diff --check` | Passed; `V/diff-check.log` |
| `.venv/bin/python V/reference.py` | 34 matching source/resource hashes; 17 passing newly evaluated comparisons, 3 density comparisons blocked; `V/reference.json`, `reference.log` |
| `.venv/bin/python examples/survey_forecast.py` | Passed; `V/example.json` and installed example output |
| `.venv/bin/python V/range_regression.py source` | Four R1 cases passed; `V/range-source.json` |
| Exact installed wheel probe, outside checkout with `-I` | Passed; `V/wheel-probe.json`, `wheel-probe.log` |

The scalar and independent two-bin oracles use rtol=5e-13, atol=0. Raw off-grid
quadratic/linear checks, signed selected-subset covariance and batching/model
comparisons use at most rtol=5e-12, atol=0. The one explicit absolute null-product
allowance is 1e-11 for cancellation of a large matrix with a known null vector;
rank is independently checked. The exponential-model refinement check verifies
approximately 1/4 error ratios at successive step halvings, testing the actual
second-order truncation behavior. No accepted tolerance was changed.

Tests cover raw masks/strict thresholds/no normalization, shuffled grids,
independent reflected Gaussian sums, exposure/pixel scaling, domain/range/value
failures, source/evaluation conversion and reader/direct weight/noise equality.
Runner checks cover shared versus distinct nuisances, priors once, signed field
matrices/selection closure, built-in versus explicit external models, analytic
and numerical/tied routes, all/1/non-divisor/oversized batches, explicit field,
pair, parameter and bin permutations, auxiliary nodes outside forecast cuts,
perturbation-only template failures, indefinite N under large signal, singular
C, invalid registry/bin/response/step structure, no repeated fixed work or file
reads, and immutable arrays. The active optional-import blocker runs normalized
multi-bin and forest orchestration and verifies actionable errors for both readers.

Initial local development tests exposed reader label aliasing, overly exact
floating-sum assertions and invalid synthetic fixture choices; these were fixed
before the final passing checks. Build/install first failed on sandbox DNS;
approved retries succeeded. Logs also contain NERSC MUNGE environment messages;
the test runner itself exited zero. No scheduler operation was requested.

## Bounded legacy comparison: explicit limitation

The exact saved Step 10 nodes and source redshifts were reused for bin
[2.47,2.705], z_eval=2.5875, with constant supplied H=250 and D_M=5500.
All 34 original source/resource/INI hashes match. Newly evaluated QSO density
normalization/dndzdm, QSO and LBG smoothed SNR/variance, QSO local galaxy density,
and QSO downstream weights/integrals/noise agree with unchanged saved named
legacy-method outputs. Maximum relative discrepancy is **3.513e-15** (atol=0).
Full nodes, policies, inputs, settings and discrepancies are in `V/reference.json`.
Downstream equality explicitly uses the recorded legacy pixel/length/rho units;
a separate new-c density comparison records the exact 299800/299792.458 ratio.
The implementation and ordinary reader composition use the new c convention.

**The LBG and LAE density assets are not uniform in redshift.** Their nodes are
[2.38,2.60,2.83,3.07,3.29], with spacings [0.22,0.23,0.24,0.22]. Consequently the
required strict reader rejects setup for lya(lbg), lbg and lae. Those three
reader-to-density/weight/noise comparisons are blocked, not passed. The assignment
simultaneously requests rejecting nonuniform density grids and comparing these
assets through that reader. Resolving this inconsistency requires a planning/user
decision; no spacing tolerance was weakened, grid changed or alternative nodes
substituted. QSO's known negative spline row is rejected at the exact recorded
query. LBG's known negative-row query is blocked earlier by its raw grid error;
that earlier rejection and the original negative values are recorded explicitly.
Both populations' SNR inputs were evaluated successfully despite this limitation.

## Installed artifact

Built one successful wheel into the fresh `V/dist/`, then installed that exact
artifact in fresh `V/wheel-env` with templates, cosmology and survey extras.
No pre-existing FishHighz install or system-site-packages was used.

- Artifact: `fishhighz-0.1.0.dev0-py3-none-any.whl`
- SHA-256: `3d7f2b9aac75b0a0689e3212d867158e330f2db06ae4fe0010d6a1b52081951b`
- Python 3.13.15; NumPy 2.5.3; SciPy 1.18.1; Astropy 8.0.1.
- All **32** source/wheel/installed production modules match byte-for-byte, as do
  wheel/installed METADATA and WHEEL. Origins are inside the fresh environment.
  NumPy is the sole unconditional dependency; optional extras are retained.
- From `/tmp` with `-I`, all four retained examples and the new raw-input example
  pass. Generated raw fixtures use temporary directories outside the checkout.
  Scalar and two-bin Fisher oracles, Astropy geometry, the active optional-import
  blocker and source/installed R1 range cases pass.
- `build.log`/`wheel-install.log` retain initial sandbox failures;
  `build-approved.log`/`wheel-install-approved.log` retain successful commands.
  `probe.py`, `base_oracle.py` and `range_regression.py` preserve executed probes.

The synthetic example reports combined data F_AA=137.31062241474916 and
prior-assisted marginalized errors [0.0853390753, 0.4999410578, 0.4999331298]
for [A,b_low,b_high]. Each bin makes three preparation P3D calls (two auxiliary
autos plus fiducial covariance) and two P1D calls per forest (auxiliary plus
vectorized final nodes). With 12 nodes in batches of five, each bin makes
15 numerical P3D calls per derivative run; no Jacobian callback is used there.
The analytic-route test separately exercises Jacobian callbacks. These numbers
are synthetic diagnostics, not physical DESI-2 forecasts or performance claims.

## Review boundary

No real full suite was requested or run; it is not a failed/outstanding gate.
The actual blocked legacy density comparison is identified above. No commit,
push, agent dispatch, governance update or Step 12 work occurred. INI translation,
CLI, serialization, overlapping-bin treatment and disk streaming remain deferred.
Stopping for user and independent review.
