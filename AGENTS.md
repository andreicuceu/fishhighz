# FishHighz contributor guide

FishHighz is a standalone Fourier-space covariance and Fisher forecasting package
for Lyα forest and galaxy clustering. The forecast formalism is documented in
`docs/methods/`; read the relevant page before changing scientific code.

## Layout

- `fishhighz/` (direct layout, no `src/`): `public.py` (`Forecast`, results),
  `cli.py`, `survey_config.py` (INI parser and preparation), `accuracy.py`
  (prescription defaults), `forecast.py` and `survey.py` (`BinSpec`, `prepare_bin`,
  `run_forecast`), `geometry.py`, `response.py`, `weights.py`, `noise.py`,
  `covariance.py`, `derivatives.py`, `fisher.py`, `results.py`, `fields.py`,
  `parameters.py`, `grids.py`, `forest_integration.py`, `full_shape.py`,
  `bao_marginalized.py`, `cosmology.py`, `magnitude.py`, `resources.py`.
- `fishhighz/survey_tables.py`: source-density and forest S/N tables
  (`DensityTable`, `SNRTable`, `SourceDensity`, `ForestSNR`, `CellHistogram2D`,
  `sample_forest_tables`).
- `fishhighz/models/`: `protocols.py` and `external.py` (external P3D/P1D
  providers), `kaiser.py`, `templates.py`, `p1d.py`, `biases.py`.
- `fishhighz/kernels/`: array-only numerical kernels (covariance, Fisher,
  derivatives, response, Kaiser, templates, full-sum and integrated forest
  weights); the Fisher contraction and the integrated forest-weight passes are
  Numba-compiled. Numba and SciPy are required dependencies.
- `fishhighz/data/`: bundled INI, CAMB configurations, FITS templates, DESI-2
  density and S/N tables; provenance and licenses in `fishhighz/data/README.md`
  and `fishhighz/data/licenses/`.
- `tests/`, `examples/`, `docs/` (MyST Markdown, Sphinx), `scripts/`
  (`check.sh`, `check_docs.py`, `check_integrated_forest_refinement.py`).

## Invariants

- Kernels accept explicit float64 arrays, scalars and integer indices, with no
  filesystem access, logging or configuration objects in numerical loops. Keep
  their arithmetic order stable: no fastmath, no reordering of reductions, no
  parallel loops. A change that alters results at the rounding level is a
  scientific change and must be documented as such.
- Covariance, forest weights and noise are fiducial and fixed during mean
  derivatives; no covariance-derivative information is added.
- The forest moments $I_1,I_2,I_3$ are full-sample sums over magnitude; the
  cumulative magnitude-prefix recurrence of lyaforecast is not supported. Adaptive
  weights raise if convergence is not confirmed; never substitute the last
  iterate or another prescription.
- Bundled data and templates keep their provenance and license notices; do not
  alter them or add data without recording origin, version and license.
- Numerical tolerances (`rtol`, finite-difference steps, quadrature orders,
  roundoff thresholds) are not changed silently; each change states the reason
  and the effect on the forecasts below.
- Keep dependencies minimal; heavy or optional ones (CAMB, Astropy) are
  imported lazily and live in extras.

## Validating a change to scientific code

1. Run `./scripts/check.sh` (pytest, `ruff check`, `ruff format --check`) on a
   compute node with single-threaded BLAS (the script sets `OMP_NUM_THREADS=1`).
2. For changes that can affect numerical output, rerun four representative
   forecasts from the same environment and CAMB version as the saved reference,
   and compare `results.npz` and `settings.json` against the reference run:
   - DESI-2 SRD `lya_qso_lbg_lae` 15x2pt, central forest sources
     (bundled recipe), and with `forest_source_integration = integrated`;
   - a `full_shape` INI and a `bao_marginalized` INI from
     `/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/forecasts/full_shape`
     and `.../forecasts/bao`.
   Expect bit-identical arrays and effective settings. If not, report the
   maximum relative difference of each Fisher matrix and BAO error, and the
   reason, before accepting the change.
3. Add focused regression tests for bug fixes, with synthetic inputs.

## Forecast selection rule

For every forecast, exclude bin-1 correlations involving LBG or LAE galaxies.
Bin-1 correlations involving lya(lbg) are permitted: the six lya(qso), QSO and
lya(lbg) auto- and cross-spectra may be retained in bin 1, including their joint
covariance. Bins 2-6 retain all 15 spectra. This is a forecast selection, not a
plot cut; keep existing result files with their original selection.

## Code style

- Four-space indentation, `snake_case` functions, modules and variables,
  `PascalCase` classes, NumPy-style docstrings for public objects, pytest, and
  the configured Ruff rules (line length 88).
- Use explicit float64 for scientific calculations and inspect numerical
  conditioning. Do not regularize singular Fisher matrices or covariance blocks;
  report rank and null directions instead.
- Separate file and configuration handling and validation from numerical
  functions. Prefer small functional interfaces and batched, vectorized
  evaluation; optimize measured bottlenecks only. No abstraction for
  hypothetical future use.
- Use scientific names (`k_grid`, `mu_grid`, `redshift_mask`), state units and
  array shapes in docstrings, and comment non-obvious numerical choices and
  physical assumptions.
- Tests are self-contained, deterministic and independent of NERSC paths; use
  analytic limits, finite-difference convergence and synthetic arrays.

## NERSC execution

Run pytest, forecasts and the integrated-source weight solve (about 1e7
pixel-magnitude nodes per forest field and bin) on a compute node, not on a
login node, with `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`.
Slurm submission requires the user's explicit approval. Write outputs to
`$SCRATCH`. The CAMB solve (about a minute) is cached in
`$FISHHIGHZ_CACHE_DIR/camb`; `FISHHIGHZ_CAMB_CACHE=0` disables the cache. Treat
CFS reference forecasts as read-only. Build the documentation with
`make -C docs html`.
