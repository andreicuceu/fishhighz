# S5.2 example and public-API audit; S5.3 focused checks

## Scientific result

The new standalone example produces a finite rank-2 forecast for the two target
parameters `ap` and `at`. Its synthetic joint uncertainties are
`(0.0468765060314, 0.0225581179462)`. The forest auto, forest–galaxy cross and
galaxy auto are selected together, so their Gaussian covariance is retained.
The joint data Fisher is not the sum of three isolated-spectrum Fishers: the
maximum absolute matrix difference divided by the largest absolute joint-Fisher
element is `1.35625090428`. This distinction is saved explicitly rather than
inferred from differing marginalized errors.

The result is an API demonstration, not a DESI-2 forecast and not new evidence
for the selected research prescription. It uses generated arrays only: no
lyaforecast, Vega, CAMB, captured bundle or input file is read. Its only optional
runtime requirement is `fishhighz[templates]`, because SciPy prepares the
in-memory synthetic spline; subsequent model and Fisher evaluation use NumPy.

## Implemented example and API documentation

- `examples/research_bao_forecast.py` exposes `run()` and is directly executable.
  It uses the public template, Kaiser, survey preparation and forecast APIs.
  The registry contains only `ap` and `at`; forest bias, forest beta, galaxy bias
  and the EdS growth rate are fixed. No prior is applied.
- The smooth template scaling is fixed at identity and `ap`/`at` act on the full
  wiggle component. Fixed per-field damping widths enter `KaiserModel`, whose
  cross width squared is the arithmetic mean of the two auto width squares.
- The forest response converts a physical resolving power defined by wavelength
  FWHM to Gaussian sigma. Early-lyaforecast weights use the explicit auxiliary
  angular/velocity mode `(2.4 deg^-1, 0.00035 s/km)`. Preparation confirms
  convergence after 10 updates from candidate 5 with `rtol=1e-4`, minimum 3
  updates, 3 stable transitions and cap 96.
- The returned report includes joint and individual AP errors, joint and summed
  individual Fisher matrices, AP correlation/rank, weighting status and stopping
  controls, selected/required pair counts, noise convention, Fourier-node and
  derivative-batch counts, and damping-width conventions.
- Public docstrings now describe the available forecasting package, auxiliary
  queries for all three compatible weighting methods, adaptive failure behavior,
  prepared weight/diagnostic metadata, and `ForecastRun` result ownership. The
  full-sum solver documentation now correctly refers to two variants. Two stale
  validation error messages now enumerate all accepted weighting methods. These
  are documentation/message changes only; no numerical implementation changed.

`tests/test_research_bao_example.py` is a scientific contract rather than an
output snapshot. It requires the full three-spectrum selection, adaptive early
weights and controls, rank-2 finite AP results, and a covariance-coupled joint
Fisher distinct from the sum of isolated Fishers.

## Checks actually run

All numerical commands used `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` and
`MKL_NUM_THREADS=1` on the login node.

- Source execution with `.validation/s1-s3-env/bin/python
  examples/research_bao_forecast.py`: passed. Exact JSON is saved in
  `.validation/s5-example/output.json` for the installed-wheel comparison.
- `pytest -q tests/test_research_bao_example.py tests/test_full_sum_weights.py
  tests/test_forecast.py tests/test_results.py`: **84 passed in 1.92 s**.
- `pytest -q tests/test_weight_range.py tests/test_weighting_w12.py
  tests/test_forecast_imports.py tests/test_import.py`: **17 passed in 5.33 s**.
- Ruff check and format check passed for the example, test and affected source
  files; `git diff --check` passed for the same bounded set.

No real survey forecast, new scientific prescription, weighting default,
negative-input policy, compatibility calculation, D1/D2 sensitivity, Slurm job,
environment installation, commit or push was performed. The parent coordinator
owns the isolated S5 wheel build, installed example comparison and final S5
governance record.
