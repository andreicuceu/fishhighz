# Step 4 independent review — public forecast interface

Status: implementation review passed. This does not constitute scientific
acceptance of the DESI-2 result or authorization to execute the real forecast.

## Result

The public boundary implements the requested lazy flow:

```python
from fishhighz import Forecast

forecast = Forecast("desi2_accuracy.ini")
prepared = forecast.prepare()
result = forecast.run()
result.save("accuracy-desi2")
```

Construction performs strict INI parsing only. The literal shipped recipe name
falls back to the bundled INI when no same-named file exists in the working
directory; an existing user file retains precedence. Background, template,
readers, survey records and covariance factors are created by `prepare()`.

Each selected individual spectrum is reconstructed with a one-spectrum
`PairSelection`, its required power closure, its own known-noise columns and a
new covariance factorization. It is therefore not extracted from a block of the
joint inverse covariance. Joint constraints continue to use the full covariance
among all selected spectra in the bin. Packed full-noise inputs are remapped by
required-pair identity, preserving the low-level external-provider path.

The result retains `FisherResult` objects, pair/bin identities, explicit
available, unavailable and excluded states, full forest convergence records,
and resolved run controls. JSON/NPZ output records the constraint order and
Fisher matrices and refuses to overwrite an existing directory. The console
entry point and `python -m fishhighz.cli` use the same implementation.

## Independent checks

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q \
  tests/test_public_forecast.py tests/test_forecast.py tests/test_results.py
79 passed in 1.42s

.venv/bin/ruff check fishhighz/public.py fishhighz/cli.py fishhighz/__init__.py \
  tests/test_public_forecast.py examples/desi2_accuracy.py
All checks passed!

.venv/bin/ruff format --check fishhighz/public.py fishhighz/cli.py \
  fishhighz/__init__.py tests/test_public_forecast.py examples/desi2_accuracy.py
5 files already formatted

.venv/bin/python -m fishhighz.cli --help
passed

git diff --check
passed
```

The synthetic run also verifies one joint three-spectrum preparation followed
by three independent one-spectrum preparations. A separate full-noise control
checks the closure mapping for a cross-spectrum.

## Remaining boundary

No CAMB solve or six-bin DESI-2 forecast was run. Installed-wheel independence
and numerical comparisons with saved cosmology, survey, weight, `A`,
`P_pixel`, and 78-individual/six-joint evidence remain Step 5 work. The final
real-data calculation remains explicitly authorization-gated.
