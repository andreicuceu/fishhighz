# Step 4 — public forecast interface

Status: bounded implementation complete and ready for independent/user review.
This handoff is implementation evidence, not scientific acceptance and does not
authorize Step 5 or a real DESI forecast.

## Implemented interface

FishHighz now exposes the requested compact entry point:

```python
from fishhighz import Forecast

result = Forecast("desi2_accuracy.ini").run()
result.save("accuracy-desi2")
```

`Forecast` parses and validates the native INI in its constructor without
invoking CAMB, loading FITS, reading survey tables, or constructing covariance
state. `prepare()` lazily resolves the configured CAMB/template resources (or
caller-injected background, template and reader factories), assembles the
reviewed `PreparedSurvey`, and exposes fixed `PreparedBin` records. Relative
external CAMB/template/input paths remain anchored to the defining INI. The
literal shipped usage name `Forecast("desi2_accuracy.ini")` resolves the bundled
recipe when no same-named user file exists, while an existing relative file is
always honored.

`SurveyResult` retains the existing `FisherResult` objects for six joint bin
constraints and each selected individual spectrum. Joint constraints use the
full selected-spectrum covariance through `run_forecast`. Each individual
spectrum is rebuilt as a separate `BinSpec` with its own covariance closure,
noise and Cholesky factors before its Fisher contraction; no joint inverse
covariance block is sliced. The result records bin/pair identities, AP
constraints, rank-derived `available`/`unavailable` statuses, explicit
bin-selection `excluded` statuses, convergence metadata, resolved settings and
the combined Fisher result. `save()` writes JSON identities/statuses and NPZ
Fisher matrices without overwriting an existing output directory.

The `fishhighz-forecast` console entry point and `python -m fishhighz.cli`
boundary call the same `Forecast(...).run().save(...)` implementation. The
former long accuracy example is reduced to the public API usage example. The
existing low-level `prepare_bin`, `run_bin`, `run_forecast`, external model
providers and result classes remain available unchanged.

Resolved settings retain the actual run controls (`batch_size`, resolved
derivative `step_scale`, and `numerical`). Forest convergence provenance keeps
the complete solver record, including status, reason, update/state counts,
candidate, last-step metrics, confirmation, forward residual and retained
states. Packed `BinSpec.full_noise` is remapped from the original required-pair
order to each one-spectrum covariance closure.

## Checks

```text
.venv/bin/pytest -q tests/test_public_forecast.py
6 passed

.venv/bin/pytest -q tests/test_public_forecast.py tests/test_forecast.py tests/test_results.py
79 passed

.venv/bin/pytest -q tests/test_public_forecast.py tests/test_import.py \
  tests/test_forecast_imports.py tests/test_survey_config.py
25 passed

.venv/bin/ruff check fishhighz/public.py fishhighz/cli.py fishhighz/__init__.py \
  tests/test_public_forecast.py examples/desi2_accuracy.py
All checks passed!

.venv/bin/ruff format --check fishhighz/public.py fishhighz/cli.py \
  fishhighz/__init__.py tests/test_public_forecast.py examples/desi2_accuracy.py
5 files already formatted

git diff --check
passed
```

The synthetic tests use injected background, template and density readers. They
check lazy construction, bundled-name resolution and preservation of an
existing relative INI, prepared-state exposure, the three-spectrum joint
versus three one-spectrum preparation calls, packed full-noise remapping,
complete convergence/run-control serialization, retained Fisher results and
CLI delegation. No CAMB solve, DESI input forecast, Slurm action, commit or push
was performed.

## Evidence boundary and limitations

The new facade has not been used for the six-bin DESI calculation. Numerical
reproduction against the saved 78-individual/6-joint evidence, installed-wheel
independence, and comparisons of growth, damping widths, densities, weights,
`A` and `P_pixel` remain Step 5 evidence. The implementation does not claim
that the native recipe is scientifically accepted. The existing public result
serialization is intentionally an interchange summary rather than a loader;
the in-memory `FisherResult` objects remain the authoritative detailed result.

Stop for independent/user review.
