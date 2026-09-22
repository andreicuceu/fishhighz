# Standalone configuration corrections — 2026-09-22

The simplified native DESI Run-2 INI expands to exactly the same cosmology,
survey, model, input policies, numerical controls, observed fields and selected
pairs as the original expanded native INI. No scientific data files, model
prescription, default integration settings or forecast selections changed.

## Changes

- `[prescription] name = accuracy` supplies the current explicit native model,
  input-policy and integration defaults, and the constant survey conventions.
  Provenance records its `early-lyaforecast-2026-09-18` revision and all expanded
  effective settings. An optional revision must match that supported revision.
- Native schema version 1 remains supported without a named prescription.
  Supported explicit overrides retain their existing validation; unknown options
  and unsupported scientific labels are rejected. Named prescriptions derive
  `num_z_bins` from `z_edges` when omitted. `selected = all` expands in declared
  field order; the bundled bin-1 three-spectrum selection remains explicit.
- Weight tolerance may be supplied in either of its existing sections. A lone
  supplied value synchronizes the other representation; explicit conflicts fail.
  The native default remains `1e-5`, rather than the historical stopping
  dictionary's separate `1e-4` value.
- Magnitude integration uses survey-wide limits. Field-level limits are optional;
  legacy matching values remain accepted and differing values are rejected as an
  unsupported per-field selection. Density normalization with the `survey`
  policy inherits the global limits; the `none` policy retains no normalization
  bounds. No new sample-selection or normalization physics was introduced.
- Configuration provenance records every scientific `FieldConfig` attribute,
  including densities, forest exposure/pixel/wavelength values and bias arrays.
  The original INI bytes are hashed when parsed; their immutable identity remains
  available after the source changes or is removed.

## Files and validation

Changed `fishhighz/accuracy.py`, `fishhighz/survey_config.py`,
`fishhighz/data/desi2_accuracy.ini`, `tests/test_survey_config.py`, and the native
INI inventory entry in `tests/test_resources.py`. The original expanded INI is
preserved as `tests/data/desi2_accuracy_expanded.ini` for direct regression tests.

Commands run from the package root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest tests/test_survey_config.py tests/test_resources.py -q
.venv/bin/ruff check fishhighz/survey_config.py fishhighz/accuracy.py tests/test_survey_config.py tests/test_resources.py
.venv/bin/ruff format --check fishhighz/survey_config.py fishhighz/accuracy.py tests/test_survey_config.py tests/test_resources.py
```

Results: **37 passed in 1.62 s**, Ruff passed, all four Python files formatted.
Coverage includes exact expanded/compact configuration equality, identical small
synthetic prepared geometry, quadrature, galaxy densities and forest inputs,
explicit quadrature changes, tolerance overrides/conflicts, unknown choices,
magnitude-selection rejection/inheritance and field metadata completeness.
The initial synthetic comparison incorrectly compared an array-containing
geometry record with scalar equality; replacing the test assertion with exact
array/scalar comparisons resolved that test-only failure.

No real CAMB calculation, real forecast, Slurm action, commit or governance edit
was performed by this assignment. Independent coordinator validation and the
configuration review are reported separately.
