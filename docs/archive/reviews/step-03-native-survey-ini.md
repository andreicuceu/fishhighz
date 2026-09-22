# Step 3 — native survey INI and preparation handoff

Status: bounded implementation complete and ready for independent/user review.
This is implementation evidence, not scientific acceptance and does not
authorize Step 4.

## Scientific result

FishHighz now has a strict versioned native INI schema (`fishhighz-native-survey`,
version 1) and a small public boundary:

```python
from fishhighz import parse_survey_ini, prepare_survey
config = parse_survey_ini("survey.ini")
prepared = prepare_survey(config, background=background, template=template,
                          readers=readers)
```

Parsing does not invoke CAMB. Preparation requires caller-prepared background
and template objects, and accepts deterministic injected density/SNR readers;
it therefore assembles `BinSpec` objects without requiring a live cosmology
solve. If readers are omitted, the bundled Step-1 resources are materialized
through context-managed `importlib.resources` paths and passed through the
strict readers plus the explicitly named accuracy adapters.

The review correction is included in the recipe and construction path: the
forest fields (`lya(qso)` and `lya(lbg)`) use `z_norm_min=2.15` and the survey
magnitude bounds `(16.1, 26.75)`, QSO uses `z_norm_min=2.15` without a density
magnitude mask, and LBG/LAE use neither.  The two numerical weighting reference
coordinates are required INI values and are passed directly to `ForestInput`;
there is no production fallback to a module constant.

The final fail-closed correction requires injected backgrounds to expose finite
`template_growth_redshift` and `damping_reference_redshift` values matching the
INI, plus a positive finite `sigma8_damping_reference`.  The damping width uses
that required value directly; it no longer falls back to `background.sigma8` or
a literal normalization.

The bundled `fishhighz/data/desi2_accuracy.ini` records the actual S2–S4
integration state: 5000 deg2, six equal bins from 2.0 to 3.41, geometric
`1+z` evaluation redshifts, 128 × order-4 k intervals on [0.01, 0.5], mu and
volume orders 32, magnitude order 16, AP step `2.5e-4`, early-lyaforecast
weighting with final `rtol=1e-5`, minimum updates 3, stable transitions 3 and
cap 96, the `(2.4 deg^-1, 0.00035 s/km)` auxiliary reference, mixed squared
damping widths, wiggle-only `ap_at` parameters, forest beta 1.45, damping
amplitude 3.26, and the adopted reader policies. Bin 1 selects exactly
lya(qso) auto, lya(qso)×QSO and QSO auto;
bins 2–6 select all 15 spectra.

The parser preserves the observed identities `lya(qso)` and `lya(lbg)` while
using shared Ly-alpha physics. Each bin receives its own AP bindings into one
shared `ParameterRegistry`; `PairSelection` retains the complete covariance
closure required by its selected means. Preparation provenance retains input
semantics, normalization and redshift-width policy, interpolation and negative
policy, SNR smoothing/clamp/sentinel, magnitude boundaries/partition, weighting
status, selected and required pairs, and bin bounds/evaluation redshifts.

Existing lyaforecast INIs fail closed with `UnsupportedSchemaError` and are not
translated. Production preparation modules contain no imports of validation,
lyaforecast or Vega. The reusable magnitude partition/quadrature and adopted
accuracy constants now live in `fishhighz.magnitude` and `fishhighz.accuracy`;
the validation modules and accuracy example re-export/use these implementations
for historical compatibility.

## Files in this bounded step

- `fishhighz/survey_config.py`, `fishhighz/accuracy.py`, and
  `fishhighz/magnitude.py`: schema dataclasses, parser, preparation assembly,
  production controls and exact reusable magnitude logic.
- `fishhighz/data/desi2_accuracy.ini`: package-data native recipe.
- `fishhighz/resources.py` and `fishhighz/__init__.py`: context-managed
  per-file bundled resource materialization (including Python 3.11-safe SNR
  directory handling) and public parse/prepare exports.
- `fishhighz/validation/accuracy.py`,
  `fishhighz/validation/profile_definitions.py`, and
  `examples/desi2_accuracy.py`: minimal compatibility import changes.
- `tests/test_survey_config.py` and the package-resource inventory update in
  `tests/test_resources.py`: deterministic schema, policy, resource, identity,
  closure, quadrature, registry and synthetic assembly controls.

The pre-existing dirty Step-1/Step-2 files and evidence were preserved; no
forecast, CAMB solve, Slurm action, commit, or push was performed for Step 3.

## Checks

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q \
  tests/test_survey_config.py tests/test_resources.py tests/test_import.py \
  tests/test_survey_imports.py tests/test_survey_composition.py \
  tests/test_cosmology.py tests/test_biases.py
35 passed in 2.14s

.venv/bin/ruff check <Step-3 Python sources, example and focused tests>
All checks passed!
.venv/bin/ruff format --check <Step-3 Python sources and focused tests>
passed
git diff --check
passed
```

The correction checks additionally assert the actual `DensityReader` provenance
for all five configured fields, reject missing sections, partial numeric lists,
unsupported fixed labels and invalid domains, validate template/background
redshift and `h_fid` consistency, and exercise a zip-style SNR directory whose
individual files remain materialized for the lifetime of the context.
They also reject missing, nonfinite, nonpositive, and mismatched background
normalization metadata before any bin assembly.

The focused synthetic assembly used a small injected template, background and
reader set. It verified six `BinSpec` objects, exact three-versus-15 selection,
required-pair closure, distinct observed field IDs and shared registry identity.
The bundled reader construction was checked separately for all five fields and
both 12-table forest populations. No real forecast or cosmology preparation was
run.

## Limits and authorization

The schema and preparation boundary are ready for independent review. Numerical
equivalence to the full DESI forecast remains the preserved S2–S4 evidence; this
step adds no new CAMB or forecast-level result. The Step-4 `Forecast` facade,
CLI, serialization and result container remain intentionally unimplemented.

Stop for independent/user review.
