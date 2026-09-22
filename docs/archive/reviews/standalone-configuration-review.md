# Independent review of native configuration simplification

The compact accuracy INI preserves the previous effective scientific inputs.
No consequential correction is requested. This review checks configuration and
synthetic preparation; it does not establish new real-survey forecast evidence.

Before implementation, I captured the original bundled parsed configuration in
session memory. Comparing every previously recorded provenance value with the
new configuration found **zero differences**, including all survey, model,
numerical and input-policy settings, field identities, bins and pair selections.
New provenance adds the named prescription, parsed INI identity and the omitted
field quantities. Separately, the regression fixture compares the complete
parsed fields and settings against the original expanded INI, and synthetic
preparation compares geometry, quadrature, galaxy densities and forest-weight
options. The native tolerance remains 1e-5, distinct from the historical
`accuracy.STOPPING` default of 1e-4.

Code inspection and independent short controls confirm:

- An explicit weighting tolerance in either supported location determines both
  resolved entries; incompatible explicit values fail validation. Preparation
  reads the resolved numerical tolerance.
- Field magnitude bounds can be omitted or match the common survey interval.
  Unequal bounds fail explicitly for both forest and galaxy fields, including
  one-sided changes. Density-normalization masking remains a separate policy.
- `selected = all` expands in configured field order. Reordering all five fields
  still produces the canonical 15-pair ordering and complete covariance closure.
  The bin-1 selection retains only its three required spectra.
- All `FieldConfig` scientific attributes enter configuration provenance, and
  model, survey, numerical and policy defaults are expanded before recording.
  Saved-file identity and serialization receive a separate coordinator review.
- README documentation states the prescription, override and magnitude-selection
  rules and describes the implemented CLI and serialization interfaces.

Validation used one numerical thread:

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest tests/test_survey_config.py -q
28 passed in 1.41s
.venv/bin/ruff check fishhighz/accuracy.py fishhighz/survey_config.py tests/test_survey_config.py
All checks passed!
.venv/bin/ruff format --check fishhighz/accuracy.py fishhighz/survey_config.py tests/test_survey_config.py
3 files already formatted
```

The additional independent controls used temporary INIs, parsing and
`PairSelection` only; no CAMB solve or forecast was run. Reviewed SHA256 values:

```text
131b7d864226521daaca9ef72b4aff129f5973dc1ad55d7e661ae885617280b4  fishhighz/accuracy.py
b361eed975a754864736ee66e16f6dfd09d1aa049c7a7f2fb978115d90288fce  fishhighz/survey_config.py
c2053fdfe3b9c89fb62bdf527fb9f053b24ee17416a6350629f80db8c6eb94e6  fishhighz/data/desi2_accuracy.ini
a2a340a27300e02e879e8e20d315422c03d9195d26f65a338399587b8529d54a  tests/test_survey_config.py
6eab9c97ea4c77c39db17029d6f42eb996aea139e383270f7de9006ccdc5893d  tests/data/desi2_accuracy_expanded.ini
```

Only this review report was written during the review; implementation files were
not modified by the reviewer.
