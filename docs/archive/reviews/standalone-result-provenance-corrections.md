# Saved survey settings and input identities

The forecast's saved JSON now contains the complete effective field settings and
named-prescription expansion provided by `SurveyConfig.provenance`, together with
`resolved_settings.inputs`. This additive metadata does not change forecast
arithmetic or the version-1 result matrix layout.

`SurveyConfig` captures the input INI identity and SHA256 from the bytes parsed.
`Forecast.prepare()` captures the CAMB configuration SHA256 before preparation,
uses the template loader's captured digest, and obtains density and every SNR-file
digest from the readers that supplied the samples. Package resources have stable
`package:` identifiers rather than temporary materialization paths. The resulting
`PreparedForecast.input_identity` is recursively immutable. Saving does not reopen
or hash input files.

Injected backgrounds, templates and readers are labelled `injected` or
`injected_factory`; supplied instances take precedence over factories. Existing
reader/template source digests are retained where available. Unavailable hashes
are explicit null values, and configured paths that were not used are not opened.
This records input identity, not a serialization of arbitrary injected callables or
a guarantee that a result can reconstruct them. External CAMB configuration files
must remain unchanged during preparation; this implementation does not snapshot
CAMB's configuration parser or its auxiliary files.

## Validation

- `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest tests/test_public_forecast.py -q`: **12 passed**.
- `.venv/bin/ruff check fishhighz/public.py tests/test_public_forecast.py`: passed.
- `.venv/bin/ruff format --check fishhighz/public.py tests/test_public_forecast.py`: passed.

Tests exercise saved target density, forest exposure count and tabulated bias;
deleted INI after construction; retained reader digest after file/reader mutation;
actual bundled density, SNR, template and CAMB-INI checksums; injected sources and
instance/factory precedence; and the existing own-covariance forecast checks.
All forecasts use small synthetic backgrounds/templates/readers. CAMB execution,
real-survey forecasts, Slurm operations and historical evidence changes were not
performed. The successful pytest process emitted environment MUNGE socket errors
after its summary; these did not affect its zero exit status.
