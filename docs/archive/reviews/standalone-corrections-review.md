# Standalone interface corrections and coordinator review

The three reported interface issues are corrected, and the native accuracy INI
is simplified without changing the effective DESI Run-2 prescription. The
coordinator dispatched three Astra agents at low reasoning effort, reviewed the
source and results, and ran the combined checks. A separate Astra review of the
configuration changes found no consequential issue.

## Scientific behavior and interface

- Saved results retain every effective field input, including target density,
  exposure count, pixel width, forest wavelength limits and tabulated bias. The
  original INI is hashed from the bytes parsed. Background, template and reader
  identities are captured during preparation, including every SNR-file digest.
  Package identities remain stable when resources are materialized temporarily;
  injected objects are labelled separately. Saving does not reread source files.
- Survey magnitude bounds define the common integration interval. A field may
  use that interval for density normalization or apply no normalization mask.
  Optional legacy field limits must equal the survey limits; different limits
  now raise an explicit error. Independent per-field selections are not added.
- CAMB preparation permits the observer limit `D_M(0)=0` for template-growth or
  damping normalization. Negative/nonfinite distances and zero distance at
  positive redshift remain invalid; survey geometry requirements are unchanged.
- `[prescription] name = accuracy` expands the adopted model, input policies,
  numerical controls and constant survey conventions. Supported overrides are
  validated and fully recorded. One weighting-tolerance override updates its
  equivalent setting; conflicting explicit values are rejected. `num_z_bins`
  can be inferred from edges, and `selected = all` expands canonical pairs.
  Bin 1 still retains only the three QSO/forest observables.

The bundled INI decreased from 4354 to 2054 bytes. The original fully explicit
INI is retained in `tests/data/desi2_accuracy_expanded.ini`. Both forms give
identical effective fields, bins, scientific settings and small synthetic
prepared inputs. The original explicit schema remains supported. The README
documents the compact interface, overrides and magnitude conventions; stale
CLI/serialization deferral statements were corrected. The source distribution
includes the original-INI regression fixture.

Implementation and independent review reports:

- [Configuration correction](standalone-configuration-corrections.md)
- [Independent configuration review](standalone-configuration-review.md)
- [Saved-input correction](standalone-result-provenance-corrections.md)
- [Zero-redshift correction](standalone-zero-redshift-correction.md)

## Coordinator validation

Evidence is retained under
`.validation/standalone-corrections-20260922T193357Z/`. All numerical commands
used `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`.

The combined command from the package root was:

```bash
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q -p no:cacheprovider \
  tests/test_desi2_accuracy_example.py tests/test_public_forecast.py \
  tests/test_survey_config.py tests/test_resources.py tests/test_import.py \
  tests/test_forecast_imports.py tests/test_survey_imports.py \
  tests/test_survey_composition.py tests/test_cosmology.py tests/test_biases.py \
  tests/test_forecast.py tests/test_results.py
```

Result: **150 passed in 4.14 s**, exit status zero. The final log is
`focused-tests-final.log`; it retains environment MUNGE messages emitted after
the successful test summary. Ruff lint and format checks passed for the four
modified implementation modules and four affected test modules;
`git diff --check` passed.

The existing read-only intermediate replay was run with the same thread limits:

```bash
.venv/bin/python .validation/step05-installation/check_saved_intermediates.py
```

All six magnitude partitions, nodes and quadrature weights were identical to
the saved accuracy evidence. All eleven forest-weight vectors, `A` and
`P_pixel` matched exactly, and every update count matched. The saved bias and
template-growth checks also matched exactly; the existing damping consistency
check differed by at most `1.12e-16`. See `saved-intermediates.log`. This uses
saved auxiliary powers and background quantities; it is not a fresh CAMB or
Fisher calculation.

Fresh wheel and source distributions were built offline from an isolated copy
of the final source using the existing `vega_test` interpreter and
`python -m build --no-isolation`. The build location and logs are retained.
The final wheel was installed with `pip --no-index --no-deps --target` into
`installed-site/`. From `/tmp`, Python with `-S` and explicit installed/dependency
paths ran the self-contained `check_installed.py` fixture. Neither lyaforecast
nor Vega was importable. The check loaded the compact bundled INI and all
resources, ran a small injected forecast (three individual results, one joint
result, joint rank two), and saved the expanded settings successfully. See
`installed-check.log`.

Every file under `fishhighz/` in the wheel was byte-compared with final package
source/data. The original explicit-INI fixture was also verified in the final
source archive. `final-sha256.json` binds the modified files and distributions:

```text
wheel: 3dfdb4c6483868ed80aec1ac6130104e3ffcdc37d248a367eb9c988c318a5180
sdist: 5fac2bd06934aaf138bff451b309ef26342b4f6ae677b50f818e1ced19942548
```

## Limits and disposition

These changes and checks establish interface correctness and preservation of
the tested preparation quantities. The earlier six-bin numerical reproduction
evidence is retained unchanged; no new real-survey forecast or CAMB solve was
run. The broader ordinary suite was not repeated; its previously reported
unrelated fixture/optional-plotting failures are outside this correction.

Input identities do not serialize arbitrary injected callables. External CAMB
configuration files must remain stable during preparation; arbitrary auxiliary
CAMB files are not snapshotted. The pre-existing data redistribution limitation
remains unchanged. No scientific prescription adoption, roadmap advancement,
Slurm action, commit or push was performed. Ready for user review.
