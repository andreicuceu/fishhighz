# Step 5 — installation and bounded reproduction

Status: the non-forecast Step-5 checks passed. This report records installation,
synthetic public-interface execution, and comparisons against the saved
accuracy-example intermediates. It is not scientific acceptance and does not
authorize the real DESI calculation.

## Focused synthetic checks

The existing Step-1--4 tests were rerun with one numerical thread:

```text
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q \
  tests/test_public_forecast.py tests/test_survey_config.py tests/test_resources.py \
  tests/test_import.py tests/test_forecast_imports.py tests/test_survey_imports.py \
  tests/test_survey_composition.py tests/test_cosmology.py tests/test_biases.py
43 passed in 8.04s
```

These controls cover strict native-schema errors and unsupported legacy INIs,
relative and package resource resolution, the five observed-field identities,
the three-versus-15 pair selection, covariance closure, shared AP registry,
background/template metadata mismatches, and a one-bin injected `Forecast`
calculation. No source correction was required by these checks.

After replacing four obsolete tests that asserted removed accuracy-tutorial
internals, the new example contract and Step-1--4 focused suite were run with:

```text
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q \
  tests/test_desi2_accuracy_example.py tests/test_public_forecast.py \
  tests/test_survey_config.py tests/test_resources.py tests/test_import.py \
  tests/test_forecast_imports.py tests/test_survey_imports.py \
  tests/test_survey_composition.py tests/test_cosmology.py tests/test_biases.py \
  tests/test_forecast.py tests/test_results.py
118 passed in 2.27s
```

The replacement tests verify delegation to
`Forecast("desi2_accuracy.ini").run()`, saving to `accuracy-desi2`, returning
the result, and the absence of legacy imports or low-level forecast execution
in the short example.

The full suite before this correction reported 7 failures, 1485 passed and 25
skipped. After the correction, the exact rerun was:

```text
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q
3 failed, 1487 passed, 25 skipped in 176.60s
```

The remaining three failures are the already-known unrelated cases: the two
`tests/test_step12_performance.py::test_independent_accuracy_check` fixtures
construct `prepared.p3d=None`, and
`tests/test_three_profile_plots.py::test_render_saved_weights_and_coefficients`
requires the unavailable optional Matplotlib dependency. No Step-4-induced
accuracy-example failure remains.

## Saved-intermediate comparisons

The reference was `.validation/desi2-examples/accuracy/settings.json`, which
contains the saved model values, magnitude grids, forest convergence states,
weights, `A`, `P_pixel`, `P` and `B`. The following bounded checks used only
the packaged input tables and the recorded intermediates; they did not invoke
CAMB or run a Fisher forecast.

The replay was run exactly with:

```text
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step05-installation/check_saved_intermediates.py
```

which printed:

```text
partition_lengths (163, 163, 163, 165, 169, 166)
partition_max_abs_difference 0.0
magnitude_node_max_abs_difference 0.0
magnitude_weight_max_abs_difference 0.0
weight_count 11
weight_max_abs_difference 0.0
A_max_relative_difference 0.0
P_pixel_max_relative_difference 0.0
all_update_counts_match True
template_growth_max_relative_difference 0.0
bias_max_relative_difference 0.0
damping_reference_min 0.30969156196323777
damping_reference_max 0.30969156196323777
damping_reference_spread 0.0
damping_reference_relative_spread 0.0
damping_width_max_relative_difference 1.1215717662161686e-16
```

* The native `breakpoints` and composite quadrature reproduced all six saved
  magnitude partitions exactly. The partition lengths were
  `(163, 163, 163, 165, 169, 166)` and the maximum absolute node difference
  was `0.0`.
* The native `prepare_forest_weights` recurrence was replayed for all 11 saved
  forest preparations using the saved `P` and `B` auxiliary values and the
  packaged density/SNR readers. Every saved weight array matched exactly;
  the maximum absolute weight difference was `0.0`. The relative differences
  in every saved `A` and `P_pixel` were `0.0`, and all 11 update counts matched.
* The native analytic and tabulated bias utilities reproduce every saved bias
  entry exactly in the recorded decimal representation. The saved template
  growth values satisfy `(sigma8/sigma8_template)**2` exactly. For every saved
  bin and field, the replay first derives
  `sigma8_damping_reference = 3.26*sigma8/(Sigma_perp*sqrt(r))` from the
  saved transverse width and reconstruction factor, then fails closed unless
  all derived values agree within `1e-12` relative spread. The 30 derived
  values have minimum and maximum `0.30969156196323777`, spread `0.0`, and
  relative spread `0.0`. Using their common mean, the damping widths satisfy
  the stated per-field reconstruction and
  `Sigma_parallel=(1+f) Sigma_perp` prescription to a maximum relative
  difference of `1.12e-16`. The common normalization is therefore a checked
  consistency inference, not a selected single record.
* The saved JSON does not contain the underlying density arrays, so a direct
  elementwise density comparison is unavailable. The packaged raw tables and
  their reader policy/partition behavior were checked instead; the prior
  Step-1 byte-identity evidence remains the applicable asset comparison.

These checks establish exact recurrence and preparation agreement for the
quantities actually saved. They do not independently recompute CAMB
`H(z)`, `D_M(z)`, `sigma8(z)` or `f(z)`: a live CAMB solve was deliberately
not run under the login-node and Step-5 scope constraints.

## Installed-wheel check

The project `.venv` lacks `setuptools.build_meta`, so its no-isolation build
failed before package compilation with `Backend 'setuptools.build_meta' is not
available`. Using the already available NERSC `vega_test` interpreter's local
setuptools (without downloading dependencies), the actual wheel build passed:

```text
/global/homes/a/acuceu/.conda/envs/vega_test/bin/python -m build \
  --wheel --no-isolation --outdir .validation/step05-installation/dist
Successfully built fishhighz-0.1.0.dev0-py3-none-any.whl
sha256: 3bd13b1a16f7fba99987bd222e9766ab081a0e7c71f63b1d1c8e48c3152689d1
```

The wheel contains 34 `fishhighz/data/` files: all 29 scientific assets, the
bundled native INI, the data README and the three license/provenance files.
It was installed with `pip --no-index --no-deps` into
`.validation/step05-installation/env-system`. From `/tmp`, a Python-3.13
interpreter was then run with `-S` and only the installed environment plus the
NERSC NumPy/SciPy site-packages directory added explicitly. This prevents the
environment's editable neighboring checkouts from participating;
`importlib.util.find_spec("lyaforecast")` and `find_spec("vega")` were both
`None`. The installed package materialized the native INI, all 12 QSO and 12
LBG SNR files, and the 28,800-byte FITS template through
`importlib.resources`.

The installation and outside-checkout calculation were run exactly with:

```text
.validation/step05-installation/env-system/bin/python -m pip install \
  --no-index --no-deps --force-reinstall \
  .validation/step05-installation/dist/fishhighz-0.1.0.dev0-py3-none-any.whl
cd /tmp && PYTHONDONTWRITEBYTECODE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/.validation/step05-installation/env-system/bin/python \
  -S /global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/fishhighz/.validation/step05-installation/check_installed.py
```

The same outside-checkout process ran a self-contained injected public API
calculation. Its minimal native INI, analytic background, density reader,
template and assertions are defined in
`.validation/step05-installation/check_installed.py`; it imports no Python
definitions from the checkout's tests or package source:

```text
fishhighz loaded from the installed wheel
{'counts': [3, 1, 0], 'rank': 2, 'neighbors': [None, None]}
```

No `lyaforecast` or Vega import was needed by this calculation.

After the README change, the wheel metadata was checked directly:

```text
unzip -p fishhighz-0.1.0.dev0-py3-none-any.whl \
  fishhighz-0.1.0.dev0.dist-info/METADATA
  contains the native Forecast usage, [camb,templates,survey] installation,
  bundled-resource description, CLI/module commands, and the explicit
  lyaforecast/Vega independence statement
  while retaining the separate compatibility-example lyaforecast requirement
git diff --check -- README.md
passed
```

## Explicitly not run

The selected six-bin DESI forecast and its comparison of 78 individual plus
six joint results were **not run**; that calculation awaits the separate
explicit authorization required by the plan. No live CAMB solve, Slurm action,
commit, push, or neighboring-checkout modification was performed.

Stop for independent and user review.
