# Step 5 independent review — installation and bounded reproduction

Status: the non-forecast Step-5 evidence passed independent review. This is not
scientific acceptance and does not authorize the real six-bin DESI-2 run.

## Reproduction result

The native preparation replay against the preserved accuracy-example settings
reproduces all saved quantities that can be reconstructed without CAMB or a
forecast:

- all six magnitude partitions and quadrature nodes/weights match exactly;
- all 11 saved forest-weight vectors, update counts, `A`, and `P_pixel` match
  exactly using the packaged density/SNR inputs and saved auxiliary `P`/`B`;
- analytic/tabulated biases and the template-growth identity match exactly;
- 30 independently inferred damping-reference sigma8 values are identical at
  `0.30969156196323777`, and reconstructed damping widths agree within
  `1.1215717662161686e-16` relative.

The saved settings do not contain sampled density arrays, so an elementwise
density comparison is unavailable. The applicable evidence is instead the
byte-identical packaged tables, exact reader policies, and exact magnitude and
weight replay. No live CAMB comparison was made.

## Installation result

The rebuilt wheel SHA-256 is
`3bd13b1a16f7fba99987bd222e9766ab081a0e7c71f63b1d1c8e48c3152689d1`.
Its 34 data files comprise the 29 scientific assets, native INI, data README,
and three provenance/license files. Wheel metadata contains the standalone
public usage and optional dependencies.

From `/tmp`, Python 3.13 with `-S` loaded FishHighz from the isolated wheel
installation. Neither `lyaforecast` nor `vega` was importable. The self-contained
injected calculation returned three individual and one joint result, with joint
rank 2, while package-resource access found the native INI, both 12-file SNR
sets, and the bundled FITS template.

## Independent checks

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python .validation/step05-installation/check_saved_intermediates.py
partition/node/quadrature/weights/A/P_pixel/growth/bias differences: 0.0
damping reference relative spread: 0.0
damping width maximum relative difference: 1.1215717662161686e-16

<isolated Python 3.13> -S \
  .validation/step05-installation/check_installed.py
counts=[3, 1, 0], rank=2, neighbors=[None, None]

.venv/bin/python -m pytest -q <Step-1--4 focused set>
118 passed

.venv/bin/python -m pytest -q
3 failed, 1487 passed, 25 skipped

git diff --check
passed
```

The three full-suite failures are outside the standalone implementation: two
historical `test_step12_performance.py` fixtures construct
`prepared.p3d=None`, and the plotting test lacks optional Matplotlib. Four
obsolete tests for the removed long accuracy tutorial were correctly replaced
with delegation and import-boundary tests; restoring removed tutorial internals
would contradict Step 4.

## Authorization boundary

The selected six-bin DESI-2 forecast, live CAMB reproduction, and comparison of
all 78 individual plus six joint results at relative tolerance `1e-8` were not
run. They remain explicitly authorization-gated by the plan. No Slurm action,
commit, or push was performed.
