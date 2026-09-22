# Step 3 independent review — native survey INI and preparation

Status: implementation review passed. This is not scientific acceptance of the
DESI-2 forecast and does not authorize a real forecast calculation.

## Scientific result

The native version-1 INI and preparation boundary reproduce the configured
six-bin survey construction without importing lyaforecast, Vega, or historical
validation modules. The observed samples retain their five distinct identities,
while the two forest samples share Ly-alpha physics. Bin 1 selects the three
Ly-alpha(QSO)/QSO spectra and bins 2--6 select all 15 spectra; the covariance
closure is constructed by `PairSelection` rather than inferred downstream.

The corrected density policies agree with the established accuracy example:
both forests use `z_norm_min=2.15` and the survey magnitude interval, QSO uses
the redshift normalization without a magnitude mask, and LBG/LAE use neither.
The magnitude partition and composite quadrature are now production utilities.
The early-lyaforecast weighting coordinates `(2.4 deg^-1, 0.00035 s/km)` are
required INI values passed explicitly to each `ForestInput`.

Template-growth and damping-reference normalizations remain separate. Native
preparation requires finite background metadata for both reference redshifts
and a positive finite `sigma8_damping_reference`, checks the redshifts against
the INI, and uses that sigma8 directly. There is no fallback to an unrelated
background attribute or unit normalization.

## Independent checks

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q \
  tests/test_survey_config.py tests/test_resources.py tests/test_import.py \
  tests/test_survey_imports.py tests/test_survey_composition.py \
  tests/test_cosmology.py tests/test_biases.py
35 passed in 1.59s

.venv/bin/ruff check <Step-3 sources and focused tests>
All checks passed!

.venv/bin/ruff format --check <Step-3 sources and focused tests>
11 files already formatted

git diff --check
passed
```

Source inspection also confirmed that package directories are materialized one
file at a time, so the SNR tables do not depend on directory support in
`importlib.resources.as_file` under Python 3.11. Relative user paths remain
anchored to the defining INI.

## Remaining boundary

This review covers parsing and assembly of immutable `BinSpec` inputs only. It
does not establish CAMB numerical equivalence, execute forest-weight iteration,
construct Fisher results, or compare the 78 individual and six joint DESI-2
outputs. Those remain later-step evidence, and the real forecast remains gated
by explicit authorization.
