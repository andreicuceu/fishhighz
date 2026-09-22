# Zero-redshift cosmology normalization correction

The CAMB preparation boundary now permits the physical observer distance
`D_M(0) = 0` when a template-growth or damping normalization is specified at
redshift zero. Previously, the universal strictly positive distance condition
rejected these otherwise valid normalization choices.

`fishhighz/cosmology.py` still requires positive Hubble parameters and sigma8,
finite nonnegative distances, and strictly positive distances at positive
redshift. The survey-geometry implementation is unchanged. The correction does
not modify the calculation for the default positive reference redshifts.

Deterministic fake-CAMB tests cover each reference redshift individually at
zero, both at zero, angular and comoving distance result surfaces, exact growth
array mapping after CAMB ordering, both named sigma8 normalizations, negative
distances, nonfinite distances and zero distance at positive redshift.

Validation from the package root:

```text
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest tests/test_cosmology.py -q
16 passed in 0.44s (initial run; final rerun also passed all 16)
.venv/bin/ruff check fishhighz/cosmology.py tests/test_cosmology.py
All checks passed!
.venv/bin/ruff format --check fishhighz/cosmology.py tests/test_cosmology.py
2 files already formatted
```

No real CAMB solve, survey forecast, Slurm action or commit was performed.
