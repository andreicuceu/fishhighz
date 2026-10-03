# FishHighz

Fourier-space Gaussian covariances and Fisher forecasts for Lyman-alpha forest
and galaxy clustering. FishHighz retains the full covariance between selected
auto- and cross-power spectra of the tracers in a redshift bin, with survey
geometry, instrumental response, aliasing and pixel noise for the forest, and
Kaiser/BAO or external P3D signal models. The documentation is built from `docs/`
(see below).

## Installation

Python 3.11 or newer. NumPy, SciPy and Numba are required and installed
automatically. From this directory:

```bash
python -m pip install -e .
python -m pip install -e '.[camb,templates]'   # bundled DESI-2 recipe
```

| Extra | Purpose |
| --- | --- |
| `camb` | CAMB background, sigma8(z) and f(z) |
| `templates` | Astropy, to read the FITS power templates |
| `cosmology` | Astropy cosmology as an alternative background (Python API) |
| `dev` | pytest, Ruff, build |
| `docs` | Sphinx documentation build |

The CAMB transfer solve is cached in `$FISHHIGHZ_CACHE_DIR/camb` (default
`~/.cache/fishhighz/camb`); set `FISHHIGHZ_CAMB_CACHE=0` to disable the cache.
Neither Vega nor lyaforecast is needed.

## Quick start

Command line (omit the INI to use the bundled DESI-2 recipe):

```bash
fishhighz-forecast my-survey.ini --output my-forecast
fishhighz-forecast --output accuracy-desi2 --joint-only
```

Python:

```python
from fishhighz import Forecast

result = Forecast("my-survey.ini").run()  # Forecast() loads the bundled recipe
result.save("accuracy-desi2")  # new directory: settings.json + results.npz
```

`result.joint` and `result.individual` give the constraints per bin on the two BAO
parameters, with the joint covariance of the selected spectra and with each
spectrum alone. Start a custom survey from
[`examples/desi2_accuracy_annotated.ini`](examples/desi2_accuracy_annotated.ini).
Long calculations belong on compute nodes, with `OMP_NUM_THREADS=1`.

## Forecast modes

`[model] mode` in the INI selects the parameters:

- `bao` (default): per-bin AP dilation parameters `ap` and `at`
  (parallel and transverse) of the wiggle component; biases and weights fixed.
- `bao_marginalized`: $(\alpha_{\rm iso},\phi)$ with the tracer biases and the Lyα
  beta marginalized.
- `full_shape`: independent dilation of the wiggle and smooth components and
  $f\sigma_8$, with the nuisances marginalized, over $k_{\min}<k<k_{\max}$.

The prescription `[prescription] name = accuracy` has revisions
`early-lyaforecast-2026-10-01` (default; source density piecewise constant in the
table cells, forest sources at one representative redshift per bin) and
`early-lyaforecast-integrated-2026-10-02`
(`[input policies] forest_source_integration = integrated`; integration over
source redshift and forest pixel inside the bin).

## Forest weighting

Forest weights follow the adaptive `early_lyaforecast` recurrence with full-sample
magnitude moments, iterated to a confirmed convergence ($10^{-5}$ relative in the
INI default). The McDonald & Eisenstein alternative is available through the
Python API (`weights.prepare_forest_weights(..., method="mcdonald")`). See
[`docs/methods/forest-weighting.md`](docs/methods/forest-weighting.md).

## External P3D models

Any package can supply the intrinsic P3D, an independent P1D and an optional
Jacobian as plain callables (`fishhighz.models.protocols`), wrapped in
`models.external.P3DProvider`. Covariance, mean derivatives (analytic, numerical
or mixed) and the Fisher matrix are then computed with the same machinery as the
built-in models. See [`examples/external_forecast.py`](examples/external_forecast.py)
and [`examples/builtin_forecast.py`](examples/builtin_forecast.py).

## Development

```bash
python -m pip install -e '.[dev,docs,templates,camb,cosmology]'
./scripts/check.sh                 # pytest, ruff check, ruff format --check
make -C docs html                  # strict Sphinx build, then scripts/check_docs.py
```

Contributor notes are in [`AGENTS.md`](AGENTS.md) and
`docs/development/contributing.md`; the HTML manual is written to
`docs/_build/html/index.html`.

## License

GNU General Public License, version 3 ([LICENSE](LICENSE)). The provenance and
licenses of the bundled survey tables and templates are recorded in
[`fishhighz/data/README.md`](fishhighz/data/README.md).
