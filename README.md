# FishHighz

FishHighz computes Fourier-space Gaussian covariances and Fisher forecasts for
high-redshift galaxy and Lyman-alpha forest surveys. It retains covariance
between selected auto- and cross-power spectra and supports external models,
Kaiser/BAO signals, survey geometry, instrumental response, and sampling noise.

Two interfaces serve different calculations:

- The native INI interface prepares and runs the implemented BAO accuracy
  prescription from survey tables and bundled or user-supplied inputs.
- The scientific Python API exposes fields, parameters, models, derivatives,
  integration grids, covariance, and Fisher assembly for more general studies.

[Documentation](docs/index.md) · [Getting started](docs/getting-started.md) ·
[INI reference](docs/user/ini.md) · [API reference](docs/api/index.md)

## Installation

Python 3.11 or newer is required. From this repository:

```bash
python -m pip install .
```

NumPy is the only unconditional runtime dependency. Install optional preparation
dependencies for the bundled standalone DESI Run-2 example with:

```bash
python -m pip install '.[camb,templates,survey]'
```

Optional extras are selected by the operation being performed:

| Extra | Purpose |
| --- | --- |
| `templates` | FITS template input and spline preparation (Astropy/SciPy) |
| `cosmology` | Optional Astropy cosmology preparation |
| `camb` | CAMB background and growth preparation |
| `survey` | Survey-table interpolation with SciPy |
| `compiled` | Optional Numba Fisher contractions |
| `docs` | Sphinx documentation builds |
| `dev` | Tests, Ruff checks, and distribution builds |

The native interface uses packaged inputs and does not require Vega or
lyaforecast installations. Historical compatibility examples have separate
reference-package requirements, described in the documentation.

## Minimal forecast

The bundled recipe is accessible from an installed package:

```python
from fishhighz import Forecast

result = Forecast("desi2_accuracy.ini").run()
result.save("accuracy-desi2")
```

The equivalent command is:

```bash
fishhighz-forecast desi2_accuracy.ini --output accuracy-desi2
```

Choose a new output directory. Saved results include individual and joint
constraints, effective configuration, and input provenance.

To define another survey, start from the
[annotated INI](examples/desi2_accuracy_annotated.ini) and consult the
[configuration guide](docs/user/ini.md). Numerical overrides or changed inputs
require their own assessment; a prescription name does not qualify a modified
calculation automatically.

The [Python guide](docs/user/python.md) explains preparation, supplied models
and readers, execution controls, and results. The
[synthetic research example](examples/research_bao_forecast.py) demonstrates
explicit scientific API choices without running a real-survey study.

## Scientific scope and active research

Calculations use independent redshift bins, a common volume per bin, and
fiducial covariance and weights held fixed during mean derivatives. Selecting
only some spectra does not remove spectra needed to construct their covariance.

The [research baseline](docs/research/RESEARCH_BASELINE.md) defines the adopted
S2–S4 prescription and its numerical qualification. The
[weighting decision](docs/research/FOREST_WEIGHTING_DECISION.md) records the
accepted early-lyaforecast recurrence and retained alternative; the
[deferred tests](docs/research/DEFERRED_SCIENTIFIC_TESTS.md) describe remaining
scientific questions. Finite-refinement and iteration checks do not establish
physical accuracy of reconstruction assumptions or a continuum limit.

The bundled DESI recipe retains only the QSO/lya(qso) spectra in bin 1 and all
15 spectra in bins 2–6. Historical studies retain their original selections
and conclusions in the [development and validation archives](docs/archive/index.md).

## Development and documentation

Contributor setup, quick checks, and distribution instructions are in the
[development guide](docs/development/contributing.md). Build the HTML manual with:

```bash
python -m pip install -e '.[docs]'
python -m sphinx -n -W --keep-going -b html docs docs/_build/html
```

Read the result at `docs/_build/html/index.html`. Read the Docs configuration is
included; hosting must be connected separately. Real-survey validation and
scheduler actions are separate from ordinary documentation and synthetic checks.

## License

See [LICENSE](LICENSE) for the GNU General Public License, version 3. Bundled
scientific input provenance and accompanying license notices are recorded in
[the data README](fishhighz/data/README.md) and its linked notices.
