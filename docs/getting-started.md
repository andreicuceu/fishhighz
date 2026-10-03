# Getting started

FishHighz requires Python 3.11 or newer. NumPy, SciPy and Numba are installed
with the package. From the package checkout:

```bash
python -m pip install -e '.[camb,templates]'
```

The `camb` extra computes the background and growth, and `templates` reads the
signed FITS power template. The `cosmology` extra provides an optional Astropy
background for the Python API. The native calculation needs neither lyaforecast
nor Vega.

Constructing a forecast parses and validates the configuration:

```python
from fishhighz import Forecast

forecast = Forecast()  # bundled DESI-2 accuracy INI
config = forecast.config
```

To calculate, prepare the inputs and run:

```python
prepared = forecast.prepare()
result = forecast.run()
result.save("accuracy-desi2")  # destination must not already exist
```

`prepare()` runs CAMB, reads the template and survey tables, computes the forest
weights and noise, and factors the fixed covariance. `run()` prepares
automatically when needed. Both perform real numerical calculations; run them on
a compute node at NERSC.

The bundled recipe selects three QSO-forest spectra in bin 1 and all 15 spectra in
bins 2-6. Its forest weights are the adaptive `early_lyaforecast` prescription
([note](methods/forest-weighting.md)). For a custom survey, copy the
{download}`annotated INI <../examples/desi2_accuracy_annotated.ini>`, edit it,
and pass its path to `Forecast("my-survey.ini")`. Continue with the
[INI reference](user/ini.md), [Python lifecycle](user/python.md),
[results](user/results.md) or [command line](user/cli.md).
