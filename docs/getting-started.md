# Getting started

FishHighz computes Fourier-space Gaussian covariances and Fisher constraints for
high-redshift surveys. Python 3.11 or newer is required. From the package checkout,
install the dependencies for the native survey interface:

```bash
python -m pip install -e '.[camb,templates,survey]'
```

NumPy is the only unconditional runtime dependency. The `camb` extra prepares the
background, `templates` reads and interpolates the signed FITS power template,
and `survey` reads the density and SNR tables. The bundled native calculation
requires neither lyaforecast nor Vega. The separate `cosmology` extra supports
the optional Astropy background adapter; `compiled` enables compiled Fisher
support. These extras serve distinct preparation or execution paths.

Constructing a forecast parses and validates its configuration:

```python
from fishhighz import Forecast

forecast = Forecast()  # bundled DESI Run-2 accuracy INI
config = forecast.config
```

When ready to perform the calculation, prepare the inputs and run:

```python
prepared = forecast.prepare()
result = forecast.run()
result.save("accuracy-desi2")  # destination must not already exist
```

`prepare()` performs CAMB, template, survey and fixed-covariance preparation;
`run()` also prepares automatically when needed. These are real numerical
calculations, not configuration checks. On shared computing systems use an
execution environment appropriate to the calculation and local policy.

The bundled recipe selects three QSO/forest spectra in bin 1 and all 15 spectra
in bins 2–6. Its [research prescription](research/RESEARCH_BASELINE.md) specifies
the scientific assumptions and qualification limits. Changing survey inputs or
numerical settings requires assessing the resulting calculation separately.

For a custom survey, download the [annotated INI](user/ini.md), edit a local
copy, and pass its path to `Forecast("my-survey.ini")`. Continue with the
[INI reference](user/ini.md), [Python lifecycle](user/python.md),
[result interpretation](user/results.md), or [command line](user/cli.md).
