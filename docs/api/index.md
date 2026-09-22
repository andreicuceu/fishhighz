# Python API

Entries are selected explicitly and grouped by scientific operation. The native
forecast interface is exported from `fishhighz`; lower-level operations are
imported from the modules shown here. Internal kernels and historical validation
studies are not stable user interfaces. Units, array ordering and ownership are
shared with the [scientific methods](../methods/index.md).

Native preparation helpers in `fishhighz.accuracy` and `fishhighz.magnitude`,
and internal survey freezing/snapshot helpers, are implementation details
rather than additional supported configuration or modelling interfaces.

```{toctree}
:maxdepth: 1

forecast
arrays
models
survey
inference
adapters
```

## Frequently used interfaces

```{autosummary}
:nosignatures:

fishhighz.public.Forecast
fishhighz.survey_config.parse_survey_ini
fishhighz.forecast.prepare_bin
fishhighz.forecast.run_forecast
fishhighz.fields.PairSelection
fishhighz.models.kaiser.KaiserModel
fishhighz.covariance.gaussian_covariance
fishhighz.fisher.fisher_matrix
fishhighz.results.FisherResult
```
