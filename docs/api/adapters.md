# Prepared bins and survey adapters

The [survey tutorial](../methods/survey.md) documents preparation, batching,
raw count normalization and independent redshift-bin combination. `prepare_bin`
returns immutable fiducial arrays; `run_bin` returns a zero-prior `BinRun`;
`run_forecast` combines bins and adds the supplied prior once. Bins share one
`ParameterRegistry` instance, preserve caller order, and must have nonoverlapping
redshift interiors. Batch size changes transient derivative storage, not the
fixed full covariance-factor storage.

`DensityReader` and `SNRReader` require the optional SciPy `survey` extra. Strict
readers reject extrapolation and invalid interpolation values. Explicit
`LegacyDensity`/`LegacySNR` policies retain their provenance; their floors and
clamps are not the strict-reader defaults. The lyaforecast intrinsic adapter is
an optional callable bridge and does not make lyaforecast a base dependency.

## `fishhighz.survey`

```{eval-rst}
.. autoclass:: fishhighz.survey.ForestInput
```

```{eval-rst}
.. autoclass:: fishhighz.survey.BinSpec
```

```{eval-rst}
.. autoclass:: fishhighz.survey.PreparedBin
```

## `fishhighz.forecast`

```{eval-rst}
.. autofunction:: fishhighz.forecast.prepare_bin
```

```{eval-rst}
.. autoclass:: fishhighz.forecast.BinRun
```

```{eval-rst}
.. autoclass:: fishhighz.forecast.ForecastRun
```

```{eval-rst}
.. autofunction:: fishhighz.forecast.run_bin
```

```{eval-rst}
.. autofunction:: fishhighz.forecast.run_forecast
```

## `fishhighz.adapters.legacy_inputs`

```{eval-rst}
.. autoclass:: fishhighz.adapters.legacy_inputs.DensityReader
   :members: query, local_galaxy_density
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.adapters.legacy_inputs.SNRReader
   :members: query, variance
   :undoc-members:
```

```{eval-rst}
.. autofunction:: fishhighz.adapters.legacy_inputs.sample_forest_readers
```

## `fishhighz.adapters.legacy_compat`

```{eval-rst}
.. autoclass:: fishhighz.adapters.legacy_compat.LegacyDensity
   :members: sample
   :undoc-members:
```

```{eval-rst}
.. autoclass:: fishhighz.adapters.legacy_compat.LegacySNR
   :members: sample
   :undoc-members:
```

```{eval-rst}
.. autofunction:: fishhighz.adapters.legacy_compat.sample_legacy_forest
```

## `fishhighz.adapters.lyaforecast`

```{eval-rst}
.. autoclass:: fishhighz.adapters.lyaforecast.IntrinsicP3D
   :special-members: __call__
```

