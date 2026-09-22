# External models, templates, and Kaiser signal

External providers return intrinsic clustering `(node, required_pair)` before
instrument response and known sampling noise. Their Jacobians have axes
`(node, pair, local_parameter)`. P1D is an independent callable with velocity
wavenumber in s/km and power in km/s. Explicit output validators below are part
of the callable boundary; internal validation-study modules are not public API.
See [external models](../methods/external.md), [templates](../methods/templates.md),
and [Kaiser/BAO](../methods/kaiser.md) for ownership and complete examples.

Templates preserve `PKSB` as smooth power and `PK-PKSB` as wiggles. Preparation
requires SciPy; FITS loading additionally requires Astropy (`templates` extra).
Evaluation uses prepared arrays. Invalid domains or mismatched metadata fail;
there is no automatic extrapolation or new wiggle decomposition.

## `fishhighz.models.protocols`

```{autoclass} fishhighz.models.protocols.P3D
:special-members: __call__
```

```{autoclass} fishhighz.models.protocols.P1D
:special-members: __call__
```

```{autoclass} fishhighz.models.protocols.P3DJacobian
:special-members: __call__
```

```{autofunction} fishhighz.models.protocols.validate_p3d
```

```{autofunction} fishhighz.models.protocols.validate_p1d
```

```{autofunction} fishhighz.models.protocols.validate_p3d_jacobian
```

## `fishhighz.models.external`

```{autoclass} fishhighz.models.external.BoundParameters
```

```{autoclass} fishhighz.models.external.P3DProvider
```

```{autoclass} fishhighz.models.external.PreparedP3D
```

```{autofunction} fishhighz.models.external.evaluate_p3d
```

```{autofunction} fishhighz.models.external.evaluate_p1d
```

## `fishhighz.models.templates`

```{autoclass} fishhighz.models.templates.PowerTemplate
:members: domain, evaluate
:undoc-members:
```

```{autofunction} fishhighz.models.templates.prepare_template
```

```{autofunction} fishhighz.models.templates.load_template
```

## `fishhighz.models.kaiser`

```{autoclass} fishhighz.models.kaiser.Scaling
```

```{autoclass} fishhighz.models.kaiser.KaiserModel
:special-members: __call__
```

## `fishhighz.models.p1d`

```{autofunction} fishhighz.models.p1d.p1d_floor
```

```{autofunction} fishhighz.models.p1d.default_p1d
```

## `fishhighz.models.biases`

```{autoclass} fishhighz.models.biases.LinearTabulatedBias
:special-members: __call__
```

```{autofunction} fishhighz.models.biases.linear_tabulated_bias
```

```{autofunction} fishhighz.models.biases.analytic_density_bias
```

```{autofunction} fishhighz.models.biases.analytic_beta_rsd
```

```{autoclass} fishhighz.models.biases.AnalyticBias
:members: set_density_bias_func, density_bias, beta_rsd, compute_bias
:undoc-members:
```

