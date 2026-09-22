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

```{eval-rst}
.. autoclass:: fishhighz.models.protocols.P3D
   :special-members: __call__
```

```{eval-rst}
.. autoclass:: fishhighz.models.protocols.P1D
   :special-members: __call__
```

```{eval-rst}
.. autoclass:: fishhighz.models.protocols.P3DJacobian
   :special-members: __call__
```

```{eval-rst}
.. autofunction:: fishhighz.models.protocols.validate_p3d
```

```{eval-rst}
.. autofunction:: fishhighz.models.protocols.validate_p1d
```

```{eval-rst}
.. autofunction:: fishhighz.models.protocols.validate_p3d_jacobian
```

## `fishhighz.models.external`

```{eval-rst}
.. autoclass:: fishhighz.models.external.BoundParameters
```

```{eval-rst}
.. autoclass:: fishhighz.models.external.P3DProvider
```

```{eval-rst}
.. autoclass:: fishhighz.models.external.PreparedP3D
```

```{eval-rst}
.. autofunction:: fishhighz.models.external.evaluate_p3d
```

```{eval-rst}
.. autofunction:: fishhighz.models.external.evaluate_p1d
```

## `fishhighz.models.templates`

```{eval-rst}
.. autoclass:: fishhighz.models.templates.PowerTemplate
   :members: domain, evaluate
   :undoc-members:
```

```{eval-rst}
.. autofunction:: fishhighz.models.templates.prepare_template
```

```{eval-rst}
.. autofunction:: fishhighz.models.templates.load_template
```

## `fishhighz.models.kaiser`

```{eval-rst}
.. autoclass:: fishhighz.models.kaiser.Scaling
```

```{eval-rst}
.. autoclass:: fishhighz.models.kaiser.KaiserModel
   :special-members: __call__
```

## `fishhighz.models.p1d`

```{eval-rst}
.. autofunction:: fishhighz.models.p1d.p1d_floor
```

```{eval-rst}
.. autofunction:: fishhighz.models.p1d.default_p1d
```

## `fishhighz.models.biases`

```{eval-rst}
.. autoclass:: fishhighz.models.biases.LinearTabulatedBias
   :special-members: __call__
```

```{eval-rst}
.. autofunction:: fishhighz.models.biases.linear_tabulated_bias
```

```{eval-rst}
.. autofunction:: fishhighz.models.biases.analytic_density_bias
```

```{eval-rst}
.. autofunction:: fishhighz.models.biases.analytic_beta_rsd
```

```{eval-rst}
.. autoclass:: fishhighz.models.biases.AnalyticBias
   :members: set_density_bias_func, density_bias, beta_rsd, compute_bias
   :undoc-members:
```

