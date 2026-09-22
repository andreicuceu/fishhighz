# Fields, pairs, parameters, and grids

The [array tutorial](../methods/arrays.md) gives a complete synthetic example.
Field and selected-pair order follow the caller; each pair is canonicalized to
`i <= j`. Required pairs are lexicographic in field indices. The covariance
requires powers beyond the selected observables. Parameter bindings explicitly
sum tied local derivatives into global registry columns.

All scientific arrays use float64. Flattening `(n_k, n_mu)` uses C order with
mu fastest. Wavenumber is in `h_fid/Mpc`, volume and P3D in `(Mpc/h_fid)^3`.
Invalid shapes, domains, identities and quadrature normalizations raise
`ValueError`; inputs are not sorted, clipped, or normalized automatically.

## `fishhighz.fields`

```{eval-rst}
.. autoclass:: fishhighz.fields.ObservedField
```

```{eval-rst}
.. autoclass:: fishhighz.fields.PairSelection
```

## `fishhighz.parameters`

```{eval-rst}
.. autoclass:: fishhighz.parameters.Parameter
```

```{eval-rst}
.. autoclass:: fishhighz.parameters.ParameterRegistry
```

```{eval-rst}
.. autoclass:: fishhighz.parameters.ParameterBinding
```

```{eval-rst}
.. autofunction:: fishhighz.parameters.gather_local
```

```{eval-rst}
.. autofunction:: fishhighz.parameters.map_jacobian
```

## `fishhighz.grids`

```{eval-rst}
.. autoclass:: fishhighz.grids.IntegrationGrid
```

```{eval-rst}
.. autofunction:: fishhighz.grids.gauss_legendre_grid
```

