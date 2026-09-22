# Covariance, derivatives, and Fisher results

For selected spectra A=(i,j), B=(m,n), Gaussian covariance is

$$C_{AB}=\frac{P^{\mathrm{tot}}_{im}P^{\mathrm{tot}}_{jn}
+P^{\mathrm{tot}}_{in}P^{\mathrm{tot}}_{jm}}{N_{\mathrm{modes}}}.$$

Power inputs have axes `(node, required_pair)`, covariance and Cholesky factors
`(node, selected_pair, selected_pair)`, and mean derivatives
`(node, selected_pair, global_parameter)`. Fisher matrices have axes
`(global_parameter, global_parameter)` in registry order. Fixed fiducial
covariance and weights enter $F_{ab}=\sum D_a^T C^{-1}D_b$; no covariance
derivative information is added.

See [covariance](../methods/covariance.md), [derivatives](../methods/external.md),
and [Fisher inference](../methods/fisher.md) for validation, stencils, priors,
singularity diagnostics and optional compilation. Singular Fisher information
is retained without a pseudoinverse or automatic regularization. Invalid
covariance blocks fail with cell context rather than receiving a diagonal floor.

## `fishhighz.covariance`

```{autofunction} fishhighz.covariance.combine_observed_power
```

```{autofunction} fishhighz.covariance.gaussian_covariance
```

## `fishhighz.derivatives`

```{autoclass} fishhighz.derivatives.Stencil
```

```{autoclass} fishhighz.derivatives.ColumnDiagnostic
```

```{autoclass} fishhighz.derivatives.CallCount
```

```{autoclass} fishhighz.derivatives.DerivativeResult
```

```{autofunction} fishhighz.derivatives.evaluate_derivatives
```

```{autoclass} fishhighz.derivatives.ConvergenceColumn
```

```{autoclass} fishhighz.derivatives.ConvergenceResult
```

```{autofunction} fishhighz.derivatives.check_convergence
```

## `fishhighz.fisher`

```{autofunction} fishhighz.fisher.factor_covariance
```

```{autofunction} fishhighz.fisher.fisher_from_factors
```

```{autofunction} fishhighz.fisher.fisher_matrix
```

## `fishhighz.results`

```{autoclass} fishhighz.results.FisherDiagnostics
```

```{autoclass} fishhighz.results.FisherResult
:members: conditional_errors, marginalized_covariance, marginalized_errors, correlations, fix_except
:undoc-members:
```

```{autofunction} fishhighz.results.diagonal_prior
```

```{autofunction} fishhighz.results.combine_results
```

