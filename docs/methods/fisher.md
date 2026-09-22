# Fisher information and constraints

This example forecasts a supplied analytic mean response at fixed fiducial
covariance. It is not an external-model derivative engine or a validated survey
forecast. The covariance already includes mode counts and quadrature weights;
Fisher assembly adds no further normalization or covariance-derivative term.

```python
import numpy as np
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.grids import gauss_legendre_grid
from fishhighz.covariance import gaussian_covariance
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.results import FisherResult, diagonal_prior

selection = PairSelection([ObservedField("sample", "galaxy", "synthetic")])
grid = gauss_legendre_grid([0.1, 0.3], k_order=2, mu_order=3, h_fid=0.7)
n_node = grid.n_k * grid.n_mu
covariance = gaussian_covariance(
    np.full((n_node, 1), 4.0), 1000.0 * grid.q_mode, selection
)
# Supplied derivatives of the selected mean, in global parameter order.
# Known subtracted noise belongs in covariance, not automatically in this J.
jacobian = np.stack([np.ones(n_node), grid.mu_flat**2], axis=-1)[:, None, :]
factors = factor_covariance(covariance)
data_information = fisher_from_factors(jacobian, factors)
registry = ParameterRegistry(
    [
        Parameter("amplitude", 1.0, "target"),
        Parameter("shape_nuisance", 0.0, "nuisance"),
    ]
)
result = FisherResult(
    registry,
    data_information,
    prior_fisher=diagonal_prior(registry, {"shape_nuisance": 2.0}),
)
conditional = result.conditional_errors()  # fix every other parameter
marginalized = result.marginalized_errors()  # all registry parameters free
amplitude_error = result.marginalized_errors(["amplitude"])
fixed_shape_error = result.fix_except(["amplitude"]).marginalized_errors()
assert result.diagnostics.rank == 2
assert np.all(marginalized >= conditional)
```

`factor_covariance(C)` returns reusable owned float64 lower Cholesky factors
`(node,selected,selected)`. `fisher_from_factors(J,L)` solves each node by forward
substitution across all parameter columns, then accumulates `Y.T @ Y` into
`(global,global)`. `fisher_matrix(J,C)` composes the two operations. Subsequent
Jacobians can reuse L; the from-factors path never refactorizes it. Direct factors
must be finite, exactly lower triangular, and have positive diagonals; their
fiducial provenance is the caller's responsibility. Covariance factorization
checks numerical solvability in a normalized basis and rejects singular blocks,
even though Step 04 correctly allows them. It judges the selected covariance,
not the rank of its parent field-power matrix. Existing array inputs need only
one node's conversion/solve workspace, plus stored factors and output information;
node slices can be evaluated independently and their Fisher matrices summed.

`FisherResult` retains the exact registry metadata and separate owned read-only
float64 `data_fisher`, `prior_fisher`, and `total_fisher`. Priors are precision
matrices in the same units and fiducial basis, not covariance matrices. Correlated
PSD priors and rank-deficient information are allowed. Bounds and derivative
steps are metadata, never priors. `diagonal_prior(registry, {id: sigma})` assigns
only explicitly requested finite positive widths.

`combine_results([data_only_result1, data_only_result2], prior_fisher=prior)`
requires identical ordered parameter metadata and sums unmarginalized data before
adding the one final prior. It rejects inputs already carrying nonzero priors.
Do not marginalize shared nuisances separately by bin before combining data;
individual contributions may be singular while their sum is informative.

Conditional errors are `1/sqrt(F_total[ii])`, with infinity for an exact zero
diagonal. Marginalized covariance is obtained through factorization and identity
solves of the entire retained information. Requesting a subset afterward still
marginalizes over all other free parameters. `fix_except(ids)` instead fixes the
complement and retains a principal information submatrix with its metadata and
separate priors. IDs must be unique and explicit; target/nuisance roles never
silently select or fix parameters. `correlations(ids=None)` derives correlations
from the marginalized covariance.

For covariance and information validation, use R=M/d/d, where d=sqrt(diag(M));
exact zero information rows use d=1 in the declared units. Zero diagonals require
exact zero rows/columns. Symmetry is tested elementwise against
`64*eps64*n*max(1,abs(Rij),abs(Rji))`; roundoff asymmetry is averaged only in owned
copies. Eigenvalue tolerance is `t=64*eps64*n*max(1,max(abs(eigenvalues(R))))`.
PSD allows eigenvalues >= -t, numerical rank counts eigenvalues > t, and
factorization requires full rank with positive variances. No eigenvalue clipping,
jitter, observable projection, or pseudoinverse is performed.

`result.diagnostics` gives ordered `ids`, numerical `rank`, `condition` (infinite
if singular), `tolerance`, `scales`, `eigenvalues`, and `null_directions`. Each
null/near-null row has unit Euclidean length in coordinates
`x = scales * delta_theta`; physical displacements are
`null_directions / scales[None,:]`, not the unconverted eigenvectors. Signs and
bases within degenerate eigenspaces are not unique. Singular results remain
inspectable, but marginalized requests fail with rank/null context. Explicit
fixing or a finite prior may resolve the degeneracy. Nonfinite arithmetic fails
rather than producing apparently reliable errors.


## Optional Fisher contraction compilation

NumPy remains the default and the only required runtime dependency. For repeated
Fisher contractions, install `fishhighz[compiled]` (Numba and SciPy) and explicitly
set `FISHHIGHZ_FISHER_BACKEND=numba`. The default is
`FISHHIGHZ_FISHER_BACKEND=numpy`. Numba is imported only on an explicitly requested
float64 contraction. If it is unavailable, or inputs require dtype conversion,
the NumPy reference executes. Explicit NumPy underflow warning/error policies also
use the reference, since compiled loops do not implement `numpy.seterr`. Invalid
backend names raise `ValueError`.

The optional kernel uses one Fourier cell of solve workspace, float64, no
fastmath and no parallel loops. Numerical failures are replayed through the
NumPy reference for the original cell/parameter diagnostic. Cold compilation is
paid separately for each encountered array signature; no compilation cache is
written. Keep `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` for bounded
login-node checks. External model callables retain their Python interface.

Step 12 revision 4 timings and limitations are recorded in
[the performance report](../archive/reviews/step-12-performance-r4.md). These synthetic
measurements do not establish real-forecast timing or forest convergence.

