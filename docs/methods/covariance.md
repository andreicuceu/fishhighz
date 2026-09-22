# Gaussian covariance

This synthetic example computes covariance, not a Fisher forecast or a validated
survey model. The supplied signal already includes any field response, and the
noise columns explicitly include cross noise. Known noise belongs in the total
power even when subtracted from the mean.

```python
import numpy as np
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.grids import gauss_legendre_grid
from fishhighz.covariance import combine_observed_power, gaussian_covariance

selection = PairSelection(
    [ObservedField("A", "galaxy", "a"), ObservedField("B", "galaxy", "b")],
    selected=[("A", "B")],
)
grid = gauss_legendre_grid([0.1, 0.3], k_order=2, mu_order=3, h_fid=0.7)
# Required order is AA, AB, BB even though only AB is selected.
signal = np.tile([3.0, -1.0, 8.0], (grid.n_k * grid.n_mu, 1))
noise = np.tile([1.0, -0.2, 1.0], (grid.n_k * grid.n_mu, 1))
total = combine_observed_power(signal, noise)
volume_fid = 1000.0  # (Mpc/h_fid)^3, supplied and fixed
mode_counts = volume_fid * grid.q_mode
covariance = gaussian_covariance(total, mode_counts, selection)
assert covariance.shape == (6, 1, 1)
np.testing.assert_allclose(covariance[:, 0, 0], (4.0 * 9.0 + 1.2**2) / mode_counts)
```

`combine_observed_power` adds exact matching `(node,required_pair)` arrays only;
it supplies no response, shot noise, smoothing, or implicit off-diagonal terms.
`gaussian_covariance` accepts total power directly in `selection.required_pairs`
order and returns owned C-contiguous float64 `(node,selected_pair,selected_pair)`
blocks. Inputs are unchanged. No matrix couples different nodes or redshift bins.
Counts can be fractional; the kernel divides by them once, without further
volume, quadrature, or conjugate-mode factors. Covariance units are
`(Mpc/h_fid)^6`. All arrays describe the fiducial model with fixed survey weights;
there is no theta argument or covariance-derivative information term.

Validation checks the full active field matrix one node at a time. Auto powers
must be nonnegative; an exactly zero auto requires an exactly zero row/column.
For positive autos, correlations R remove field-amplitude scales before checking
PSD. The allowed negative eigenvalue tolerance is
`64*eps64*n_positive*max(1,max(abs(eigvalsh(R))))`; each absolute correlation
must also be at most `1 + 64*eps64*n_positive`. Diagnostics identify invalid
nodes/fields. No jitter, clipping, noise floor, or regularization is applied.
Valid singular total powers and covariance blocks are returned without an
invertibility claim. Fisher factorization (below) assesses selected-block solvability.

The private `_gaussian_covariance_kernel` in `fishhighz.kernels.covariance`
accepts validated numeric power/count arrays, the four lookup tables, and a
preallocated output. It fills one triangle and mirrors it, using one extra
node-vector. The public wrapper remains the safe default. Consecutive node
slices may be computed independently and concatenated. Nonfinite addition or
covariance arithmetic raises an error; amplitudes or units are never silently
adjusted to avoid overflow.

