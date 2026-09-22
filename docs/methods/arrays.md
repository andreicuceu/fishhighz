# Fields, parameters, and Fourier quadrature

The example below prepares and validates arrays; it does not compute covariance
or a Fisher forecast. IDs are opaque strings. Shared physical labels never merge
samples or tie parameters. Select pairs using two-ID tuples or integer pairs.

```python
import numpy as np
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.grids import gauss_legendre_grid
from fishhighz.parameters import (
    Parameter,
    ParameterRegistry,
    ParameterBinding,
    gather_local,
    map_jacobian,
)
from fishhighz.models.protocols import validate_p3d, validate_p3d_jacobian

pairs = PairSelection(
    [
        ObservedField("A", "galaxy", "galaxy_model"),
        ObservedField("forest(qso)", "forest", "lya_model", background="qso"),
    ],
    selected=[("A", "forest(qso)")],
)
# Cross-only selection still requires AA, AB, BB for covariance.
assert pairs.required_pairs.tolist() == [[0, 0], [0, 1], [1, 1]]
grid = gauss_legendre_grid(
    [0.02, 0.1, 0.3],
    k_order=2,
    mu_order=3,
    h_fid=0.67,
)
registry = ParameterRegistry([Parameter("shared_amplitude", 1.0, "target")])
binding = ParameterBinding(
    registry,
    ["x", "y"],
    {"x": "shared_amplitude", "y": "shared_amplitude"},
)
local = gather_local(registry.fiducials, binding.local_to_global)


# Arbitrary plain functions may instead wrap an external model.
def p3d(theta_local, z, k, mu, pairs):
    sign = np.where(pairs[:, 0] == pairs[:, 1], 1.0, -1.0)
    shape = (1 + k[:, None] + mu[:, None] ** 2) * sign[None, :]
    return (2 * theta_local[0] + 3 * theta_local[1]) * shape


def jacobian(theta_local, z, k, mu, pairs):
    base = p3d([0.5, 0.0], z, k, mu, pairs)
    return base[:, :, None] * np.array([2.0, 3.0])


args = (local, 2.4, grid.k_flat, grid.mu_flat, pairs.required_pairs)
n_node = grid.n_k * grid.n_mu
power = validate_p3d(p3d(*args), n_node, len(pairs.required_pairs))
local_jac = validate_p3d_jacobian(
    jacobian(*args),
    n_node,
    len(pairs.required_pairs),
    len(local),
)
global_jac = map_jacobian(local_jac, binding.local_to_global, len(registry.ids))
mean = power[:, pairs.selected_to_required]
assert mean.shape == (12, 1)
assert global_jac.shape == (12, 3, 1)  # tied derivatives sum: 2 + 3 = 5
```

`PairSelection` preserves field and selected-pair order, canonicalizing each
pair to i <= j. `required_pairs` follows increasing i, then j.
`im`, `jn`, `in_`, and `jm` have shape `(n_selected,n_selected)` and index
`required_pairs` for A=(i,j), B=(m,n). A provider must supply every required
spectrum; missing cross powers are never reconstructed or treated as zero.
Mappings use owned, read-only, C-contiguous int64 arrays. Registries and grids
likewise copy inputs and protect their stored float64 arrays from accidental writes.
Fixed inputs, such as BAO widths and cuts, stay outside the free parameter registry.

Grid axes `k`, `mu` have integration weights `w_k` (dk) and `w_mu` (dmu).
Flattening `(n_k,n_mu)` uses C order, with mu fastest. `weights` is the flattened
product rule; `q_mode = k_flat**2 * weights / (2*pi**2)` gives
`N_modes = V_fid*q_mode`, counting both conjugate hemispheres for real, even
spectra with mu in [0,1]. k uses h_fid/Mpc; P3D and volume use (Mpc/h_fid)^3.
`h_fid` is fixed metadata, never an automatic coordinate rescaling.

For custom positive rules, use
`IntegrationGrid(k, w_k, mu, w_mu, k_min=..., k_max=..., h_fid=...)`
from `fishhighz.grids`. Axes must be strictly increasing and inside the cuts;
endpoints are allowed. Weight sums must equal k_max-k_min and 1 within relative
1e-12 (zero absolute tolerance). No sorting, normalization, clipping, or inferred
weights occur. Polynomial exactness does not establish survey-forecast convergence;
choose quadrature orders through later scientific convergence checks. Nodes are
not bin-averaged bandpowers. Padded model domains and arbitrary provider sampling
are separate from the fixed integration grid and add no forecast modes.

`P3D`, `P1D`, and `P3DJacobian` in `fishhighz.models.protocols` are structural
callable protocols. Providers evaluate one scalar z and local parameter vector;
P3D uses paired k/mu nodes and explicit pair indices. Output validators return
owned C-contiguous float64 arrays with exact shapes `(node,pair)`, `(node,)`,
and `(node,pair,local)`, respectively. Signed cross powers and derivatives are
preserved; complex, bool, object/string, nonfinite, and mismatched outputs fail.

P3D returns clustering before instrument response and known sampling noise.
The provider/adapter owns AP and cosmological physics and handles any preexisting
response/noise explicitly to prevent double counting. P1D is independent:
`p1d(theta_local, z, k_parallel_velocity)` uses s/km inputs and km/s outputs,
validated by `validate_p1d`. Custom P3D never implies or requires custom P1D.
Adapters own unit conversions, keeping h_fid fixed when physical h varies.

The Fisher convention keeps covariance and survey weights fixed at
the fiducial model and differentiate only the predicted mean. Parameter changes
do not move integration nodes, cuts, weights, or mode counts. This preparation example does not compute covariance or derivatives.

