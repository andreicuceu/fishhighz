# FishHighz

FishHighz is a Python package for high-redshift Fourier-space Gaussian
covariances and Fisher forecasts. It provides field/pair identities, parameter
bindings, fixed integration grids, external callable models, signed linear-power
templates, intrinsic Kaiser/BAO models, survey geometry and response, the default
P1D prescription, forest weighting and sampling noise, fixed-covariance Fisher
assembly, named results, and Python survey orchestration. NumPy is the only
unconditional runtime dependency. Template, cosmology, survey-reader and compiled
Fisher support are explicit optional extras.

The recommended research prescription tested in S2--S4 is documented in
[RESEARCH_BASELINE.md](RESEARCH_BASELINE.md). It uses converged
`method="early_lyaforecast"` forest weights at the representative mode
`(k_transverse, k_parallel) = (2.4 deg^-1, 0.00035 s/km)`, physical FWHM
resolution, per-field crossed squared-width damping, and full wiggle inverse-AP
derivatives with fixed fiducial weights and covariance. The compact executable
[research_bao_forecast.py](examples/research_bao_forecast.py) exposes these
public API choices with synthetic inputs; it requires `fishhighz[templates]` and
can also be imported as `run()`.

Three [DESI Run-2 examples](examples/) read the real survey inputs:
[full compatibility](examples/desi2_full_compatibility.py),
[fixed compatibility](examples/desi2_fixed_compatibility.py), and the commented
[accuracy tutorial](examples/desi2_accuracy.py). Each computes individual and
joint BAO constraints in six independent redshift bins, retaining only the
Lyα(QSO) and QSO auto/cross spectra in bin 1 and all 15 spectra in bins 2–6.
The compatibility examples use the optional installed `lyaforecast` reference
package; the tutorial explicitly constructs the research-baseline FishHighz API
objects. See each script's docstring and `--help` for dependencies, input paths
and output files. No captured validation bundle is required.

The package API is more general than that tested prescription. Callers choose
fields, models, parameters, priors, selections and numerical controls explicitly.
Validation profiles under `fishhighz.validation` reproduce named historical
studies; they are not a general public profile factory.

## Development

Use Python 3.11 or newer. From this directory, create an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pip install -e '.[dev,templates]'
./scripts/check.sh
```

Source lives directly in `fishhighz/`. Package discovery includes only `fishhighz`
and its package descendants, with implicit namespaces disabled. Reinstall the
editable package after changing its layout. Imports from the checkout root alone
do not verify installation; the smoke test imports in an isolated subprocess
outside the checkout.

The runner uses the active environment and works when invoked by absolute or
relative path from another directory. It stops on a failed check and does not
install tools or change formatting. Its individual commands are:

```bash
python -m pytest
python -m ruff check .
python -m ruff format --check .
```

Ruff targets Python 3.11 with an 88-character line length. Enabled rules are E4
(imports), E7 (statement errors), E9 (runtime/syntax errors), F (Pyflakes), and I
(import sorting). Pytest, Ruff, and the build frontend are declared in the `dev`
extra. Setuptools >=68 is the isolated build requirement. The version is defined
only in `pyproject.toml`; installed versions can be read with
`importlib.metadata.version("fishhighz")`.

Ruff excludes generated `.validation/` evidence even before Git is initialized.

## Distributions

```bash
python -m build
```

This creates a wheel and source distribution in `dist/`. The retained local
acceptance helper builds from a clean source snapshot into a new output directory
and creates two fresh installation environments on every invocation:

```bash
python .validation/validate.py
```

The helper is local review evidence, not a distributed package tool. It records
exact wheel paths and hashes, rejects ambiguous artifact selection, and compares
installed package bytes with each selected wheel. It repeats those checks for a
wheel rebuilt from a separately extracted source archive. Run directories and
helper/probe snapshots remain under `.validation/r4-*`; no cleanup is needed
between runs. Isolated builds may download the declared build requirements.

The [Step 01 report](reviews/step-01.md) records the exact isolated interpreter,
commands, artifact hashes, and independent installation checks used for verification.
Workspace contributors can consult the [design](../../FISHHIGHZ_DESIGN.md) and
[roadmap](../../FISHHIGHZ_IMPLEMENTATION_PLAN.md); neither is needed to install,
import, or test the package.

## Lyaforecast reference capture

Step 02's maintained tool defaults to a quick scientific check: one unchanged,
full-resolution 15x2pt forecast compared with a previously accepted full bundle.
Supply the baseline, reference checkout, and scientific Python explicitly. The
destination is created exclusively and must not already exist:

```bash
python scripts/lyaforecast_baseline.py capture \
  --suite quick \
  --baseline .validation/baseline/<existing-full-bundle> \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python
```

`--suite quick` is the default, but is shown for clarity. It requires a compatible
full baseline and never falls back to a full run. A complete capture remains
available only as explicit opt-in; it runs all seven authoritative DESI-2 cases
and an independent 15x2pt repeat:

```bash
python scripts/lyaforecast_baseline.py capture \
  --suite full \
  --reference-checkout ../lyaforecast \
  --python ../lyaforecast/.validation/dev-env/bin/python
```

The default destination is a unique, ignored
`.validation/baseline/<UTC-timestamp>-<suffix>/` directory. Each case runs
serially in a fresh isolated subprocess with single-thread settings. Bundles
contain original/effective INIs, configuration/input/source snapshots and
inventories, Git and environment provenance, complete pickle and typed-JSON
results, logs, timings, and checksums. Quick bundles copy the full reference
result, compatibility fields, and numerical comparison so offline checks survive
relocation. Existing destinations are never resumed or overwritten.

Validate a bundle without rerunning forecasts, NumPy, CAMB, or the original
checkout:

```bash
python scripts/lyaforecast_baseline.py check \
  .validation/baseline/<UTC-timestamp>-<suffix>
python scripts/lyaforecast_baseline.py check \
  .validation/baseline/<UTC-timestamp>-<suffix> --require-suite full
```

The checker supports the accepted schema-v1 full evidence and identifies
inventory guarantees unavailable in that legacy schema. Real reference forecasts
are explicit validation work and are not run by `scripts/check.sh`; ordinary
pytest coverage uses synthetic inputs and stub workers only. Run the real full
suite only when it is specifically requested.

## Preparing synthetic model arrays (Step 03)

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

The initial Fisher convention will keep covariance and survey weights fixed at
the fiducial model and differentiate only the predicted mean. Parameter changes
do not move integration nodes, cuts, weights, or mode counts. Step 03 introduces
no covariance computation, derivative engine, or weight-recomputation mechanism.

## Covariance from supplied observed powers (Step 04)

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

## Fisher information from supplied mean derivatives (Step 05)

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

## External callables and derivatives (Step 06)

Run `python examples/external_forecast.py` for a standalone synthetic galaxy/forest
forecast. It explicitly applies supplied fixed responses and noise, factors the
covariance once, and compares analytic and numerical Fisher matrices and named
marginalized errors at three steps. Its independent P1D example uses velocity
units. Production survey preparation, built-in models and real external-package
validation remain later work.

Use `BoundParameters(registry, local_names, bindings)` from
`fishhighz.models.external` to create an existing `ParameterBinding` together
with its exact registry identity. `bindings` explicitly maps each local name to
a global ID. Multiple names may share one ID; no sharing is inferred from fields.
A binding from a separately constructed registry is rejected, even if its indices
look compatible. Empty local dependencies are allowed.

Create `P3DProvider(label, model, bound_parameters, pairs, jacobian=None,
analytic_ids=None)` and collect providers with
`PreparedP3D(registry, selection, providers)`. Labels must be unique. Each required
covariance pair must have exactly one owner, including unselected spectra.
Pair declarations use field IDs or original field indices, with reversed aliases
canonicalized. Duplicate, overlapping, missing, unknown and unused declared routes
are rejected. Wholly unused fields require no provider. A single provider can own
the entire closure. Providers receive only their owned required pairs, ordered by
increasing original field indices; output is scattered to `required_pairs` order.

The unchanged callable signatures are:

```python
model(theta_local, z, k, mu, pairs)  # -> (node, pair)
jacobian(theta_local, z, k, mu, pairs)  # -> (node, pair, local)
p1d(theta_local, z, k_parallel_velocity)  # -> (node,)
```

`evaluate_p3d(prepared, theta_global, z, k, mu)` returns intrinsic clustering power
in `(Mpc/h_fid)^3`. k is in `h_fid/Mpc`, at fixed h_fid. k and mu are paired,
nonempty 1D nodes (including slices), with k>0 and 0<=mu<=1; z is finite and
nonnegative. There is no sorting, Cartesian expansion or coordinate conversion.
Providers/wrappers own cosmology, AP, physical smoothing, and explicitly modeled
stochastic terms. FishHighz applies none implicitly. Apply supplied fixed response
products once to the mean Jacobian; known noise contributes to fiducial covariance.

Call `evaluate_derivatives(prepared, theta_global, z, k, mu, steps=None,
step_scale=1.0, numerical=False)` from `fishhighz.derivatives`. The returned record
contains owned `power` `(node,required_pair)` and `jacobian`
`(node,required_pair,global)`, plus per-provider/global `columns` diagnostics and
actual `calls`. Gather axis 1 using `selection.selected_to_required` before
`fisher_from_factors`; the global axis already matches the registry.

With a supplied Jacobian, all dependencies are analytic by default. Set
`analytic_ids` to an explicit subset of dependent global IDs for mixed methods,
or to `[]` for fully numerical treatment. The supplied callable must still return
all finite local columns; only selected global columns are consumed. All local
columns tied to an analytic ID sum through `map_jacobian`. A numerical ID perturbs
the global vector before gathering, moving every tied local entry together.
Different providers may use different methods for a shared ID. Missing Jacobians
mean numerical dependencies; unused global columns are exactly zero.

Every numerical ID requires a positive finite absolute `Parameter.step` or an
ID-to-step override. `step_scale` multiplies these explicit steps. Analytic-only
and unused parameters need no steps. Schedules are validated before any provider
calls. Prefer central x±h; near bounds use forward x+h,x+2h, then backward
x-h,x-2h if feasible. Inclusive bounds, distinct finite float64 points and their
ordering are checked. No clipping or step reduction occurs; unresolved steps
fail. Actual offsets define scaled second-order coefficients, accounting for
unequal floating-point spacing. Diagnostics retain requested/actual offsets,
step and stencil family. A provider incurs one fiducial call plus two calls per
unique numerical dependency, and one Jacobian call if any columns are analytic.
Only affected providers are perturbed. No P1D, noise, covariance, grid or weight
preparation occurs in this path.

`check_convergence(..., atol=..., rtol=..., refinements=1)` explicitly compares
numerical columns at h and h/2; higher `refinements` add successive halvings.
It uses the infinity norm over each provider's nodes/pairs and accepts
`change <= atol + rtol*max(norm(coarse), norm(fine))`. Relative change is zero
when both norms are zero. Reports include both stencils, absolute/relative
changes, pass flags and stencil-family changes. Failure does not alter methods;
agreement is not proof of accuracy. Before any provider or Jacobian call, the
entire study rejects consecutive requested steps that round to identical actual
float64 point sets, including at later refinements. Such repeated stencils cannot
assess convergence: use a better-resolved starting step or fewer refinements.
Exact point identities are compared; a refinement moving just one point is valid.
Analytic columns are not convergence tests.
Use `numerical=True` with explicit steps to override analytic strategies for
analytic-versus-FD verification, both here and in `evaluate_derivatives`.

`evaluate_p1d(callable, independent_bound_parameters, independent_theta_global,
z, k_parallel_velocity, label="p1d")` is explicitly separate. It accepts exact
nonnegative velocity nodes in s/km (including zero) and returns km/s power.
No custom P1D means no evaluation, integration from P3D, or implicit default.
There is no assumed common P1D/P3D parameter basis or comoving conversion.

Callable inputs are read-only owned snapshots. Outputs are immediately validated
and copied before another call, permitting providers to reuse mutable buffers.
Signed powers/derivatives are retained; wrong shapes, nonfinite, complex, bool
and object outputs fail. Provider exceptions retain their cause and identify the
label, redshift, parameter/stencil state and requested pairs (where applicable).
Model-specific domain violations fail without retries or substituted values.
Providers must be deterministic for explicit arguments and own hidden-cache
invalidation. Prepared routes cache structure only. Evaluations keep current
buffers and outputs, with no persistent value cache or all-stencil/node history;
pointwise deterministic providers support full-versus-sliced equivalence.

## Linear-power template input (Step 07)

External-model users can install `fishhighz` with NumPy alone. To prepare local
Vega-format FITS templates, install `fishhighz[templates]`; Astropy reads FITS
and SciPy prepares coefficients. Development acceptance uses `.[dev,templates]`
and executes every template test (missing extras are not silently skipped).
Neither extra is imported by package/module import or repeated evaluation.
Missing preparation dependencies give an installation message when requested.

Obtain a decomposed template externally, for example a file produced by Vega's
`bin/make_template.py`. FishHighz reads that format without running Vega, CAMB or
the generator; it bundles no cosmological template assets. Use a local file whose
provenance and scale coverage are appropriate for your application:

```python
import numpy as np
from fishhighz.models.templates import load_template

# Replace this path with your externally prepared Vega-format file.
template = load_template("/path/to/linear_power.fits", h_fid=0.67)
k_min, k_max = template.domain  # converted h_fid/Mpc; not forecast cuts
k = np.geomspace(k_min, k_max, 25)
components = template.evaluate(k)  # z_ref, G=1; (query,2)
smooth, wiggle = components.T  # template.component_names == ("smooth", "wiggle")
dp_dk = template.evaluate(k, derivative=1)
# Illustrative caller-supplied growth ratio 0.9, squared exactly once:
at_other_z = template.evaluate(k, z=template.z_ref + 0.2, growth=0.9**2)
```

The loader selects the unique named `PK` binary table, wherever it occurs in the
file. Scalar real columns `K`, `PK`, `PKSB` must have at least four finite samples,
with strictly positive increasing K. Extra columns/extensions are ignored.
`PKSB` is the supplied smooth component and `PK-PKSB` is the wiggle component;
negative values in **either** are preserved. No positivity constraint, logarithm
of amplitudes, sideband fit, smoothing or new decomposition is applied.

The table header's `ZREF` and `H0/100` specify z_ref and h_template. If absent,
supply `z_ref=...` and/or `h_template=...` explicitly to `load_template`.
Primary-header values do not substitute for missing table metadata. Supplied and
header values must agree within `64*eps64*max(1,abs(header),abs(supplied))`, comparing
h after H0/100; consistent table values take precedence. ZREF must be finite and
nonnegative and H0/h values finite and positive. Original scalar header metadata
(including H0, F_ZREF, SIGMA8 and cosmology) remains immutable provenance; it never
implies a growth function or a forecast parameter. COMMENT/HISTORY cards are not
retained. `source_path` and `source_sha256` identify the exact single byte snapshot
parsed, even if the file changes later. Prepared objects retain no file handles.

Absent TUNIT cards mean the documented Vega convention: k in h_template/Mpc,
power in (Mpc/h_template)^3. The bounded unit-spelling table accepts the following
case-sensitive forms, ignoring spaces and treating `**` like `^`:

| Column | Accepted explicit units |
| --- | --- |
| K | `h/Mpc`, `h Mpc^-1`, `h Mpc-1` |
| PK, PKSB | `(Mpc/h)^3`, `(Mpc/h)3`, `Mpc^3/h^3`, `Mpc3/h3`, `Mpc^3 h^-3`, `Mpc3 h-3` |

Other spellings, including physical `1/Mpc` or `Mpc^3`, are rejected. This is
not a general unit conversion library. Explicit h_fid is always required, and
preparation applies `k_fid=k_file*h_template/h_fid` and
`P_fid=P_file*(h_fid/h_template)^3` once. The same factor scales smooth and wiggle.
Source `k_file`, `pk_file`, `pksb_file` and both h values remain available to audit
this conversion. It changes units for a fixed cosmology; it is not AP scaling.
Unrepresentable conversions or collapsed logarithmic knots fail.

For in-memory input, use
`prepare_template(k, pk, pksb, z_ref=..., h_template=..., h_fid=..., metadata=None)`.
It follows the same source units and numerical validation and only needs SciPy.
Its optional metadata is copied scalar/string provenance; required numerical
inputs are explicit. Both factories return `PowerTemplate`. All source/converted
samples and prepared arrays are owned read-only float64. Component order is
`("smooth", "wiggle")`; converted `full` retains the original full-power samples.

Interpolation is a **not-a-knot cubic spline in ln(k_fid)** for both signed linear
amplitudes. `log_k` stores knots; `coefficients` stores descending cubic powers
of `ln(k)-log_k[interval]`, shaped `(4,interval,component)`. SciPy constructs these
once. NumPy interval lookup and Horner kernels evaluate them; no SciPy object is
retained or invoked. Interpolating both components and summing reproduces full
power interpolation up to floating-point arithmetic. `derivative=1` includes the
1/k chain factor and returns dP/dk in `(Mpc/h_fid)^4`, not dP/dln(k).

Queries are positive finite nonempty 1D arrays. Order/repetitions, read-only and
noncontiguous arrays and consecutive slices are supported; outputs are owned
C-order float64 `(query,2)`. Exact converted endpoints are accepted, with the
upper endpoint in the last interval. Even a one-ULP excursion fails before
logarithms with requested/available domain context. No extrapolation, clipping
or added modes occurs; later transformed queries require caller-provided coverage.

At exactly z_ref, G defaults to 1. At another nonnegative finite z, supply the
positive finite **power factor** `growth=G=[D(z)/D(z_ref)]^2`. Both powers and
k derivatives are multiplied by G once. Explicit G at z_ref must equal 1 within
`64*eps64*max(1,abs(G))`; arbitrary free amplitude changes are not this operation.
Reference metadata/coefficients stay fixed. This layer supplies no growth solver,
Kaiser factors, BAO/AP transformations, damping, survey physics or forecast cuts.
Analytic toy sampling convergence and local format checks do not establish
physical cosmology or survey-forecast accuracy.

## Built-in intrinsic Kaiser/BAO signal (Step 08)

Run `python examples/builtin_forecast.py` after installing `.[dev,templates]`
(or use `fishhighz[templates]` for template preparation). This standalone example
constructs a synthetic oscillatory template and compares BAO-only dilation and
common AP plus galaxy f forecasts. It supplies fixed responses and diagonal
noise explicitly, checks grid/template h_fid equality, and factors covariance
once per mode. Nuisance priors are explicit toy external constraints. No survey
noise preparation or physical forecast-accuracy claim is implied.

`KaiserModel` and `Scaling` live in `fishhighz.models.kaiser`. Construct one model
per redshift bin using a prepared `PowerTemplate` and the original ordered
`ObservedField` sequence. All listed fields are supported by this model; every
one needs a bias and both fixed widths. The signature is:

```python
model = KaiserModel(
    template,
    fields,
    biases={"F": "bf", "g": "bg"},
    betas={"F": "beta"},  # exactly the forest IDs; no galaxy beta
    widths={"F": (4.0, 2.0), "g": (3.0, 1.0)},  # parallel, transverse
    f="rate",  # one shared galaxy f; omit for forest-only models
    local_names=("bf", "bg", "beta", "rate", "w_ap", "w_at"),
    wiggle=Scaling("ap_at", ap="w_ap", at="w_at"),
    # omitted smooth is fixed identity: BAO mode
)
```

Every bias, forest beta, shared galaxy f, or scaling coordinate accepts a fixed
finite number or an explicitly named local slot string. `local_names` defines
exact theta_local order. Each declared slot must occur exactly once in these
settings: unknown, missing, unused and duplicate assignments fail. Use distinct
local slots mapped to one global ID through existing `BoundParameters` to tie
fields or smooth/wiggle blocks; neither field labels nor names imply sharing.
There is no additional registry or expression language. In particular, a width
or G cannot be a free slot. Explicitly choose all finite biases/betas/f or caller
registry bounds; no positivity prior is inferred for them.

Wrap the callable with `P3DProvider(label, model, bound_parameters, owned_pairs)`
and `PreparedP3D` exactly as for an external provider. Pairs use original integer
field indices, with requested order preserved, including reversed pairs and
covariance-required unselected spectra. The callable signature is unchanged:
`model(theta_local, z, k, mu, pairs) -> (node,pair)`. Use the Step 06 numerical
derivative path with explicit steps; a production analytic Jacobian is not added.

Each scaling block chooses exactly one basis and its two coordinates:

| Public basis | Coordinates | a_parallel | a_perp | Vega name |
| --- | --- | --- | --- | --- |
| `ap_at` | `ap`, `at` | ap | at | ap_at |
| `alpha_phi` | `alpha`, `phi` | alpha/sqrt(phi) | alpha*sqrt(phi) | phi_alpha |
| `alpha_iso_epsilon` | `alpha_iso`, `epsilon` | alpha_iso*(1+epsilon)^2 | alpha_iso/(1+epsilon) | aiso_epsilon |

Positive scales/alpha/phi and epsilon>-1 are required, with finite representable
derived factors and Q. `aiso_aap` is unsupported. Smooth and wiggle blocks may
use separate bases, fixed values, partial freedom or independent local slots.
To implement common AP scaling, tie distinct slots in both same-basis blocks to
the same global IDs. No nonlinear cross-basis ties are inferred. Dilations are
template-coordinate parameters; interpreting them as BAO distance/sound-horizon
ratios or geometric AP ratios belongs to the caller, without inferred cosmology.

For each component independently, Fourier mapping uses
`k_parallel=k*mu/ap`, `k_perp=k*sqrt(1-mu²)/at`, their radial magnitude and
`mu_component=k_parallel/k_component`, plus `Q=1/(ap*at²)`. Both the template and
Kaiser angular factor use these component coordinates. Exact identity scaling
preserves the supplied radial coordinates, including template endpoints.

Forest factors are `b*(1+beta*mu_component²)`; galaxy factors are
`b+f*mu_component²`, including b=0 without division. Signed pair products preserve
negative forest crosses. Wiggle damping is
`exp[-(k_parallel,w²*Sigma_parallel,ij²+k_perp,w²*Sigma_perp,ij²)/2]`, with each
pair width squared the mean of its two field widths squared. Autos recover their
widths and cross damping is the geometric mean of the two auto factors. Widths
are explicit nonnegative constants in Mpc/h_fid; zero gives unit damping and a
large finite exponent may underflow to zero. Smooth power is never damped.
Changing f does not change widths or G; varying wiggle coordinates changes D.

The output is `G*(Q_s*B_i,s*B_j,s*P_s + Q_w*B_i,w*B_j,w*D_ij*P_w)` in
`(Mpc/h_fid)^3`, before instrument response or noise. G defaults to 1 only at
z_ref; otherwise prepare with explicit fixed `z=...` and `growth=G`, following
Step 07's squared-growth-ratio convention. Calls at a different z fail. The
model evaluates the reference template and multiplies by G once; retained F_ZREF
and SIGMA8 metadata do not set f or infer evolution. Observed-total PSD validation
remains in covariance preparation; separate scalings/signed components need not
produce a physical total for arbitrary parameters, and are never repaired here.

Provide template coverage for **both mapped components at every fiducial and
finite-difference state**. Domain failures report component, mapped range,
available domain and scales; the derivative engine adds provider/stencil context.
No clipping, extrapolation, dropped nodes, changed steps or automatic regridding
occurs. The model has no second set of forecast cuts or modes. Observed nodes,
weights, volume and k_min/k_max remain fixed, and callers must use the template's
h_fid convention without silent coordinate conversion during derivatives.

Preparation owns fixed settings/widths and integer free-slot maps. Repeated
calls use only NumPy, evaluate each template/factor batch once per component,
reuse results across pairs, and hold no mutable parameter-dependent cache.
Read-only/noncontiguous paired nodes and consecutive slices are supported.
Lazy optional-import checks, independent scalar/product-rule tests and bounded
individual Vega-method comparisons cover the implementation; they do not
validate a complete Vega pipeline, survey model or JIT performance.

## Survey geometry, response, and intrinsic P1D

These composable primitives use supplied backgrounds and synthetic noise; they
are not a survey loader or realistic survey forecast. Run the small three-field
example with `python examples/survey_primitives.py`. The retained external and
built-in examples remain available.

`fishhighz.geometry.prepare_geometry(z_min, z_max, *, z_eval, area_deg2, h_fid,
hubble, transverse_distance, z_order)` prepares one common bin. Supply explicit
`0 <= z_min < z_max`, positive `z_eval` inside the bin, area in square degrees
(up to the full sky), and positive fixed `h_fid`. Ordinary callables receive
immutable 1D redshift arrays and return exactly matching finite positive arrays:
H in km/s/Mpc and transverse comoving D_M in Mpc. D_M is neither angular-diameter
distance nor an already h-scaled distance. Both backgrounds are sampled in
batches at interior quadrature nodes and at `z_eval`, during preparation only.

The immutable `BinGeometry` retains bin/evaluation redshifts, area/solid angle,
`h_fid`, `z_order`, `z_nodes`, dz weights `w_z`, sampled `hubble_nodes` and
`transverse_distance_nodes`, evaluation values, and conversion scalars. It
integrates `volume = Omega*h_fid**3*sum(w_z*c*D_M**2/H)` in `(Mpc/h_fid)^3`.
The explicit Gauss–Legendre order controls volume accuracy; check its convergence
separately from Fourier-grid convergence. No bin-centre approximation, flatness,
background consistency, or interpolation is inferred. `speed_light_kms` records
299792.458. The record contains no cosmology object or background callable.

`prepare_astropy_geometry(cosmology, z_min, z_max, *, z_eval, area_deg2, h_fid,
z_order)` accepts a caller-created Astropy FLRW and delegates to the same
integrator, using unit-converted H and **transverse** comoving distance even for
curved cosmologies. It never substitutes `cosmology.h` for `h_fid` or selects a
default cosmology. Install `pip install 'fishhighz[cosmology]'` for this optional
Astropy/SciPy preparation path. Development acceptance uses
`pip install -e '.[dev,templates,cosmology]'`. Both extras remain optional;
NumPy alone supports supplied backgrounds, response, P1D, and external models.

At `z_eval`, `a_v = H/((1+z_eval)*h_fid)` is in `(km/s)/(Mpc/h_fid)` and
`d_deg = h_fid*D_M*pi/180` is in `(Mpc/h_fid)/degree`. `a_v` is not the scale
factor. Explicit geometry helpers take keyword `a_v`:

- `wavenumber_comoving_to_velocity(k_parallel_comoving, *, a_v)` divides by a_v;
  `wavenumber_velocity_to_comoving` is its inverse.
- `p1d_velocity_to_comoving(power_velocity, *, a_v)` divides km/s power by a_v;
  `p1d_comoving_to_velocity` is its inverse.
- `width_velocity_to_comoving(width_velocity, *, a_v)` divides km/s widths by a_v;
  `width_comoving_to_velocity` is its inverse. Neither changes full-width/sigma
  meaning. These conversion helpers accept finite real arrays or scalars.
- `mode_counts(geometry, grid)` returns positive `volume*grid.q_mode`, shape
  `(node,)`, and rejects unequal geometry/grid h_fid. It does not change the
  grid, cuts, weights, node order, or conjugate-mode normalization.

In `fishhighz.response`, `InstrumentResponse(pixel_width_velocity,
gaussian_sigma_velocity)` requires nonnegative full top-hat pixel width and
Gaussian **one-sigma** in km/s. Use explicit `InstrumentResponse(0, 0)` for
identity, including galaxies. `prepare_response(fields, k, mu, *, a_v, settings)`
requires a mapping from every field ID to its own settings; shared physical
labels imply no sharing. It returns immutable `W(node,field)` in field order.
Inputs are arbitrary paired 1D observed k (h_fid/Mpc) and mu in [0,1], with
slicing and zero k supported. No AP transformation is applied. With q=k*mu/a_v,
`W = np.sinc(q*pixel_width_velocity/(2*pi))*exp(-(q*sigma_velocity)**2/2)`.
Zero q gives exactly one; negative sinc lobes are preserved. Gaussian
attenuation may underflow to zero; unrepresentable arguments fail.
`velocity_response(q, *, pixel_width_velocity, gaussian_sigma_velocity)` returns
the corresponding `(node,)` transfer directly for nonnegative s/km q.

`pair_response(W, selection)` returns `W_i*W_j` in original required-pair order.
Multiply intrinsic P3D by these products, and intrinsic Jacobians by
`products[:, :, None]`, then gather `selection.selected_to_required` for means
and mean derivatives. Supplied noise receives **no automatic response factor**.
The explicit P1D consumer applies W² and the 1/a_v conversion once. Keep geometry,
W, modes, noise, and covariance factors fixed during mean perturbations.

Instrument conversion helpers have deliberately distinct meanings:

| Helper in `fishhighz.response` | Conversion to km/s |
| --- | --- |
| `pixel_width_angstrom_to_velocity(width, *, lambda_obs_angstrom)` | c*full_pixel_width/lambda_obs |
| `gaussian_sigma_angstrom_to_velocity(sigma, *, lambda_obs_angstrom)` | c*sigma_lambda/lambda_obs |
| `gaussian_fwhm_velocity_to_sigma(fwhm_velocity)` | FWHM/(2*sqrt(2*log(2))) |
| `resolving_power_fwhm_to_sigma(resolving_power_fwhm)` | c/(R_FWHM*2*sqrt(2*log(2))) |
| `legacy_resolving_power_to_sigma(resolving_power_legacy)` | c/R_legacy, explicitly interpreted as sigma |

Wavelength inputs are Angstrom and use the local narrow-width approximation.
Supply observed wavelength explicitly; the Lyα convention is
`lambda_obs = 1215.67*(1+z_eval)` Angstrom. FWHM converts once. The labeled legacy
helper uses the new c=299792.458; lyaforecast used 299800, so their values differ
by the known ratio. Direct velocity sigma has no light-speed conversion.

`fishhighz.models.p1d.default_p1d(theta_local, z, k_parallel_velocity)` is a
plain-compatible intrinsic PD2013 callable. Its local vector must be empty;
`BoundParameters(registry, (), {})` connects it to `evaluate_p1d`. The registry
remains nonempty. Direct calls accept finite scalar z>-1; the existing evaluator
retains its nonnegative-redshift contract. k is nonempty finite nonnegative
1D in s/km, including zero. Output is positive `(node,)` km/s power, with no
response, noise, G, BAO, or template normalization. `p1d_floor(z)` exposes the
original redshift-dependent stationary floor. Below it the spectrum is flat;
its first k derivative joins continuously at zero, while the second derivative
need not. Mathematical validity does not establish empirical accuracy outside
the fit's calibration range. Constants and floor are adapted from lyaforecast's
`analytic_p1d_PD2013.py` under its GPLv3 terms; see LICENSE. External P1D uses its
own explicit callable/binding and the same downstream response/conversion.
Changing P3D never selects, integrates, or modifies P1D.

## Fixed forest weighting and sampling noise

`fishhighz.weights.prepare_forest_weights(field, geometry, response, *,
z_source, magnitudes, quadrature, rho, variance, length_velocity, method,
weights=None, iterations=None, signal=None, alias=None, auxiliary=None,
rtol=1e-4, min_updates=3, stable_steps=3, max_updates=96)` prepares
one immutable forest sample per field/bin. `field` is an `ObservedField`,
`geometry` a `BinGeometry`, and `response` the same `InstrumentResponse` used for
signal and final noise. Its full pixel width must be positive; Gaussian sigma
may be zero. Source redshift must exceed `geometry.z_eval`: densities/variances
are sampled at the source, while clustering, P1D and distance conversions use
`z_eval`. Field identity, redshifts, h convention, geometry conversion factors,
response and input arrays are retained. A fixed-count `iterations` value and an
adaptive outcome/counts/residuals are retained, but adaptive input controls must
be preserved separately in caller provenance. Noise evaluation rejects a
different field, evaluation geometry or response. Changing area alone does not
change this local noise preparation.

All magnitude arrays have shape `(n_magnitude,)`, including a one-node sample.
Magnitudes are strictly increasing, quadrature weights are explicitly positive
in magnitudes, and `rho` is already selected/normalized in
`deg^-2 (km/s)^-1 mag^-1`. `variance` is dimensionless delta-flux pixel **variance**,
not RMS. Forest length and response widths are km/s. No file counts, target
density, magnitude spacing, survey area, or source selection are inferred.
`density_per_velocity(dndzdm, *, z_source)` converts a normalized
`dN/(dz_source dm deg^2)` row by `(1+z_source)/299792.458` without normalization.
Legacy rectangular sums require constant `dm` for **every** node, including
both endpoints; their sum is not the node-span. Nonuniform positive quadrature
weights are equally valid.

Choose explicitly:

- `method="supplied", weights=w`: same-shaped nonnegative weights with positive
  weighted support. No auxiliary model calls or legacy settings are used.
  Uniform rescaling leaves final noise unchanged.
- `method="legacy", iterations=3, signal=S, alias=B`: initialize
  `w=B/(B+Delta_v*variance)` and apply exactly three simultaneous updates
  `w_new=S/(S+variance/(cumsum(rho*quadrature*w)*L_v/Delta_v))`.
  `iterations=0` returns initialization; any explicit nonnegative integer is
  supported. Zero-density samples get zero weight; positive-density samples
  with zero variance get one. This is the legacy cumulative recipe, not a
  convergence algorithm or a claim of global optimality. `S` and `B` must be
  strictly positive finite representable scalars.
- `method="early_lyaforecast"`: use the full-sample `sum_historical`
  recurrence. Its update signal is `S + B/(I1*L_v)`, with `I1` recomputed from
  the current weights. This is the recommended research method. With
  `iterations=None`, the adaptive solver requires at least three updates, three
  stable transitions, doubled-count confirmation, positive finite `rtol`
  (default `1e-4`) and a finite cap (default 96). Failure raises and no last
  iterate is substituted. An explicit nonnegative `iterations` value requests
  a labeled fixed-count diagnostic instead. The result retains the convergence
  outcome, counts and residuals; callers must retain nondefault input controls
  in their own provenance.
- `method="mcdonald"`: use the optional full-sample `sum_aliasing` recurrence,
  with update signal `S + B*I2/(I1**2*L_v)`. Adaptive and explicit fixed-count
  behavior is identical to `early_lyaforecast`, including explicit failure.
  This retained alternative is not the selected research baseline.
- `method="inverse_variance", alias=B_star`: compute
  `w=B_star/(B_star+Delta_v*variance)` once on positive-density nodes; store
  zero on zero-density nodes and exactly one for zero variance on support.
  Supply a strictly positive finite scalar `B_star` in km/s, already including
  the P1D response squared at your chosen fixed weighting mode. There is no
  default mode, automatic P1D sampling or 3D signal dependency. `weights`,
  `iterations` (even zero), `signal` and `auxiliary` must all be `None`.
  The returned method and alias record this choice, with no iteration history.
  It equals the legacy zero-iteration seed algebraically. For a nonnegative
  measure, common forest length/pixel width and independent sightline noise,
  it minimizes `Q=A*B_star+P_pixel` under the diagonal-covariance approximation
  at that weighting mode; this does not establish a multi-mode BAO optimum.

For the explicit inverse-variance option, pass these keywords directly in
`ForestInput.weight_options` and omit `auxiliary_coordinates`. Final noise
still uses its own mode-dependent **intrinsic** P1D and applies that mode's
response squared once to aliasing only. `B_star` does not replace this P1D;
pixel noise remains unsmoothed. Weights and coefficients, as well as final
fiducial noise, remain fixed during mean-model differentiation. Legacy
compatibility and all profile defaults are unchanged.

For actual callable integration, `sample_auxiliary(field, geometry, response,
prepared_p3d, theta_p3d, *, p1d_model, p1d_parameters, theta_p1d, k_t_deg,
k_p_velocity)` returns immutable samples usable as `auxiliary=sample` in place
of `signal/alias`. It routes only the matching forest auto through the original
field indices and binding, and calls independently bound P1D once. It computes
`k_parallel=a_v*k_p_velocity`, `k_transverse=k_t_deg/d_deg`,
`S=P3D_auto*W^2*a_v/d_deg^2` (deg² km/s) and `B=P1D*W^2` (km/s).
`S` is the caller's full-auto P3D prediction after the field response and
angular/velocity unit conversion shown above; it is not an unsmoothed intrinsic
scalar. A legacy method name containing “smooth” does not require a smooth-only
external provider. Explicit historical
auxiliary examples are `(2.4, 0.00035)` for BAO and `(7, 0.001)` for P1D, in
`deg^-1` and `s/km`. There is no default. Auxiliary modes may lie outside
forecast cuts, add no modes, and must remain inside provider/template domains.
The retained fiducial states have no provider object or mutable model cache.

The result exposes `weights`, `I1/I2/I3` prefix integrals, `A` in deg² and
`P_pixel` in deg² km/s, plus maximum absolute successive `weight_changes`.
For masses `r=rho*quadrature`, these are prefixes of `r*w`, `r*w²`, and
`r*w²*variance`, with `A=I2[-1]/(I1[-1]²*L_v)` and
`P_pixel=I3[-1]*Delta_v/(I1[-1]²*L_v)`. Zero prefixes are allowed without
optional divisions; unrepresentable arithmetic and empty weighted support fail.
Individual integral-product underflows are rejected before prefix accumulation,
even when another positive term would leave a nonzero total. Exact-zero products
remain valid. Extreme weight normalizations can therefore be rejected while an
equivalent representable normalization succeeds; exposed inputs are never rescaled.

`fishhighz.noise.forest_noise(weights, field, geometry, response, k, mu, p1d)`
accepts paired observed `(n_node,)` arrays and intrinsic nonnegative P1D samples
in km/s at **q=k*mu/a_v**. Obtain these through `evaluate_p1d` with its own
callable/binding/state; replacing P3D never chooses P1D. The immutable result has
`aliasing=A*P1D*W(q)^2*d_deg^2/a_v`,
`pixel=P_pixel*d_deg^2/a_v`, and `total`, all `(Mpc/h_fid)^3`.
Pixel noise is unsmoothed; signed sinc lobes are squared for aliasing.
`galaxy_noise(n_bar)` gives `1/n_bar` for explicit positive comoving number
density `(h_fid/Mpc)^3`. `local_galaxy_density(dndzdm, quadrature, geometry)`
uses `sum(quadrature*dndzdm)*(1+z_eval)/c*a_v/d_deg²`; this is a local-density
approximation, not a count integrated over the bin. Supply the chosen `n_bar`
directly for a count/volume convention.

`prepare_noise(selection, n_node, *, diagonal=None,
independent_sampling=None, full=None)` returns immutable
`(n_node, n_required)` noise in `selection.required_pairs` order. Generated
noise requires `diagonal={active_field_id: node_array, ...}` with exact active
coverage and explicit `independent_sampling=True`. Unused fields need no noise
or model calls. Full packed `full=N` is a **replacement**, not an addition;
combining it with generated settings is rejected. Signed crosses, zero and
singular matrices are supported, with the existing normalized PSD tolerance
checked independently of signal. No overlap is inferred from physical labels.

Prepare each field's weights/noise once and share them across pairs. Add noise
once using `combine_observed_power(W_i*W_j*P_ij, N)` for covariance; the default
mean excludes known noise. Freeze weights, response, P1D/noise, modes and factors
while differentiating intrinsic P3D. A different P1D with fixed supplied weights
changes aliasing only; explicitly regenerating legacy weights can also change
pixel noise. No automatic P1D/noise or covariance derivatives are added.

Run `python examples/weighted_noise.py` for two distinct forests, a galaxy,
negative signal crosses, explicit independent sampling, an independent matrix
covariance/amplitude Fisher oracle and fixed-state checks. Inputs are synthetic;
this is not a physical DESI-2 forecast. Raw readers and overlap-derived noise
remain outside these array APIs. The cumulative formulas retain
lyaforecast GPLv3/McDonald & Eisenstein (2007) provenance in the source.

## Python survey orchestration and raw inputs

`fishhighz.survey.BinSpec(id, geometry, grid, p3d, responses, *, ...)` supplies
explicit existing objects. `p3d` owns the only pair selection and field order.
Use `forests={id: ForestInput(...)}`, `galaxies={id: n_bar}` and
`independent_sampling=True`, or `full_noise=N` as a complete replacement with
none of those generated inputs. Every field requires an `InstrumentResponse`,
including identity settings for galaxies. Only active fields require noise data.

```python
from fishhighz.forecast import prepare_bin, run_bin, run_forecast

prepared = [prepare_bin(spec) for spec in bin_specs]
one = run_bin(prepared[0], batch_size=64, step_scale=1.0)
combined = run_forecast(prepared, batch_size=64, prior_fisher=prior)
refined = run_forecast(prepared, step_scale=0.5, prior_fisher=other_prior)
```

Bins share the **same** `ParameterRegistry` object. Explicit bindings tie shared
physics or give distinct nuisance IDs; unused global columns stay zero. Bin IDs
must be unique, redshift interiors must not overlap, and h_fid and reused field
identities must agree. Unequal widths, gaps and non-midpoint z_eval values are
supported: each geometry integrates its own bounds and area. Raw source-table
cells never define forecast-bin bounds or evaluation redshifts. Caller bin order is preserved. `BinRun.result` is a
zero-prior `FisherResult`; `ForecastRun.combined` adds the supplied prior once.
Singular bin information is valid; no automatic inversion or marginalization
occurs. Inspect `result.diagnostics` or explicitly request marginalized errors.

Preparation freezes intrinsic fiducial powers, response products, known noise,
mode counts and Cholesky factors. Prepared arrays are immutable owned snapshots:
`k`, `mu`, `modes` are `(node,)`; `power`, `noise`, `products`, `total` are
`(node, required_pair)`; `response` is `(node, field)`; `factors` are
`(node, selected_pair, selected_pair)`. Noise is absent from the default mean.
Repeated runs differentiate that mean and reuse every fixed survey quantity.
Numerical step overrides, `step_scale`, and `numerical` follow the existing
`evaluate_derivatives` contract. Errors retain bin, provider/field, global node
slice and stencil context. No cuts or steps are adjusted automatically.

`batch_size=None` uses all nodes; a positive integer bounds transient Jacobians
with consecutive C-order slices. Full factors remain in memory, requiring
O(n_node*n_selected²) storage. This is not disk streaming or a performance claim.
Small accumulation roundoff can differ with batch size. `BinRun.calls` records
actual derivative provider calls aggregated across batches, including repeated
fiducial calls; `columns` records methods/stencils. `PreparedBin.diagnostics`
records preparation calls, units, IDs, volume, response/noise conventions and
source provenance; `weights` retains per-field weighting diagnostics. Model
callables must be deterministic and must not be externally mutated between
preparation and runs. Changing the fiducial/model requires fresh preparation;
callable internals are not deep-copied, hashed, or serialized.

Install `fishhighz[survey]` for the lazy SciPy readers in
`fishhighz.adapters.legacy_inputs`; normalized arrays and external models need
only NumPy. Standalone reader signatures are:

- `DensityReader(path, *, semantics='cell_count_per_deg2', target_density=...,
  z_norm_min=..., magnitude_bounds=None, redshift_widths=None, width_policy=None,
  label='density')`. The semantics argument
  is required. Set target and threshold explicitly to `None` for no normalization
  and whole-grid selection. Inclusive raw magnitude masks precede the strict
  `z > z_norm_min` raw count sum. Divide each row by Delta_z[i]*dm exactly once;
  later magnitude quadrature is separate. Supply `redshift_widths` as a positive
  (n_z,) vector aligned with reconstructed sorted redshifts for physical cell
  measures, including irregular centres. Widths are owned and immutable; no
  physical widths are inferred from centres or forecast bins. Without widths,
  the default (`width_policy=None` or `'uniform'`) requires uniform redshifts.
  Explicit `width_policy='legacy_first_spacing'` uses z[1]-z[0] for every row,
  including irregular tables, solely as the labeled legacy conversion. It is
  never selected automatically. Widths and compatibility policies conflict.
  Provenance retains axes, effective widths, their ordering and policy label. `query(z, magnitudes)` returns ordered differential
  density. `local_galaxy_density(geometry, magnitudes, quadrature)` explicitly
  chooses the local z_eval approximation, without an area multiplier.
- `SNRReader(paths, *, smoothing, label='SNR')` parses header magnitudes and
  metadata, sorts files, and requires identical grids. `smoothing='legacy'` uses
  sigma=10 **sample indices**, reflect boundaries and truncate=4 on wavelength
  only; `'none'` explicitly disables smoothing. Linear interpolation follows.
  `query(*, z_source, magnitudes, wavelength)` returns SNR per Angstrom.
  `variance(..., pixel_width_angstrom, exposure_count, exposure_time=None)` returns
  dimensionless delta-flux variance `1/(SNR²*Delta_lambda*N_exp/N_exp_file)`.
  Explicit exposure time must match file EXPTIME.
- `sample_forest_readers(density, snr, geometry, response, *, z_source,
  magnitudes, pixel_width_angstrom, exposure_count, exposure_time=None)` returns
  immutable `rho`, `variance`, and plain `provenance`. Both readers use z_source;
  observed wavelength uses z_eval. It verifies the pixel velocity width using
  c=299792.458 km/s. Pass these arrays in `ForestInput.weight_options` with explicit
  quadrature, forest velocity length, and the Step 10 weighting settings.
  `ForestInput` also requires an independent P1D callable, binding and state.
  Optional `auxiliary_coordinates=(k_t_deg,k_p_velocity)` selects one auto/P1D
  query for `legacy`, `early_lyaforecast`, or `mcdonald` weights; supplied and
  inverse-variance weights bypass auxiliary calls.

Readers retain resolved paths, SHA-256 hashes, owned tables, domains, units,
normalization/exposure settings and interpolation provenance. They reject
extrapolation, negative spline overshoot, irregular magnitude grids, inconsistent
SNR grids and unrepresentable arithmetic. Uniform-axis checks allow only
64*eps64*max(1,max(abs(axis))) absolute coordinate roundoff. There are no floors,
bright caps, population averaging or automatic resource lookup. Nonuniform raw
redshifts without explicit widths or legacy_first_spacing fail with an actionable
error. Normalization remains a raw count sum before cell-width division; its
measure is distinct from the later magnitude integral.
Live readers never enter the prepared forecast or derivative loops. Adapted
interpolation conventions retain lyaforecast GPLv3/source provenance in the module.

Run `python examples/survey_forecast.py` for generated raw fixtures, two distinct
forests plus a galaxy, two unequal independent bins with explicit non-midpoint
evaluation redshifts, shared target/distinct nuisances,
repeated runs and explicit prior/error diagnostics. Its outputs are synthetic,
not DESI-2 forecasts. INI translation, CLI and serialization remain deferred.

Input snapshots preserve non-object ndarray dtypes and values until scientific
validation. Complex, Boolean, string and object scientific arrays are rejected
consistently with the direct weighting API; no imaginary parts are discarded.
Accepted prepared scientific arrays are immutable float64. Plain metadata arrays
retain their dtype (including indices and Boolean flags). Object arrays cannot be
frozen as immutable byte-backed values and are rejected; use typed arrays or
plain nested metadata instead.

## Historical DESI-2 validation and explicit legacy compatibility (Step 12)

This section records the earlier Step-12 comparison and remains available for
reproduction. Its fixed-reference `inverse_variance` accuracy profile is not the
current research baseline. The current three-profile identity and selected
S2--S4 prescription are given in [RESEARCH_BASELINE.md](RESEARCH_BASELINE.md).

`python examples/desi2_synthetic.py` runs all seven explicit selections on tiny
synthetic arrays, checking their correlated Fisher matrices against an independent
five-field covariance. It requires only NumPy. `fishhighz.validation.cases`
provides `CASE_IDS`, `recipe(case)`, `selection(case)`, `bins(case)` and
`verify_inventory(directory)`. Recipes preserve original tracer indices and
reference upper-triangle pair order, including covariance-required unselected
spectra. Every original option is checked; this is not a general INI translator.

Strict `DensityReader` and `SNRReader` remain the default. Separately construct
`LegacyDensity(reader, negative_policy="reject" | "floor_negative")` and
`LegacySNR(reader)` from validated tables. Their `sample` methods return immutable
`values`, `raw` and `provenance` records. `sample_legacy_forest(...)` converts
these into normalized arrays and plain provenance for `ForestInput`:

- Density strictly outside the raw magnitude domain becomes **1e-20**; closed
  endpoints remain inside. Redshift uses the reference spline domain extension,
  with original and effective coordinates recorded. This is not validated
  physical extrapolation. Exact zeros and positive densities below the floor
  remain unchanged. Negative in-domain interpolants fail with `reject`;
  `floor_negative` replaces only negatives by **1e-20**, an explicit extension
  beyond legacy behavior. Compatibility preserves the reference normalization
  reduction order, including zero-masked forest rows.
- SNR clamps bright magnitudes to the first node. Faint magnitudes or source
  redshift/wavelength outside closed domains return variance **1e20**, independently
  of pixel/exposure scaling. Inside, the **1e-10 SNR floor** applies after
  `sqrt(pixel_width_angstrom)*sqrt(exposure_count/file_exposure_count)` scaling.
  Exposure time must match the file. Arithmetic failures and invalid inputs fail.
- Diagnostics include source hashes, normalization/units, original/effective
  queries, raw values, reason masks/counts and constants. Out-of-range SNR raw
  placeholders are zero and explicitly marked unevaluated. Diagnostics are
  snapshots, with no counters updated during derivatives. Each forest owns its
  own source, SNR, response, weights and independent P1D.

Revision 2 historically provided two explicit validation profiles for the seven
original cases (39 bins, 78 primary records). They share observed cuts and
selected spectra, but represent different estimators:

- **Maximum compatibility** captures actual upstream mean/total powers and BAO
  derivatives, then independently constructs Wick covariance and Fisher
  information in FishHighz. It preserves legacy endpoint sums, polynomial peak
  extraction, backward derivative, pair-specific noise and redshift conventions.
  Literal field-matrix PSD failures are reported separately; selected covariance
  must still be positive definite. This path is validation-only.
- **Historical maximum accuracy** uses the accepted `prepare_bin`/`run_bin` and `KaiserModel`
  engines: geometric evaluation redshift, integrated volume, Gauss–Legendre k/mu
  and composite magnitude quadrature, CAMB growth at actual redshifts, physical
  FWHM resolution, accepted auto/cross damping and full wiggle-mapping derivatives.
  Each forest now uses fixed inverse-variance weights at the declared
  `q_star=0.00035 s/km`: `B_star=P1D(z_eval,q_star)*W_field(q_star)^2` and
  `nu=B_star/(B_star+l_p*v)`. The intrinsic P1D, field response and `B_star` are
  recorded per forest; no P3D signal or auxiliary coordinates enter this choice.
  Noise and weights are fixed during differentiation. An explicitly selected
  `accuracy_method="legacy"` retains cumulative compatibility/replay behavior.
  “Maximum” is bounded by this existing model and the prescribed convergence
  tests; the fixed mode is a convention, not a fitted or optimized scale.

The default is the 15x2pt bin [2.47,2.705]. Full execution must be explicitly
selected and authorized. The interpreter must already provide the requested
reference checkout, CAMB, and FishHighz optional extras. The scripts never modify
reference code or install dependencies automatically. Use fresh output paths:

```bash
python scripts/desi2_validation.py preflight --reference /path/to/lyaforecast
# Separately authorized fresh reference capture:
python scripts/capture_desi2_reference.py --suite full \
  --reference /path/to/lyaforecast --output /path/to/new-reference
# Both profiles; omit --suite full for the one-bin quick default:
python scripts/desi2_profiles.py run --suite full \
  --reference /path/to/lyaforecast --template /path/to/template.fits \
  --reference-bundle /path/to/new-reference --wheel /path/to/exact.whl \
  --output /path/to/new-profiles
python scripts/desi2_profiles.py check --output /path/to/new-profiles
# Strict full-assignment inventory and scientific gate:
python scripts/check_desi2_full.py /path/to/new-profiles
# Offline tables/plots, including failed convergence (Matplotlib optional):
python scripts/plot_desi2.py --bundle /path/to/new-profiles \
  --output /path/to/new-plots
# Sensitivity panels with explicit gaps and labels for failed variants:
python scripts/plot_desi2_sensitivities.py --plots /path/to/new-plots \
  --output /path/to/new-sensitivity-plots
python scripts/audit_desi2.py --bundle /path/to/new-profiles \
  --output /path/to/new-audit.json
python scripts/combine_desi2_bins.py --bundle /path/to/new-profiles \
  --output /path/to/new-independent-bin-check
python scripts/attribute_desi2.py --bundle /path/to/new-profiles \
  --template /path/to/template.fits --output /path/to/new-attribution
# Verify saved attribution matrices and redraw without model evaluation:
python scripts/attribute_desi2.py --replot /path/to/new-attribution \
  --output /path/to/redrawn-attribution
python scripts/check_desi2_plots.py --plots /path/to/new-plots \
  --reference /path/to/new-reference
```

`--profiles`, `--cases` and `--bins` select explicit subsets; subsets never count
as the full 78-record comparison. The Python entry point is
`fishhighz.validation.profiles.run(...)`. `wheel_identity` verifies every imported
module against the exact supplied artifact; use an installed copy of that wheel
or explicitly place that wheel on the interpreter's import path.

Schema 2 binds case/profile/bin/field/pair/parameter identities and reconstructs
nodes, Wick covariance, Fisher information and rank-aware errors offline. It
recomputes recorded convergence metrics and requires boolean scientific success.
Execution finished and scientific acceptance are distinct. Real reports live in
individually hashed JSON files, keeping the manifest bounded as diagnostics grow.
After an interrupted validation run, `--reuse-completed /path/to/bundle` can reuse
immutable numerical inputs only when every scientific/schema module and the
actual reference interpreter, source and resources match. Report-I/O/controller
changes are the only permitted code differences. C/F are reassembled with the
current exact wheel; both the original producer and reuse provenance are retained. Unconstrained errors
use explicit availability masks; plots/tables use nulls. Historical schema-1
bundles can only be inspected with `evidence.inspect_legacy`, which returns a
limited result and cannot pass this gate. No pickle or external model imports
are needed to read evidence.

Fixed-reference accuracy studies use trial-contract version 2 and vary k, mu,
volume, magnitude and derivative step independently, then test their combined
lower controls. Iterations and a `weights` convergence metric are explicitly
inapplicable. Historical version-1 evidence retains the six-control cumulative
schedule; missing method metadata means legacy, never fixed-reference. Cache
reuse checks the requested method, and `three_weights` remains a legacy-only
diagnostic. The existing cumulative rule may underflow at higher counts; such
historical failures and unresolved single-pair convergence remain findings,
without altered tolerances or a new optimizer. Floors/clamps and
`legacy_first_spacing` remain labeled input policies.
Floor changes, support removal and artificial ±10% cell widths are diagnostics,
not calibrated systematic uncertainties. Additional one-setting physical swaps
and common-node array chains record order dependence and unresolved grouped
causes. Significant unexplained differences prevent an unqualified pass.

The historical r1 one-bin study used 107 rectangular magnitude nodes, three weight
updates, EdS power growth and legacy c/R resolution. It is preserved separately
and is not the revision-2 accuracy result. All six examples and the explicit
legacy adapters remain available.

The callable `IntrinsicP3D(provider, *, routes, n_fields, h_source, h_fid,
k_domain, z_domain)` wraps a caller-prepared object's **intrinsic**
`compute_p3d_hmpc(z,k,mu,corr)`. It validates all paired queries and explicit routes
before calling a provider that might clamp. It preserves signed crosses and
arbitrary pair order, maps `k_source=k_fid*h_fid/h_source` and
`P_fid=P_source*(h_fid/h_source)**3`, and adds no response/noise. It has zero
local parameters; the real validation uses a separate `A*P_external` wrapper
and independent matrix amplitude-Fisher oracle. That amplitude example is not
the primary BAO forecast; P1D remains independently chosen.

Run the bounded actual external-model example in the identified reference
environment with the exact FishHighz wheel imported:

```bash
python scripts/external_desi2.py --reference /path/to/lyaforecast \
  --template /path/to/template.fits --wheel /path/to/exact.whl \
  --output /path/to/new-external-amplitude
```

This caller constructs the actual reference `PowerSpectrum` once and saves a
separately typed, independently reconstructed amplitude record. It does not run
`NewForecast` or replace the full BAO comparison.

General INI translation and production CLI/serialization remain deferred. Step 13
now provides a validation-only amplitude/shape diagnosis of the unchanged
cumulative weighting recurrence; it is not a production weighting mode. See
`reviews/step-13.md` for its bounded scientific conclusion and
`reviews/step-12.md` for the preceding forecast evidence.

### Optional Fisher contraction compilation

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
[the performance report](reviews/step-12-performance-r4.md). These synthetic
measurements do not establish real-forecast timing or forest convergence.

### Step 12 revision 5: six-bin 15x2pt diagnosis

Revision 5 retains the revision-4 NumPy default and optional compiled contraction.
New evidence uses schema 3: every accuracy metric references distinct actual
refinement trials, including their controls, combined and individual-spectrum
Fisher matrices, volumes and failed outcomes. The offline checker replays the
bounded control sequence from saved summaries without evaluating a model.
Schema-2 bundles remain historical, with explicitly limited convergence checks;
they cannot satisfy the new 15x2pt acceptance gate.

Run only the explicitly selected case for this assignment, using the exact newly
built wheel in an isolated environment and the preserved upstream reference:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export FISHHIGHZ_FISHER_BACKEND=numba # optional fishhighz[compiled]
python scripts/desi2_profiles.py run --suite full \
  --cases lya_qso_lbg_lae_15x2pt --bins 0 1 2 3 4 5 \
  --reference ../lyaforecast \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --reference-bundle .validation/step12-r2-implementation/reference-01 \
  --wheel /absolute/path/to/exact.whl --output /new/exclusive/bundle
python scripts/check_desi2_15x2pt.py /new/exclusive/bundle
python scripts/plot_desi2.py --bundle /new/exclusive/bundle --output /new/plots
python scripts/check_desi2_plots.py --plots /new/plots \
  --reference .validation/step12-r2-implementation/reference-01 \
  --cases lya_qso_lbg_lae_15x2pt
```

The scoped gate requires 12 primary records, 180 individually selected-spectrum
results and 72 diagnostic requests. The historical seven-case gate is separate.
Plots retain failed convergence and unavailable diagnostics. Compatibility
reconstructs the legacy estimator from saved upstream powers and Jacobians;
accuracy uses the existing physical-response, growth and full-wiggle-derivative
prescription. Neither profile introduces sample-overlap noise or calibrated
coverage beyond the retained input policies. Cumulative weighting stabilization
remains a scientific limitation; a changed prescription requires user review.

The cumulative-update diagnosis can be reproduced from the saved six-bin
samples without a CAMB solve:

```bash
python scripts/recheck_desi2_15x2pt.py --source /historical/schema2/bundle \
  --output /new/derived/evidence
python scripts/diagnose_desi2_weights.py --bundle /six-bin/bundle \
  --reference ../lyaforecast \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --output /new/weight-diagnosis
python scripts/plot_desi2_weight_diagnosis.py --source /new/weight-diagnosis \
  --output /new/weight-figures
```

This diagnostic holds the initial weight function fixed while refining magnitude
integration, then compares the retained finite cumulative updates. It does not
replace the primary weights. Decimal products in the diagnosis are evaluated to
80 significant digits to identify unrepresentable operations; they are not used
in forecast calculations. Historical failure controls reconstructed by exact
controller replay are explicitly distinguished from newly recorded attempts.

For the scoped independent-bin check and controlled comparisons:

```bash
python scripts/combine_desi2_bins.py --bundle /six-bin/bundle \
  --output /new/independent-bins --cases lya_qso_lbg_lae_15x2pt
python scripts/attribute_desi2.py --bundle /six-bin/bundle \
  --template ../vega/vega/models/Planck18/Planck18_z_2.406.fits \
  --output /new/attribution
python scripts/plot_desi2_sensitivities.py --plots /new/plots \
  --output /new/sensitivity-figures
```

After resuming a run that reassembles primary Fisher matrices, bind saved
sensitivity comparisons to those final matrices before plotting:

```bash
python scripts/bind_desi2_sensitivities.py --source /resumed/bundle \
  --output /new/checked-comparisons
```

This derives a new bundle, preserving numerical arrays and prior reports. It
recomputes comparison metrics rather than relaxing checks for nearly zero effects.

### Step 12 revision 6: weak-spectrum evidence checks

Revision 6 compares each individual-spectrum Fisher matrix and derived summary
on its own scale, including primary/final-trial and replayed operands. Scaling
before subtraction and norm evaluation avoids range loss for very weak or strong
information. The checker reconstructs convergence metrics and the verdict from
saved study arrays; detached operands cannot certify convergence. Evidence
consistency remains 5e-12, separately from the unchanged scientific thresholds.

This repair uses synthetic regressions and offline inspection only. The preceding
real-run commands document revision 5 and require explicit authorization for any
new execution. Saved revision-5 compatibility passes in six bins; all six accuracy
records remain unconverged, with 12 of 72 diagnostics unavailable. Scientific
weighting decisions and Step 12 acceptance remain open. See the
[current handoff](reviews/step-12.md) and [preserved revision-5 handoff](reviews/step-12-r5.md).

### Step 13 revision 2: cumulative-weight limit diagnosis

`fishhighz.validation.weight_limit` evolves the exact nonlinear recurrence in a
range-safe log-amplitude/log-shape representation. It records ordinary-value
availability, the fixed-grid linearized spectrum, the normalized sightline
distribution, tail measure and coefficient concentration without changing
`prepare_forest_weights` or its underflow guards. The explicit controller
`scripts/diagnose_weight_limit.py` reads only preserved Step 12 arrays and binds
every population, quadrature order and attempted checkpoint to those immutable
inputs. `--check-finalized` independently reconstructs all final tables and
figures before accepting a saved bundle.

The five available magnitude orders show strong, consistently positive growth of
the fixed-grid asymptotic pixel-noise coefficient, but they are finite and
non-nested. They do not establish the required non-atomic continuum measure,
uniform nonlinear remainder control or an exchange of iteration and refinement
limits. Revision 2 therefore classifies all twelve saved forest populations as
unresolved in the continuum limit. No forecast was reassembled and no production
prescription was adopted; Step 12 scientific acceptance remains open.
