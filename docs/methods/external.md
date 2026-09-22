# External models and mean derivatives

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

