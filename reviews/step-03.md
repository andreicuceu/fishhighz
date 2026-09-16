# Step 03 handoff — revision 1

Implemented against `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-11).
Ready for user review; acceptance and progression are not claimed. Step 04 has
not started. Base commit: `0d69786a06d5d676564a51fad14a7156951c2228`.

## Delivered contracts

| Module | Public names and behavior |
| --- | --- |
| `fields.py` | `ObservedField(id, kind, physical_model, background=None)` and `PairSelection(fields, selected=None)`. Ordered opaque identities; forest backgrounds need no observed galaxy counterpart. Sharing a physical label never ties parameters or merges samples. |
| `parameters.py` | `Parameter(id, fiducial, role, bounds=None, step=None)`, `ParameterRegistry(parameters)`, `ParameterBinding(registry, local_names, bindings)`, `gather_local(theta_global, local_to_global)`, `map_jacobian(local_jacobian, local_to_global, n_global)`. Explicit equality bindings; repeated local columns sum into the shared global column. Fixed inputs are outside the registry. |
| `grids.py` | `IntegrationGrid(k, w_k, mu, w_mu, *, k_min, k_max, h_fid)` and `gauss_legendre_grid(k_edges, *, k_order, mu_order, h_fid)`. Fixed tensor coordinates; custom positive normalized rules and configurable per-bin Gauss–Legendre quadrature. |
| `models/protocols.py` | Structural callables `P3D`, `P1D`, `P3DJacobian`; separate `validate_p3d(value, n_node, n_pair)`, `validate_p1d(value, n_node)`, `validate_p3d_jacobian(value, n_node, n_pair, n_local)`. No provider inheritance or autodiff requirement. |

`_arrays.py` is a small private boundary-helper module for copying, dtype,
finite-value, scalar, label, and index validation. No routing framework or
scientific model was introduced. Package-level import remains unchanged and
loads neither NumPy nor neighboring packages.

Pair arrays use `(n_pair,2)` canonical integer indices, i <= j. Available and
required pairs sort by i then j; selection order is retained.
`selected_to_required` gathers means in selection order. `im`, `jn`, `in_`,
`jm` each have shape `(n_selected,n_selected)` and index the required-pair array
for selected A=(i,j), B=(m,n). Missing required spectra are errors for later
provider integration, never implicitly zero or reconstructed. All stored indices
are owned, read-only, C-contiguous int64 arrays.

Registry `fiducials` and every stored grid array are owned, read-only,
C-contiguous float64. Bindings retain immutable local names and an owned
`local_to_global` vector. Optional finite bounds are inclusive, may use None for
an open endpoint, and require lower < upper when both endpoints exist. Negative
and zero fiducials are supported; positive steps need not fit within bounds.

Grid flattening is C-order `(n_k,n_mu)`, mu fastest:
`node = i_k*n_mu + i_mu`. Exposed arrays include `k`, `w_k`, `mu`, `w_mu`,
`k_flat`, `mu_flat`, `w_k_flat`, `w_mu_flat`, `weights`, and `q_mode`.
`weights` integrates dk*dmu; `q_mode = k_flat**2*weights/(2*pi**2)` gives
`N_modes = V_fid*q_mode`. It counts both conjugate hemispheres on mu in [0,1]
without an extra half factor. k uses h_fid/Mpc; P3D and volume use
(Mpc/h_fid)^3. Positive `h_fid` metadata does not rescale coordinates. Custom
weight sums use rtol=1e-12, atol=0; rules are never sorted, clipped, or normalized.
Numerically unrepresentable interior k nodes are rejected by the default factory.

Providers evaluate one scalar z. P3D returns `(n_node,n_pair)`; supplied local
Jacobians return `(n_node,n_pair,n_local)` and map to
`(n_node,n_pair,n_global)`. P1D takes velocity wavenumbers in s/km and returns
`(n_node,)` power in km/s independently of P3D. Output validators copy to
C-contiguous float64, preserving signed powers/derivatives and rejecting wrong
shapes, broadcasting, complex/bool/object/string values, and nonfinite values.
Providers may evaluate arbitrary paired points, independently of tensor grids.
Adapters own AP/cosmological physics, unit conversion, and explicit handling of
preexisting response/noise. P3D's contract precedes response and known sampling
noise. Future Fisher evaluation holds covariance and survey weights fiducial;
this step implements no covariance or weight recomputation.

## Acceptance evidence

All commands below ran from the package root unless explicitly stated. Local
checks used Python 3.13.15, NumPy 2.5.3, pytest 9.1.1, Ruff 0.16.7, build 1.6.1,
and pip 26.2.1. NumPy is the sole runtime dependency, without an unnecessary
version pin; the APIs used are longstanding and the resolver respects supported
Python versions. Validation here exercised Python 3.13, not a Python-version
matrix. Dependency records: `.validation/step03/dev-versions.json` and
`.validation/step03/wheel-probe.json`.

- `.venv/bin/python -m pip install -e '.[dev]'`: passed after an approved
  network-enabled retry; `.validation/step03-install-retry.log`.
- `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh`: final result **146 passed
  in 23.52 s**, Ruff lint passed, 27 files already formatted. Log:
  `.validation/step03/checks-final.log`. The runner sets OMP_NUM_THREADS,
  OPENBLAS_NUM_THREADS, and MKL_NUM_THREADS to 1.
- `git diff --check`: passed.
- `.venv/bin/python -m build --wheel --outdir .validation/step03/final-dist`:
  passed with isolated build dependencies; `.validation/step03/build-final.log`.
- `.venv/bin/python -m venv .validation/step03/wheel-env`: created a new isolated
  environment for this assignment, before any FishHighz installation there.
- `.validation/step03/wheel-env/bin/python -m pip install .validation/step03/final-dist/fishhighz-0.1.0.dev0-py3-none-any.whl`:
  passed, installing FishHighz and declared NumPy 2.5.3 together;
  `.validation/step03/wheel-install.log`. No same-version install was reused.

The isolated wheel probe ran from `/tmp`, with all three numerical thread limits
set to 1. Its exact invocation used the following absolute-path structure,
where ROOT was `/global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ROOT/.validation/step03/wheel-env/bin/python" -I \
  "$ROOT/.validation/step03/probe.py" \
  "$ROOT/.validation/step03/final-dist/fishhighz-0.1.0.dev0-py3-none-any.whl" \
  "$ROOT/.validation/step03/readme-example.py"
```

Output was saved to `.validation/step03/wheel-probe.json`. All seven installed
Python module files matched the exact wheel bytes. The probe checked quiet,
minimal import before importing NumPy; executed the README integration example;
and exercised custom-grid and P1D behavior. It resolved imports to the fresh
environment's site-packages, not the checkout. The home-path workspace resolves
to the CFS path shown in the probe; all writes stayed within the authorized
FishHighz package.

Final wheel SHA-256:
`373bce8b93bb7af97213348ca04c2bac14aaae28bada9377dd36e1c6d08427f7`.

Coverage includes opaque identities, 15 pairs from five fields, cross-only and
auto-only covariance dependencies, unused fields, explicitly named products for
all four lookup tables with permuted selection, invalid pairs/metadata, distinct
forest samples, local ordering, empty dependencies, explicit global/bin/field/pair
scope IDs, repeated derivative contributions (2+3=5) checked analytically and by
a tiny test-only finite difference, and input/stored-array ownership. Quadrature
tests check unequal bins, analytic k and mu moments, custom rules, flattening,
known mode/volume normalization, fixed h conventions, and malformed inputs.
Synthetic providers exercise unequal output dimensions, signed cross powers,
shared Jacobians, unchanged grids under theta changes, and arbitrary model points.
The retained 52 scaffold/reference-tool tests are included; none run real forecasts.

## Resolved check failures and limits

The first dependency installation and isolated build attempts failed because the
sandbox blocked PyPI DNS. Approved retries succeeded. An initial no-isolation
build also found no setuptools backend in the development environment; the final
build used the declared isolated build requirements. Logs are retained under
`.validation/step03/` (and the two install logs in `.validation/`). An intermediate
check passed 145 tests but flagged README code-block formatting; formatting was
fixed, and the added narrow-bin boundary test brings the final total to 146.
An intermediate wheel in `.validation/step03/dist` is superseded by `final-dist`;
it was never used as final installation evidence.

No acceptance check remains failing. No real lyaforecast quick capture or full
suite was run: neither is required here, and full execution was not requested.
Accepted baseline evidence and reference tooling were preserved. No covariance,
Fisher engine, finite-difference implementation, survey integration, or scientific
forecast convergence is claimed. Polynomial quadrature tests establish only the
specified integration identities. Array write protection prevents accidental
mutation, not deliberate re-enabling of NumPy write flags.

## Git state and source identity

`AGENTS.md` and `IMPLEMENTATION_STEP.md` were already modified when work began;
they were read and preserved. This assignment changes README.md and
pyproject.toml, adds six implementation files (including the model initializer),
four focused test files, and this report. Existing source, tests, reference
scripts, plans, and historical reports were not edited. No commit, push, or
roadmap advancement was performed.

The source hashes below identify the final implementation and tests. The separate
`.validation/step03/source-sha256.txt` also records this report, unchanged package
initializer, and current governance inputs; `.validation/step03/git-status.txt`
and `tracked.diff` retain the review state. The report itself is omitted from
its inline hash list to avoid a self-referential hash.

```text
3e7e9aececcc0bdd54023d6746dc4b03d0223152cdbe42659d64b6157b8c16a5  README.md
0659b35674927445b95d450d39b87664f246ca8282c1557470b3bf5c55f2bafb  pyproject.toml
113df72b69cac65b64e3c8bd73765398270bce5a9856a07d8038ba7b362ad0a6  fishhighz/_arrays.py
27c3d17b587b3f6a126723fd015b5bb8db1f59420887314d850163d0b02fb1bd  fishhighz/fields.py
824ff40892a2414a07cb0959d32a28a97889eccc61d65463b2b858f1a88b3345  fishhighz/parameters.py
8fcbd6f012d9a0b0cf5c25912c66c54b51f93581d13b1540b1d3d21b2960dd6e  fishhighz/grids.py
f1801f533e9f06d84a671f6296c76c2b70bf5b9033e9cf083c514d9f391378d0  fishhighz/models/__init__.py
3b145b2b6b3f59e77d2fa7b0a07d5ab487e6acef475b3f90809deb707e1053ab  fishhighz/models/protocols.py
d6de4fb44b3e72956b35e36c962807052fa4a4b2678fb3ddc4313aad0df11e18  tests/test_fields.py
86af59bfaf3680e0f656b18d41c8061b27d30904a9271ee50522ec2f558642bd  tests/test_parameters.py
ffe2597dbb9e50e49dac859d8eb7dc7be5821d198a4dd285585df5815af3304c  tests/test_grids.py
d9ed28e2fe08f1c3246bda838fd8e5fc4c8749083890a310fa87ecc0c18f1442  tests/test_models.py
```
