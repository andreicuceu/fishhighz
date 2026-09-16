# Step 05 handoff — revision 1

Implemented `IMPLEMENTATION_STEP.md`, revision 1 (2026-09-11). Ready for user
review; acceptance and progression are not claimed. Step 06 has not started.
Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228`. Accepted Steps
03–04 and the current planning documents were already uncommitted at dispatch
and have been preserved.

## Interfaces and normalization

Added `fishhighz/fisher.py`, `fishhighz/kernels/fisher.py`,
`fishhighz/results.py`, and private `fishhighz/_information.py`, plus
`tests/test_fisher.py` and `tests/test_results.py`. README.md now documents the
supplied-array Fisher calculation. `_information.py` centralizes the shared
matrix-boundary checks for this step; it does not change Step 04's validation.
No existing Python source/test, dependency, reference tool, or governance file
was edited. NumPy remains the sole runtime scientific dependency.

| Interface | Contract |
| --- | --- |
| `factor_covariance(covariance)` | Returns owned C-contiguous float64 lower factors `(node,selected,selected)`. Normalized covariance validation rejects numerically unresolved blocks with node/rank context. |
| `fisher_from_factors(jacobian, factors)` | Accepts `(node,selected,global)` supplied mean derivatives and reusable factors; returns owned float64 `(global,global)` data information. |
| `fisher_matrix(jacobian, covariance)` | Convenience composition of the preceding operations. |
| `FisherResult(registry, data_fisher, prior_fisher=None)` | Stores exact registry metadata and separate owned, read-only float64 data, prior, and total information. `prior_fisher` is keyword-only. |
| `diagonal_prior(registry, sigmas)` | Constructs precision from an explicit ID-to-sigma mapping, leaving other IDs unconstrained. |
| `combine_results(results, prior_fisher=None)` | Sums data-only contributions in identical ordered metadata, then adds the one optional keyword-only prior. |

The numerical kernels are `_cholesky(matrix)`, returning one lower factor;
`_forward_substitute(lower, rhs, out)`, overwriting one node's solution across
all RHS columns and returning None; and `_accumulate_fisher(solved, out)`, adding
`solved.T @ solved` and returning None. They accept numeric arrays only and
contain no metadata, boundary validation, model calls, logging, or filesystem
work. Forward substitution uses no general LU or explicit covariance inverse.

At each node, `C=L@L.T`, `Y=solve(L,J)`, and data information accumulates
`Y.T@Y`. Step 04 covariance already contains mode counts, volume, and quadrature.
No extra normalization or covariance-derivative information is added. Supplied
Jacobians must follow selected mean-spectrum order and the global registry basis;
known subtracted noise belongs in fiducial covariance, while a noise parameter
requires an explicit supplied mean response. There is no theta argument or
covariance/weight recomputation. Zero derivative columns and zero data Fisher
matrices are legal.

Factors are explicit reusable artifacts, not a cache. Reusing the from-factors
path performs no Cholesky or covariance eigensolve. Direct supplied factors are
validated for finiteness, exact lower-triangular structure, and positive diagonals;
their scientific provenance and conditioning are the caller's responsibility.
Use `factor_covariance` for the normalized solvability check. Its decision is
based on selected covariance, not the rank of the parent field-power matrix.

Existing array inputs retain views at the outer boundary. Float64 conversion,
finiteness checks, and solve buffers are limited to one node. Python array-like
inputs necessarily materialize an initial array. Working numerical memory is
one node's factor/derivative/solution and contraction buffers, plus stored factors
and output information; there are no all-node Fisher intermediates. Kernel
inspection and whole-versus-sliced tests support this contract. Inputs, including
read-only/noncontiguous views, remain unchanged.

## Numerical thresholds and result behavior

For a square matrix M of size n, require nonnegative diagonals exactly. An exact
zero diagonal requires an exact zero row and column. Set d=sqrt(diag(M)), using
1 in the declared units for zero information rows, and form R=M/d/d by sequential
axis division. This removes disparate observable/parameter amplitudes from rank
decisions. Covariance factorization additionally requires positive variances.

The exact symmetry tolerance is elementwise
`64*eps64*n*max(1,abs(Rij),abs(Rji))`. Substantive asymmetry fails; accepted
roundoff asymmetry is averaged in local/owned copies only. Exactly symmetric
physical entries are left untouched, including subnormal values. With
`eps64 = np.finfo(np.float64).eps`, the exact spectral threshold is
`t = 64*eps64*n*max(1,max(abs(eigenvalues(R))))`:

- PSD information permits eigenvalues >= -t; no eigenvalue is modified.
- Numerical rank counts eigenvalues > t.
- Covariance factorization and marginalized uncertainties require full numerical
  rank; no jitter, clipping, projection, node removal, or pseudoinverse is used.
- The normalized condition diagnostic is max/min eigenvalue for full rank,
  otherwise infinity. Null/near-null directions include eigenvalues <= t.

Covariance failures identify the node and numerical issue; unresolved rank reports
include the normalized minimum eigenvalue and threshold. Solve failures include
node/global-column context; arithmetic overflow fails explicitly rather than
altering physical units or returning nonfinite information. Result-level errors
include the registry IDs or rank/null context. Finite prior widths producing
unrepresentable zero/infinite precision are rejected.

`result.diagnostics` stores ordered `ids`, `rank`, `condition`, `tolerance`,
`scales`, `eigenvalues`, and `null_directions`. Directions are unit Euclidean
row vectors in coordinates `x = scales * delta_theta`, with scales d defined
above. Physical displacements are `null_directions / scales[None,:]`; the stored
eigenvectors are never labeled physical-coordinate combinations. Signs and
bases within a degenerate eigenspace are nonunique. Diagnostic arrays are owned
and read-only. Diagnostics assess total data-plus-prior information; singular
data or prior information alone is allowed.

`conditional_errors(ids=None)` returns ordered `1/sqrt(F_total[ii])`, fixing all
other parameters, with infinity for exactly zero diagonals.
`marginalized_covariance(ids=None)` factors the full retained total information,
solves against identity, and forms the covariance before gathering a requested
subset. `marginalized_errors` and `correlations` derive from that covariance.
Other parameters remain free during these requests. A singular full result
raises even when a requested subset might be partially identifiable.
`fix_except(ids)` instead explicitly fixes the complement, taking principal
data/prior submatrices and preserving corresponding parameter metadata. No role
implicitly removes nuisance parameters. Unknown/duplicate/empty subsets fail.
Construction, combination, and diagnostics perform no eager uncertainty solve.

Combination requires identical ordered parameter records, including fiducials,
roles, bounds, and steps. It rejects any nonzero prior already attached to an
input. A shared finite diagonal or correlated PSD precision is applied exactly
once, after data summation; bounds and steps never create information. Individual
bin contributions may be singular. Combining before marginalizing is essential
for shared nuisances, as the complementary theta+eta/theta-eta tests demonstrate.

## Acceptance evidence

All commands below ran from the package root unless noted. Development versions:
Python 3.13.15, NumPy 2.5.3, pytest 9.1.1, Ruff 0.16.7, build 1.6.1, pip 26.2.1.
Full records are in `.validation/step05/dev-versions.json`. Numerical checks used
OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, and MKL_NUM_THREADS=1.

| Command | Actual result and evidence under `.validation/step05/` |
| --- | --- |
| `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest tests/test_fisher.py tests/test_results.py` | Initial 72 tests passed in 3.26 s; `focused.log`. |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` | Final **261 tests passed in 25.47 s**, including 75 new tests and all 186 accepted tests; Ruff lint passed and 41 files already formatted. `checks-final.log`. |
| `git diff --check` | Passed. |
| `.venv/bin/python -m build --wheel --outdir .validation/step05/dist` | One exact wheel built in a fresh output directory with isolated dependencies; `build.log`. |
| `.venv/bin/python -m venv .validation/step05/wheel-env` | Created a new isolated environment. |
| `.validation/step05/wheel-env/bin/python -m pip install .validation/step05/dist/fishhighz-0.1.0.dev0-py3-none-any.whl` | Installed the exact artifact and declared NumPy 2.5.3; `wheel-install.log`. |

Build/install used approved network-enabled execution because PyPI access had
required leaving the sandbox in earlier steps. No test failure remains. An
intermediate full quick run passed 260 tests (`checks.log`); code inspection
then identified the exact-symmetry/subnormal averaging issue. That was corrected
and a regression test added before the final suite and the single wheel build.

Coverage includes the one-spectrum analytic normalization, fractional counts,
grid/volume scaling, dense multi-node solve oracles, five fields/all 15 spectra,
permuted and cross-only selection, signed derivatives, direct triangular-solve
comparison, factor reuse with refactorization prohibited, batching, and array
ownership. Selected covariance singularity is distinct from parent-field rank.
Near-rank-threshold cases and substantive asymmetry/indefiniteness are tested
under widely disparate units.

Results tests verify the analytic `[[4,1],[1,1]]` conditional/marginalized/fixed
answers; explicit bin-local bindings; complementary shared nuisance constraints;
priors exactly once; correlated/rank-deficient priors; metadata rejection; zero,
disconnected, and null information; diagnostic null residuals; parameter-unit
transformation of information/covariance; no eager solves; roundoff thresholds;
and invalid/overflow inputs. The wheel probe also assembles singular information
from supplied dependent derivative columns and combines complementary responses.

### Fresh installed-wheel probe

Artifact: `.validation/step05/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`.
SHA-256: `72cdba93af0a2c7d84c4c73fcc393f2137b3be7f4cc5308bcd16718d20e732e8`.
No earlier same-version installation was reused.

From `/tmp`, the probe used `-I` and all three thread limits set to one. The exact
absolute-path invocation was equivalent to the following, with ROOT equal to
`/global/homes/a/acuceu/desi_acuceu/vega_dev/lib/fishhighz`:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ROOT/.validation/step05/wheel-env/bin/python" -I \
  "$ROOT/.validation/step05/probe.py" \
  "$ROOT/.validation/step05/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" \
  "$ROOT/.validation/step05/readme-example.py" "$ROOT"
```

`wheel-probe.json` records absolute module origins and results. Assembly,
results, and kernel modules resolved under the fresh environment's
`lib/python3.13/site-packages/fishhighz/`. All 14 installed Python modules matched
the selected wheel, and every wheel module matched current source bytes. Quiet
package import without NumPy or neighboring packages passed. The README example,
assembly/factor reuse, analytic named covariance, conditional/fixed uncertainties,
singular-result rejection, and prior-once combination all passed. Probe Python
and NumPy matched development versions. Home-path workspace and resolved CFS
artifact paths identify the same authorized package location.

## Git state and review limits

`before.json` snapshots preexisting source, tests, reports, scripts, and governance
inputs; `preservation.json` confirms README.md was the only changed preexisting
file. `git-status.txt` and `tracked.diff` preserve the full final Git state,
including prior uncommitted work. Hashes below identify this step's source/tests
and README; `source-sha256.txt` also records this handoff and governance inputs.
The handoff's own hash is external to avoid self-reference.

No real lyaforecast quick capture or full suite ran; neither is required here and
full execution was not requested. Existing reference evidence was preserved.
These tests validate supplied-array mean-response Fisher calculations on Python
3.13, not external-model derivative generation, survey preparation/convergence,
a Python-version matrix, or JIT performance. Marginalization of partially
identifiable combinations inside a singular joint result remains unsupported by
design. No commit, push, planning revision, or Step 06 work was performed. Stop
for user review and independently requested review before progression.

```text
518dec5fcc27800d5b6a40f8d7c83cb4d30ee1c43b9aa2d446d2a228a0e70a8f  README.md
acbbf98c4107ae6c6d6255b2d3f0fdbe656faee690e0a31812e6bd047bf6bbfb  fishhighz/_information.py
cbd888e9cf8724e34cae3ce92da73c753923e11a85fe0c0d7f1a0bcc0f9e4d21  fishhighz/fisher.py
8fbf8977005493980fbaec56abe984368602574777cab5cfcd005cb20ade7aee  fishhighz/kernels/fisher.py
5cacef349f3de6c2c584b0e557afaa44540c452a73680666b4fe1f0bbf78a75e  fishhighz/results.py
31e47f10cd3f62e67c4d417428937d42fa9ae2bc699f8d19c2174380d20e72bd  tests/test_fisher.py
36839cb24f42227bc96671b1b957dd94affcc3b762d690aa8b97aacf60502bdb  tests/test_results.py
```
