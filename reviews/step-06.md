# Step 06 handoff — revision 2

Implemented `IMPLEMENTATION_STEP.md`, revision 2 (2026-09-12), addressing finding
R1 in [the revision-1 review](step-06-review-r1.md). Ready for user and independent
review; acceptance and progression are not claimed. Step 07 has not started.
Base commit remains `0d69786a06d5d676564a51fad14a7156951c2228` with prior
uncommitted Steps 03–06 preserved.

## Repair and unchanged contracts

`check_convergence` now compares consecutive actual perturbation-point sets for
every numerical provider/global column throughout the full-study preflight.
Sorted two-point tuples use exact float64 equality, independent of point order;
there is no approximate comparison or inference from derivative agreement.
Unchanged points raise `ValueError` before any model or analytic Jacobian call,
even when the duplicate first occurs at a later refinement or other providers
have valid schedules. Moving only one point remains an effective refinement.

The error identifies provider label, global parameter ID, one-based refinement
level, both requested steps and the repeated actual points. It explains that no
effective refinement occurred and recommends a better-resolved starting step or
fewer refinements. Steps and strategies are never changed automatically.

Preflight retains only the previous and current structural schedules, so memory
does not grow with refinement count. It creates no model-value arrays or cache.
Ordinary `evaluate_derivatives`, actual-offset coefficients, mixed methods,
tied global perturbations, output ownership and legitimate stencil-family
changes retain their reviewed behavior. The public signatures and numerical
kernels are unchanged. The convergence docstring and README explain the new
failure condition and remedy.

The retained interfaces remain `BoundParameters`, `P3DProvider`, `PreparedP3D`,
`evaluate_p3d`, independent `evaluate_p1d`, `evaluate_derivatives`, and
`check_convergence`. Exclusive required-pair routes use original field indices;
selected means/Jacobians use `selected_to_required`. Supplied local derivatives
sum tied columns through `map_jacobian`; numerical changes perturb global entries
before gathering locals. Numerical steps are explicit and bounds-aware. Each
provider uses one fiducial call plus two calls per numerical global dependency,
and one Jacobian call when analytic columns are used. `_combine_three` and
`_scatter_column` remain array-only kernels. No covariance, weights, P1D or noise
is differentiated or recomputed by this repair.

## Tests and outcomes

Evidence directory: `.validation/step06-r2-20260912T170452/` (`$RUN` below).

- Focused derivative suite: **36 passed in 0.22 s**, including five new cases.
- Complete ordinary quick runner, executed **once** for this revision:
  **332 passed in 17.60 s**, Ruff lint passed, 51 files passed formatting.
- Fresh-wheel probe from `/tmp` with `-I`: passed source/wheel/installed module
  and metadata identity, duplicate-stencil rejection, quiet package import,
  tied external derivative/Fisher calculation and the retained working example.
- Final handoff formatting and `git diff --check`: passed.

New regressions reproduce exp(x), x=1.5, h=1.2*spacing(x), and the later-level
repeat at h=2.4*spacing(x), refinements=2. Both reject before any calls. Mixed
analytic/numerical columns and a separate valid provider also reject with zero
model/Jacobian calls, in both provider orders. The single-step derivative path
still executes normally. A successful x=1, h=1.2*spacing(x) shifted-linear test
has one identical point and one changed point, exercising exact rather than
approximate comparison. Existing zero/canceling, well-resolved nonlinear,
second-order error decrease and forward-to-central cases still pass.

The installed probe reports these diagnostics at levels 1 and 2 respectively:
requested steps `2.6645352591003756e-16` and `1.3322676295501878e-16`, repeated
points `(1.4999999999999998, 1.5000000000000002)`, provider `exp`, parameter `x`,
and **zero calls**. Neither study returns a successful convergence assessment.
See `wheel-probe.json` for the complete messages.

The unchanged installed example reuses one fixed covariance factorization.
At step scales 1, 1/2 and 1/4, maximum analytic/FD Fisher relative differences
remain `8.6415264e-6`, `2.1603755e-6`, `5.4009350e-7`. Analytic marginalized
errors are `[0.7354307832838902, 0.34736532355330246]`. Explicit convergence
passes, and independent P1D returns `[4, 3.6363636363636362, 2]` km/s.

## Commands and exact installation

Commands below ran from the component root unless otherwise noted. Set `$RUN`
to the evidence directory above. Lightweight checks used
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`; no Slurm work ran.

```bash
.venv/bin/python -m pytest tests/test_derivatives.py -q > "$RUN/focused.log" 2>&1
PATH="$PWD/.venv/bin:$PATH" scripts/check.sh > "$RUN/checks.log" 2>&1
.venv/bin/python -m build --wheel --outdir "$RUN/dist" > "$RUN/build.log" 2>&1
.venv/bin/python -m venv "$RUN/wheel-env"
"$RUN/wheel-env/bin/python" -m pip install numpy==2.5.3 > "$RUN/numpy-install.log" 2>&1
"$RUN/wheel-env/bin/python" -I - <<'PY'
import importlib.util
assert importlib.util.find_spec('fishhighz') is None
PY
"$RUN/wheel-env/bin/python" -m pip install --no-index --no-deps \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" > "$RUN/wheel-install.log" 2>&1
```

The new environment had no preexisting FishHighz (`fresh-environment.log`).
One wheel was built in the fresh directory and installed by its exact filename,
without index/dependency resolution. Isolated build dependencies and NumPy used
approved network access; no shared environment changed. No build/test failure
occurred during this revision. The focused subprocess emitted NERSC MUNGE socket
messages after pytest success; its exit status was zero and no scheduler action
was requested.

With absolute `$RUN` and `$ROOT` paths, the extended probe ran from `/tmp`:

```bash
"$RUN/wheel-env/bin/python" -I "$RUN/probe.py" \
  "$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl" "$ROOT" \
  "$RUN/external_forecast.py" > "$RUN/wheel-probe.json"
```

The probe executes the copied standalone example using the installed package.
All **17 package Python modules** match source/wheel/installed bytes, and wheel
METADATA/WHEEL match installed metadata. The example copy matches source bytes.
Module origins resolve to the fresh environment's
`lib/python3.13/site-packages/fishhighz/`; the probe checks that the working
directory is outside the source checkout. The home workspace and resolved CFS
paths refer to the same authorized package location.

Exact wheel: `$RUN/dist/fishhighz-0.1.0.dev0-py3-none-any.whl`.
SHA-256: `54ea2f32b1613cf97ee1c1fbb9cb649a37245f9a0d643ac0c46d1684a3db86b8`.
Development/probe Python: 3.13.15; NumPy: 2.5.3. Development tools: pytest 9.1.1,
Ruff 0.16.7, build 1.6.1, pip 26.2.1. Package version remains 0.1.0.dev0 and
NumPy is still its only runtime dependency. See `dev-versions.json`.

## Preservation, source identity and limits

This revision changes only `fishhighz/derivatives.py`, `tests/test_derivatives.py`,
README.md and this handoff. `before.json` and `preservation.json` verify all other
preexisting maintained files are unchanged, including planning/governance files
and the independent revision-1 review. No source contracts or other scientific
code were modified. `tracked.diff`, `git-status.txt`, `base-commit.txt` and
`source-sha256.txt` record the final dirty state; the hash manifest includes this
report without self-reference. `revision.diff` isolates this revision's changes
from the pre-dispatch snapshot.

The revision-1 handoff is preserved verbatim in `$RUN/step-06-r1.md`.
Revision-1 logs/wheel/probe remain in `.validation/step06-20260912T163621/`;
the independent review and reproducer remain in
`.validation/step06-review-r1/` and `reviews/step-06-review-r1.md`. The old
reproducer asserts the formerly broken behavior, so it was not rerun or edited.
No reference evidence was recaptured or modified.

No real lyaforecast quick capture or full suite was run or required. Synthetic
checks establish this preflight repair and retained toy behavior, not external
cosmology accuracy, survey convergence, autodiff, JIT speedup or a Python-version
matrix. Distinct stencils alone still do not prove derivative accuracy.
No commit, push or Step 07 work was performed. Stop for user and independent
review before progression.

Changed source/documentation hashes (handoff hash is in the external manifest):

```text
1ef979411db49033d23095a287bcac1bc5ada8f277160dd8ad8e36b29349bb97  fishhighz/derivatives.py
0e86ff3aac5b0b1b41aa33b993a21105cd7a4ac0b4062f85f34e554dfa0f8b30  tests/test_derivatives.py
435cf9cea12d6ad6b9034d705a2c45c65999bff6cf221b95c281c23b7885332f  README.md
```
