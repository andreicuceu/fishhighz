# Step 12, revision 6: convergence evidence for weak individual spectra

Status: independent review passed for the bounded revision-6 evidence repair.
R2 is resolved for the reviewed contract; no further code correction was found.
R3 remains unresolved and requires the user's scientific decision. Step 12 is
not scientifically accepted. Await the user's next request; no new revision,
real run, implementation dispatch or Step 13 is authorized. See
revision-6 review (local-only path: `reviews/step-12-review-r6.md`). The exact reviewed instructions
are preserved in `.validation/step12-review-r6/IMPLEMENTATION_STEP.md`.

Read AGENTS.md (local-only path: `AGENTS.md`), design section 0 (local-only path: `../../FISHHIGHZ_DESIGN.md`), the
roadmap (local-only path: `../../FISHHIGHZ_IMPLEMENTATION_PLAN.md`), the
r5 handoff (local-only path: `reviews/step-12-r5.md`), and r5 review (local-only path: `reviews/step-12-review-r5.md`).
The exact reviewed assignment is preserved in
step-12-instructions-r5.md (local-only path: `reviews/step-12-instructions-r5.md`). The original
revision-3 archive (local-only path: `reviews/step-12-instructions-r3.md`) remains unchanged.

## 1. Bounded assignment

Repair the semantic binding of convergence information to actual saved trials,
including weak individual spectra, and demonstrate rejection through the writer
and offline checker. Preserve the successful revision-5 trial inventory,
controller replay, metadata correction, diagnostics and performance changes.
This revision needs only self-contained tests and offline inspection of the saved
15x2pt evidence. No new real forecast, CAMB preparation, sensitivity run,
NewForecast recapture, or runtime benchmark is authorized by this repair plan.
Do not run the other six DESI-2 cases. Do not implement another weighting rule.

The saved scientific scope remains `lya_qso_lbg_lae_15x2pt`: six redshift bins,
12 primary records across both profiles, 180 individually selected spectrum
results, and 72 diagnostic requests. Revision 5 has six passing compatibility
records, six unconverged accuracy records, 60 completed diagnostics and 12
unavailable diagnostics. Preserve these outcomes; they do not become scientific
passes because the checker is repaired. Any later real run requires an explicit
user request. Any consequential weighting or acceptance change requires the
user's scientific decision before dependent work.

## 2. R2 finding and required correction

The review reproducer is
`.validation/step12-review-r5/weak_pair_probe.py`, with its recorded result in
`weak-pair-probe.json`. It uses a synthetic five-field, 15-spectrum payload.
One spectrum has a Jacobian amplitude of 1e-8 relative to the others and hence
Fisher information of order 1e-16 relative to the dominant spectra. Its last
weight refinements change its marginalized errors by 0.01941932430908011,
exceeding the required 1e-3. Replacing that spectrum's lower metric operand with
the upper one, recalculating the reported metrics and asserting `passed=True`
leaves all actual study arrays unchanged. The validator, writer and offline
checker nevertheless accept it and report completion.

`validation/trials.py` compares the entire pair-Fisher stack with a single
relative norm via `schema._close`. Dominant spectra conceal the weak spectrum's
change. The final replay check uses the same aggregate comparison and does not
compare the replayed numerical metrics or convergence verdict. The existing
whole-stack mutation tests do not cover this dynamic range.

Required behavior:

1. Bind every individual-spectrum Fisher block to its referenced trial using
   that block's own numerical scale. Apply the same requirement to the primary
   versus final-trial comparison and replayed operands. Do not normalize one
   spectrum or trial by another. Retain finite-value, shape, rank/null and
   exact-zero semantics; make scale calculations robust to floating-point range.
2. Independently derive convergence metrics from the referenced successful
   study arrays, and require the reported metrics and convergence verdict to
   agree. A controller replay that disagrees about convergence cannot certify
   the payload. Detached metric operands alone cannot establish a pass.
3. Preserve failed-attempt inventories and scientific-failure inspection:
   `require_pass=False` may admit a truthful unconverged record, but must still
   reject inconsistent numerical evidence or a false pass. Keep the existing
   failed-trial behavior, including unrepresentable arithmetic and absent operands.
4. Inspect related pair-summary consistency checks for the same aggregate-scale
   weakness and correct those needed to enforce this contract. Keep changes
   local to validation/evidence and its tests. Do not change general Fisher,
   covariance, physical weighting, model or response code for this repair.
5. Keep valid identical results at distinct controls, finite roundoff differences,
   and valid null or partially constrained information admissible. Do not require
   artificial nonzero changes or use bitwise equality as a universal substitute.

## 3. Scientific and numerical invariants

Retain the archived r5 scientific conventions: fixed fiducial covariance and
weights in derivatives, independent bins, common volume within each bin,
observed-coordinate k cuts and mode counts, response applied exactly once,
independent P3D/P1D, complete covariance closure and original pair ordering.
Retain wiggle-only ap/at derivatives, signed PK-PKSB and unchanged smooth PKSB,
CAMB growth and explicit geometry/units, retained input-policy labels, strict
reader defaults and the existing cumulative-weight arithmetic guards.

Convergence thresholds remain 1e-3 for relative joint Fisher change, joint
marginalized-error change and the maximum individual-spectrum marginalized-error
change; volume remains 1e-6. Preserve the existing definitions, required
refinement families and bounded trial levels. The 5e-12 evidence consistency
scale is separate from these scientific convergence thresholds. Do not loosen
either to pass the reproducer or saved unconverged results.

Preserve revision-4 batched matrix checks/factorization, optional compiled
contraction, NumPy default and direct NumPy oracle, factor/Jacobian reuse and
bounded three-payload study caching. Retain their regression tests. The NumPy
matrix speed target remains unmet; no new performance work is requested.

## 4. Acceptance tests closing R2

Use tiny deterministic synthetic studies; a real 15x2pt forecast is unnecessary.
Test the semantic validator, writer and offline reader, with consistent rewritten
hashes/inventories where appropriate. Record rejection reasons separately.

- Reproduce the review's 1.9419% weak-spectrum error change at information ratios
  including 1, 1e-8 and 1e-16 relative to the other spectra. Copy only its lower
  operand from its upper trial and recalculate metrics/pass. Every inconsistent
  case must reject, including when all saved study arrays remain unchanged.
- Permute the weak spectrum's location and include both auto and cross-spectrum
  selections. Test isolated-family and combined-refinement operands separately.
- Mutate a weak primary/final pair block and a weak replayed pair block. Verify
  that each mismatch rejects at its own scale even when aggregate changes are
  smaller than 5e-12. Assert the actual comparison causing rejection.
- Include a truthful unconverged control admitted with `require_pass=False`, and
  a converged weak-spectrum control admitted normally. Distinct controls giving
  identical valid matrices must pass. Harmless within-tolerance roundoff must
  pass at both large and small scales. Retain exact-zero, null/unconstrained,
  finite-range and nonfinite rejection tests with explicit expected behavior.
- A changed verdict or metric inconsistent with the saved-trial reconstruction
  must reject even if other report fields remain consistent. Preserve missing,
  duplicate, reordered, unrelated and failed-trial rejection coverage from r5.
- Retain all five R1 and original R2 repeated-operand regressions. Exercise the
  exact 12-primary/72-diagnostic inventory gate with a self-contained fixture;
  reject omissions, extras, duplicate bins and mislabeled payload kinds. A
  one-record writer probe does not alone demonstrate full scoped acceptance.
- Inspect the saved r5 primaries offline with per-spectrum normalization. Their
  primary C/F/error values and truthful six accuracy failures must be preserved.
  Keep all 12 diagnostic failures with reasons. The scientific gate must remain
  false; do not treat this expected outcome as an implementation-test failure.

Avoid adding a large real-shaped allocation to every parametrization. The saved
review reproducer is a concrete starting point, not a required test-array size.
Use independent analytic two-parameter inverses/errors where possible, rather
than duplicating the repaired validator in the tests.

## 5. R3: decision retained, no scientific alternative selected

Revision 5 establishes a useful diagnosis. Holding the initial weight function
fixed gives magnitude-order 32/64 individual-spectrum errors consistent at about
1e-14. The cumulative update changes the weight function as the quadrature is
refined; its permitted refinements still fail the 1e-3 requirement. The saved
zero-weight Jacobians have spectral radii below one (maximum 0.4072476439), and
some intermediate products are below float64 range.

The unnormalized weights tending to zero does not prove that the scale-invariant
noise coefficients lack a finite limit. No normalized asymptotic limit or its
quadrature independence has yet been established. Preserve the handoff's three
unselected alternatives: an explicitly finite-iteration convention with a user
review of its acceptance requirement; investigation of a scale-normalized limit
of the same rule; or an independently justified different estimator. Ask the user
before choosing or implementing any of them. Do not waive convergence, alter
floors/support, or silently normalize weights as part of R2. The evidence repair
can be completed without resolving this scientific choice.

## 6. Quick checks, preservation and handoff

Run the ordinary checks with one-thread settings:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PATH="$PWD/.venv/bin:$PATH" scripts/check.sh
```

Retain all 1208 current parametrized tests (r5 independent review: 1183 passes,
25 optional-compiler skips), then add the focused regression cases. Run affected
validation/performance tests with both available backends. Preserve the disclosed
r4 inherited optional-dependency failure when SciPy is deliberately blocked while
the compiled backend is forced; report it separately from default-NumPy checks.
Do not silently count that combination as a pass or expand scope to repair it.

Build one exact wheel, verify source/wheel/installed module and METADATA/WHEEL
identity, and run the affected regressions and six examples outside the checkout
with isolated imports. Record interpreter/backend/versions and any skipped or
blocked checks. Use immutable r5 inputs for offline checks; do not overwrite
`.validation/step12-r5-20260914T191855Z/`, `.validation/step12-review-r5/`, the r3/r5
instruction archives, r2 handoff, r4 handoff/review, or earlier evidence.
Before replacing the active r5 handoff, archive it byte-for-byte under a distinct
revision-qualified name, verifying rather than overwriting any existing target.

Write new evidence under a unique `.validation/step12-r6-.../` directory. The
handoff must identify the exact plan/source/wheel, list changed files and commands,
show before/after results for every mutation and valid control, and distinguish
new quick checks from historical real forecasts. State whether R2 is repaired,
retain R3 and the six accuracy failures, and link the preserved scientific figures.
Do not claim scientific acceptance. Stop for independent/user review; no agent
dispatch, next step, commit, push, Slurm action or new real forecast.
