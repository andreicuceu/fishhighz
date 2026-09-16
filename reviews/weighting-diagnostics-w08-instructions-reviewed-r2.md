# Forest weighting W08, revision 2: close the positive-product range gap

Status: revision 2 is implemented and passes independent scientific review.
The positive `l_p*v` underflow finding W08-R1 is closed; exact-zero and exactly
representable subnormal controls pass, as do the ordinary analytic, noise,
survey and saved-reference checks. No consequential correction is needed.
Awaiting user review and acceptance; this pass does not adopt a profile default
or authorize progression. The selected fixed-reference option, cumulative
compatibility, package Step 13 and forecast acceptance remain unchanged.
Created: 2026-09-16.

## 1. Scientific objective and bounded change

Expose the tested fixed reference as an explicit weight-preparation option:
`method="inverse_variance"`. It prepares nu=B_star/(B_star+l_p v) once from
supplied arrays and a positive scalar B_star, with no cumulative updates and
no dependence on a 3D weighting signal. The existing `legacy` and `supplied`
methods retain their numerical behavior and contracts. No default is changed.

W05 identifies nu with the legacy initialization and establishes its minimum
of Q=A B_star+P_pixel for a nonnegative measure under the stated covariance
approximation. W06 verifies its finite magnitude-refinement stability. W07
finds that, for one accuracy-profile QSO-forest auto-spectrum/bin, cumulative
three-update radial/transverse BAO errors exceed the reference by
7.63057%/10.32563%; six-update errors exceed it by 8.61956%/11.70862%.
These are conditional results, not a survey-wide or multi-mode optimum.

W08 is a small package implementation step within this separately reviewed
sequence. It replaces no estimator automatically. Revision 1 added the explicit
method using the existing ForestWeights and noise APIs. Its independent review
reproduces the ordinary rational coefficients, supplied/legacy equivalence,
response/noise conversion, fixed differentiation semantics and saved W06
coefficients. The only required revision-2 change is to reject an inexact
underflow in the positive instrumental-power product before it can alter a
weight. No new public class, mode axis, provider protocol, survey configuration
layer or production optimizer is needed.

## 2. Read and verify live contracts

Read workspace/package AGENTS.md, design section 0, main roadmap handover/
progress register, the [diagnostic roadmap](../../FISHHIGHZ_WEIGHTING_DIAGNOSTICS_PLAN.md),
[W07 report](reviews/weighting-diagnostics-w07-r1.md),
[W07 review](reviews/weighting-diagnostics-w07-review-r1.md), and the reference
conventions in [W05](reviews/weighting-diagnostics-w05-r1.md) and
[W06](reviews/weighting-diagnostics-w06-r1.md). The
[reviewed W07 assignment](reviews/weighting-diagnostics-w07-instructions-reviewed-r1.md)
is archived. Also read the [W08 r1 handoff](reviews/weighting-diagnostics-w08-r1.md),
the [W08 r1 review](reviews/weighting-diagnostics-w08-review-r1.md), and the
[exact r1 assignment](reviews/weighting-diagnostics-w08-instructions-r1.md).
`IMPLEMENTATION_STEP.md` remains the separate package Step 13.

Inspect `fishhighz/weights.py::prepare_forest_weights`,
`fishhighz/kernels/weights.py`, `fishhighz/noise.py::forest_noise`,
`ForestInput` in `survey.py`, and forest preparation in `forecast.py`.
Existing tests cover supplied weights, explicit legacy counts (including zero),
arithmetic guards, immutable snapshots, and model-call counts. Preserve all
uncommitted/untracked work and historical evidence; do not reset or refactor it.

The live API already permits the reference through supplied weights or legacy
zero iterations. W08 gives it a scientifically explicit name and removes the
unnecessary S/iteration arguments for this option. Do not interpret that naming
change as a different estimator from the reviewed nu.

## 3. Public API and scientific conventions

Retain this method in the existing keyword interface:

```python
prepare_forest_weights(
    field, geometry, response,
    z_source=z_source, magnitudes=m, quadrature=q,
    rho=rho_v, variance=v, length_velocity=L,
    method="inverse_variance", alias=B_star,
)
```

- Reuse `alias` as the scalar B_star argument; do not add a synonymous parameter.
  It is strictly positive, finite, and in km/s. It is the fiducial 1D forest
  power including the instrument response squared at a caller-selected weighting
  mode. There is no default mode or B_star. The caller retains responsibility
  for supplying the appropriate fixed P1D value and response convention.
- Require `weights`, `iterations`, `signal` and `auxiliary` to be None for this
  method. Even `iterations=0` is a conflicting setting: it describes the legacy
  interface, not this iteration-free method. Give clear ValueError messages.
- W08 supports direct scalar B_star only. Do not add automatic P1D sampling,
  change `sample_auxiliary`, or allow `ForestInput.auxiliary_coordinates` for
  this method. That existing joint S/B route remains a legacy facility. A
  caller can use a separately prepared B_star in `weight_options` with no
  auxiliary coordinates. This keeps all provider calls outside the kernel and
  avoids introducing an unnecessary P3D dependency.
- Apply all existing shape, ordering, nonnegative-density/variance, strictly
  positive quadrature, redshift, length, response and weighted-support checks.
  No density floor, support extension, sorting or interpolation is added.

For positive-density nodes, set

```text
r_i = rho_i q_i
w_i = B_star / (B_star + l_p v_i).
```

At zero-density nodes store w_i=0, matching the legacy seed's canonical support
convention. With v_i=0 on positive support, w_i=1 exactly. Do not normalize the
weights or evolve them. Use a small numeric helper only if it simplifies the
host/kernel separation; do not call `_iterate` with a dummy signal or alter
its legacy initialization/update arithmetic.

Use the existing `_integrals` and final coefficient checks unchanged:

```text
I1 = sum r_i w_i; I2 = sum r_i w_i^2; I3 = sum r_i w_i^2 v_i
A = I2/(L I1^2); P_pixel = l_p I3/(L I1^2).
```

Retain the common underflow/overflow and representability policy. Positive
support must not silently acquire zero/nonfinite weights through arithmetic;
nonrepresentable integral products still fail, including mixed contributions
hidden by a positive total. Exact zeros remain valid. Do not broaden the
extreme-input domain, add logarithmic/rescaled coefficients, floor results or
relax guards in this step. Report a consequential obstruction for review.

Revision 2 must also reject an inexact underflow of the positive product
`l_p*v_i` before adding it to B_star. The r1 implementation evaluates this
product inside an `under="ignore"` block. With `l_p=B_star=nextafter(0,1)` and
`v_i=1/2`, it becomes zero and the accepted weight is one, although the positive
instrumental contribution is not representable. Use a local raised-underflow
operation or an equally direct check only for this product. Preserve exact
zero variance and zero-density behavior. Do not change the legacy calculation,
`_iterate`, `_integrals`, ordinary arithmetic order, or the accepted numerical
domain beyond rejecting this already unrepresentable case.

The returned existing ForestWeights object records method="inverse_variance",
alias=B_star, signal=None, iterations=None, auxiliary=None and an empty
weight_changes array. Keep owned immutable arrays and context validation.
Expose the ordinary I1/I2/I3 prefixes as before; prefixes in the output do not
imply prefix feedback in the preparation.

`forest_noise` remains unchanged: at a forecast mode, it uses the intrinsic
P1D supplied there and its response once,

```text
N_F(k,mu) = [A P1D(q) W(q)^2 + P_pixel] d_deg^2/a_v,
q = k mu/a_v.
```

B_star is used to prepare weights only; do not substitute it for every mode's
P1D, apply its response twice or smooth the pixel noise. Weights, A and P_pixel
remain fiducial and fixed while mean-model parameters are differentiated.
No claim of per-mode reoptimization is made. Units, common L/l_p assumptions,
independent P3D/P1D and sampling-independence rules remain unchanged.

## 4. Minimal acceptance checks

Use focused self-contained tests, with one bounded saved-array check outside
the ordinary test suite. Do not reproduce all W01–W07 diagnostics.

1. **Independent analytic coefficient test.** Use q=(1/2,2), rho=(2,1/2),
   v=(1,4), L=l_p=B_star=1, so r=(1,1). Require w=(1/2,1/5), A=29/49,
   P_pixel=41/49 and A B_star+P_pixel=10/7 from exact rational arithmetic.
   Check nonuniform measure enters once, and compare with `supplied` using
   independently formed nu and with unchanged legacy zero iterations. Use a
   tiny zero-density/zero-variance variant to check the stated support convention.
   Confirm the new branch never invokes `_iterate`.
2. **Method contract and retained safety.** Check the explicit conflicting
   keywords above, nonpositive/nonfinite B_star, immutable result/metadata,
   and one ordinary positive example for each QSO/LBG field identity. Existing
   common invalid-input and weight-range tests must still pass. Add only the
   new-branch assertions needed to demonstrate that zero/nonfinite computed
   weights and unrepresentable products cannot bypass the existing checks.
   No exhaustive validator or hypothetical mutation inventory is required.
   Add the concrete r1 review regression with
   `tiny=np.nextafter(0.0,1.0)`, `l_p=B_star=tiny`, positive density and
   `v=(0,1/2)`. Use masses for which the current final A and P_pixel remain
   representable (for example q=(1,1), rho=(1/8,1/8), L=1). Revision 1
   incorrectly accepts w=(1,1), A=4 and P_pixel=tiny; revision 2 must raise a
   contextual ValueError for the underflowed positive product. Pair it with an
   exact-zero-variance control that still returns one on positive support.
3. **Noise and survey integration.** On a tiny synthetic field/bin, compare
   direct new-method preparation with independently supplied nu: weights,
   coefficients and noise at at least two different q values must agree.
   Check the response/conversion formula above with a scalar calculation.
   Exercise `ForestInput` -> `prepare_bin` -> the existing small `run_bin`
   fixture: no auxiliary S/B query, no extra provider call for weight generation,
   and no weight/noise/P1D reevaluation during mean finite differences. Compare
   the new and supplied-weight results. Extend `test_generated_spies` or add
   one equally small test. Synthetic execution is allowed; no real forecast.
4. **One saved-reference coefficient check.** Load only the accuracy bin-0
   order-32 QSO sample from the report identified below. Form rho_v from its
   saved normalized density with `(1+z_source)/c`, use its own quadrature,
   variance, L=44147.09373976799, l_p=63.328211161339375 and
   B_star=16.356748967852234. Prepare with the new public method and compare
   against W06's reference:

   | A [deg^2] | P_pixel [deg^2 km/s] | Q [deg^2 km/s] |
   | ---: | ---: | ---: |
   | 0.0294198456633004 | 0.314354922180624 | 0.795567952368186 |

   Independently reconstruct the three moments and coefficients using scalar
   80-digit Decimal arithmetic from the exact saved binary inputs, forming nu
   independently. Do not use the implementation's helpers for the oracle.
   A valid synthetic test geometry with the saved z_eval=2.1152848986890427
   and constant supplied backgrounds is sufficient for this coefficient-only
   check: A/P_pixel depend on L/l_p and the arrays, not a_v/d_deg or volume.
   Label that context honestly; do not claim a geometry or forecast replay.
   Match the field/response context used for preparation. No raw reader,
   interpolation, cosmology/P3D/P1D call or model evaluation is needed.

For ordinary nonzero quantities use rtol=5e-12, atol=0, exact reference-zero
handling and finite checks; do not loosen tolerances to accommodate a mismatch.
The scalar control can use tighter exact/roundoff comparisons. Preserve actual
errors and stop for review on a scientifically relevant discrepancy.

The one saved report is
`.validation/step12-r5-20260914T191855Z/profiles-checked/records-001.report.json`,
SHA-256 `9b64dd91adfb84d3c7a5ecb5fb1db0b4206abcd4f8c1ebfe1124698a97a28666`.
Verify accuracy/bin/field/order identity, and select
`settings.samples["lya(qso)"]`. Its nonnegative density/SNR policies are already
applied; neither reapply them nor endorse them physically. The normal pytest
suite must stay self-contained and must not depend on this local file.

## 5. Edit scope and quick validation

Expected revision-2 production edit: `fishhighz/weights.py`; a small additional
numeric helper in `fishhighz/kernels/weights.py` is allowed only if necessary.
Retain the r1 method, metadata, documentation and focused tests. Preserve
`_iterate` and `_integrals` numerically. Survey/forecast/noise implementations
should accept this through their existing direct-options route; do not modify
them or broaden callable preparation in this step. A relevant docstring can
be clarified without changing behavior. If a production integration change
beyond this boundary proves necessary, document the concrete obstruction for
review rather than silently extending the architecture.

Extend `tests/test_inverse_variance_weights.py` only with the new regression and
exact-zero control; no new test module is needed. The README's existing r1
description of the direct-scalar option, response ownership, fixed-fiducial
semantics and conditional optimality claim remains applicable. Do not change it
unless the error description needs one minimal clarification. Do not change
profile recipes, legacy adapters, CLI/INI defaults, public result classes or
existing examples except a minimal relevant documentation snippet.
Compatibility remains legacy.

Run the new focused tests, existing `test_weights.py`, `test_weight_range.py`,
`test_weight_models.py`, `test_noise.py`, and the affected tiny survey spy/
equivalence tests (not the complete forecast suite). Run Ruff lint/format on
changed Python and `git diff --check`. Use the existing interpreter, one process
and numerical thread limits of one. The saved-array check has a 30-second cap;
the focused pytest invocation has a 60-second cap. Retain partial evidence on
any cap, report it and do not automatically extend the scope or run a full suite.
No wheel, new environment, broad performance exercise or forecast controller.

## 6. Handoff and independent review

Preserve `reviews/weighting-diagnostics-w08-r1.md`. Write the revision-2 handoff
as `reviews/weighting-diagnostics-w08-r2.md`. Any new check output goes in a
new `.validation/forest-weight-diagnostics/w08-r2-<UTC>/`; do not overwrite r1
or review evidence. Lead with closure of the one range finding and what the
unchanged ordinary checks establish.
Include exact commands, changes, source/input identity, oracle agreement,
legacy regression results, provider-call evidence, limitations and failures.
State clearly that W07 BAO improvements remain historical evidence; W08 has
not measured a new survey improvement or changed any profile's results.
Do not edit planning/governance files during implementation.

The reviewer checks only the new method, affected contracts and decisive
independent arithmetic. Verify that signal/iterations are unnecessary, legacy
behavior and safety checks are preserved, noise response is counted once, and
the survey route freezes weights during differentiation. Write
`reviews/weighting-diagnostics-w08-review-r2.md`, beginning with **Scientific
conclusion and impact of proposed changes**. Propose corrections only if likely
to affect the scientific conclusion or its justified scope, including a real
scientific regression of compatibility/noise behavior. Give the possible impact
and smallest resolving check; explicitly state when none are needed.

Stop for independent/user review. No profile-default adoption, changed physical
input policy, mode-dependent weighting, new real forecast, full 15x2pt suite,
package Step 13 repair, subsequent step, Slurm, agent dispatch, commit or push
is authorized. The user retains acceptance and progression.
