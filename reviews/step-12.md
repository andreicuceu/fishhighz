# Step 12 revision 6: weak-spectrum convergence evidence

**R2 is repaired in this implementation and ready for independent/user review.
R3 remains unresolved. Step 12 is not scientifically accepted.** This revision
contains synthetic tests and offline inspection only. No real forecast, CAMB
preparation, sensitivity calculation, runtime benchmark or weighting change was
performed. The preceding handoff's R2 closure claim was disproved by independent
review; that document is preserved byte-for-byte as [step-12-r5.md](step-12-r5.md).

New evidence: [step12-r6-20260914T222237Z](../.validation/step12-r6-20260914T222237Z/).
The exact dispatched plan, initial source/tests and Git status are snapshotted
there. `preservation.json` records the plan hash, changed-file hashes, all 35
unchanged non-validation modules and the historical array/report checksums.
The r3/r5 instruction archives, r2 handoff, r4 performance handoff/review and all
previous scientific evidence remain unchanged. No governance document was edited.

## Evidence correction

The previous whole-stack norm could hide a weak spectrum's changed Fisher matrix
behind stronger spectra. `validation/schema.py` now compares each spectrum's
Fisher, covariance, errors, correlation, rank and constrained mask separately.
`validation/trials.py` applies the same matrix-local normalization to primary
versus final-trial information, each metric operand and every replayed operand.
No spectrum, refinement family or trial supplies another's normalization.

`validation/numerics.py` scales both operands before subtraction and squaring,
retaining the Frobenius relative discrepancy while avoiding norm overflow and
underflow. An exact-zero reference still requires exact zero. Tests cover
information scales from 1e-300 to 1e300, nonfinite rejection, null and partially
constrained information, identical matrices at distinct controls and harmless
roundoff. The consistency tolerance remains 5e-12; scientific convergence limits
remain 1e-3 for Fisher/joint-error/individual-error changes and 1e-6 for volume.

The saved-trial controller replay now supplies an independent reconstruction of
the reported metrics and convergence verdict. These must agree even under
`require_pass=False`; truthful unconverged records remain inspectable. Failed
attempts and their errors remain in the replay, including the failed-combined
trial placeholder. Detached metric operands cannot establish convergence.
No covariance, general Fisher, response, model, weighting or compiled kernel was
changed. The existing three-payload cache and all earlier regressions remain.

## Mutation and valid-control evidence

The tiny reproducer retains five fields, all 15 auto/cross-spectra and the actual
bounded controller, while substituting four deterministic quadrature nodes.
The weak Jacobian gives information ratios 1, 1e-8 and 1e-16. An analytic inverse
of a scalar multiple of a two-parameter Fisher matrix gives the final weight
error change **1.9419324309%**, independently of that amplitude.

For the review's complete lower-operand substitution at ratio 1e-16, all actual
study arrays remain unchanged. The old wheel accepts the false zero-change pass
through validator, writer and offline reader. The new wheel rejects all three
at `metric/trial operand weights 0 metric_pair_fisher block (0,)`.
See `weak-probe-before.json` and `weak-probe-after.json`; these deliberately
synthetic probes disable source-path verification and are not real forecasts.

The new 60 parametrized cases test isolated weights, combined refinements and
complete substitution separately; positions 0, 1, 8 and 14 include auto and
cross spectra. They also test weak primary/final and failed-trial replay blocks,
detached metrics, both inconsistent verdict directions, pair summaries, valid
weak/roundoff/null controls, range failures and exact scoped inventory.
`mutation-results.md` and `mutation-results.json` list every before/after result;
`after-numpy.xml` records validator, writer and offline rejection reasons
separately. All saved corruptions have consistent hashes/inventories, so rejection
reaches numerical semantics.

The old wheel passes 30 of these tests and fails 30. These are **test outcomes**,
not 29 falsely certified forecasts: eight isolated weak mutations were already
rejected because another metric still failed, several assertions require the
new precise rejection location, and seven exercise the new comparison helper.
The full substitution and replay/summary holes are distinguished in the table.
All 60 new tests pass with the repaired wheel.

The complete synthetic scoped fixture passes the actual 12-primary/72-diagnostic
gate, retaining 180 selected-spectrum results. Omissions, extras, duplicate bins,
mislabeled kinds and missing diagnostics reject. This fixture uses synthetic
arrays and explicit fixture provenance; it demonstrates inventory and evidence
semantics, not physical adequacy or real scientific acceptance.

## Checks and exact installation

- Ordinary one-thread suite: **1243 passed, 25 skipped in 185.17s (0:03:05)** (`quick-complete.log`).
- Installed explicit-Numba affected regressions: **251 passed in 120.46s (0:02:00)** (`after-numba.log`).
- Installed NumPy validation/performance regressions: **435 passed in 158.57s (0:02:38)** (`after-numpy.log`).
- Ruff lint/format and `git diff --check`: passed. All 1208 previous parametrized cases retained; 60 added.

The exact wheel SHA256 is `7286a2c34014807c2d3309f9b9f140594aeb64fd8fe415932ba234cae399697c`. `installed-identity.json` verifies all
52 source/wheel/installed modules, METADATA/WHEEL bytes, import origin and isolated
execution outside the checkout. All six examples pass there. The fresh environment
uses Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, Numba
0.67.0, llvmlite 0.49.0 and Astropy 8.0.1. No reference
package or historical environment was installed into or modified.

`commands.md` records commands and backend selection; logs retain the sandboxed
build/install DNS failures and the successful approved retries. The first candidate wheel is retained in `dist-initial/`: offline inspection exposed exact compatibility-metric dictionary equality after the range-safe norm change. The corrected final wheel is in `dist/`; all final installed checks use it. The compatibility comparison now uses the existing metric consistency rule (rtol 5e-12, atol 1e-15); two additional controls accept roundoff and reject a false metric. `comparison-roundoff.json` records the original differences. No scientific threshold changed. The inherited explicit-Numba/SciPy-blocked subprocess test still fails
with `scipy 0.16+ is required for linear algebra`; its default-NumPy counterpart
passes in the installed suite. This disclosed combination was not repaired or
counted as a compiled pass. The r4 NumPy matrix speed target remains unmet; this
revision makes no new timing claim.

## Saved scientific results and remaining decision

The installed repaired checker inspects all saved r5 primaries with
`require_pass=False`. The independent saved-array Wick/direct-solve oracle and
analytic two-parameter inverse check all 12 primary and 180 individual-spectrum
results. Maximum relative discrepancies are covariance 0, fisher 1.34e-15, pair_fisher 2.2e-15, errors 1.33e-15, pair_errors 2.11e-15. Every successful metric operand is also checked at
its own spectrum's scale. No primary C/F/error values are changed.

Compatibility remains **6/6 passing**; accuracy remains **0/6 passing**.
All 72 diagnostic requests remain: 60 completed and 12 unavailable, with original
reasons retained in `preservation.json` (floor_low in every bin, remove_bright in
bins 0–4, width_plus in bin 2). The scoped scientific gate still returns
`partial or scientifically failed validation`, as expected. This is not an
implementation-test failure and does not authorize a further real run.

The preserved [joint comparison](../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt.png),
[15 selected-spectrum pages](../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt-pairs.pdf),
[attribution](../.validation/step12-r5-20260914T191855Z/attribution/lya_qso_lbg_lae_15x2pt.png),
[sensitivities](../.validation/step12-r5-20260914T191855Z/sensitivities/lya_qso_lbg_lae_15x2pt.png)
and [magnitude/weight diagnosis](../.validation/step12-r5-20260914T191855Z/weight-figures/magnitude-convergence.png)
are historical r5 figures, not regenerated results.

The finite cumulative weighting refinements remain unconverged. Vanishing
unnormalized weights do not exclude a finite scale-invariant noise limit, whose
existence and quadrature independence remain unestablished. Three alternatives
remain unselected: an explicit finite-iteration convention with user review of
its convergence requirement; investigation of a normalized limit of the same
rule; or an independently justified estimator. No floor, support, tolerance or
weight normalization was changed. That scientific choice remains with the user.

Stopped for independent/user review. No Step 13, agent dispatch, Slurm action,
commit, push or new real forecast occurred.
