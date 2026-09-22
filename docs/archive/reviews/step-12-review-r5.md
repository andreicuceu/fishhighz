# Step 12 revision 5 independent review

**Changes required: R2 is not fully repaired. R3 remains unresolved and requires
a user scientific decision. Step 12 is not scientifically accepted.**

Reviewed the [r5 handoff](step-12.md), live code/tests, exact installed wheel and
saved numerical evidence against the
[archived r5 assignment](step-12-instructions-r5.md). Fresh review artifacts are
in .validation/step12-review-r5 (local-only path: `../.validation/step12-review-r5/`). No real
forecast, new installation, production fix, agent dispatch, Slurm action, commit
or push was performed during this review.

## Required finding: R2 weak-spectrum convergence can still be falsely certified

**P2 — Bind each individual spectrum at its own scale and check the replayed
metrics/verdict.** At trials.py:170 (local-only path: `../fishhighz/validation/trials.py#L170`),
the metric/trial comparison passes the entire 15-spectrum Fisher stack to
schema.py:203 (local-only path: `../fishhighz/validation/schema.py#L203`). Its global relative norm
allows dominant spectra to hide inconsistencies in a weak spectrum. The
replay check (local-only path: `../fishhighz/validation/trials.py#L217`) omits the replayed numerical
metrics and `passed` verdict, then applies the same aggregate norm to all operands.
The primary/final pair comparison at trials.py:124 uses the same pattern.

The fresh synthetic reproducer
weak_pair_probe.py (local-only path: `../.validation/step12-review-r5/weak_pair_probe.py`) gives one
of 15 spectra a Jacobian amplitude of 1e-8 relative to the others. Its last
weight refinements change marginalized errors by **0.01941932430908011 (1.9419%)**,
above the 0.1% requirement. Replacing only that spectrum's lower metric operands
with its upper operands, recalculating metrics and setting `passed=True` reduces
the claimed change to zero. All saved study Fisher, pair-Fisher and volume arrays
remain unchanged. The semantic validator accepts, the writer reports
`complete=True`, and the offline checker also reports `complete=True`.

weak-pair-probe.json (local-only path: `../.validation/step12-review-r5/weak-pair-probe.json`) records
both metric sets and these outcomes. This is a deliberately synthetic, real-shaped
per-record evidence test, with fixture provenance and source-path verification
disabled in the offline call. It is not a real scientific result or a demonstrated
bypass of the full 12-primary/72-diagnostic scope inventory gate.

The original repeated-whole-stack mutation is rejected by r5, but that does not
close R2. Require per-spectrum primary/trial/operand consistency and independent
metrics/verdict reconstruction from actual trial arrays, with weak-spectrum
mutations and valid small-scale controls through both writer and offline reader.
The [revision-6 instructions](../notes/IMPLEMENTATION_STEP.md) specify these checks.
No change to physical weighting or a new real forecast is needed for this repair.

**No corresponding inconsistency was found in the submitted numerical arrays.**
All six actual accuracy reports bind their successful metric operands correctly,
including when checked separately for each spectrum. They truthfully report
failed convergence. The r5 handoff's opening claim “R2 is repaired” is therefore
too strong; its R3 qualification is correct. Preserve the submitted handoff as
historical evidence and correct the claim in the next handoff.

## Scientific assessment of the submitted 15x2pt results

Both profiles retain all six bins, with 12 primary results and 180 individually
selected spectrum forecasts. Compatibility passes 6/6 and reproduces the saved
actual lyaforecast reference. Accuracy passes 0/6. All 72 diagnostic requests are
retained: 60 complete, 12 unavailable (floor_low in every bin, remove_bright in
bins 0–4, width_plus in bin 2). The scoped scientific gate correctly rejects this
bundle. The other six real benchmark cases were not rerun in this revision.

The saved joint accuracy uncertainties differ from lyaforecast by -3.31% to
+2.36% in alpha_parallel and -0.61% to +3.32% in alpha_perpendicular. These are
finite-setting comparisons, not converged accuracy forecasts. The handoff's
ordered model/noise/geometry comparisons retain their interpolation residual;
they are not independent additive physical effects or calibrated systematics.

The weight diagnosis is supported by independent saved-array calculations:

- With the initial weight function held fixed, magnitude orders 32 and 64 change
  the individual-spectrum errors by at most **1.09e-14** in our analytic 2x2
  inversion. Thus the tested fixed-weight magnitude integral is well converged.
- For the cumulative rule, with r_i the quadrature-weighted sightline density,
  alpha=L/pixel and a_i=variance_i/S,
  `w_i(next) = alpha*sum(j<=i, r_j*w_j)/(alpha*sum(j<=i, r_j*w_j)+a_i)`.
  Its zero-weight Jacobian is lower triangular with diagonal alpha*r_i/a_i.
  All 60 tested bin/field/order combinations have spectral radius below one;
  the maximum is **0.4072476438698275**. The first-cell positive fixed point is
  absent in these samples. This supports decay of the unnormalized weights and
  shows why quadrature and finite iterations cannot be treated independently.
- Forty matching coupled-weight trial matrices reproduce the saved study results
  with maximum individual pair-Fisher element discrepancy **9.99e-15**. All
  **56** recorded unrepresentable multiplication products reproduce with an
  independent 80-digit Decimal calculation. The bin-0/order-16/12-update example
  is about 8.12e-566, below float64 range; the arithmetic guard is appropriate.

The cumulative refinements still change joint errors by up to 5.017% and
individual-spectrum errors by 38.50% in the saved reports. Increasing the number
of iterations is not an established remedy. Vanishing weights alone do not
exclude a finite limit for A and P_pixel, since these coefficients are invariant
under common weight rescaling. A normalized asymptotic limit, its quadrature
independence and any optimality claim remain unproved. The handoff correctly
leaves three alternatives to the user: a specified finite-iteration convention
with reconsidered acceptance criteria; investigation of a normalized limit of
the same rule; or a separately justified estimator. None is selected by this
review. No floor/support change or tolerance waiver is implied.

The joint comparison (local-only path: `../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt.png`)
and magnitude-convergence figure (local-only path: `../.validation/step12-r5-20260914T191855Z/weight-figures/magnitude-convergence.png`)
were visually inspected. The joint plot explicitly marks all six accuracy bins
unresolved; the magnitude plot distinguishes fixed weights from cumulative
updates and records unavailable orders. The
15 individual-spectrum pages (local-only path: `../.validation/step12-r5-20260914T191855Z/plots/lya_qso_lbg_lae_15x2pt-pairs.pdf`),
attribution (local-only path: `../.validation/step12-r5-20260914T191855Z/attribution/lya_qso_lbg_lae_15x2pt.png`)
and sensitivity figure (local-only path: `../.validation/step12-r5-20260914T191855Z/sensitivities/lya_qso_lbg_lae_15x2pt.png`)
remain linked in the handoff. No plots or real forecasts were regenerated here.

## Fresh validation and preservation

| New review check | Result |
| --- | --- |
| `PATH="$PWD/.venv/bin:$PATH" scripts/check.sh` (script sets the three thread limits) | 1183 passed, 25 optional-compiler skips, 118.94 s; Ruff lint and format pass. |
| Existing r5 installed interpreter, `-I -m pytest -q tests/test_step12_revision5.py tests/test_step12_revision2.py tests/test_step12_performance.py`, copied fixtures under `/tmp/fishhighz-step12-review-r5-tests` | 191 passed, 95.38 s, outside checkout; includes parametrized NumPy/compiled checks. No fresh installation. |
| `numerical_review.py`, installed interpreter, one thread | Reconstructed Wick C from five-field total powers, solved C directly and accumulated joint F plus all individual-spectrum F; 12 records/180 spectra pass. Max C discrepancy 0; joint F 1.34e-15; maximum per-spectrum F 2.20e-15; errors 2.11e-15. Six truthful failed convergence reports pass semantic inspection with `require_pass=False`. |
| `weight_review.py`, installed interpreter, one thread | All six bins, 60 linearizations, 40 coupled-trial comparisons, fixed-weight coefficients/error refinements and 56 Decimal products pass. |
| `check_desi2_plots.py --plots …/plots --reference …/reference-01 --cases lya_qso_lbg_lae_15x2pt` | 240 plotted series, six bins and 12 profiles pass values/reference/ratio/artist checks. |
| `check_desi2_15x2pt.py …/profiles-checked` | Expected exit 1, `partial or scientifically failed validation`; not scientific acceptance. |
| `weak_pair_probe.py` | Reproduces the required R2 correction described above. |

The full commands, logs, review scripts and JSON results are under the separate
review directory. The quick log contains nonfatal site MUNGE messages; the plot
check used a temporary Matplotlib cache after the default cache was unwritable.
Neither caused a failed check. The implementation's disclosed SciPy-blocked,
explicit-Numba subprocess failure is inherited from r4; it was not rerun or
silently counted as a new pass here. Its default-NumPy coverage is included in
the affected installed checks. Historical r5 counts (including 373 installed
checks and six examples) remain implementation evidence, distinct from this
review's 191 installed checks.

Fresh identity checks match all 128 implementation snapshot entries, all 52
live/wheel/installed Python modules and installed METADATA/WHEEL bytes. Exact
wheel SHA256 is
`acb7ae177d2d20bdc270c088c940bb80c6b42fa05c8c26a911e8a438f76e7e69`.
All 35 non-validation modules including `__init__.py` match r4 (the handoff's
34-module count excludes that initializer). All 78 historical r2 primary
NPZ/report hashes, the original r3 instruction hash and the r2 handoff archive
match. No production source or tests were edited during review.

The exact r5 instructions and handoff were snapshotted before review. The active
instructions now contain only the proposed Step 12 revision-6 evidence repair;
design, roadmap and package guidance are synchronized. Implementation dispatch,
scientific weighting decisions, acceptance and progression remain with the user.
