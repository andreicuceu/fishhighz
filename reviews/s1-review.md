# Independent scientific review: S1 profile definitions and baseline weights

## Scientific verdict

**PASS.** I find no scientifically consequential correction required before the
coordinator proceeds to the authorized S2 calculation. The implementation preserves
the accepted early-lyaforecast recurrence and fixed-compatibility evidence, gives
the revised accuracy profile its own fiducial per-field inputs, and applies the new
bin-1 observable selection before covariance construction and inversion.

This verdict establishes implementation consistency and saved-input reproduction.
It does not qualify the revised accuracy calculation numerically, accept a new
forecast scientifically, or establish magnitude-grid convergence. Those are S2/S3
outputs and remain subject to explicit numerical status and user review.

## Independent scientific checks

I independently evaluated the current frozen source with one numerical thread.
The review script `.validation/s1-review/check_science.py` passed and verifies:

- new full-compatibility, fixed-compatibility and revised-accuracy requests contain
  exactly lya(qso) auto, QSO auto and their cross in bin 1, while bins 2--6 retain
  all 15 spectra; a historical accuracy request remains readable with its original
  15-spectrum bin 1;
- direct two-cell early and McDonald updates agree with independently written scalar
  equations; rescaling the weights changes the nonlinear update while leaving final
  `A` and `P_pixel` invariant;
- the strict public preparation rejects signed density, while the literal signed
  compatibility calculation remains confined to the validation path;
- the adaptive solver tests amplitude, signed normalized shape, full vector, `A`
  and `P_pixel`, and returns the actually computed doubled-count state;
- the selected bin-1 covariance has shape 3x3 at each Fourier cell and its Fisher
  matrix agrees with an independent direct NumPy solve; it differs from both the
  old 15-spectrum joint result and the sum of the three individual Fishers;
- full-compatibility retains the selected legacy powers, Jacobians and individual
  Fishers exactly; fixed-compatibility changes only the required QSO-forest auto
  total-power column and its three individual Fishers agree with the accepted
  Stage-4 `sum_historical`, `rtol=1e-4` evidence.

The six-bin saved-input replay also passed. Early weights converge in every required
forest auto, including the selected QSO forest in bin 1 at update 18. The known
McDonald LBG-forest cap in user bin 2 remains explicit. Across all bins, independent
selected-covariance solves agree with the saved joint Fishers to at most
`8.57e-15` relative. Bins 2--6 reproduce the accepted fixed-compatibility joint
Fishers bitwise; the new bin-1 fixed-compatibility errors are
`(0.0251434661515, 0.0183510893061)`.

## Physical ownership and numerical status

The accuracy path queries each active forest auto from its own fiducial P3D provider
and independent P1D model at `(k_t_deg, k_p_velocity) = (2.4, 0.00035)`. The squared
field response enters representative P and B once. The final noise remains
`A P1D_resp + P_pixel`, with the angular/velocity-to-comoving conversion applied
once and the pixel term unsmoothed. QSO and LBG forests retain separate density,
noise, response and weights; their prepared states are reused across every relevant
auto/cross covariance and remain frozen through BAO finite differences. Focused
unit, response, field-reuse and derivative tests passed, including an explicit
angular/comoving equivalence check.

The evidence identities distinguish recipe revision, method, representative mode,
quadrature and stopping controls. Historical trial contracts 1 and 2, including
the W12 fixed-reference option, remain intact; adaptive trials use contract 3.
Passing full/fixed records are labelled reproducibility controls, passing accuracy
records require finite-refinement qualification, and failures remain unresolved.
The comparison tables retain 78 selected individual results and six selected joint
results per profile, plus 12 explicit bin-1 exclusion records; plotting cuts do not
remove values from the tables or joint contractions. Four comparison-output tests,
including rendering, passed in the coordinator environment.

Final focused checks were `10 passed` for `tests/test_full_sum_weights.py`; the
implementation handoff records the broader 87-test weighting/API set, 54-test
W12/trial/historical set and 98-test forecast/W12 regression set passing. Ruff,
formatting, `git diff --check`, the saved-input replay, wheel build and installed
61-module wheel identity check also passed. I ran no real accuracy forecast, full
suite, Slurm action, commit or push.

S1 is ready for the authorized S2 numerical qualification. Scientific acceptance
remains with the user.
