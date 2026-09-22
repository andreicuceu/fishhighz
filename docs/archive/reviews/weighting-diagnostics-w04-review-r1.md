# Scientific conclusion and impact of proposed changes

**W04 revision 1 passes independent scientific review within its assigned
finite-count scope. No scientifically consequential changes are needed.**
Additional iterations alone reduce the fixed-intrinsic full-sample drift:
both branches satisfy the prescribed 0.1% coefficient-stability criterion over
t=6/12/24. Explicitly refreshed aliasing is therefore not necessary for this
bounded stability. It reduces the measured iteration drift and changes the
coefficient levels: at t=24, A is 5.08747% lower and P_pixel is 12.60455% higher
than with fixed P0. The implementation report states these conclusions and
their limits correctly. No repair revision is required.

This is one fixed signed magnitude grid, not evidence of a continuum limit,
asymptotic convergence, physical admissibility, optimality, or a production
prescription. The experiment does not resolve the historical signal-refresh
convention. **Awaiting user review and acceptance; this pass does not authorize
progression.** Package Step 13 and its acceptance status remain unchanged.

Reviewed on 2026-09-15: [W04 handoff](weighting-diagnostics-w04-r1.md),
the [exact r1 assignment](weighting-diagnostics-w04-instructions-r1.md), and
`scripts/compare_forest_weight_signals.py`. No substantive feedback was supplied;
the request retained its feedback placeholder.

## Controlled scientific comparison

I read workspace/package AGENTS.md, design section 0, the main handover/progress
register, diagnostic roadmap, W03 report/review, and W02's signal assessment.
The relevant W03 initialization, full-sample update, moment products and scalar
calculation were inspected alongside W04. W02's historical source assessment
is retained as context; no historical algorithm or literature search was reopened.

W04 starts two separate arrays from the same seed (B/lp)/(B/lp+v). Both use
the full-sample I1, with N=(L/lp)I1 common to all magnitude nodes. W04 retains
W03's float64 moment product order and `cumsum(...)[-1]` reduction. The only
branch-dependent operation is S=P0 versus S=P0+A[w]B. Its `update` function
computes moments from the current branch's weights before simultaneously
updating them. There is no borrowed, frozen, prefix, or post-update aliasing
coefficient, and no inter-update normalization. Both branches extract
A=I2/(L I1²) and P_pixel=lp I3/(L I1²) from their own weights.

The authoritative saved report has SHA-256
`3e5200d3e5b1bab3e7570c61fc02651caa907ae30816036b0b1b323098a87ed6`:
`.validation/step12-r5-20260914T191855Z/profiles-checked/records-000.report.json`.
I independently checked compatibility bin 0, z=[2.0,2.235], and selected
`settings.pair_inputs["lya(qso)_lya(qso)"]`. All 107 magnitudes, signed densities
and variances match W03 and W04 exactly and are finite. Nine densities are
negative; none is floored. The seed and all held scalars also match W03.

| Held quantity | Value |
| --- | ---: |
| Magnitude endpoints | 16.1, 26.75 mag |
| Rectangular dm | 0.100471698113207 mag |
| L | 44148.20436604321 km/s |
| lp | 63.32980433473595 km/s |
| P0 | 1.7189455896517238 deg² km/s |
| B | 16.333147249645222 km/s |

I1/I2/I3 have units deg⁻² (km/s)⁻¹ and A has units deg², so A B, P0,
v/N and P_pixel all have units deg² km/s. W02's inspected response ownership
places instrumental smoothing in the supplied P0/B. W04 uses those exact saved
scalars without another transfer factor, conversion or physical-input evaluation.
P_pixel is absent from S; coefficient extraction adds no second aliasing term.

## Independent scalar evidence

The review-only calculation uses 80-digit Decimal arithmetic with exact
binary-float conversion. It forms direct scalar sums and evolves each branch
independently through 24 updates. It uses the algebraically equivalent form
w_new=S N/(S N+v), imports neither implementation's real-input recurrence,
and never calls their moment or Decimal helpers. It checks the seed, all eight
checkpoint weight vectors, output-only shapes, moments and coefficients, all
48 pre-update signals, matched-count contrasts, and each stability classification.
The maximum relative discrepancy across these comparisons, including the
fractional contrasts, is **2.360e-14**, below rtol=5e-12 with atol=0.
Finite operands and exact reference-zero handling are explicit. No precision
increase was needed.

Fixed-P0 t=3/6 weights and all five moment/coefficient entries independently
reproduce the saved W03 full-sample arrays. Their coefficients reproduce the
assignment's values, approximately (0.0220213796793063, 0.529854638024920)
and (0.0218188652536176, 0.540542582564187). W03 supplied no t=12/24 result;
those are W04 measurements.

Independently computed matched-count fractional shifts, refreshed minus fixed:

| t | ΔA/A_fixed | ΔP_pixel/P_pixel,fixed |
| ---: | ---: | ---: |
| 3 | −5.398858% | +13.264440% |
| 6 | −5.104704% | +12.646487% |
| 12 | −5.087495% | +12.604612% |
| 24 | −5.087471% | +12.604552% |

The independent maximum absolute fractional change over both coefficients and
all three comparisons in each checkpoint triple is:

| Branch | (3,6,12) | Classification | (6,12,24) | Classification |
| --- | ---: | --- | ---: | --- |
| Fixed P0 | 0.0208308023 | Not stable at 1e-3 | 0.000646928289 | Stable at 1e-3 |
| Refreshed aliasing | 0.0148839004 | Not stable at 1e-3 | 0.000274414011 | Stable at 1e-3 |

Thus continuing at t=12 and stopping at t=24 is consistent with the assignment.
The saved attempt records contain exactly updates 0→1 through 23→24 for each
branch. Every independent state has finite positive moments, coefficients,
weights, N and S; no arithmetic or diagnostic-domain obstruction occurs.
This positivity does not turn the signed input density into a physical measure.
The independently checked refreshed A B/P0 is about 0.196714 at t=24;
that final-state ratio does not imply an update beyond the cap.

For the one-bin analytic control, A=1 and N=w give fixed f(w)=w/(1+w)
and refreshed f(w)=2w/(1+2w). Starting at 1/2 yields 1/3 and 1/2 respectively.
For two bins, exact rational arithmetic gives A0=29/49, S0=78/49 and
w1=(39/74,39/179); recomputation gives A1=37517/64009,
S1=101526/64009 and w2=(1979757/3655376,1979757/8682233).
These independently derived values agree with the live float `update` function.
Its second signal differs from the first, verifying the refresh timing.

## Reproducibility and disposition

The original W04 artifacts remain in
`.validation/forest-weight-diagnostics/w04-r1-20260915T225726Z/`.
All seven source/input hashes recorded in its summary matched live files
before status synchronization. The exact assignment is archived byte for byte
with SHA-256 `2dd75313f84c404a159bb16872b567d3fece57207583698dd61ae0516c55eb02`.
The reviewed W04 script hash is
`1adfa3354fc887d5008c3a37596b494ba8585c774a9866aaa366bdf8a462b788`.

New review evidence is in
`.validation/forest-weight-diagnostics/w04-review-r1-20260915T230147Z/`:
`scalar_check.py` and `result.json`. From the package root:
The check now reads the exact instruction archive for the recorded assignment
hash, preserving reproducibility after the current file's status update; that
path substitution was checked against the identical archived bytes and does not
change the numerical calculation.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 30 .venv/bin/python -B .validation/forest-weight-diagnostics/w04-review-r1-20260915T230147Z/scalar_check.py
.venv/bin/ruff check scripts/compare_forest_weight_signals.py .validation/forest-weight-diagnostics/w04-review-r1-20260915T230147Z/scalar_check.py
.venv/bin/ruff format --check scripts/compare_forest_weight_signals.py
```

The single numerical invocation passed in 0.1573 s internally (0.638 s command
wall time), on login32, with one process, numerical thread limits of one, and
a 30-second cap. Ruff passed. Only the review evidence, this report, exact
instruction archive and diagnostic-status text were added or changed. The
scientific assignment is unchanged; no repair or subsequent step is prepared.
All production code, historical evidence and unrelated dirty/untracked work
remain preserved. Package `IMPLEMENTATION_STEP.md` remains at SHA-256
`3675cd861fc7909d464d9f7e60f7725f1a9e2f5cb3fa8a9f9849f3a30552c310`.

No broad regression suite, W01 hardening, general validator, prefix extension,
other population/bin, magnitude refinement, new physical inputs, forecast,
Slurm action, agent dispatch, commit or push was performed. Stop for user review.
