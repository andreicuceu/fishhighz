# W10 revision 1: two saved forest samples

Both assigned samples pass numerical equivalence and the 0.1% bounded
magnitude-refinement criteria. The historical fixed-reference BAO errors are
smaller than the cumulative t=3 and t=6 errors in both samples. No reversal,
numerical obstruction or reference-stability failure occurred. This is an
implementation handoff for independent/user review, not acceptance or W11
progression.

## Evidence and calculation

The single numerical invocation used
scripts/check_inverse_variance_samples.py (local-only path: `../scripts/check_inverse_variance_samples.py`)
and wrote summary.json (local-only path: `../.validation/forest-weight-diagnostics/w10-r1-20260916T024517442610Z/summary.json`).
The summary contains all five verified input SHA-256 identities, decisive source
hashes, complete scalar moments/residuals, refinement operands, eight saved 2x2
matrices, determinants, normalized minimum eigenvalues, errors and contrast
denominators. Source and input identities were unchanged across calculation.
Input root is `.validation/step12-r5-20260914T191855Z/`; no original evidence was
modified. Live HEAD was `0d69786a06d5d676564a51fad14a7156951c2228`, with extensive
pre-existing modified/untracked package work. That state was preserved.

Read workspace/package guidance, design section 0, the main-roadmap handover
and progress register, separate package Step 13 context, the diagnostic roadmap,
W10 r1 assignment and W06/W07/W09 reviews. The documents still call W10 proposed;
the user's explicit approval supersedes that status without a plan edit.

Selection was by report/diagnosis identity, in this order:

| Sample | Bounds | Field / auto-pair indices | Parameters | Native report order | Nodes at 16/32/64 |
| --- | --- | --- | --- | --- | --- |
| A: `lya(lbg)`, bin 0 | [2, 2.235] | 4 / 14 | ap_0, at_0 | 32 | 2592/5184/10368 |
| B: `lya(qso)`, bin 5 | [3.175, 3.410] | 0 / 0 | ap_5, at_5 | 64 | 2640/5280/10560 |

Both are accuracy-profile forest fields with physical label `lya` and the
respective `lbg`/`qso` background. At the actual native order, magnitudes,
quadrature and variance match the report exactly, as do masses reconstructed
from density times `(1+z_source)/299792.458` times quadrature. At other orders,
`rho=masses/measure` reconstructs masses exactly without losing support. All
saved masses are positive; quadrature is positive and sums to each partition
width; variance is finite and nonnegative.

## Held inputs and inherited policies

| Quantity | A | B |
| --- | ---: | ---: |
| z_eval | 2.1152848986890427 | 3.2908915157575356 |
| z_source | 2.3830099582078716 | 3.6596472552176555 |
| L [km/s] | 44147.09373976799 | 44147.09373976799 |
| l_p [km/s] | 63.328211161339375 | 45.9777226171801 |
| B_star [km/s] | 16.356748967852234 | 64.49126195731874 |
| a_v | 102.61326549228777 | 117.92808174108923 |
| d_deg | 64.41180749124948 | 79.67758690215015 |

Both hold h_fid=0.6736, resolution R=2500 with FWHM convention,
sigma_velocity=50.92405402826606 km/s, and 0.8 Angstrom pixels. Saved wavelengths
are 3787.158392789309 and 5216.308088960964 Angstrom. The historical recipe hash
matches each report. W09's constant-background adapter recovers the saved
geometry conversions exactly; its artificial volume is unused. B_star is the
same positive diagnosis alias at all three orders within each sample, already
response-smoothed at k_parallel=0.00035 s/km. No response is applied to it again.

Inherited policies are `floor_negative` density with floor 1e-20 (explicit
reference extension), `legacy_floor` magnitude support,
`legacy_spline_extension` redshift support, `legacy_floor_clamp` SNR with floor
1e-10 and sentinel variance 1e20, four exposures, and
`legacy_first_spacing; unknown physical cells`. Density normalization measures
are 380.0 and 3789.61349941401, using the saved legacy reduction convention.
These assumptions are retained, not physically validated. Only saved magnitude
quadrature varies within each sample.

## Newly prepared coefficients

Each order calls public `method="inverse_variance", alias=B_star` once, with
no signal, iterations, auxiliary sampling or model evaluation. Independent
Python scalar products and `math.fsum` compute I1, I2, I3 directly from saved
masses and nu=B_star/(B_star+l_p*v). They give A=I2/(L*I1^2),
P_pixel=l_p*I3/(L*I1^2), and Q_star=A*B_star+P_pixel=B_star/(L*I1).
A is in deg^2; P_pixel and Q_star are in deg^2 km/s.

| Sample / order | Public A | Public P_pixel | Public Q_star |
| --- | ---: | ---: | ---: |
| A / 16 | 0.15489596981717563 | 42.0349589673603 | 44.568553461791865 |
| A / 32 | 0.15489596981717518 | 42.034958967360154 | 44.56855346179171 |
| A / 64 | 0.15489596981717602 | 42.034958967360446 | 44.56855346179201 |
| B / 16 | 0.2783040496036248 | 2.560857491877908 | 20.509036858647903 |
| B / 32 | 0.2783040496036282 | 2.5608574918779383 | 20.509036858648155 |
| B / 64 | 0.278304049603625 | 2.560857491877902 | 20.509036858647907 |

Largest public/scalar relative residuals are 4.74e-15 (A) and 8.16e-15 (B).
Moment residuals are below 3.43e-15; historical fixed-coefficient residuals
are below 8.46e-16. The independent Q identity agrees within 3.47e-16.
All satisfy rtol=5e-12, atol=0, with exact handling of expected zeros.

Signed fractional refinement changes, in order (A, P_pixel, Q_star):

| Sample | 32/16 minus 1 | 64/32 minus 1 | 64/16 minus 1 |
| --- | --- | --- | --- |
| A | (-2.89, -3.55, -3.55)e-15 | (+5.33, +6.88, +6.66)e-15 | (+2.44, +3.33, +3.11)e-15 |
| B | (+1.22, +1.18, +1.22)e-14 | (-1.15, -1.42, -1.21)e-14 | (+0.666, -2.44, +0.222)e-15 |

Every absolute change is below 0.001. Decimal division of the full binary
operands independently checks contrasts, with absolute residuals below 8.67e-17.
This is bounded stability, not a continuum theorem.

## Historical Fisher matrices, newly inverted

Only after each sample's coefficient checks passed were its unique selected
auto-spectrum matrices read. These are historical `pair_fisher` matrices for an
individual observable, not blocks of a joint inverse and not new mode sums.
Analytic determinants and marginalized errors sqrt(d/det), sqrt(a/det) agree
with ordinary 2x2 solves and saved `pair_errors` within 2.91e-16 relative.
All matrices are exactly symmetric, positive definite and adequately conditioned
for these checks. No prior, regularization or pseudoinverse enters.

| Sample / order / weights | sigma_parallel | sigma_perp |
| --- | ---: | ---: |
| A / 32 / reference | 0.7288907443946463 | 0.9889518040663758 |
| A / 64 / reference | 0.7288907443946513 | 0.9889518040663827 |
| A / 32 / t=3 | 1.5733923238865393 | 2.1510891378372503 |
| A / 32 / t=6 | 2.55166452158379 | 3.497722105304831 |
| B / 32 / reference | 0.12640403641329562 | 0.17544181456742458 |
| B / 64 / reference | 0.12640403641329434 | 0.17544181456742267 |
| B / 32 / t=3 | 0.13413019378759533 | 0.18725382396522464 |
| B / 32 / t=6 | 0.1366919681615717 | 0.19103589602370208 |

Reference 64/32 fractional error changes are (+6.8834e-15, +6.8834e-15)
for A and (-1.0103e-14, -1.0880e-14) for B. Both pass 0.1% before cumulative
comparisons proceed. Each percentage below uses its own **order-32 reference
error** in the denominator: `100*(sigma_cumulative/sigma_reference-1)`.

| Sample / numerator | Radial contrast | Transverse contrast |
| --- | ---: | ---: |
| A / t=3 | +115.86120224% | +117.51202930% |
| A / t=6 | +250.07503404% | +253.67973352% |
| B / t=3 | +6.11227109% | +6.73272186% |
| B / t=6 | +8.13892660% | +8.88846339% |

Independent Decimal divisions agree within 2.23e-16 in fractional units,
well below 5e-12. No contrast sign or magnitude was a pass criterion.

## Execution and limitations

From the package root, the sole numerical invocation was:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 \
  .venv/bin/python -B scripts/check_inverse_variance_samples.py
```

It performed exactly six preparations and eight saved-matrix inversions in one
process, with all three thread limits equal to one. Internal elapsed time was
0.920416214 s; shell-tool wall time was 3.72 s, below the 60-second cap. No early
stopping was required. Ruff lint, Ruff format checking on the new script, and
`git diff --check` passed. No broad tests ran.

After calculation, one output-only edit omitted the array-valued `scaled_snr`
policy diagnostic from serialization. The numerical calculations were not
rerun. The 70 kB summary records both executed and final script hashes and the
exact text replacement that reconstructs the executed source; that reconstruction
was hash-verified. Numerical statements/results and scalar policies are unchanged.

Relative to W06–W09, this extends the stable reference and historical error
ranking to one LBG-background sample and one higher-redshift QSO-background
sample. It does not repeat W09's public mode sum for these samples. The large
LBG Fisher errors are local Gaussian curvature estimates, not evidence that a
Gaussian likelihood accurately describes such weak BAO constraints. These
samples establish no general multi-mode optimum, population/redshift survey,
physical policy validation or acceptance of the failed Step 12 forecasts.

The historical bin-0 order-64/t=6 diagnostic remains unavailable and was neither
requested nor regenerated; no individual-spectrum failure was inferred from it.
No new modes, derivatives, P3D/P1D evaluations, production fixes, profile changes,
plan/manuscript edits, installations, forecasts, Slurm actions, agent dispatch,
commits or pushes occurred. Only the new script, this handoff and fresh small
outputs were written. Stop for review and acceptance; W11 remains unstarted.
