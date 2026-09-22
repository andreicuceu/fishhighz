# Forest-weighting baseline: early lyaforecast weights

Decision recorded: 2026-09-18.

## Decision and scope

**Early lyaforecast weights are the chosen FishHighz baseline. McDonald weights
will remain available as an explicit alternative.** Both prescriptions use
full-sample scalar magnitude integrals and include an aliasing contribution in
the signal used to determine the weights. They differ in that contribution's
dependence on the weight moments. Both use the same moment-based aliasing and
pixel-noise expressions in the final covariance.

The choice follows the five-stage compatibility weighting study and the
additional magnitude-grid-to-BAO calculation performed during the assessment in
Step 5. Early lyaforecast weights achieved finite, nonzero iteration convergence
for both forest populations in all six tested redshift bins on all three
magnitude grids. McDonald weights did not meet this criterion in the two lowest
redshift LBG-forest bins. For early lyaforecast weights, the final grid refinement
changed every individual BAO uncertainty in the subsequently tested bins 2–6 by
less than 0.1%, and joint uncertainties by at most 0.00562%.

This is the user's accepted scientific choice, not a claim that the production
default has already changed. The tested variants and adaptive stopping currently
reside in the validation modules. This report records the intended baseline,
the retained alternative, and the evidence for the decision; it does not change
code or complete the remaining accuracy-profile reassessment and package plan.
Older status statements in the planning documents and `AGENTS.md` predate the
Step-4 calculation and this decision. The historical Step 13 and W12 assignments
remain closed.

Throughout this report, bins are numbered **1–6**, with mean redshifts
2.1175, 2.3525, 2.5875, 2.8225, 3.0575 and 3.2925. Saved files and some earlier
reports instead use internal indices 0–5.

## 1. Shared definitions and final covariance

Consider one observed forest population at a fixed redshift, with a selected
magnitude range. Let

- $n(m)$ be its source density per magnitude;
- $q_i$ be the magnitude quadrature weights and $r_i=n(m_i)q_i$;
- $L$ be the common forest length and $l_p$ the pixel width;
- $v_i=\sigma_{N,i}^2$ be the dimensionless pixel-noise variance;
- $P$ be the intrinsic 3D forest power at the representative weighting mode;
- $B$ be the 1D forest power at that mode, with consistent instrumental response.

In comoving units, $P$ has volume units, while $B$ and $l_pv_i$ have
length units. The same equations apply to the saved angular/velocity quantities
after their explicit coordinate conversions. Define the full-sample moments

$$
I_1=\sum_i r_iw_i,\qquad
I_2=\sum_i r_iw_i^2,\qquad
I_3=\sum_i r_iv_iw_i^2.
$$

Each is a single scalar for the selected population. Every source in that
population uses the same $I_1$ during an update. A magnitude-prefix array
$I_{1,i}=\sum_{j\leq i}r_jw_j$ would define a different recurrence and is not
the chosen normalization. Using `cumsum(...)[-1]` to evaluate a scalar full sum
is acceptable; the diagnostic implementation uses this reduction to preserve
the compatibility arithmetic order.

The common initial weights are

$$
w_i^{(0)}=\frac{B}{B+l_pv_i}.
$$

For either prescription, the noise assigned to a source in the update is

$$
P_{N,i}^{(t)}=\frac{l_pv_i}{L I_1^{(t)}},\qquad
w_i^{(t+1)}=\frac{S^{(t)}}{S^{(t)}+P_{N,i}^{(t)}}.
$$

All weights are updated simultaneously from the previous iterate. Three
iterations means three updates after the seed. The source-dependent
$P_{N,i}$ above is distinct from the population-averaged final pixel power.

After choosing the weights, **both prescriptions** use

$$
A=\frac{I_2}{L I_1^2},\qquad
P_{\rm pixel}=\frac{l_pI_3}{L I_1^2},
$$

and, in a common coordinate system,

$$
P_{\rm obs}(\mathbf{k})=
P_{F,\rm resp}^{3D}(\mathbf{k})
+A P_{F,\rm resp}^{1D}(k_\parallel)+P_{\rm pixel}.
$$

Here “resp” indicates that the appropriate instrumental response has already
been applied. The 1D contribution receives the squared amplitude response once;
the additive pixel-noise term is not smoothed again. The scalar $B$ used to
choose the weights is evaluated at the representative mode, whereas the final
1D power is evaluated across the forecast modes.

Thus choosing early lyaforecast weights does **not** replace the final aliasing
coefficient by $1/(L I_1)$. The approximation concerns the weight update only.
Common rescaling of the final weights cancels in $A$ and $P_{\rm pixel}$,
but generally changes the next nonlinear update. Arbitrary renormalization
during iteration is therefore not part of either tested prescription.

## 2. Chosen baseline: early lyaforecast weights

The baseline uses the early Python expression

$$
S_{\rm early}^{(t)}=P+\frac{B}{L I_1^{(t)}}.
$$

Equivalently, its recurrence is

$$
\boxed{
w_i^{(t+1)}=
\frac{P L I_1^{(t)}+B}
     {P L I_1^{(t)}+B+l_pv_i}.
}
$$

Only $I_1$ enters the update; $I_2$ and $I_3$ are still required for
the final covariance and convergence diagnostics. The validation identifier is
`sum_historical`.

The [early Python source](https://github.com/igmhub/lyaforecast/blob/8f2ddf1787f63eb17f1f48a59456b7ddcd716520/py/forecast.py)
also expresses this in terms of an effective pixel density and a variance-like
signal:

$$
n_{p,\rm eff}=\frac{L I_1}{l_p},\qquad
V_S=P n_{p,\rm eff}+\frac{B}{l_p},\qquad
w_i=\frac{V_S}{V_S+v_i}.
$$

These are algebraically the same update. The $B/l_p$ term represents intrinsic
forest fluctuations at the chosen mode in this approximation; it is not the
pixel variance obtained by integrating the entire 1D power spectrum. Historical
comments associate the alternative formulation with McDonald's C++ forecasts,
but that C++ implementation has not been independently verified here.

This prescription uses $L I_1$ as an effective sightline density for its
weighting signal. When the intrinsic 3D term is small, it tends towards
$B/(B+l_pv_i)$: noisy spectra are weighted according to instrumental noise
relative to intrinsic forest fluctuations. In particular, for positive $B$,
the numerator retains a positive contribution when $I_1$ becomes small.
This helps explain the numerical robustness found in the low-redshift
LBG-forest cases. It is a motivation and a limiting-case interpretation, not a
proof of convergence for arbitrary inputs or of optimal multi-mode BAO weights.

## 3. Retained alternative: McDonald weights

The optional prescription uses the current moments in the representative
aliasing contribution:

$$
S_{\rm McD}^{(t)}=P+
\frac{B I_2^{(t)}}{L[I_1^{(t)}]^2},
$$

or

$$
\boxed{
w_i^{(t+1)}=
\frac{P L I_1^{(t)}+B I_2^{(t)}/I_1^{(t)}}
     {P L I_1^{(t)}+B I_2^{(t)}/I_1^{(t)}+l_pv_i}.
}
$$

The validation identifier is `sum_aliasing`. Both moments and the aliasing
contribution are recomputed from the previous iterate at every update.

This is the full-moment, aliasing-inclusive interpretation of
[McDonald & Eisenstein (2007), section II.B](https://arxiv.org/pdf/astro-ph/0607122).
Their equations 14–19 and accompanying discussion give full-sample moments,
the FKP-like weight, and an aliasing-inclusive representative signal, using
$k=0.07\,h\,\mathrm{Mpc}^{-1}$ and $\mu=0.5$.
[Font-Ribera et al. (2014), equations 23–29](https://arxiv.org/pdf/1308.4164)
give the corresponding moment and weighting expressions. The papers do not
uniquely specify the initialization, refresh schedule or stopping rule. Here,
“McDonald weights” names the explicit recurrence above, rather than asserting
exact reproduction of a historical C++ algorithm or a globally optimal estimator.

For nonnegative density and $0\leq w_i\leq1$, $I_2/I_1\leq1$.
At the same input weights, early lyaforecast therefore uses a larger or equal
aliasing contribution in the update and relatively less strongly suppresses
noisy spectra. These statements need the stated positivity assumptions; some
compatibility interpolants have signed density and do not inherit this argument.

A homogeneous positive-density example clarifies the differing fixed points.
For common weight $w$, total source density $n$, noise $N=l_pv>0$, and
$s=P L n$,

$$
w_{\rm McD}^{(t+1)}=
\frac{(s+B)w^{(t)}}{(s+B)w^{(t)}+N},\qquad
w_{\rm early}^{(t+1)}=
\frac{s w^{(t)}+B}{s w^{(t)}+B+N}.
$$

The McDonald example has a positive fixed point only if $s+B>N$; otherwise
positive iterates approach zero. The early update tends to $B/(B+N)>0$ as
its input weight approaches zero. This example explains a possible amplitude
collapse despite stable normalized covariance coefficients. It does not prove
the asymptotic behavior of the measured magnitude-dependent populations.

## 4. Numerical prescription accompanying the choice

Use converged weights for the baseline, with the tested adaptive prescription
as the starting point for production integration:

1. Start from the common seed, using full-sample scalar moments.
2. Require at least three updates and three consecutive stable transitions.
   Test relative weight amplitude, normalized signed shape, the complete weight
   vector, $A$, and $P_{\rm pixel}$, each against its own scale, with
   `rtol=1e-4` and no absolute floor that could mistake decay for convergence.
3. For a candidate at update $t$, require stability through an actual doubled
   count $2t$, including every intervening transition and every comparison
   with that candidate. Return the state reached at $2t$.
4. Cap the calculation at 96 updates. Report convergence, cap, arithmetic
   failure or ineligibility explicitly; never label the last attempted state
   as converged or silently substitute another weighting prescription.

The validation API currently permits adaptive stopping only for full-sum,
positive-signal forest-auto preparations. Weights belong to an observed forest
field and must be used consistently in every covariance involving that field.
The extra historical pair-specific preparations are diagnostic contexts, not
independent physical weights for each cross-spectrum.

The 425-node magnitude grid is a supported working resolution for the tested
endpoints, interpolants and survey configurations in bins 2–6. Resolution should
remain configurable: the evidence does not establish a universal node count or
a continuum error bound. The original rectangular quadrature, including full
endpoint weights, was retained in these refinements. No accuracy-profile
quadrature or density regularization was introduced.

The McDonald option should use the same explicit convergence reporting and
retain the known unavailable cases. Fixed three-update calculations remain
useful for reproducibility and comparisons, with their iteration count stated.
Neither fixed-count results nor a switch to early weights should be a hidden
fallback for a failed McDonald solve. Public option names and production
integration remain implementation work; the identifiers above describe the
existing validation calculations.

## 5. Evidence from the five steps

### Step 1: isolate the five variants and reproduce compatibility

The study implemented the following controlled variations, with the same seed,
input arrays, arithmetic conventions and final covariance expressions:

| Validation identifier | Moments in update | Representative signal |
| --- | --- | --- |
| `prefix_intrinsic` | Magnitude prefixes | $P$; literal compatibility control |
| `sum_intrinsic` | Full-sample sums | $P$ |
| `prefix_aliasing` | Magnitude prefixes | $P+B I_{2,i}/(L I_{1,i}^2)$ |
| `sum_aliasing` | Full-sample sums | $P+B I_2/(L I_1^2)$; McDonald |
| `sum_historical` | Full-sample sum | $P+B/(L I_1)$; early lyaforecast |

There were 54 saved preparations: nine forest-related pair contexts in each
of six bins. Only 12 forest-auto preparations supply additive forest noise in
the saved covariance. The remaining 42 were retained to expose the behavior of
the historical pair preparation, including signed auxiliary quantities.

All 54 three-update control weight vectors and coefficient pairs reproduced
exactly. Individual Fisher matrices and errors also reproduced exactly; maximum
relative joint Fisher and marginalized-error discrepancies were respectively
$4.67\times10^{-17}$ and $9.11\times10^{-17}$.
Independent analytic tests checked the response and coordinate conversion of
the forest noise, so agreement was not based solely on subtracting and adding
the same cached term. Eighteen focused tests passed.

The calculation retained the existing distinction between the arithmetic mean
redshift used in the mean/Jacobian and the geometric covariance redshift. It did
not adjust intrinsic spectra to force those two quantities to coincide.

### Step 2: iteration convergence and magnitude-grid sensitivity

The study attempted 810 trajectories: 54 contexts, five variants and three
endpoint-inclusive magnitude grids of 107, 213 and 425 nodes. Every finite
iterate was saved through update 96, or until arithmetic failure. The raw
density/SNR reconstruction reproduced all 54 original input preparations
bitwise. Refinements retained the original interpolation and input policies.

At relative tolerance $10^{-4}$, confirmed finite, nonzero forest-auto
convergence was:

| Variant | 107 nodes, out of 12 | 213 nodes, out of 12 | 425 nodes, out of 12 | Confirmation updates where converged |
| --- | ---: | ---: | ---: | --- |
| `prefix_intrinsic` | 0 | 0 | 0 | None |
| `sum_intrinsic` | 10 | 10 | 10 | 14–40 |
| `prefix_aliasing` | 11 | 11 | 9 | Grid-dependent; up to 94 |
| McDonald | 10 | 10 | 10 | 12–34 |
| Early lyaforecast | 12 | 12 | 12 | 10–30 |

Convergence required amplitude and shape stability as well as coefficient
stability, doubled-count confirmation, and agreement with every later saved
state through update 96. This last check rejected transient plateaus in the
prefix-aliasing trajectories; an early doubled-count comparison in one mixed
context passed before a later drift of approximately 82%.

The intrinsic-prefix trajectories encountered overflow or zero-prefix divisions
at updates 15–17. This is a failure of the literal floating-point recurrence,
not a mathematical proof that every normalized quantity diverges. Full-sum
intrinsic and McDonald weights failed the finite-nonzero convergence criterion
for LBG forests in bins 1 and 2. Early lyaforecast passed for all forest autos,
but did not converge in every one of the additional signed pair contexts.

The two McDonald failures have different character. On the original grid, bin 1
reached a maximum weight amplitude of only $6.758\times10^{-92}$ at update 96,
while normalized quantities were stable; it was classified as amplitude decay
with stable normalized quantities. Bin 2 retained amplitude 0.866836 at update
96 but was still classified as continued decay, without a confirmed stable
candidate. The latter is an unresolved result within the tested cap, not proof
that a nonzero fixed point cannot exist at a larger iteration count.

For converged autos, the largest absolute coefficient changes under the final
213-to-425-node refinement were approximately:

| Prescription | Aliasing coefficient $A$ | Pixel power $P_{\rm pixel}$ | Population entering comparison |
| --- | ---: | ---: | --- |
| McDonald | 0.00289% | 0.00966% | Ten converged autos; excludes LBG bins 1–2 |
| Early lyaforecast | 0.194% | 0.117% | All twelve autos |

The early-scheme maxima came from low-redshift LBG forests: bin 1 for $A$,
bin 2 for pixel power. The two rows therefore compare different sets of
converged cases. Small coefficient changes for McDonald do not resolve its
missing low-redshift solutions, and coefficient changes alone do not determine
the BAO uncertainty. Twenty-one focused tests passed.

### Step 3: validate adaptive stopping against the saved trajectories

The adaptive prescription was tested in 216 eligible combinations of forest
auto, full-sum variant, magnitude grid and tolerance. There were 192 converged
results and 24 caps. The caps were precisely the two lowest-redshift LBG autos
for intrinsic-only and McDonald sums, on every grid at both tolerances.

All returned weights and coefficients matched their saved trajectory states
bitwise. Every later state through update 96 agreed with each returned converged
state within the specified tolerance. The largest stopped-to-later discrepancy
across the five monitored quantities at `rtol=1e-4` was
$2.00\times10^{-9}$ for early lyaforecast and $6.46\times10^{-9}$ for
McDonald, among their converged cases.

The checks also covered 1,404 ineligible combinations and verified that all 810
exactly-three-update results were unchanged. Thirty focused tests passed.
These results support the bounded stopping rule for the tested autos, rather
than an asymptotic theorem for arbitrary signed measures.

### Step 4: propagate the weighting choices into 15×2pt BAO forecasts

The five observed fields were Lyα through QSO spectra, Lyα through LBG spectra,
QSO positions, LBG positions and LAE positions. All 15 auto/cross spectra were
forecast individually and jointly in all six bins, on the original 107-node
magnitude grid. The calculation varied only the two forest-auto noise
contributions. Intrinsic powers, cross powers, galaxy noise, responses,
Jacobians, Fourier modes, volumes, redshift conventions, parameter definitions
and input policies remained fixed. The two BAO parameters had no additional
priors, and weights were held fixed when differentiating the mean.

Individual constraints used each spectrum's own variance. Joint constraints
used the full inter-spectrum covariance, not the sum of the individual Fisher
matrices. All five variants supplied 90 individual and six joint forecasts at
three updates. At convergence with `rtol=1e-4`, availability was:

| Variant | Individual results / 90 | Joint results / 6 | Missing converged forest autos |
| --- | ---: | ---: | --- |
| `prefix_intrinsic` | 36 | 0 | Both populations in every bin |
| `sum_intrinsic` | 80 | 4 | LBG bins 1–2 |
| `prefix_aliasing` | 85 | 5 | LBG bin 3 |
| McDonald | 80 | 4 | LBG bins 1–2 |
| Early lyaforecast | 90 | 6 | None |

The 36 prefix-intrinsic entries are galaxy-only constraints, which do not depend
on forest weights. All galaxy-only individual Fisher matrices and uncertainties
remained bitwise unchanged across variants and stopping choices.

Converged joint uncertainties differed from the three-update compatibility
control by approximately **−7.15% to +0.896%**, pooling available variants, bins
and components. Within a given variant, changing from three updates to
converged weights changed joint uncertainties by approximately **−0.0371% to
+0.945%**. The individual LBG-forest auto was more sensitive: in bin 3, its
transverse uncertainty increased by **6.318% for McDonald** and **4.539% for
early lyaforecast** relative to the respective three-update results. Thus a
fixed count of three can miss a scientifically relevant individual-spectrum
change even when the joint change is modest.

Tightening the stopping tolerance from $10^{-3}$ to $10^{-4}$ changed the
largest available BAO uncertainty by only **0.0000866922%**. This comparison
covered 780 components, including valid points omitted from plots, and passed
the assigned 0.1% criterion. It establishes iteration stability at fixed grid.

Some explicit LBG-forest auto uncertainties illustrate the distinction between
nonconvergence and a weak but finite forecast. These are the original-grid
results; each entry is $(\sigma_{\alpha_\parallel},\sigma_{\alpha_\perp})$.

| Bin | Mean redshift | Early lyaforecast, converged | McDonald, converged |
| --- | ---: | --- | --- |
| 1 | 2.1175 | (0.75859, 1.04154) | Unavailable |
| 2 | 2.3525 | (0.058370, 0.067444) | Unavailable |
| 3 | 2.5875 | (0.031229, 0.031662) | (0.031153, 0.031560) |
| 4 | 2.8225 | (0.028055, 0.028252) | (0.028010, 0.028191) |
| 5 | 3.0575 | (0.047482, 0.056153) | (0.046536, 0.054864) |
| 6 | 3.2925 | (0.039499, 0.045742) | (0.039076, 0.045162) |

Early lyaforecast gives usable formal LBG-forest BAO constraints in bins 2–6;
bin 1 has very little information and its order-unity Fisher errors should not
be interpreted as a useful BAO detection. McDonald gives sensible formal
constraints in bins 3–6, but has no demonstrated converged result in bins 1–2
under the tested stopping rule. Finite three-update values in those bins do
not resolve that distinction. These statements concern forests measured through
LBG spectra; the LBG galaxy-only auto spectrum is unaffected by forest weighting.

The user-selected plotting cut omits both components for a correlation/bin when
either uncertainty exceeds 0.2, or the result is unavailable/nonfinite. Ratios
require both operands to survive. The cut leaves gaps and affects figures only;
it neither removes valid spectra from the joint Fisher calculation nor removes
valid large errors from the saved tables and numerical checks. Eleven figures
were produced. Independent scalar Wick covariances agreed exactly in the tested
cells, and independent Fisher contractions agreed to relative
$1.89\times10^{-15}$. Thirty-three focused tests passed.

### Step 5: assess the choice and resolve the relevant grid-to-BAO sensitivity

The assessment compared the two full-sum aliasing prescriptions and identified
the early scheme's remaining low-redshift coefficient sensitivity as the
relevant question before adoption. The additional authorized check propagated
the saved converged early-scheme coefficients on all three magnitude grids into
BAO forecasts, using **only bins 2–6**. Bin 1 was explicitly excluded.

This calculation retained all other compatibility quantities and did not
repeat raw interpolation or weight iteration. It used the previously validated
stopped states at `rtol=1e-4`. On every grid, confirmed counts in bin order were
16, 16, 14, 12, 10 for QSO forests and 28, 30, 26, 24, 22 for LBG forests.
The 30 selected forest-auto states produced 240 forecasts and 480 component
comparisons across the three grid pairs.

For $\Delta\sigma/\sigma=\sigma_{\rm finer}/\sigma_{\rm coarser}-1$, the
signed changes at the global maximum absolute differences were:

| Magnitude-grid comparison | Largest individual change | Largest joint change |
| --- | ---: | ---: |
| 107 → 213 | −0.201877% | −0.0158272% |
| 213 → 425 | −0.0843233% | −0.00562013% |
| 107 → 425 | −0.286030% | −0.0214464% |

All these extrema are transverse, in bin 2; the individual extrema are the
LBG-forest auto. Every individual and joint component passes the 0.1% diagnostic
benchmark for the final 213-to-425 refinement. The largest absolute individual
changes in that refinement for bins 3–6 are only 0.00391%, 0.00331%, 0.00254%
and 0.00186%, respectively.

The original 107-node grid is nevertheless not uniformly within 0.1% of the
finest grid. In bin 2, **both components of all five spectra involving the LBG
forest** exceed 0.1% in the 107-to-425 comparison. The four cross spectra change
by approximately −0.112% to −0.117% parallel and −0.140% to −0.143% transverse.
The LBG-forest auto changes by −0.241641% and −0.286030%, respectively. All other
individual and joint components remain below 0.1% in that comparison.

For this most sensitive auto, the actual uncertainties are:

| Magnitude nodes | $\sigma_{\alpha_\parallel}$ | $\sigma_{\alpha_\perp}$ |
| --- | ---: | ---: |
| 107 | 0.05836980079 | 0.06744376343 |
| 213 | 0.05827028233 | 0.06730761020 |
| 425 | 0.05822875526 | 0.06725085418 |

The corresponding joint transverse error changes from 0.01429039692 to
0.01428733213 between the coarsest and finest grids. The decline in the final
refinement increment supports using the refined calculation for this sample;
three finite grids do not bound the remaining continuum error. Coefficient
shifts alone overstate or understate BAO effects because aliasing and pixel
noise can change in opposite directions and enter differently across modes.

All original-grid results reproduced Step 4 bitwise. Galaxy-only results were
bitwise independent of magnitude grid. An independent analytic LBG-auto Fisher
calculation on the finest grid agreed to relative $1.22\times10^{-15}$.
Both refinement figures retain all 15 individual spectra and the joint result;
none of the selected-bin points is hidden by the 0.2 plotting cut. The bounded
calculation passed independent scientific review and the focused checks.

## 6. Interpretation, remaining limits and retained alternative

The case for the early baseline is numerical robustness across the tested
forest populations, a plausible treatment of intrinsic sightline fluctuations,
and stable downstream BAO uncertainties after the tested grid refinement.
It is not based on obtaining the smallest Fisher errors: McDonald often gives
slightly tighter constraints where both schemes converge. A smaller formal
uncertainty does not establish a better physical approximation.

McDonald remains scientifically useful because it directly uses the published
moment expression for the aliasing contribution in the weighting signal. Its
failure to establish finite nonzero weights in two LBG bins is a limitation of
this explicit recurrence and stopping criterion, not a demonstration that the
published observed-power formula is wrong. Stable normalized coefficients
under vanishing weight amplitude also deserve a different mathematical
interpretation from convergence of the weights themselves; the present choice
does not redefine that criterion to accept such states.

The following limits remain explicit:

- The extra grid-to-BAO test excludes bin 1. Early weights converge there at
  fixed grid, but the $A$ refinement sensitivity of approximately 0.194% has
  not been propagated into a corresponding grid-accuracy statement for BAO.
- McDonald has no qualified converged LBG-auto solution in bins 1–2 within
  96 updates. The extra refinement calculation tested early weights only.
- The evidence concerns the retained compatibility interpolants, signed
  measures, density/SNR treatment, response and representative mode. It does
  not establish their physical accuracy or adopt signed inputs as a production
  policy. Strict nonnegative public preparation remains a separate concern.
- The 0.1% criteria test stopping sensitivity and finite grid increments;
  they are not bounds on total forecast error or evidence of global BAO
  optimality. No new survey configurations or other reference forecast cases
  were tested in this study.
- Step 5 has supplied the weighting assessment, extra grid diagnostic and
  accepted choice recorded here. Reassessment of the accuracy profile and
  the new overall package-completion plan remain to be done. The choice does
  not silently replace the separate W12 implementation or reopen closed work.

## 7. Evidence and reproducibility

The five-step sequence is defined in the
[compatibility weighting plan](../archive/planning/FISHHIGHZ_COMPATIBILITY_WEIGHTING_PLAN.md).
The [historical findings](../archive/notes/LYAFORECAST_WEIGHTING_FINDINGS.md) document the
scalar-to-prefix change, removal of aliasing from the weighting signal, and
the distinction between early and moment-based aliasing formulas. Their
pre-calculation status statements should be read at their original dates.

| Calculation | Scientific report | Independent review | Saved numerical evidence |
| --- | --- | --- | --- |
| Step 1 | [Variant/control report](../archive/reviews/compatibility-weighting-stage1.md) | [Review](../archive/reviews/compatibility-weighting-stage1-review.md) | Baseline (local-only path: `.validation/compatibility-weighting-stage1/baseline.json`) |
| Step 2 | [Trajectory/grid report](../archive/reviews/compatibility-weighting-stage2.md) | [Review](../archive/reviews/compatibility-weighting-stage2-review.md) | Summary and per-bin trajectories (local-only path: `.validation/compatibility-weighting-stage2/`) |
| Step 3 | [Adaptive stopping report](../archive/reviews/compatibility-weighting-stage3.md) | [Review](../archive/reviews/compatibility-weighting-stage3-review.md) | Summary and per-bin checks (local-only path: `.validation/compatibility-weighting-stage3/`) |
| Step 4 | [15×2pt BAO report](../archive/reviews/compatibility-weighting-stage4.md) | [Review](../archive/reviews/compatibility-weighting-stage4-review.md) | Forecasts, changes and eleven figures (local-only path: `.validation/compatibility-weighting-stage4/`) |
| Step 5 additional check | [Early-weight grid-to-BAO report](../archive/reviews/early-lyaforecast-grid-bao.md) | [Review](../archive/reviews/early-lyaforecast-grid-bao-review.md) | Forecasts, changes and two figures (local-only path: `.validation/early-lyaforecast-grid-bao/`) |

The implementation calculations used Astra (Light) agents and the independent
scientific reviews used Sol (High) agents, with coordinator checks between
assignments. Each linked scientific review passed within its stated scope.
The reports retain exact commands and numerical checks; the saved artifacts
retain counts, status, coefficients, Fisher matrices and source provenance.
The baseline source records are in the
Step-12 revision-5 compatibility bundle (local-only path: `.validation/step12-r5-20260914T191855Z/profiles-checked/`).

Relevant implementations are the
weight variants (local-only path: `fishhighz/validation/compatibility_weights.py`),
trajectory diagnostics (local-only path: `fishhighz/validation/weight_convergence.py`),
adaptive stopping (local-only path: `fishhighz/validation/adaptive_weights.py`),
BAO reassembly (local-only path: `fishhighz/validation/compatibility_bao.py`), and
grid-to-BAO calculation (local-only path: `scripts/check_early_lyaforecast_grid_bao.py`).
This report summarizes the existing evidence; writing it did not rerun forecasts
or alter numerical outputs, production defaults, or historical records.
