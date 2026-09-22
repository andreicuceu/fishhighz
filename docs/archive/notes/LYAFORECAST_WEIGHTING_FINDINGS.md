# Historical changes to the lyaforecast forest weights

Date: 2026-09-17. Scope: mathematical interpretation, Python source history, and
implications for FishHighz. This note records findings; it introduces no change
to either package's scientific prescription or implementation.

Two independent changes distinguish the current lyaforecast weighting iteration
from the full-sample, aliasing-inclusive interpretation of McDonald & Eisenstein:

1. A February 2025 plotting change replaced scalar magnitude sums with cumulative
   sums. Those arrays also entered the weight update, changing its normalization
   from one full-sample integral to a different magnitude prefix for each source.
2. An August 2025 change commented out the aliasing contribution to the signal
   power used in that update. Aliasing remains in the final covariance.

The first change has strong evidence of an unintended scientific consequence of
a plotting feature. The second is an explicit departure from the signal choice
stated by McDonald & Eisenstein, but its motivation is undocumented in the
identified diff. Neither the magnitude of their separate BAO impacts nor the
accuracy of a repaired iteration has been established by this investigation.

## 1. Sources and conventions

The source inspection used clean checkouts at:

| Repository | Inspected HEAD |
| --- | --- |
| lyaforecast | `5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81` |
| fishhighz, before adding this note | `0cfcf2a0e7ae25c424a05604b64835f2c1bcde40` |

The papers were inspected through their arXiv PDFs. Equation numbers below refer
to those PDFs, not a separately verified journal typesetting:

- [McDonald & Eisenstein (ME07), section II.B, equations 10–19](https://arxiv.org/pdf/astro-ph/0607122).
- [Font-Ribera et al. (FR14), section IV.B.1, equations 23–29](https://arxiv.org/pdf/1308.4164).

Existing local evidence includes the [W02 equation audit](../reviews/weighting-diagnostics-w02-r1.md)
and its [independent review](../reviews/weighting-diagnostics-w02-review-r1.md).
The Git history below adds a concrete Python origin for both discrepancies. It
does not establish the behavior of the historical C++ forecasting code.

For the algebra, fix a forest redshift and one selected source population. Let
`n(m) = dn_q/dm` be its differential comoving density, `L` a common forest length,
`lp` the pixel width, and `v(m) = sigma_N(m)^2` the fractional-flux pixel variance.
Let `P` and `B` denote the intrinsic 3D and 1D forest powers at the representative
weighting mode, with consistent instrumental response. Thus `P` has volume units,
`B` and `lp*v` have length units, and weights are dimensionless. The same algebra
holds in lyaforecast's angular/velocity coordinates after the stated conversions.
The noise assigned to a magnitude in the update, `PN(m)`, is distinct from the
final population-averaged pixel-noise power, `P_pixel`.

## 2. What the published integrals normalize

ME07 constructs a weighted forest field

$$
\delta_{\rm obs}(\mathbf x)
=\frac{W(\mathbf x)}{\overline W}
 [\delta_F(\mathbf x)+\delta_N(\mathbf x)].
$$

The sampling weight vanishes outside sightlines. Their approximation treats it
as uncorrelated transversely and correlated along a spectrum, with constant noise
within each spectrum. For a small transverse cell of area `a`, equation 14 gives
`Wbar = a L I1`, integrated over the full magnitude-selected population.
[ME07, equations 10–15](https://arxiv.org/pdf/astro-ph/0607122).

FR14 writes the moments and observed power as

$$
I_1=\int_{\rm sample}n(m)w(m)\,dm,\qquad
I_2=\int_{\rm sample}n(m)w(m)^2\,dm,
$$

$$
I_3=\int_{\rm sample}n(m)v(m)w(m)^2\,dm,
$$

$$
A=\frac{I_2}{L I_1^2},\qquad
P_{\rm pixel}=\frac{l_p I_3}{L I_1^2},\qquad
P_{\rm obs}(\mathbf k)=P_F^{3D}(\mathbf k)
 +A P_F^{1D}(k_\parallel)+P_{\rm pixel}.
$$

These are FR14 equations 23–28, corresponding to ME07 equations 13–18. The
integrals characterize one selected population; FR14 does not specify a running
upper bound tied to the magnitude being weighted.
[FR14, section IV.B.1](https://arxiv.org/pdf/1308.4164).

Interpretation: `I1` normalizes the mean sampled field; `I2` and `I3` enter its
second moments. The denominator is squared because a two-point function contains
two field normalizations. This construction provides one normalization per
selected sample, not one normalization per source magnitude. Common rescaling
`w -> c*w` cancels in `A` and `P_pixel`.

The choice of weights is an additional approximation. ME07 adopts

$$
w(m)=\frac{P_S}{P_S+P_N(m)},\qquad
P_N(m)=\frac{v(m)l_p}{L I_1}.
$$

After equation 19, ME07 explicitly includes aliasing in the representative
`P_S`, evaluated at `k = 0.07 h/Mpc`, `mu = 0.5`, while acknowledging ambiguity
in the noise assignment. FR14 equation 29 uses the same weight and noise formulas
and representative mode, referring back to ME07, but does not explicitly repeat
the inclusion of aliasing. Its wording alone does not establish a deliberate
switch to intrinsic 3D power.
[ME07, discussion following equation 19](https://arxiv.org/pdf/astro-ph/0607122);
[FR14, discussion following equation 29](https://arxiv.org/pdf/1308.4164).

The paper discussions do not fully specify how an aliasing-inclusive signal is
initialized or refreshed during iteration. Recomputing it from the current
moments is a concrete implementation choice; it should not be presented as a
uniquely specified published algorithm. Nor is the FKP-like choice a proof of
globally optimal BAO weighting.

## 3. February 2025: scalar sums became magnitude prefixes

### Concrete change and its context

Commit [1eeb28ee62ab058b026d208763ae5f3f2d7d713e](https://github.com/igmhub/lyaforecast/commit/1eeb28ee62ab058b026d208763ae5f3f2d7d713e),
dated 2025-02-13, is titled **Add options for BAO as a function of magnitude**.
In `lyaforecast/covariance.py::_compute_int_1`, it replaced

```python
I1 = np.sum(dndm_degkms*weights)*dm
```

with

```python
integrand = dndm_degkms * weights
# weighted density of quasars
#I1 = np.sum(dndm_degkms * weights) * dm
#move to using cumsum so we can plot as a function of magnitude
int_1 = np.cumsum(integrand) * dm
```

`I2` and `I3` also changed from sums to cumulative sums. The commit introduced
`per_mag` and selected `power_variance[-1]` when only the final magnitude limit
was requested. It also changed the magnitude grid to `self.survey.maglist`.
Consequently, a comparison across the entire commit would conflate the recurrence
change with grid and other changes.

The scientifically consequential call chain was already shared:

```text
_compute_weights -> _weights1 -> _compute_noise_power_m
                                -> _get_np_eff -> _compute_int_1
```

`_get_np_eff` multiplied the returned integral by `L/lp`; the noise routine divided
the magnitude-dependent pixel variance by that result. Returning an array thus
changed the update through elementwise array arithmetic. Selecting the final
forecast element afterward did not restore a scalar normalization during the
preceding iterations.

The current equivalents are weights.py (local-only path: `../lyaforecast/lyaforecast/weights.py`),
lines 135–156, 184–206, 255–273, and 339–358. The present forest auto-power takes
the final coefficient elements in
covariance.py (local-only path: `../lyaforecast/lyaforecast/covariance.py`), lines 335–365.

### Different equations, independent of quadrature accuracy

Let `r_j = n(m_j) q_j`, with magnitude quadrature weights `q_j`. A full-sample
update uses

$$
I_1^{(t)}=\sum_{j=1}^{N}r_jw_j^{(t)},\qquad
P_{N,i}^{(t)}=\frac{v_i l_p}{L I_1^{(t)}}.
$$

The prefix update uses

$$
I_{1,i}^{(t)}=\sum_{j\le i}r_jw_j^{(t)},\qquad
P_{N,i}^{(t)}=\frac{v_i l_p}{L I_{1,i}^{(t)}}.
$$

Bright sources therefore use a density excluding fainter sources in the same
catalogue. This remains a different equation with a perfectly known smooth
luminosity function. Increasing the iteration count, or replacing rectangular
quadrature by Gauss–Legendre quadrature, does not convert one equation into the
other.

A cumulative integral is appropriate for plotting successive magnitude cuts with
fixed weights. For self-consistent weights at a limiting magnitude `M`, however,
the calculation must solve for `w(m; M)` with one common `I1(M)` for every source
in that selected sample. A single vector updated using its elementwise prefixes
does not solve that family of separate sample problems.

An elementary physical check is a positive-density population whose spectra have
identical noise, length, and response. The scalar update gives equal weights.
The current intrinsic-signal prefix update gives unequal weights solely because
their magnitude ranks differ. With `n_tot = integral n(m) dm`, Cauchy–Schwarz gives

$$
\frac{I_2}{I_1^2}\ge\frac{1}{n_{\rm tot}},
$$

with equality for constant weights. For identical noise, both aliasing and pixel
noise are minimized by those constant weights. The introduced magnitude
dependence has no advantage in this example.

**Assessment:** the diff and call chain strongly support an unintended change to
the weighting equation during a plotting feature. Author intent is inferred,
not proven. Arbitrary finite prefix-derived weights can still define a weighted
estimator with a consistently normalized covariance; the objection is that they
do not implement the stated full-sample prescription and have no demonstrated
optimality.

## 4. August 2025: aliasing was removed from the weighting signal

### Concrete removal and retained covariance term

Commit [6d079143a1bbe6d4056d79280c974afd7ed551a6](https://github.com/igmhub/lyaforecast/commit/6d079143a1bbe6d4056d79280c974afd7ed551a6),
dated 2025-08-06, is titled **Initial commit wiht new modules** (original spelling).
Its change to `lyaforecast/weights.py` is:

```diff
 aliasing = int_2 / (int_1**2 * self._forest_length)
 # weights include aliasing as signal
-signal_power = self._p3d_w + self._p1d_w * aliasing
+signal_power = self._p3d_w #+ self._p1d_w * aliasing
```

Thus the pre-change expression was `P + B*A`, with `A = I2/(L*I1**2)`; afterward
it was just `P`. At this date the integrals were already cumulative, so restoring
this line alone would retain the prefix issue.

The current weight update (local-only path: `../lyaforecast/lyaforecast/weights.py`), lines 149–154,
still has the aliasing comment followed by `signal_power = self._p3d_w`.
Tracing power_spectrum.py (local-only path: `../lyaforecast/lyaforecast/power_spectrum.py`), lines
109–183, shows that `_p3d_w` contains intrinsic clustering and instrumental
smoothing, not sampling aliasing. The final covariance (local-only path: `../lyaforecast/lyaforecast/covariance.py`),
lines 354–361, still adds `p3d + aliasing + noise`.

This changes the weights used to form the covariance, not the presence of
aliasing in its final formula. The identified commit supplies no scientific
explanation for the removal.

### Earlier history: there was more than one aliasing formula

The immediately preceding August implementation should not be conflated with the
earlier Python approximation:

| Revision | Weighting signal and context |
| --- | --- |
| February scalar-to-prefix change, `1eeb28e` | Retained `P + B/(L*I1)` in `_weights1`; the change made `I1` an array. |
| [6fe1083](https://github.com/igmhub/lyaforecast/commit/6fe10831d3865fa561e5ca0410fda97879574428) | Introduced `I2` in the weighting aliasing term, but wrote `int_2 / (int_1**2 + forest_length)`. That addition is dimensionally inconsistent. Separately corrected the final covariance coefficient to `I2/(I1**2*L)`. |
| [45f61c8](https://github.com/igmhub/lyaforecast/commit/45f61c8e3ad2a7a3d311c3dbec789312d91fc53a) | Corrected the weighting denominator to `int_1**2 * forest_length`, giving `P + B*I2/(L*I1**2)`. |
| [4a1b5d7](https://github.com/igmhub/lyaforecast/commit/4a1b5d7e970711a78fedc5f51c54241cc2f1b1b0) | Moved the weighting calculation into `weights.py`. |
| August removal, `6d07914` | Commented out `B*I2/(L*I1**2)` in the weight update. |

This sequence documents relevant expressions, not a validation of every historical
revision. In particular, `B/(L*I1)` and `B*I2/(L*I1**2)` are not interchangeable
for general weights. Their ratio involves `I2/I1`, and their response to a common
weight rescaling differs.

### Mathematical effect on weights and iteration

To separate signal choice from the prefix issue, the following comparisons use
full-sample scalar moments, positive density, and positive noise. Define

$$
C=P L I_1,\qquad D=B\frac{I_2}{I_1},\qquad N_i=l_p v_i.
$$

For a given input weight vector, the three update rules are

| Signal choice | Updated weight |
| --- | --- |
| Intrinsic 3D only | `C / (C + N_i)` |
| Recomputed aliasing from the moments | `(C + D) / (C + D + N_i)` |
| Earlier Python approximation `P + B/(L*I1)` | `(C + B) / (C + B + N_i)` |

Including positive aliasing raises the weights at fixed input moments and gives
relatively more weight to noisier spectra. For example, with `D = 9*C`, spectra
with `N_i = C` and `N_i = 10*C` receive `(0.5, 1/11)` without aliasing and
`(10/11, 0.5)` with it. This is one algebraic update, not a measured survey result.

The physical distinction is between instrumental noise and the intrinsic forest
fluctuations sampled by each sightline. Even a noiseless spectrum retains the
latter. Their contribution limits the gain from improving one spectrum compared
with adding independent sightlines. Including aliasing in an auxiliary weighting
signal does not, by itself, add new BAO information to the mean derivatives.

Restoring a scalar sum does not guarantee a nonzero fixed point. Consider a
homogeneous population of density `n`, common weight `w`, and common 1D noise
power `N > 0`. Set `s = P*L*n > 0`. The intrinsic-only update is

$$
w^{(t+1)}=\frac{s w^{(t)}}{s w^{(t)}+N}.
$$

Its positive fixed point is `1 - N/s` only if `s > N`; otherwise positive
iterates approach zero. With the aliasing coefficient recomputed from the
moments, the common-weight example instead gives

$$
w^{(t+1)}=\frac{(s+B)w^{(t)}}{(s+B)w^{(t)}+N}.
$$

The threshold becomes `s+B > N`, so aliasing does not guarantee a nonzero limit
either. The earlier approximation gives

$$
w^{(t+1)}=\frac{s w^{(t)}+B}{s w^{(t)}+B+N},
$$

which remains positive as the previous weight approaches zero when `B > 0`.
These are distinct recurrences. None of these scalar examples proves convergence
of the actual magnitude-dependent prefix calculation.

Also distinguish zero amplitude from divergent physical coefficients: for common
positive weight, `I2/I1**2 = 1/n` remains finite as that amplitude tends to zero.
The exact all-zero vector gives undefined moment ratios in direct arithmetic.
Amplitude therefore matters in the nonlinear update even when it cancels in the
limiting normalized covariance.

An independent check, under fixed-mode, common-response, nonnegative-density
assumptions, is to minimize the weight-dependent covariance contribution

$$
Q[w]=\frac{\int n(m)[B+l_pv(m)]w(m)^2\,dm}
 {L[\int n(m)w(m)\,dm]^2}.
$$

Cauchy–Schwarz gives `w(m) proportional to 1/[B+lp*v(m)]`. Its normalized shape
is the legacy seed `B/[B+lp*v(m)]`, before subsequent updates. This restricted
variational result explains the role of intrinsic forest variance without
claiming that either iterative prescription is the global multi-mode BAO optimum.
It is also documented in the [W05 review](../reviews/weighting-diagnostics-w05-review-r1.md).

## 5. Downstream consequences and limits of the evidence

The direct propagation is

```text
normalization and representative signal in the update
    -> relative weights versus magnitude
    -> I1, I2, I3
    -> aliasing coefficient A and pixel-noise power P_pixel
    -> observed auto powers and Gaussian covariance blocks
    -> Fisher matrix and marginalized BAO uncertainties.
```

For an auto-spectrum, the Gaussian variance contains `2*P_obs**2/N_modes`.
Cross-spectrum covariance contains products of the relevant observed auto powers,
so a forest weight change can affect joint constraints even when intrinsic
cross-power predictions and BAO derivatives are held fixed. Agreement of the
final moment formulas does not establish agreement of their numerical values.

The following conclusions have different evidential status:

- **Established by source history:** the scalar-to-array change and the later
  explicit removal of the weighting aliasing term occurred in separate commits.
  Their effects enter the executable weight update, and final aliasing remains.
- **Established by algebra:** the resulting equations differ; scalar and prefix
  normalization are not interchangeable; the signal variants have different
  relative weights and fixed-point behavior. Global weight rescaling cancels in
  final coefficients but generally changes the next update.
- **Supported interpretation:** the cumulative feedback is a likely regression
  associated with magnitude-limit plotting. The aliasing removal departs from
  ME07's explicit signal choice, but its intention remains unknown.
- **Not quantified here:** separate or combined changes to survey BAO errors.
  No new forecast or controlled four-way comparison was run. These findings
  establish neither a less-than-5% error bound nor a greater-than-5–10% shift.

Magnitude-grid sensitivity and iteration sensitivity must be distinguished.
For the current intrinsic-only prefix update, the first magnitude cell uses

$$
w_1^{(t+1)}=
\frac{P L r_1 w_1^{(t)}}{P L r_1 w_1^{(t)}+l_pv_1}.
$$

Refining a smooth-density grid reduces the first-cell mass `r1`, changing this
feedback directly. This identifies a mechanism for sensitivity; it does not
prove a continuum result for the normalized full population. Running more
iterations addresses a fixed-grid recurrence, not the choice of equation or the
order of the grid and iteration limits.

Finite-resolution luminosity-function measurements require an explicit
within-bin representation or interpolation in either prescription. Quadrature
can approximate that chosen representation; it cannot recover unresolved
luminosity-function structure or repair a changed normalization. Positivity and
optimality arguments in this note assume nonnegative density. They do not justify
silently clipping the signed interpolated densities present in some saved
diagnostics.

FishHighz deliberately preserves the inherited recurrence in its explicit
`method="legacy"` path: see kernels/weights.py (local-only path: `fishhighz/kernels/weights.py`),
lines 9–37, and weights.py (local-only path: `fishhighz/weights.py`). Its agreement with lyaforecast
is a compatibility statement. The [W12 independent review](../reviews/weighting-diagnostics-w12-review-r1.md)
records a separate accuracy-profile route using fixed inverse-variance weights,
with bounded saved-input validation. That route does not execute either nonlinear
iteration discussed here. Its existence does not repair, relabel, or scientifically
accept the historical cumulative forecasts.

An attribution study would need to vary full-sample versus prefix normalization
and intrinsic-only versus a precisely specified aliasing-inclusive signal
separately, holding source inputs, quadrature, response, seed, Fourier modes,
derivatives, and stopping criteria fixed. That comparison is not part of this
documentation change. Reverting a whole historical commit would also change
unrelated numerical choices and would not isolate either effect.

## 6. Reproducing the source inspection

From the workspace root, these read-only commands expose the principal evidence:

```bash
git -C lib/lyaforecast show 1eeb28ee62ab058b026d208763ae5f3f2d7d713e -- lyaforecast/covariance.py
git -C lib/lyaforecast show 1eeb28e:lyaforecast/covariance.py
git -C lib/lyaforecast show 6fe10831d3865fa561e5ca0410fda97879574428 -- lyaforecast/covariance.py
git -C lib/lyaforecast show 45f61c8e3ad2a7a3d311c3dbec789312d91fc53a -- lyaforecast/covariance.py
git -C lib/lyaforecast show 4a1b5d7e970711a78fedc5f51c54241cc2f1b1b0 -- lyaforecast/weights.py
git -C lib/lyaforecast show 6d079143a1bbe6d4056d79280c974afd7ed551a6 -- lyaforecast/weights.py
git -C lib/lyaforecast show 5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81:lyaforecast/weights.py
git -C lib/lyaforecast show 5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81:lyaforecast/covariance.py
```

The findings were checked against the paper text, these diffs, the current
call chains, and the linked diagnostic reports. The examples above are analytic
controls, not new numerical forecast evidence. No production code, historical
review, saved result, or implementation plan was changed for this note.
