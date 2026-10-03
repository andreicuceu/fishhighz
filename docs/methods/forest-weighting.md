# Forest weighting

Each Lyα-forest sample (forest field, redshift bin) carries a scalar weight
function $w_i$ over source magnitude $m_i$. The weights set the aliasing and
pixel-noise power that enter the covariance. This note gives the adopted
prescription, its stopping criteria and the alternatives. Call signatures are in
[Forest weighting and sampling noise](weights.md).

## Definitions

For one forest population at fixed redshift let $n(m)$ be the source density per
magnitude, $q_i$ the magnitude quadrature weights, $r_i=n(m_i)q_i$, $L$ the forest
length, $l_p$ the pixel width, $v_i$ the dimensionless pixel-noise variance, $P$
the response-smoothed 3D forest power and $B$ the response-smoothed 1D power at
the representative weighting mode (for BAO,
$(k_{t},k_{p})=(2.4\ \mathrm{deg^{-1}},\,3.5\times10^{-4}\ \mathrm{s\,km^{-1}})$;
supplied through `auxiliary` samples or explicit `signal` and `alias`). The
full-sample moments are scalars for the whole population,

$$
I_1=\sum_i r_i w_i,\qquad I_2=\sum_i r_i w_i^2,\qquad I_3=\sum_i r_i v_i w_i^2 .
$$

The weights are seeded with $w_i^{(0)}=B/(B+l_p v_i)$ and updated simultaneously,
$w_i^{(t+1)}=S^{(t)}/(S^{(t)}+l_p v_i/(L I_1^{(t)}))$, where $S^{(t)}$ is the
representative signal of the prescription. After the weights are fixed, the
noise coefficients of the final covariance are

$$
A=\frac{I_2}{L I_1^{2}},\qquad P_{\rm pixel}=\frac{l_p I_3}{L I_1^{2}},\qquad
P_{\rm obs}=P^{\rm 3D}_{F}+A\,P^{\rm 1D}_{F}(k_\parallel)+P_{\rm pixel},
$$

with the instrumental response applied once to the aliasing term and not to
$P_{\rm pixel}$. The prescription affects only the weights; the final aliasing
coefficient is always $A$ as above. Weights, $A$ and $P_{\rm pixel}$ are held
fixed when differentiating the mean model.

## Adopted prescription: `early_lyaforecast`

$$
S^{(t)}=P+\frac{B}{L I_1^{(t)}}
\;\Rightarrow\;
w_i^{(t+1)}=\frac{P L I_1^{(t)}+B}{P L I_1^{(t)}+B+l_p v_i}.
$$

Only $I_1$ enters the update. The term $B/(L I_1)$ keeps the numerator positive
as $I_1\to0$ and, for small $P$, the weights tend to $B/(B+l_pv_i)$. This
prescription follows the early `lyaforecast` implementation. It is the only
weighting selectable from the INI interface.

**Stopping.** The iteration is adaptive. At each update five relative changes
are formed, each in the infinity norm against the previous state and without an
absolute floor: the weight amplitude, the normalised signed weight shape, the
complete weight vector, $A$, and $P_{\rm pixel}$. A state at update $t$ becomes a
candidate when $t\ge3$ and the last three transitions all satisfy
`rtol = 1e-4`. The candidate is accepted only if every transition up to update
$2t$, and every comparison of the state at $2t$ with the candidate, also
satisfies `rtol`; the state at $2t$ is returned. The cap is 96 updates, so
candidates with $2t>96$ are not nominated. Failure to converge, non-finite
arithmetic, or a non-positive $P$ or $B$ raises `ValueError`; the last iterate is
never substituted and no other prescription is used as a fallback. The
controls are `rtol`, `min_updates`, `stable_steps` and `max_updates`. The
convergence record (status, update count, candidate, residuals) is stored in the
prepared weights.

## Alternative: `mcdonald`

$$
S^{(t)}=P+\frac{B\,I_2^{(t)}}{L\,[I_1^{(t)}]^2},
$$

i.e. the aliasing term in the representative signal uses the current moments
(McDonald & Eisenstein 2007, sec. II.B; Font-Ribera et al. 2014, eqs. 23-29).
Since $I_2\le I_1$ for $0\le w\le1$ and non-negative density, this gives a smaller
aliasing contribution than the early prescription at the same weights and
suppresses noisy spectra more strongly. In the homogeneous limit with
$s=PLn$ and noise $N=l_pv$, the McDonald map is $w\mapsto(s+B)w/[(s+B)w+N]$,
which has a positive fixed point only for $s+B>N$, whereas the early map tends to
$B/(B+N)>0$ as $w\to0$. In the DESI-2 SRD tests the McDonald recurrence did
not meet the convergence criterion in the two lowest-redshift LBG-forest bins
(it raises there); the early prescription converged for both forest populations
in all six bins. `mcdonald` is available through the Python API
(`prepare_forest_weights(..., method="mcdonald")`, central mode only), not the
INI interface.

## Central and integrated forest sources

In the default *central* mode the density and pixel variance are sampled at one
source redshift and the quadrature is over magnitude. The *integrated* mode
(`[input policies] forest_source_integration = integrated`,
`early-lyaforecast-integrated-2026-10-02`) integrates source redshift $z_q$ and
forest pixel inside the bin with the measure
$\mu_{pj}=(1+z_{q,p})\,w_{y,p}\,c\,w_{u,p}\,\mathrm{d}n/\mathrm{d}z\,w_{m,j}/L_{\rm bin}$.
The same full-sum kernels run on the flattened $(p,j)$ arrays with
$\rho\to\mu$, $q\to1$, $L\to1$, so that $L I_1\to N_1=\sum\mu w$,
$N_2=\sum\mu w^2$, $N_3=\sum\mu w^2v$, $S=P+B/N_1$, $A=N_2/N_1^2$ and
$P_{\rm pixel}=l_pN_3/N_1^2$. Only `early_lyaforecast` is available. When
$(1+z_q)\,\mathrm{d}n/\mathrm{d}z_q$ is constant and the variance depends on
magnitude alone, the integrated moments reduce to the central ones (relative
differences $\sim10^{-13}$). The S/N tables are per angstrom; the pixel variance
is queried at $\Delta\lambda=\Delta v\,\lambda/c$, so $v\,\Delta v$ is independent
of the angstrom width. The quadrature orders are set by `forest_zq_order`,
`forest_lambda_order`, `forest_lambda_panels` (see the
[INI reference](../user/ini.md#integrated-forest-sources));
`scripts/check_integrated_forest_refinement.py` compares successive orders.

## Magnitude cumulative recurrence is not supported

The original `lyaforecast` code evaluated the weight update with a cumulative sum
over magnitude, $I_{1,i}=\sum_{j\le i}r_jw_j$, so each source saw a different
normalisation that depends on its position in the magnitude array. This does not
correspond to the intended full-population moments of the weighting
derivation and was identified as a bug. FishHighz uses the scalar full-sum
moments throughout, and the cumulative recurrence is not implemented. Forecasts
made with it are not reproduced.

## Numerical and open points

- The magnitude quadrature is the rectangular rule with full endpoint weights.
  A 425-node grid reproduced individual BAO uncertainties in redshift bins 2-6 to
  better than 0.1% (joint uncertainties to $<6\times10^{-5}$). No universal node
  count or continuum bound is established, so the grid is configurable.
- The weighting mode is fixed in angular and velocity coordinates, so the
  corresponding comoving mode varies with redshift. The sensitivity of the BAO
  forecasts to this choice (e.g. a fixed comoving mode matched at a pivot
  redshift) has not been measured.
- The weights are optimal for neither a multi-mode BAO estimator nor a global
  figure of merit; they approximate a representative-mode FKP-type weighting.
- Signed or interpolated densities can violate the positivity assumptions behind
  the McDonald comparison above.
