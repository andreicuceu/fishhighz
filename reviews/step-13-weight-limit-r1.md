# Step 13 revision 1: cumulative-weight limit

## Result

The existing cumulative forest-weight recurrence has a well-defined
amplitude/shape evolution and a finite asymptotic noise pair on every fixed
saved quadrature grid examined. It does **not** have a suitable finite,
magnitude-grid-independent limiting pair on the saved DESI-2 samples. This
negative result covers `lya(qso)` and `lya(lbg)` in all six 15x2pt redshift
bins. It is a statement about the present recurrence and inputs; it neither
selects an alternative estimator nor changes the production implementation.

## Recurrence, domain, and units

Let (r_i=\rho_iq_i>0), (v_i\geq0),
(\alpha=L/\Delta_v), (a_i=v_i/S), and let (B) denote the auxiliary
alias term. The live increasing-magnitude, inclusive-prefix recurrence is

\[
w_i^{(0)}=\left(1+\Delta_vv_i/B\right)^{-1},\qquad
N_i^{(t)}=\alpha\sum_{j\leq i}r_jw_j^{(t)},\qquad
w_i^{(t+1)}={N_i^{(t)}\over N_i^{(t)}+a_i}.
\]

Here (r_i) has units deg\(^{-2}\)(km s\(^{-1}\))\(^{-1}\), while (v_i)
and (w_i) are dimensionless; (S) has units deg\(^2\) km s\(^{-1}\), and
(B,L,\Delta_v) are in km s\(^{-1}\). Thus (N_i) and (a_i) have the
same units. With

\[
I_1=\sum_i r_iw_i,\quad I_2=\sum_i r_iw_i^2,\quad
I_3=\sum_i r_iw_i^2v_i,
\]

the coefficients (A=I_2/(LI_1^2)) and
(P_{\rm pixel}=\Delta_vI_3/(LI_1^2)) have units deg\(^2\) and
deg\(^2\) km s\(^{-1}\), respectively.

The diagnostic requires finite, non-negative masses and variances, positive
(L,\Delta_v,S,B), strictly increasing physical magnitude coordinates, and
positive total mass. Zero-density cells are excluded exactly. Supported
zero-variance cells retain the live recurrence's exact branch; the discrete
spectral reduction below is restricted to strictly positive variance.

## Exact range-safe formulation

Write (w_i=s u_i), with (s>0) and \(\max_i u_i=1\). Defining

\[
h_i=\alpha\sum_{j\leq i}r_ju_j,\qquad
b_i={h_i\over a_i+s h_i},
\]

gives the algebraic identity

\[
s' =s\max_i b_i,\qquad u'_i={b_i\over\max_jb_j}.
\]

`weight_limit.py` evolves \(\log s\) and \(\log u_i\) separately with
`logaddexp`, so it never reconstructs an unrepresentably small common
amplitude. This distinction is essential: setting \(\max w=1\) before the next
nonlinear update changes the map. Coefficients are evaluated after removing
the common amplitude, before moment summation. If

\[
J_n=\sum_i r_i u_i^n,\qquad J_{2v}=\sum_i r_i u_i^2v_i,
\]

then

\[
\log A=\log J_2-2\log J_1-\log L,
\quad
\log P_{\rm pixel}=\log\Delta_v+\log J_{2v}-2\log J_1-\log L.
\]

Ordinary values are accompanied by explicit availability flags; mathematical
zeros remain exact. The saved distribution
(p_i=r_iu_i/\sum_jr_ju_j) obeys

\[
\sum_i {p_i^2\over r_i}=LA,
\]

so its inverse is a continuous-measure effective support. Any relative
components lost on exponentiation are bounded in distribution mass, with a
mandatory (10^{-12}) failure threshold.

## Discrete and continuum asymptotics

For positive variance, the derivative at zero amplitude is the lower
triangular matrix

\[
M_{ij}={\alpha r_j\over a_i}\;\mathbf{1}_{j\leq i},
\qquad d_i=M_{ii}={\alpha r_i\over a_i}.
\]

Its spectral radius is \(\max_i d_i\). Values below, equal to, and above unity
give locally decaying, marginal, and unstable zero-amplitude branches,
respectively; the one-cell solution makes these cases explicit. For one cell,

\[
w_t={d^t\over 1/w_0+d(1-d^t)/(1-d)}
\]

for (d\ne1), while (w_t=w_0/(1+tw_0)) for (d=1). Despite amplitude
decay for (d\leq1), (A=1/(Lr)) and
(P_{\rm pixel}=\Delta_vv/(Lr)): vanishing weights alone do not imply
divergent noise.

Repeated dominant diagonals can introduce generalized-eigenvector polynomial
factors; near degeneracy produces long transients, and the prefix matrix is
nonnormal. Accessibility also matters. For a unique dominant index (k),
positive initial support, and \(\rho(M)<1\), the asymptotic fixed-grid shape is
the positive eigenvector obtained without a dense eigensolver:

\[
x_i=0\;(i<k),\quad x_k=1,\quad
x_i={\alpha S\sum_{j<i}r_jx_j/v_i\over \lambda-d_i}\;(i>k),
\quad\lambda=d_k.
\]

The nonlinear map satisfies (F(w)=Mw+O(\lVert w\rVert^2)); componentwise,
the difference between (N/a) and (N/(N+a)) is bounded by
((N/a)^2). Once the saved trajectories enter the small-amplitude regime,
the quadratic remainder is summable under the observed contraction. All saved
grids have positive initial support, a unique dominant diagonal, and
(2.99338\times10^{-5}\leq\rho(M)\leq0.407248). The relative spectral gaps
range from (3.08112\times10^{-5}) to (5.87705\times10^{-3}), explaining
slow normalized-shape evolution on some grids.

The fixed-grid conclusion does not survive continuum refinement. For the
constant-coefficient Volterra control
(Kf(x)=c\int_0^x f(y)\,dy), (K^t1=c^tx^t/t!\). After normalization the
shape is (x^t), but for total measure (R)

\[
A_t={(t+1)^2\over(2t+1)LR},\qquad
P_{{\rm pixel},t}=\Delta_vvA_t,
\]

which diverges linearly. Thus fixed-grid iteration, continuum refinement, and
joint limits cannot be exchanged without an additional uniform-integrability
or non-concentration result.

## Saved-input experiment

The experiment used only the hash-verified Step 12 r5 arrays. It evaluated 60
independent batches: six bins, two forest populations, and quadrature orders
4/8/16/32/64. Each batch recorded 14 checkpoints through 1024 updates and
finished in at most 0.812 s. All 840 checkpoint attempts completed. At 1024
updates, \(\log s\) spans -9136.91 to -825.159, while the shape-domain
coefficients remain representable.

The exact fixed-grid asymptotes for the three finest saved orders are listed as
`A / P_pixel`; (p_P) is the order exponent inferred from orders 16 and 64.

| Bin | Population | Order 16 | Order 32 | Order 64 | \(p_P\) |
|---:|:---|---:|---:|---:|---:|
| 0 | lya(qso) | 0.0374662 / 12.6340 | 0.0572858 / 24.7726 | 0.0971897 / 49.0883 | 0.9790 |
| 0 | lya(lbg) | 2.34876 / 7671.75 | 4.39865 / 14849.5 | 8.26859 / 29171.9 | 0.9635 |
| 1 | lya(qso) | 0.0516760 / 9.66180 | 0.0755132 / 18.9483 | 0.123829 / 37.5533 | 0.9793 |
| 1 | lya(lbg) | 0.249034 / 440.670 | 0.465054 / 852.387 | 0.870604 / 1673.26 | 0.9624 |
| 2 | lya(qso) | 0.0700742 / 17.2293 | 0.101165 / 33.8050 | 0.165005 / 67.0146 | 0.9798 |
| 2 | lya(lbg) | 0.0868774 / 167.482 | 0.162333 / 324.000 | 0.304206 / 636.159 | 0.9627 |
| 3 | lya(qso) | 0.114856 / 18.9915 | 0.160518 / 37.2486 | 0.254435 / 73.8281 | 0.9794 |
| 3 | lya(lbg) | 0.122153 / 155.152 | 0.227660 / 300.048 | 0.424282 / 588.357 | 0.9615 |
| 4 | lya(qso) | 0.228518 / 24.4476 | 0.312158 / 47.9641 | 0.482309 / 95.0825 | 0.9797 |
| 4 | lya(lbg) | 0.509868 / 474.238 | 0.950585 / 917.203 | 1.77301 / 1798.96 | 0.9617 |
| 5 | lya(qso) | 0.815482 / 47.5971 | 0.931840 / 93.3102 | 1.18096 / 184.934 | 0.9790 |
| 5 | lya(lbg) | 0.479528 / 390.957 | 0.894049 / 756.113 | 1.66781 / 1483.05 | 0.9617 |

For every population, (P_{\rm pixel}) rises by more than 80% on both
16-to-32 and 32-to-64 refinement. Its measured order exponent is
0.9615--0.9798. Meanwhile the effective-measure exponent is negative in every
case (-0.9079 to -0.2671), directly diagnosing concentration onto shrinking
quadrature mass. The dominant physical magnitude moves by at most 0.00234 and
the displacement decreases on the second refinement, so the trend is not an
artifact of a migrating physical feature. The finite 1024-update trajectories
also fail the required late-doubling screen: their last-doubling changes span
0.128--0.377 in (A) and 0.342--0.528 in (P_{\rm pixel}), far above
(10^{-4}).

The prototype reproduces all 144 historically available 3/6/12/24 coefficient
pairs with maximum relative discrepancies (8.95\times10^{-14}) in (A) and
(1.22\times10^{-13}) in (P_{\rm pixel}). Of 240 historical coefficient
attempts, 96 were unavailable; 56 of 200 saved-weight states are rejected by
the unchanged production integrals. These are retained historical arithmetic
failures, not missing Step 13 trajectories. Twelve order-64 samples also match
the independently saved profile arrays at relative tolerance (5\times10^{-12}).

## Scientific implication

The unique fixed-grid asymptote concentrates as cell mass decreases, and its
pixel-noise coefficient grows approximately linearly with quadrature order for
all 12 populations. More floating-point precision or more recurrence updates
cannot turn that refinement trend into a finite grid-independent limit. Holding
the other factors fixed in

\[
N_F(k,\mu)={\left[A P_{1D}(k_\parallel)W^2(k_\parallel)+P_{\rm pixel}\right]
d_{\rm deg}^2\over a_v},
\]

the diverging additive (P_{\rm pixel}) term is sufficient to reject a finite
forest-noise limit for this recurrence. No forecast was reassembled, so this
does not establish new Fisher results or scientifically accept Step 12. The
next scientific decision is to specify an alternative forest estimator or an
explicit finite-iteration/discretization prescription; neither is selected by
this investigation.

Numerical records, decision table, and separate-population figures are under
`.validation/step13-r1-20260914T235633Z/`.
