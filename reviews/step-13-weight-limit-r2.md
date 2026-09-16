# Step 13 revision 2: cumulative-weight limit

## Result and domain

The exact finite recurrence and every fixed saved grid remain well defined in a
log-amplitude/log-shape representation. The available saved refinements do not,
however, establish either a finite continuum limit or divergence of the full
nonlinear estimator. All twelve DESI-2 forest populations are therefore
**unresolved within the stated bounds**. This corrects revision 1's promotion of
a strong finite-refinement trend to a continuum theorem.

The diagnostic domain is finite non-negative cell masses (r_i), finite
non-negative variances (v_i), positive finite (L,\Delta_v,S,B), positive
total mass, and strictly increasing physical magnitude coordinates. Positive
supported weights retain finite logarithms; exact zero mass and zero variance
use explicit branches. Ordinary coefficient values carry availability flags.
The production recurrence, guards and scientific prescription are unchanged.

## Recurrence, units and range-safe form

For (r_i=\rho_iq_i), (a_i=v_i/S), and
\(\alpha=L/\Delta_v\),

\[
w_i^{(0)}=(1+\Delta_vv_i/B)^{-1},\quad
N_i^{(t)}=\alpha\sum_{j\leq i}r_jw_j^{(t)},\quad
w_i^{(t+1)}={N_i^{(t)}\over N_i^{(t)}+a_i}.
\]

Here (r_i) has units deg\(^{-2}\)(km s\(^{-1}\))\(^{-1}\),
(v_i,w_i) are dimensionless, (S) has units deg\(^2\) km s\(^{-1}\),
and (B,L,\Delta_v) are in km s\(^{-1}\). Thus (N_i) and (a_i)
have matching units. The scale-free noise coefficients are

\[
A={\sum r_iw_i^2\over L(\sum r_iw_i)^2},\qquad
P_{\rm pixel}={\Delta_v\sum r_iw_i^2v_i\over L(\sum r_iw_i)^2}.
\]

Writing (w_i=su_i), with \(\max u=1\), gives

\[
h_i=\alpha\sum_{j\leq i}r_ju_j,\quad
b_i={h_i\over a_i+sh_i},\quad
s'=s\max b,\quad u_i'=b_i/\max b.
\]

The implementation retains \(\log s\) and \(\log u_i\). It evaluates every
scalar ratio through differences of logarithms, including
\(\log\alpha=\log L-\log\Delta_v\); it never first forms
\(L/\Delta_v\). Moment logarithms cancel the common amplitude before summation.
The required one-cell controls

\[
(r,v,L,\Delta_v)=(10^{300},10^{-300},10^{-300},10^{300})
\]

and its reciprocal-range counterpart both give (w_0=1/2), (d=1),
(w_1=1/3), and (A=P_{\rm pixel}=1). The diagnostic also agrees with an
independent Decimal recurrence after all ordinary float64 weights underflow.
Production still rejects the Step 10 mixed-product and saved
\(8.12\times10^{-566}\) cases.

## Fixed-grid result and nonlinear remainder

For strictly positive variance the derivative at zero is the lower-triangular
positive matrix

\[
M_{ij}={\alpha r_j\over a_i}\,\mathbf 1_{j\leq i},\qquad
\rho(M)=\max_i {\alpha r_i\over a_i}.
\]

Componentwise (0\leq F(w)\leq Mw\). Hence, for any fixed finite grid with
\(\rho(M)<1\), (0\leq w_t\leq M^tw_0\to0\). This conclusion does not require
immediate contraction in an ordinary norm and remains valid with nonnormal
transients or repeated eigenvalues. For a unique accessible dominant diagonal,
the normalized asymptotic shape is the corresponding positive triangular
eigenvector. Repeated dominance does not yield that unique-shape conclusion;
near degeneracy produces the observed long transients.

The nonlinear correction is exactly

\[
(Mw)_i-F_i(w)={(Mw)_i^2\over1+(Mw)_i}\leq(Mw)_i^2.
\]

Each batch records its final one-step correction in log range. No accumulated
future bound is claimed: the 1024-update state is not used as an error-bounded
approximation to the asymptote. Fixed-grid coefficients come from the exact
triangular eigenvector instead. The saved 1024-update coefficient differences
from that eigenvector remain finite-trajectory diagnostics only.

The scalar solution provides a positive control. With
\(d=\alpha r/a\),

\[
w_t={d^t\over 1/w_0+d(1-d^t)/(1-d)}\quad(d\ne1),\qquad
w_t={w_0\over1+tw_0}\quad(d=1),
\]

while (A=1/(Lr)) and (P_{\rm pixel}=\Delta_vv/(Lr)). Vanishing amplitude
therefore does not imply divergent noise.

## Continuum measure and limit ordering

The frozen recipe represents the magnitude measure

\[
dR(m)=\rho_v(m)\,dm,qquad
\rho_v(m)=\rho_z(m){1+z_{\rm source}\over c},
\]

by positive Gauss--Legendre cell masses
\(r_i=\rho_v(m_i)q_i\) on a bounded, ordered magnitude support. The saved
variance (v(m_i)) is strictly positive. Every diagnostic batch is reconstructed
from the immutable diagnosis NPZ plus the accuracy report's (L,\Delta_v,S,B).
Four profile matches use order 32 (bins 0--1, both fields) and eight use order 64
(bins 2--5); all sixty have a separate diagnosis-source comparison.

Revision 2 establishes the iteration-first limit on each fixed grid. It does not
establish the refinement-first limit at fixed iteration or a joint limit. The
five available quadratures are finite and non-nested, so they do not verify a
uniform quadrature family or justify exchange of these limits.

For the linearized continuum control

\[
(Kf)(x)=c\int_0^x f(y)\,dy,
\]

starting from (f=1), the normalized shape is (x^t) and

\[
A_t={(t+1)^2\over(2t+1)LR},\qquad
P_{{\rm pixel},t}=\Delta_vvA_t.
\]

Both diverge linearly. This is a theorem for the constant-variance linearized
Volterra control, not for the saved nonlinear recurrence; nonconstant variance
already breaks its proportional (P_{\rm pixel}/A) relation.

A sufficient iteration-first concentration argument can be written from the
backward eigenvector relation

\[
C_i=\prod_{j>i}(1-d_j/\lambda)
\leq\exp[-\sum_{j>i}d_j/\lambda].
\]

If a consistently refined non-atomic measure has
\(\max_i r_i/v_i\to0\), every fixed terminal interval retains positive
\(\int dR/v\), and the limiting probability puts a positive mass
\(p_{\rm tail}\) on a tail whose source measure (R_{\rm tail}\to0\), then

\[
LA\geq {p_{\rm tail}^2\over R_{\rm tail}}.
\]

A positive lower variance bound similarly gives a lower bound on
\(P_{\rm pixel}\). The saved finite orders do not verify these asymptotic
assumptions. In particular, the actual faintest-10% tail probabilities of the
fixed-grid distributions are only (1.87\times10^{-41}) to
\(4.98\times10^{-29}); the largest Jacobian diagonal is neither their mean nor
their concentration location. The fixed-grid means span 22.519--24.358 and the
widths 0.0478--0.497. These measured distributions support an unresolved verdict,
not the revision-1 boundary-concentration claim.

## Saved populations

The table gives fixed-grid eigenvector coefficients (A/P_{\rm pixel}) at the
three finest saved orders. Units are deg\(^2\) / deg\(^2\) km s\(^{-1}\).

| Bin | Population | Order 16 | Order 32 | Order 64 | Outcome |
|---:|:---|---:|---:|---:|:---|
| 0 | lya(qso) | 0.0374662 / 12.6340 | 0.0572858 / 24.7726 | 0.0971897 / 49.0883 | unresolved |
| 0 | lya(lbg) | 2.34876 / 7671.75 | 4.39865 / 14849.5 | 8.26859 / 29171.9 | unresolved |
| 1 | lya(qso) | 0.0516760 / 9.66180 | 0.0755132 / 18.9483 | 0.123829 / 37.5533 | unresolved |
| 1 | lya(lbg) | 0.249034 / 440.670 | 0.465054 / 852.387 | 0.870604 / 1673.26 | unresolved |
| 2 | lya(qso) | 0.0700742 / 17.2293 | 0.101165 / 33.8050 | 0.165005 / 67.0146 | unresolved |
| 2 | lya(lbg) | 0.0868774 / 167.482 | 0.162333 / 324.000 | 0.304206 / 636.159 | unresolved |
| 3 | lya(qso) | 0.114856 / 18.9915 | 0.160518 / 37.2486 | 0.254435 / 73.8281 | unresolved |
| 3 | lya(lbg) | 0.122153 / 155.152 | 0.227660 / 300.048 | 0.424282 / 588.357 | unresolved |
| 4 | lya(qso) | 0.228518 / 24.4476 | 0.312158 / 47.9641 | 0.482309 / 95.0825 | unresolved |
| 4 | lya(lbg) | 0.509868 / 474.238 | 0.950585 / 917.203 | 1.77301 / 1798.96 | unresolved |
| 5 | lya(qso) | 0.815482 / 47.5971 | 0.931840 / 93.3102 | 1.18096 / 184.934 | unresolved |
| 5 | lya(lbg) | 0.479528 / 390.957 | 0.894049 / 756.113 | 1.66781 / 1483.05 | unresolved |

All 24 signed order-doubling ratios for (P_{\rm pixel}) are positive and span
0.934--0.982. The three-order exponent is 0.9615--0.9798, while the effective-
measure exponent is -0.9079 to -0.2671. These are strong empirical refinement
trends and fail the (10^{-4}) screening target; they do not determine the
continuum theorem. A bounded sequence with the same long rising transient and a
false late plateau are explicit counterexamples to either inference from a
screen alone.

At fixed geometry, P1D and response,

\[
N_F(k,\mu)=\left[A P_{1D}(k_\parallel)W^2(k_\parallel)+P_{\rm pixel}\right]
{d_{\rm deg}^2\over a_v}.
\]

The observed positive refinement trend raises the saved fixed-grid additive
noise term, but an adopted coefficient prescription and a separately authorized
forecast are required to quantify Fisher effects. The six Step 12 accuracy bins
and twelve unavailable diagnostics remain unchanged.

## Evidence integrity

The source diagnosis JSON SHA256 is
`069a713b0a5df67556aed992374cb8a8ef5a0c9201ef40919f3137de634b1008`;
the profile manifest SHA256 is
`70edc87acf3c64229a144793c63528632a55b5e1763e7993d163083c71684474`.
All 60 batches and 840 checkpoints completed; maximum trajectory time was
0.801 s. Historical finite-step discrepancies are at most
\(9.49\times10^{-14}\) in (A) and \(1.22\times10^{-13}\) in
\(P_{\rm pixel}\).

The finalizer reconstructs each source identity, array, scalar, attempted scope,
historical comparison, fixed-grid result and summary operand. Offline validation
rebuilds all final tables and SVGs byte-for-byte. It rejects swapped fields and
orders, self-consistent replacement arrays with stale provenance, an omitted
checkpoint/cap, falsified attempt status, altered historical operands and an
altered final conclusion. Complete, explicitly failed/capped, and identical-input
controls remain inspectable.

The immutable numerical bundle and separate-population figures are under
`.validation/step13-r2-20260915T010410Z/`. No real model, reader or forecast was
run, and no production method was changed.
