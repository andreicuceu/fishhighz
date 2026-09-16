"""Array-only intrinsic Kaiser/scaling/BAO operations; no preparation libraries."""

import numpy as np


def _resolve(fixed, destinations, sources, theta):
    """Gather explicit free slots into an owned numeric settings vector."""
    values = fixed.copy()
    values[destinations] = theta[sources]
    return values


def _scales(values, basis):
    """Return ap,at,Q; basis codes 0=ap_at, 1=alpha_phi, 2=alpha_iso_epsilon."""
    a, b = values
    if basis == 0:
        ap, at = a, b
    elif basis == 1:
        ap, at = a / np.sqrt(b), a * np.sqrt(b)
    else:
        ap, at = (a * (1 + b)) * (1 + b), a / (1 + b)
    q = np.exp(-np.log(ap) - 2 * np.log(at))
    return ap, at, q


def _coordinates(k, mu, ap, at):
    """Inverse Fourier mapping with exact isotropic/identity radial coordinates."""
    parallel = k * mu / ap
    transverse = k * np.sqrt(1 - mu**2) / at
    if ap == at:
        mapped_k, mapped_mu = k / ap, mu.copy()
    else:
        mapped_k = np.hypot(parallel, transverse)
        mapped_mu = parallel / mapped_k
    return mapped_k, mapped_mu, parallel, transverse


def _field_factors(mu, biases, betas, growth_rate, forest):
    """Signed factors, without division by galaxy bias (including b=0)."""
    mu2 = mu[:, None] ** 2
    return np.where(
        forest[None, :], biases * (1 + betas * mu2), biases + growth_rate * mu2
    )


def _damping(parallel, transverse, widths_squared, pairs):
    """Fixed cross widths squared are the arithmetic mean of two auto squares."""
    pair_widths = 0.5 * widths_squared[pairs[:, 0]] + 0.5 * widths_squared[pairs[:, 1]]
    exponent = (
        (parallel[:, None] * np.sqrt(pair_widths[None, :, 0])) ** 2
        + (transverse[:, None] * np.sqrt(pair_widths[None, :, 1])) ** 2
    ) / 2
    return np.exp(-exponent), exponent


def _assemble(
    smooth,
    wiggle,
    factors_smooth,
    factors_wiggle,
    damping,
    pairs,
    q_smooth,
    q_wiggle,
    growth,
):
    """Assemble all requested pairs from two shared template/factor evaluations."""
    i, j = pairs.T
    return growth * (
        q_smooth * factors_smooth[:, i] * factors_smooth[:, j] * smooth[:, None]
        + q_wiggle
        * factors_wiggle[:, i]
        * factors_wiggle[:, j]
        * damping
        * wiggle[:, None]
    )
