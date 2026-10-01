"""Array-only intrinsic Kaiser/scaling/BAO operations; no preparation libraries."""

import numpy as np


def _resolve(fixed, destinations, sources, theta):
    """Gather free parameters into a copy of the fixed model settings.

    Parameters
    ----------
    fixed : ndarray of shape (n_setting,)
        Fiducial settings, with units determined by each model slot.
    destinations, sources : ndarray of shape (n_free_slot,)
        Integer destination settings and corresponding source parameter indices.
    theta : ndarray of shape (n_parameter,)
        Parameter values in the units of their destination settings.

    Returns
    -------
    values : ndarray of shape (n_setting,)
        Owned settings array with free slots replaced.
    """
    values = fixed.copy()
    values[destinations] = theta[sources]
    return values


def _scales(values, basis):
    """Convert dilation parameters to parallel and transverse scale factors.

    Parameters
    ----------
    values : ndarray of shape (2,)
        Dimensionless dilation parameters in the selected basis.
    basis : int
        Basis code: 0 for (ap, at), 1 for (alpha, phi), 2 for
        (alpha_iso, epsilon), and otherwise (alpha_iso, phi).

    Returns
    -------
    parallel_scale, transverse_scale : float
        Dimensionless parallel and transverse dilation factors.
    volume_factor : float
        Inverse dilation volume, 1 / (ap * at**2), evaluated logarithmically.

    Notes
    -----
    The alpha/phi basis uses ap=alpha/sqrt(phi), at=alpha*sqrt(phi).
    The alpha_iso/phi basis uses ap=alpha_iso*phi**(-2/3),
    at=alpha_iso*phi**(1/3). Input-domain validation belongs to the caller.
    """
    first_scale, second_scale = values
    if basis == 0:
        parallel_scale, transverse_scale = first_scale, second_scale
    elif basis == 1:
        parallel_scale, transverse_scale = (
            first_scale / np.sqrt(second_scale),
            first_scale * np.sqrt(second_scale),
        )
    elif basis == 2:
        parallel_scale, transverse_scale = (
            (first_scale * (1 + second_scale)) * (1 + second_scale),
            first_scale / (1 + second_scale),
        )
    else:
        parallel_scale, transverse_scale = (
            first_scale * second_scale ** (-2 / 3),
            first_scale * second_scale ** (1 / 3),
        )
    volume_factor = np.exp(-np.log(parallel_scale) - 2 * np.log(transverse_scale))
    return parallel_scale, transverse_scale, volume_factor


def _coordinates(k, mu, ap, at):
    """Map observed Fourier nodes to the dilated template coordinates.

    Parameters
    ----------
    k : ndarray of shape (n_node,)
        Observed wavenumbers in h/Mpc.
    mu : ndarray of shape (n_node,)
        Dimensionless line-of-sight direction cosines.
    ap, at : float
        Dimensionless parallel and transverse dilation factors.

    Returns
    -------
    mapped_k : ndarray of shape (n_node,)
        Dilated wavenumbers in h/Mpc.
    mapped_mu : ndarray of shape (n_node,)
        Dilated direction cosines.
    parallel, transverse : ndarray of shape (n_node,)
        Dilated parallel and transverse wavenumbers in h/Mpc.

    Notes
    -----
    An isotropic dilation preserves mu exactly and avoids roundoff from
    reconstructing k through its parallel and transverse components.
    """
    parallel = k * mu / ap
    transverse = k * np.sqrt(1 - mu**2) / at
    if ap == at:
        mapped_k, mapped_mu = k / ap, mu.copy()
    else:
        mapped_k = np.hypot(parallel, transverse)
        mapped_mu = parallel / mapped_k
    return mapped_k, mapped_mu, parallel, transverse


def _field_factors(mu, biases, betas, growth_rate, forest):
    """Evaluate signed Kaiser factors for forests and discrete tracers.

    Parameters
    ----------
    mu : ndarray of shape (n_node,)
        Dimensionless direction cosines in the relevant dilation coordinates.
    biases, betas : ndarray of shape (n_field,)
        Dimensionless tracer biases and forest redshift-space distortion factors.
    growth_rate : float
        Dimensionless logarithmic growth rate used for discrete tracers.
    forest : ndarray of bool, shape (n_field,)
        Mask selecting forest fields.

    Returns
    -------
    factors : ndarray of shape (n_node, n_field)
        Forest b*(1+beta*mu**2) or galaxy b+f*mu**2. The latter expression
        remains defined at zero galaxy bias.
    """
    mu_squared = mu[:, None] ** 2
    return np.where(
        forest[None, :],
        biases * (1 + betas * mu_squared),
        biases + growth_rate * mu_squared,
    )


def _damping(parallel, transverse, widths_squared, pairs):
    """Evaluate Gaussian BAO damping with fixed pair broadening widths.

    Parameters
    ----------
    parallel, transverse : ndarray of shape (n_node,)
        Parallel and transverse wavenumbers in h/Mpc.
    widths_squared : ndarray of shape (n_field, 2)
        Squared parallel and transverse broadening lengths in (Mpc/h)^2.
    pairs : ndarray of int, shape (n_pair, 2)
        Field indices of the requested spectra.

    Returns
    -------
    damping, exponent : ndarray of shape (n_node, n_pair)
        Dimensionless exp(-exponent) and its nonnegative exponent.

    Notes
    -----
    Each squared cross width is the arithmetic mean of the two squared auto
    widths. Widths are fixed inputs rather than functions of varied growth.
    """
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
    """Combine smooth and damped wiggle spectra for the requested pairs.

    Parameters
    ----------
    smooth, wiggle : ndarray of shape (n_node,)
        Smooth and signed wiggle template powers in (Mpc/h)^3.
    factors_smooth, factors_wiggle : ndarray of shape (n_node, n_field)
        Dimensionless signed Kaiser factors for each component.
    damping : ndarray of shape (n_node, n_pair)
        Dimensionless BAO damping, applied only to the wiggle contribution.
    pairs : ndarray of int, shape (n_pair, 2)
        Indices of the fields forming each spectrum.
    q_smooth, q_wiggle : float
        Dimensionless inverse dilation-volume factors for each component.
    growth : float
        Dimensionless power growth relative to the template redshift.

    Returns
    -------
    power : ndarray of shape (n_node, n_pair)
        Intrinsic pair power in (Mpc/h)^3, in the requested pair order.
    """
    i, j = pairs.T
    return growth * (
        q_smooth * factors_smooth[:, i] * factors_smooth[:, j] * smooth[:, None]
        + q_wiggle
        * factors_wiggle[:, i]
        * factors_wiggle[:, j]
        * damping
        * wiggle[:, None]
    )
