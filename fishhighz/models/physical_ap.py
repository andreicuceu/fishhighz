"""Physical-distance AP mapping in fixed fiducial-observed coordinates."""

import numpy as np


def physical_ap_coordinates(k_observed, mu_observed, q_perp, q_parallel):
    """Return intrinsic (k, mu) and power-volume Jacobian.

    k_observed is in fixed h_fid/Mpc, mu dimensionless, q_perp=DM/DMfid and
    q_parallel=Hfid/H. Both dilation factors are positive dimensionless
    scalars. Inputs broadcast; the output k and mu preserve that shape.
    The Jacobian is 1/(q_perp**2*q_parallel). No sound horizon enters these
    factors when the cosmological linear spectrum is regenerated.
    """
    if not all(
        np.ndim(q) == 0 and np.isfinite(q) and q > 0 for q in (q_perp, q_parallel)
    ):
        raise ValueError("AP factors must be finite positive scalars")
    k, mu = np.broadcast_arrays(
        np.asarray(k_observed, dtype=np.float64),
        np.asarray(mu_observed, dtype=np.float64),
    )
    if (
        np.any(~np.isfinite(k))
        or np.any(k <= 0)
        or np.any(~np.isfinite(mu))
        or np.any(abs(mu) > 1)
    ):
        raise ValueError("require positive finite k and mu in [-1,1]")
    parallel = k * mu / q_parallel
    transverse = k * np.sqrt(np.maximum(0.0, 1 - mu * mu)) / q_perp
    intrinsic = np.hypot(parallel, transverse)
    return intrinsic, parallel / intrinsic, 1.0 / (q_perp * q_perp * q_parallel)
