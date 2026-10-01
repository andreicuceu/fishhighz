"""Supplied-array causal swaps on common nodes, with exact chain endpoints.

These diagnostics isolate measured arrays, not independent physical parameters.
The noise/signal split of legacy pair inputs is not generally field representable;
its residual is explicitly retained instead of inventing an overlap model.
"""

import numpy as np

from .numerics import contract, relative, summaries, wick


def interpolate(values, k, mu, target_k, target_mu):
    """Linear k/mu interpolation with explicit endpoint extrapolation in mu.

    Parameters
    ----------
    values : array_like, shape (n_mu * n_k, ...)
        Values in mu-major, k-minor order.
    k : array_like, shape (n_k,)
        Increasing source wavenumbers in h/Mpc.
    mu : array_like, shape (n_mu,)
        Increasing source direction cosines.
    target_k : array_like, shape (n_target,)
        Paired target wavenumbers in h/Mpc.
    target_mu : array_like, shape (n_target,)
        Paired target direction cosines.

    Returns
    -------
    interpolated : ndarray, shape (n_target, ...)
        Interpolated values with original trailing scientific axes and units.
    """
    values = np.asarray(values)
    tail = values.shape[1:]
    grid = values.reshape(len(mu), len(k), -1)
    along = np.stack(
        [
            np.stack(
                [np.interp(target_k, k, row[:, j]) for j in range(grid.shape[2])],
                axis=1,
            )
            for row in grid
        ]
    )
    index = np.clip(np.searchsorted(mu, target_mu) - 1, 0, len(mu) - 2)
    mu_fraction = (target_mu - mu[index]) / (mu[index + 1] - mu[index])
    nodes = np.arange(len(target_k))
    return (
        (1 - mu_fraction[:, None]) * along[index, nodes]
        + mu_fraction[:, None] * along[index + 1, nodes]
    ).reshape((len(target_k),) + tail)


def chain(legacy, accuracy, task, legacy_volume, accuracy_volume):
    """Change grid, volume, total power, then J; store interaction/order caveat.

    Parameters
    ----------
    legacy : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    accuracy : dict of str to ndarray
        Numerical evidence arrays; Fourier-cell axes and pair order follow the
        declared task. Powers use (Mpc/h)^3 and volumes use (Mpc/h)^3 unless
        separately labeled.
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    legacy_volume : float
        Compatibility volume in (Mpc/h)^3.
    accuracy_volume : float
        Accuracy volume in (Mpc/h)^3.

    Returns
    -------
    arrays : dict of str to ndarray
        Total powers, derivatives, modes and Fisher summaries at every
        substitution stage.
    report : dict
        Stage descriptions, endpoint discrepancy and order-dependence
        limitation.
    """
    target_k_grid, target_mu_grid = accuracy["k"], accuracy["mu"]
    legacy_k_grid, legacy_mu_grid = np.unique(legacy["k"]), np.unique(legacy["mu"])
    total = interpolate(
        legacy["total"], legacy_k_grid, legacy_mu_grid, target_k_grid, target_mu_grid
    )
    observed_jacobian = interpolate(
        legacy["observed_j"],
        legacy_k_grid,
        legacy_mu_grid,
        target_k_grid,
        target_mu_grid,
    )
    modes = accuracy["modes"] * legacy_volume / accuracy_volume
    required = np.asarray(task["required_pairs"])
    selected = np.asarray(task["selected_pairs"])
    n_fields = len(task["fields"])
    stages = []
    arrays = {}

    def save(name, t, j, n, description):
        """Record one ordered substitution and its independent Fisher contraction.

        Parameters
        ----------
        name : str
            Quantity or record label used in diagnostics.
        t : array_like, shape (n_cell, n_required_pair)
            Signed total field-pair powers, including noise, in (Mpc/h)^3.
        j : ndarray, shape (n_cell, n_selected_pair, n_parameter)
            Observed mean-spectrum derivatives, including field responses, in power
            units per parameter unit.
        n : array_like, shape (n_cell,)
            Fourier-mode counts for this substitution.
        description : str
            Scientific interpretation of the changed operands.

        Notes
        -----
        Appends to the enclosing stage list and numerical-array mapping.
        """
        covariance = wick(t, n, required, selected, n_fields)
        joint_fisher, individual_fisher = contract(covariance, j, independent=True)
        index = len(stages)
        arrays.update(
            {
                f"total_{index}": t,
                f"jacobian_{index}": j,
                f"modes_{index}": n,
                **{
                    f"{key}_{index}": value
                    for key, value in summaries(joint_fisher, individual_fisher).items()
                },
            }
        )
        stages.append(dict(name=name, description=description, index=index))

    save(
        "compatibility",
        legacy["total"],
        legacy["observed_j"],
        legacy["modes"],
        "Captured legacy nodes and inputs; independent contraction",
    )
    save(
        "grid",
        total,
        observed_jacobian,
        modes,
        "Change only quadrature; linearly interpolate legacy T/J, extrapolate mu endpoints; interpolation is included in this diagnostic",
    )
    save(
        "volume",
        total,
        observed_jacobian,
        accuracy["modes"],
        "Change only integrated volume/modes; hold common-node T/J fixed",
    )
    save(
        "total_power",
        accuracy["total"],
        observed_jacobian,
        accuracy["modes"],
        "Replace full required total power, including model/response/pair-specific versus per-field noise; hold J/modes fixed",
    )
    save(
        "full_jacobian",
        accuracy["total"],
        accuracy["observed_j"],
        accuracy["modes"],
        "Replace observed J, including peak/full-wiggle mapping and response; hold T/modes fixed",
    )
    endpoint = relative(arrays["fisher_4"], accuracy["fisher"])
    return arrays, dict(
        context=task,
        stages=stages,
        endpoint_relative=endpoint,
        complete=bool(endpoint < 5e-12),
        limitation="Array-level sequential attribution is order dependent. Total-power and Jacobian groups combine physical changes; this chain alone does not isolate template, redshift, reconstruction and noise causes.",
        unresolved_physical_attribution=True,
    )
