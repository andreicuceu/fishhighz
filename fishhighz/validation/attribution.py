"""Supplied-array causal swaps on common nodes, with exact chain endpoints.

These diagnostics isolate measured arrays, not independent physical parameters.
The noise/signal split of legacy pair inputs is not generally field representable;
its residual is explicitly retained instead of inventing an overlap model.
"""

import numpy as np

from .numerics import contract, relative, summaries, wick


def interpolate(values, k, mu, target_k, target_mu):
    """Linear k/mu interpolation with explicit endpoint extrapolation in mu."""
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
    frac = (target_mu - mu[index]) / (mu[index + 1] - mu[index])
    nodes = np.arange(len(target_k))
    return (
        (1 - frac[:, None]) * along[index, nodes]
        + frac[:, None] * along[index + 1, nodes]
    ).reshape((len(target_k),) + tail)


def chain(legacy, accuracy, task, legacy_volume, accuracy_volume):
    """Change grid, volume, total power, then J; store interaction/order caveat."""
    tk, tm = accuracy["k"], accuracy["mu"]
    lk, lm = np.unique(legacy["k"]), np.unique(legacy["mu"])
    total = interpolate(legacy["total"], lk, lm, tk, tm)
    jac = interpolate(legacy["observed_j"], lk, lm, tk, tm)
    modes = accuracy["modes"] * legacy_volume / accuracy_volume
    required = np.asarray(task["required_pairs"])
    selected = np.asarray(task["selected_pairs"])
    nf = len(task["fields"])
    stages = []
    arrays = {}

    def save(name, t, j, n, description):
        c = wick(t, n, required, selected, nf)
        f, p = contract(c, j, independent=True)
        index = len(stages)
        arrays.update(
            {
                f"total_{index}": t,
                f"jacobian_{index}": j,
                f"modes_{index}": n,
                **{f"{key}_{index}": value for key, value in summaries(f, p).items()},
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
        jac,
        modes,
        "Change only quadrature; linearly interpolate legacy T/J, extrapolate mu endpoints; interpolation is included in this diagnostic",
    )
    save(
        "volume",
        total,
        jac,
        accuracy["modes"],
        "Change only integrated volume/modes; hold common-node T/J fixed",
    )
    save(
        "total_power",
        accuracy["total"],
        jac,
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
