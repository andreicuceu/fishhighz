"""Reproduce the idealized full-EFT Fisher study on a single local CPU.

Run with both repository roots importable. No remote execution, CAMB solve or
external survey files are required. Outputs include means, Jacobians, Fisher
matrices, rank diagnostics and dimensionless nuisance PCA modes.
"""

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from eftkernels.utils import read_template

from fishhighz.covariance import gaussian_covariance
from fishhighz.derivatives import evaluate_derivatives
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.fisher import factor_covariance, fisher_from_factors
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.eft import EFTModelAdapter
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.parameters import Parameter, ParameterRegistry

TARGETS = ("alpha_w", "phi_w", "alpha_s", "phi_s", "f")
KMAX = [0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0]
PAIR_NAMES = ("FF", "FG", "GG")


def configuration():
    d = read_template()
    d = {
        key: np.asarray(value, dtype=float) if isinstance(value, np.ndarray) else value
        for key, value in d.items()
    }
    values = dict(
        alpha_w=1.0,
        phi_w=1.0,
        alpha_s=1.0,
        phi_s=1.0,
        f=d["growth_rate"],
        bias_F=-0.14,
        bias_eta_F=-0.21,
        bias_G=2.4,
    )
    scales = {}
    for name, value in dict(
        b2=0.05,
        bG2=-0.02,
        bEta2=0.03,
        bdeltaeta=0.08,
        bKKPar=-0.02,
        bPi2Par=0.02,
        bPi3Par=0.0,
        bdeltaPi2Par=0.0,
        bKPi2Par=0.0,
        bEtaPi2Par=0.0,
        bGamma3=0.0,
    ).items():
        values[name + "_F"] = value / 0.14
        scales[name + "_F"] = 0.14
    for name, value in dict(b2=0.5, bG2=-0.3, bGamma3=0.2).items():
        values[name + "_G"] = value
        scales[name + "_G"] = 1.0
    for tracer, scale in [("F", 0.14), ("G", 2.4)]:
        for coefficient in ("cs", "cs2", "cs4"):
            values[coefficient + "_" + tracer] = 0.0
            scales[coefficient + "_" + tracer] = scale
    values["cx"] = 0.0
    scales["cx"] = 1.0
    # Residual stochastic amplitudes: no known observational noise in the mean.
    for pair, scale in [("FF", 100.0), ("FG", np.sqrt(100.0 * 5000.0)), ("GG", 5000.0)]:
        for suffix in ["0", "k2", "mu2"]:
            name = f"stoch_{pair}_{suffix}"
            values[name] = 0.0
            scales[name] = scale
    with fits.open(
        Path(__import__("eftkernels").__file__).parent / "templates/Planck18_z2.3.fits"
    ) as hd:
        sigma8 = float(hd[1].header["SIGMA8_ZREF"])
    return d, values, scales, sigma8


def identifiable_errors(design, tolerance=1e-9):
    """Give infinite coordinate errors whenever the null space contains it.

    SVD is applied to column-normalized, whitened mean derivatives. Solves on
    the identifiable subspace supply variances only for estimable coordinates;
    a null direction never becomes a finite error via a pseudoinverse.
    """
    norm = np.linalg.norm(design, axis=0)
    norm = np.where(norm > 0, norm, 1.0)
    _, s, vh = np.linalg.svd(design / norm, full_matrices=False)
    rank = int(np.count_nonzero(s > tolerance * s[0])) if s.size else 0
    null = vh[rank:] / norm[None, :]
    if len(null):
        null /= np.linalg.norm(null, axis=1)[:, None]
    errors = np.sqrt(np.sum((vh[:rank] / s[:rank, None]) ** 2, axis=0)) / norm
    affected = np.linalg.norm(vh[rank:], axis=0) > 1e-6
    errors[affected] = np.inf
    return errors, rank, null, s


def marginal_pca(design, names):
    """Marginalize targets and linear biases before dimensionless higher EFT PCA."""
    higher = [
        i
        for i, n in enumerate(names)
        if n not in TARGETS and n not in ("bias_F", "bias_eta_F", "bias_G")
    ]
    other = [i for i in range(len(names)) if i not in higher]
    target = design[:, other]
    u, s, _ = np.linalg.svd(target, full_matrices=False)
    u = u[:, s > s[0] * 1e-9]
    nuisance = design[:, higher]
    nuisance = nuisance - u @ (u.T @ nuisance)
    _, singular, vh = np.linalg.svd(nuisance, full_matrices=False)
    modes = []
    for sval, vec in zip(singular, vh):
        indices = np.argsort(np.abs(vec))[-5:][::-1]
        modes.append(
            dict(
                sigma=float(1 / sval) if sval > singular[0] * 1e-9 else None,
                terms={names[higher[i]]: float(vec[i]) for i in indices},
            )
        )
    return modes


def run(args):
    t = time.perf_counter()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    d, values, scales, sigma8 = configuration()
    names = list(values)
    grid = gauss_legendre_grid(
        [0.02] + KMAX, k_order=args.k_order, mu_order=args.mu_order, h_fid=d["H0"] / 100
    )
    fields = [
        ObservedField("F", "forest", "F", "idealized"),
        ObservedField("G", "galaxy", "G"),
    ]
    all_selection = PairSelection(fields)
    loop = EFTModelAdapter(
        d,
        names,
        values,
        scales,
        num_mu=args.eft_mu,
        nmax=args.nmax,
        k_output=np.geomspace(0.005, 1.12, args.output_k),
    )
    k, mu = grid.k_flat, grid.mu_flat
    fid = np.array([values[n] for n in names])
    signal = loop(fid, d["z_fiducial"], k, mu, all_selection.required_pairs)
    noise = np.tile([100.0, 0.0, 5000.0], (len(k), 1))
    total = signal + noise
    field_corr = total[:, 1] / np.sqrt(total[:, 0] * total[:, 2])
    if np.any(total[:, [0, 2]] <= 0) or np.any(np.abs(field_corr) >= 1):
        raise ValueError(
            "Nonpositive fiducial field covariance; no clipping is allowed"
        )
    summaries = []
    arrays = dict(k=k, mu=mu, q_mode=grid.q_mode, signal=signal, noise=noise)
    for order in ["tree", "loop"]:
        order_names = names[:8] if order == "tree" else names
        registry = ParameterRegistry(
            [
                Parameter(
                    n,
                    values[n],
                    "target" if n in TARGETS else "nuisance",
                    step=(2e-4 if n in TARGETS[:4] else 1e-3) * args.step_scale,
                )
                for n in order_names
            ]
        )
        model = (
            loop
            if order == "loop"
            else EFTModelAdapter(
                d,
                order_names,
                values,
                scales,
                order="tree",
                num_mu=args.eft_mu,
                nmax=args.nmax,
                k_output=np.geomspace(0.005, 1.12, args.output_k),
            )
        )
        binding = BoundParameters(registry, order_names, {n: n for n in order_names})
        prepared = PreparedP3D(
            registry,
            all_selection,
            [P3DProvider("eft", model, binding, all_selection.required_pairs)],
        )
        deriv = evaluate_derivatives(
            prepared, registry.fiducials, d["z_fiducial"], k, mu
        )
        arrays[order + "_power"] = deriv.power
        arrays[order + "_jacobian"] = deriv.jacobian
        print(order, "derivatives completed", time.perf_counter() - t, flush=True)
        for data, selected in [
            ("FF", [("F", "F")]),
            ("GG", [("G", "G")]),
            ("3x2", None),
        ]:
            selection = PairSelection(fields, selected=selected)
            all_pairs = list(map(tuple, all_selection.required_pairs))
            required = [all_pairs.index(tuple(p)) for p in selection.required_pairs]
            chosen = [all_pairs.index(tuple(p)) for p in selection.selected_pairs]
            active = []
            for i, n in enumerate(order_names):
                relevant = (
                    n in TARGETS
                    or data == "3x2"
                    or data == "FF"
                    and (n.endswith("_F") or n.startswith("stoch_FF_"))
                    or data == "GG"
                    and (n.endswith("_G") or n.startswith("stoch_GG_"))
                )
                if relevant:
                    active.append(i)
            local_names = [order_names[i] for i in active]
            covariance = gaussian_covariance(
                total[:, required], 20e9 * grid.q_mode, selection
            )
            factors = factor_covariance(covariance)
            jac = deriv.jacobian[:, chosen, :][:, :, active]
            white = np.linalg.solve(factors, jac)
            for cut in KMAX:
                mask = k <= cut
                design = white[mask].reshape(-1, len(active))
                information = fisher_from_factors(jac[mask], factors[mask])
                np.testing.assert_allclose(
                    information, design.T @ design, rtol=2e-12, atol=1e-7
                )
                for prior in [None] if order == "tree" else [None, 0.3, 1.0, 3.0]:
                    augmented = design
                    prior_rows = []
                    if prior is not None:
                        # Bias/counterterm priors only. Residual stochasticity stays free.
                        for j, n in enumerate(local_names):
                            if (
                                j >= 5
                                and n not in ("bias_F", "bias_eta_F", "bias_G")
                                and not n.startswith("stoch_")
                            ):
                                row = np.zeros(len(active))
                                row[j] = 1 / prior
                                prior_rows.append(row)
                        augmented = np.vstack((design, prior_rows))
                    errors, rank, null, s = identifiable_errors(augmented)
                    record = dict(
                        order=order,
                        data=data,
                        kmax=cut,
                        prior=prior,
                        names=local_names,
                        errors=[float(e) if np.isfinite(e) else None for e in errors],
                        rank=rank,
                        nparameters=len(active),
                        singular_values=s.tolist(),
                        null_directions=null.tolist(),
                        fsigma8_error=float(errors[4] * sigma8)
                        if np.isfinite(errors[4])
                        else None,
                    )
                    if prior is None:
                        record["pca"] = (
                            marginal_pca(design, local_names) if len(active) > 5 else []
                        )
                        arrays[f"{order}_{data}_{cut}_fisher"] = information
                    summaries.append(record)
        np.savez_compressed(output / "arrays.npz", **arrays)
        (output / "forecasts.json").write_text(json.dumps(summaries, indent=2))
    metadata = dict(
        volume_mpc_h3=20e9,
        z=d["z_fiducial"],
        sigma8=sigma8,
        f=d["growth_rate"],
        fsigma8=d["growth_rate"] * sigma8,
        noise_F=100.0,
        galaxy_density=2e-4,
        defaults=values,
        scales=scales,
        settings=vars(args),
        python=platform.python_version(),
        max_field_correlation=float(np.max(np.abs(field_corr))),
        seconds=time.perf_counter() - t,
    )
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2))
    print("completed", metadata["seconds"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/main")
    parser.add_argument("--k-order", type=int, default=14)
    parser.add_argument("--mu-order", type=int, default=20)
    parser.add_argument("--eft-mu", type=int, default=64)
    parser.add_argument("--output-k", type=int, default=420)
    parser.add_argument("--nmax", type=int, default=256)
    parser.add_argument("--step-scale", type=float, default=1.0)
    run(parser.parse_args())
