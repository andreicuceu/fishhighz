"""Fixed-survey preparation and independent-bin Fisher execution (NumPy only)."""

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from ._arrays import integer
from .covariance import combine_observed_power, gaussian_covariance
from .derivatives import CallCount, _schedule, evaluate_derivatives
from .fisher import factor_covariance, fisher_from_factors
from .geometry import _immutable, mode_counts, wavenumber_comoving_to_velocity
from .models.external import evaluate_p1d, evaluate_p3d
from .noise import forest_noise, galaxy_noise, prepare_noise
from .response import InstrumentResponse, pair_response, prepare_response
from .results import FisherResult, combine_results
from .survey import BinSpec, ForestInput, PreparedBin, freeze, snapshot_p3d
from .weights import prepare_forest_weights, sample_auxiliary


def _validate_spec(spec):
    if not isinstance(spec, BinSpec):
        raise ValueError("require BinSpec")
    from .models.kaiser import KaiserModel

    fields = spec.p3d.selection.fields
    for route in spec.p3d.routes:
        model = route.provider.model
        if isinstance(model, KaiserModel):
            if (
                model.fields != fields
                or model.template.h_fid != spec.geometry.h_fid
                or model.z != spec.geometry.z_eval
            ):
                raise ValueError(
                    f"{route.provider.label}: built-in field/z/h_fid context mismatch"
                )
    if (
        spec.responses is None
        or set(spec.responses) != {f.id for f in fields}
        or any(not isinstance(r, InstrumentResponse) for r in spec.responses.values())
    ):
        raise ValueError(
            "responses must cover exactly all fields with InstrumentResponse"
        )
    if spec.full_noise is not None:
        if any(
            x is not None
            for x in (spec.forests, spec.galaxies, spec.independent_sampling)
        ):
            raise ValueError("full noise conflicts with generated inputs")
    else:
        if spec.independent_sampling is not True:
            raise ValueError(
                "generated noise requires explicit independent_sampling=True"
            )
        active = [fields[i] for i in np.unique(spec.p3d.selection.selected_pairs)]
        for kind, data in [("forest", spec.forests), ("galaxy", spec.galaxies)]:
            expected = {f.id for f in active if f.kind == kind}
            if set(data or {}) != expected:
                raise ValueError(
                    f"{kind} noise inputs must cover exactly active fields {expected}"
                )
        for name, source in (spec.forests or {}).items():
            if not isinstance(source, ForestInput):
                raise ValueError(f"{name}: require ForestInput")
        for value in (spec.galaxies or {}).values():
            galaxy_noise(value)


def prepare_bin(spec):
    """Fix covariance at registry fiducials, evaluating active inputs once.

    Returns PreparedBin, with immutable (node,required_pair) N/R/P/T, (node,)
    modes/k/mu, (node,field) W and (node,selected,selected) Cholesky factors.
    No derivative schedule is assumed until run_bin's explicit options arrive.
    """
    _validate_spec(spec)
    selection, geometry = spec.p3d.selection, spec.geometry
    k, mu = spec.grid.k_flat, spec.grid.mu_flat
    theta = _immutable(spec.p3d.registry.fiducials)
    weights, provenance, p1d_counts = {}, {}, {}
    model_counts = {r.provider.label: 0 for r in spec.p3d.routes}
    try:
        modes = mode_counts(geometry, spec.grid)
        response = prepare_response(
            selection.fields, k, mu, a_v=geometry.a_v, settings=spec.responses
        )
        products = pair_response(response, selection)
        if spec.full_noise is not None:
            noise = prepare_noise(selection, len(k), full=spec.full_noise)
            convention = "full replacement"
        else:
            diagonal = {
                name: np.full(len(k), galaxy_noise(n))
                for name, n in (spec.galaxies or {}).items()
            }
            q = wavenumber_comoving_to_velocity(k * mu, a_v=geometry.a_v)
            for field in selection.fields:
                if field.id not in (spec.forests or {}):
                    continue
                source = spec.forests[field.id]
                options = dict(source.weight_options)
                p1d_counts[field.id] = 0
                try:
                    if source.auxiliary_coordinates is not None:
                        kt, kp = source.auxiliary_coordinates
                        options["auxiliary"] = sample_auxiliary(
                            field,
                            geometry,
                            spec.responses[field.id],
                            snapshot_p3d(spec.p3d),
                            theta,
                            p1d_model=source.p1d_model,
                            p1d_parameters=source.p1d_parameters,
                            theta_p1d=source.theta_p1d,
                            k_t_deg=kt,
                            k_p_velocity=kp,
                        )
                        index = selection.fields.index(field)
                        for route in spec.p3d.routes:
                            if np.any(np.all(route.pairs == (index, index), axis=1)):
                                model_counts[route.provider.label] += 1
                        p1d_counts[field.id] += 1
                    w = prepare_forest_weights(
                        field, geometry, spec.responses[field.id], **options
                    )
                    power1d = evaluate_p1d(
                        source.p1d_model,
                        source.p1d_parameters,
                        source.theta_p1d,
                        geometry.z_eval,
                        q,
                    )
                    p1d_counts[field.id] += 1
                    diagonal[field.id] = forest_noise(
                        w, field, geometry, spec.responses[field.id], k, mu, power1d
                    ).total
                except Exception as error:
                    raise ValueError(
                        f"field {field.id}, P1D/weights: {error}"
                    ) from error
                weights[field.id] = w
                provenance[field.id] = dict(
                    source=source.provenance,
                    theta_p1d=source.theta_p1d,
                    p1d_local_names=source.p1d_parameters.binding.local_names,
                    p1d_global_ids=source.p1d_parameters.registry.ids,
                    p1d_label=getattr(
                        source.p1d_model,
                        "__qualname__",
                        type(source.p1d_model).__name__,
                    ),
                )
            noise = prepare_noise(
                selection, len(k), diagonal=diagonal, independent_sampling=True
            )
            convention = "independent sampling"
        power = evaluate_p3d(spec.p3d, theta, geometry.z_eval, k, mu)
        for name in model_counts:
            model_counts[name] += 1
        total = combine_observed_power(products * power, noise)
        factors = factor_covariance(gaussian_covariance(total, modes, selection))
    except Exception as error:
        raise ValueError(
            f"bin {spec.id}, preparation nodes [0:{len(k)}): {error}"
        ) from error
    arrays = {
        name: _immutable(value)
        for name, value in dict(
            theta=theta,
            k=k,
            mu=mu,
            modes=modes,
            response=response,
            products=products,
            noise=noise,
            power=power,
            total=total,
            factors=factors,
        ).items()
    }
    diagnostics = freeze(
        dict(
            node_count=len(k),
            selected_pair_count=len(selection.selected_pairs),
            required_pair_count=len(selection.required_pairs),
            volume=geometry.volume,
            field_ids=[f.id for f in selection.fields],
            selected_pairs=selection.selected_pairs,
            required_pairs=selection.required_pairs,
            responses={
                name: (r.pixel_width_velocity, r.gaussian_sigma_velocity)
                for name, r in spec.responses.items()
            },
            noise_convention=convention,
            p3d_calls=model_counts,
            p1d_calls=p1d_counts,
            galaxies=spec.galaxies,
            forest_provenance=provenance,
            units=dict(power="(Mpc/h_fid)^3", k="h_fid/Mpc", modes="dimensionless"),
        )
    )
    return PreparedBin(
        spec.id,
        geometry,
        snapshot_p3d(spec.p3d),
        **arrays,
        weights=MappingProxyType(weights),
        diagnostics=diagnostics,
    )


@dataclass(frozen=True)
class BinRun:
    """One zero-prior FisherResult and actual batched derivative diagnostics."""

    id: str
    result: FisherResult
    columns: tuple
    calls: tuple[CallCount, ...]
    node_slices: tuple


@dataclass(frozen=True)
class ForecastRun:
    """Caller-ordered bin results and a combined result with one explicit prior."""

    bin_ids: tuple[str, ...]
    bins: tuple[BinRun, ...]
    combined: FisherResult


def _batch_size(batch_size, n):
    return n if batch_size is None else integer(batch_size, "batch_size", minimum=1)


def run_bin(prepared, *, batch_size=None, steps=None, step_scale=1.0, numerical=False):
    """Reuse fixed factors; transient required-pair Jacobians are node-bounded.

    Batches are consecutive C-order nodes, never spectra. Derivative evaluation
    may repeat fiducial P3D calls. Survey quantities and factors are never rebuilt.
    """
    if not isinstance(prepared, PreparedBin):
        raise ValueError("require PreparedBin")
    n = len(prepared.k)
    size = _batch_size(batch_size, n)
    _schedule(prepared.p3d, prepared.theta, steps, step_scale, numerical)
    data = np.zeros((len(prepared.theta), len(prepared.theta)))
    calls, slices, columns = {}, [], ()
    selected = prepared.p3d.selection.selected_to_required
    for start in range(0, n, size):
        stop = min(start + size, n)
        part = slice(start, stop)
        try:
            derivative = evaluate_derivatives(
                prepared.p3d,
                prepared.theta,
                prepared.geometry.z_eval,
                prepared.k[part],
                prepared.mu[part],
                steps=steps,
                step_scale=step_scale,
                numerical=numerical,
            )
            jac = (prepared.products[part, :, None] * derivative.jacobian)[
                :, selected, :
            ]
            data += fisher_from_factors(jac, prepared.factors[part])
        except Exception as error:
            raise ValueError(
                f"bin {prepared.id}, global nodes [{start}:{stop}) (local node + {start}): {error}"
            ) from error
        columns = derivative.columns
        slices.append((start, stop))
        for count in derivative.calls:
            old = calls.get(count.provider, (0, 0))
            calls[count.provider] = (old[0] + count.model, old[1] + count.jacobian)
    return BinRun(
        prepared.id,
        FisherResult(prepared.p3d.registry, data),
        columns,
        tuple(CallCount(name, *counts) for name, counts in calls.items()),
        tuple(slices),
    )


def _validate_bins(bins):
    bins = tuple(bins)
    if not bins or any(not isinstance(b, PreparedBin) for b in bins):
        raise ValueError("require nonempty sequence of PreparedBin")
    if len({b.id for b in bins}) != len(bins):
        raise ValueError("duplicate bin IDs")
    registry = bins[0].p3d.registry
    identities = {}
    for i, b in enumerate(bins):
        if b.p3d.registry is not registry:
            raise ValueError(f"bin {b.id}: different registry identity")
        if not np.array_equal(b.theta, bins[0].theta):
            raise ValueError("bins were prepared at different global fiducials")
        if b.geometry.h_fid != bins[0].geometry.h_fid:
            raise ValueError("bins have different h_fid")
        for field in b.p3d.selection.fields:
            if field.id in identities and field != identities[field.id]:
                raise ValueError(f"inconsistent field identity {field.id}")
            identities[field.id] = field
        for other in bins[:i]:
            if max(b.geometry.z_min, other.geometry.z_min) < min(
                b.geometry.z_max, other.geometry.z_max
            ):
                raise ValueError(f"overlapping bin interiors: {other.id}, {b.id}")
    return bins


def run_forecast(
    bins,
    *,
    prior_fisher=None,
    batch_size=None,
    steps=None,
    step_scale=1.0,
    numerical=False,
):
    """Run prepared independent bins in one registry; add prior exactly once.

    No per-bin marginalization, implicit preparation/cache refresh or partial
    success return. Changing priors reuses the same prepared bins.
    """
    bins = _validate_bins(bins)
    for b in bins:
        _batch_size(batch_size, len(b.k))
        _schedule(b.p3d, b.theta, steps, step_scale, numerical)
    # Validate the requested prior before invoking any models.
    FisherResult(
        bins[0].p3d.registry,
        np.zeros((len(bins[0].theta),) * 2),
        prior_fisher=prior_fisher,
    )
    runs = tuple(
        run_bin(
            b,
            batch_size=batch_size,
            steps=steps,
            step_scale=step_scale,
            numerical=numerical,
        )
        for b in bins
    )
    return ForecastRun(
        tuple(b.id for b in bins),
        runs,
        combine_results([r.result for r in runs], prior_fisher=prior_fisher),
    )
