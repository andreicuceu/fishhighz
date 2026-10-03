"""Fixed-survey preparation and independent-bin Fisher execution (NumPy only)."""

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from ._arrays import integer
from .covariance import (
    combine_observed_power,
    gaussian_covariance,
    gaussian_variances,
)
from .derivatives import CallCount, _schedule, evaluate_derivatives
from .fisher import factor_covariance, fisher_from_factors
from .geometry import _immutable, mode_counts, wavenumber_comoving_to_velocity
from .models.external import evaluate_p1d, evaluate_p3d
from .noise import forest_noise, galaxy_noise, prepare_noise
from .response import InstrumentResponse, pair_response, prepare_response
from .results import FisherResult, combine_results
from .survey import BinSpec, ForestInput, PreparedBin, freeze, snapshot_p3d
from .weights import (
    prepare_forest_weights,
    prepare_integrated_forest_weights,
    sample_auxiliary,
)


def _validate_spec(spec):
    """Validate bin model, response, and noise preparation contracts.

    Parameters
    ----------
    spec : BinSpec
        Fixed geometry, Fourier grid, model routes, and explicit noise inputs.

    Returns
    -------
    None
        Validate structural and physical context without evaluating model power.

    Raises
    ------
    ValueError
        If built-in model context, per-field responses, noise settings, or
        positive galaxy densities are inconsistent.
    """
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
    """Fix per-field weights, known noise, and covariance at registry fiducials.

    Parameters
    ----------
    spec : BinSpec
        Complete bin definition with fixed geometry, quadrature, model bindings,
        responses, and either generated or full supplied noise.

    Returns
    -------
    prepared : PreparedBin
        Immutable fiducial state: power/noise/total arrays (n_node,
        n_required_pair) in (Mpc/h_fid)^3, dimensionless field and pair
        responses, observed k in h_fid/Mpc and mu, mode counts, and Cholesky
        factors (n_node, n_selected, n_selected) in power units.

    Raises
    ------
    ValueError
        If bin context, weighting, provider evaluation, physical noise, or
        covariance factorization fails; errors include bin and field context.

    Notes
    -----
    Active inputs are evaluated during preparation only. Per-field forest
    weights, convergence metadata, provenance, and model call counts are
    retained. Numerical derivative schedules are checked later when run
    options are supplied. Failed preparation returns no partial bin.
    """
    _validate_spec(spec)
    selection, geometry = spec.p3d.selection, spec.geometry
    k_grid, mu_grid = spec.grid.k_flat, spec.grid.mu_flat
    theta = _immutable(spec.p3d.registry.fiducials)
    weights, provenance, p1d_counts = {}, {}, {}
    model_counts = {r.provider.label: 0 for r in spec.p3d.routes}

    try:
        modes = mode_counts(geometry, spec.grid)
        response = prepare_response(
            selection.fields, k_grid, mu_grid, a_v=geometry.a_v, settings=spec.responses
        )
        products = pair_response(response, selection)

        # Prepare noise once per field, independently of the derivative calls.
        if spec.full_noise is not None:
            noise = prepare_noise(selection, len(k_grid), full=spec.full_noise)
            convention = "full replacement"
        else:
            diagonal = {
                name: np.full(len(k_grid), galaxy_noise(n))
                for name, n in (spec.galaxies or {}).items()
            }
            velocity_wavenumber = wavenumber_comoving_to_velocity(
                k_grid * mu_grid, a_v=geometry.a_v
            )
            for field in selection.fields:
                if field.id not in (spec.forests or {}):
                    continue
                source = spec.forests[field.id]
                options = dict(source.weight_options)
                p1d_counts[field.id] = 0
                try:
                    if source.auxiliary_coordinates is not None:
                        transverse_angular_wavenumber, parallel_velocity_wavenumber = (
                            source.auxiliary_coordinates
                        )
                        options["auxiliary"] = sample_auxiliary(
                            field,
                            geometry,
                            spec.responses[field.id],
                            snapshot_p3d(spec.p3d),
                            theta,
                            p1d_model=source.p1d_model,
                            p1d_parameters=source.p1d_parameters,
                            theta_p1d=source.theta_p1d,
                            k_t_deg=transverse_angular_wavenumber,
                            k_p_velocity=parallel_velocity_wavenumber,
                        )
                        index = selection.fields.index(field)
                        for route in spec.p3d.routes:
                            if np.any(np.all(route.pairs == (index, index), axis=1)):
                                model_counts[route.provider.label] += 1
                        p1d_counts[field.id] += 1

                    if source.integrated is None:
                        forest_weights = prepare_forest_weights(
                            field, geometry, spec.responses[field.id], **options
                        )
                    else:
                        # Integrated source: same method and stopping options.
                        forest_weights = prepare_integrated_forest_weights(
                            field,
                            geometry,
                            spec.responses[field.id],
                            source.integrated,
                            **options,
                        )
                    power1d = evaluate_p1d(
                        source.p1d_model,
                        source.p1d_parameters,
                        source.theta_p1d,
                        geometry.z_eval,
                        velocity_wavenumber,
                    )
                    p1d_counts[field.id] += 1
                    diagonal[field.id] = forest_noise(
                        forest_weights,
                        field,
                        geometry,
                        spec.responses[field.id],
                        k_grid,
                        mu_grid,
                        power1d,
                    ).total
                except Exception as error:
                    raise ValueError(
                        f"field {field.id}, P1D/weights: {error}"
                    ) from error
                weights[field.id] = forest_weights
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
                selection, len(k_grid), diagonal=diagonal, independent_sampling=True
            )
            convention = "independent sampling"

        # The fiducial total power fixes covariance for the entire bin run.
        power = evaluate_p3d(spec.p3d, theta, geometry.z_eval, k_grid, mu_grid)
        for name in model_counts:
            model_counts[name] += 1
        total = combine_observed_power(products * power, noise)
        factors = factor_covariance(gaussian_covariance(total, modes, selection))
    except Exception as error:
        raise ValueError(
            f"bin {spec.id}, preparation nodes [0:{len(k_grid)}): {error}"
        ) from error

    arrays = {
        name: _immutable(value)
        for name, value in dict(
            theta=theta,
            k=k_grid,
            mu=mu_grid,
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
            node_count=len(k_grid),
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
    """One zero-prior FisherResult and actual derivative columns/calls/batches.

    ``individual`` holds, when requested, one zero-prior FisherResult per
    selected spectrum (selected-pair order) using that spectrum's own variance.
    """

    id: str
    result: FisherResult
    columns: tuple
    calls: tuple[CallCount, ...]
    node_slices: tuple
    individual: tuple | None = None


@dataclass(frozen=True)
class ForecastRun:
    """Caller-ordered BinRuns and their combined FisherResult with one prior."""

    bin_ids: tuple[str, ...]
    bins: tuple[BinRun, ...]
    combined: FisherResult


def _batch_size(batch_size, n):
    """Resolve the requested maximum number of Fourier nodes per batch.

    Parameters
    ----------
    batch_size : int or None
        Positive explicit node count, or None for all nodes.
    n : int
        Number of prepared Fourier nodes.

    Returns
    -------
    size : int
        Explicit batch size or n when batch_size is None.

    Raises
    ------
    ValueError
        If the explicit batch size is not a positive integer.
    """
    return n if batch_size is None else integer(batch_size, "batch_size", minimum=1)


def _individual_factors(prepared):
    """One-spectrum (node,1,1) factors from the diagonal of the fixed covariance.

    Parameters
    ----------
    prepared : PreparedBin
        Fixed fiducial total power, mode counts, and selected spectrum
        definitions.

    Returns
    -------
    factors : list of ndarray
        One array of shape (n_node, 1, 1) per selected spectrum, containing
        Cholesky factors in (Mpc/h_fid)^3.

    Raises
    ------
    ValueError
        If an individual variance cannot be factored; includes selected-pair
        context.

    Notes
    -----
    Each factor is exactly what factor_covariance returns for the one-spectrum
    covariance of an independently prepared single-pair bin, because forest
    weights, noise and response are per field and the variance arithmetic is
    shared with the full covariance.
    """
    selection = prepared.p3d.selection
    variances = gaussian_variances(prepared.total, prepared.modes, selection)
    factors = []
    for spectrum, pair in enumerate(selection.selected_pairs.tolist()):
        try:
            factors.append(factor_covariance(variances[:, spectrum, None, None]))
        except ValueError as error:
            raise ValueError(f"selected spectrum {tuple(pair)}: {error}") from error
    return factors


def run_bin(
    prepared,
    *,
    batch_size=None,
    steps=None,
    step_scale=1.0,
    numerical=False,
    individual=False,
):
    """Reuse fixed factors; transient required-pair Jacobians are node-bounded.

    Parameters
    ----------
    prepared : PreparedBin
        Validated immutable fiducial bin state.
    batch_size : int, optional
        Positive maximum Fourier nodes per consecutive batch. Default None
        processes all nodes at once.
    steps : mapping of str to float, optional
        Absolute finite-difference steps by global ID in parameter units.
        Default None uses registry steps.
    step_scale : float, default=1.0
        Positive dimensionless multiplier of the numerical steps.
    numerical : bool, default=False
        Force numerical evaluation of otherwise analytic derivative columns.
    individual : bool, default=False
        Also accumulate each selected spectrum using its own one-spectrum
        covariance.

    Returns
    -------
    run : BinRun
        Zero-prior global Fisher information in inverse products of parameter
        units, derivative/call diagnostics, node slices, and optional
        individual-spectrum results.

    Raises
    ------
    ValueError
        If the prepared type, derivative schedule, provider evaluation, or
        contraction fails; includes bin and node context.

    Notes
    -----
    Batches are consecutive C-order nodes, never spectra. Derivative evaluation
    may repeat fiducial P3D calls. Survey quantities and factors are never rebuilt.
    Only the scheduled (nonzero) global columns enter the contraction; the other
    rows and columns of the returned information are exactly zero. With
    ``individual=True`` the same Jacobian batches also give each selected
    spectrum's Fisher matrix under its own one-spectrum covariance.
    """
    if not isinstance(prepared, PreparedBin):
        raise ValueError("require PreparedBin")
    n_node = len(prepared.k)
    size = _batch_size(batch_size, n_node)
    schedules = _schedule(prepared.p3d, prepared.theta, steps, step_scale, numerical)
    active = np.array(
        sorted({index for schedule in schedules for index, _ in schedule}),
        dtype=np.intp,
    )
    if not len(active):
        active = np.arange(len(prepared.theta))
    block = np.ix_(active, active)
    data = np.zeros((len(prepared.theta), len(prepared.theta)))
    calls, slices, columns = {}, [], ()
    selected = prepared.p3d.selection.selected_to_required
    if individual:
        own_factors = _individual_factors(prepared)
        own_data = [np.zeros_like(data) for _ in own_factors]
    for start in range(0, n_node, size):
        stop = min(start + size, n_node)
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
            observed_jacobian = (
                prepared.products[part, :, None] * derivative.jacobian
            )[:, selected, :][:, :, active]
            data[block] += fisher_from_factors(
                observed_jacobian, prepared.factors[part]
            )
            if individual:
                for spectrum, factor in enumerate(own_factors):
                    own_data[spectrum][block] += fisher_from_factors(
                        observed_jacobian[:, spectrum : spectrum + 1], factor[part]
                    )
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
        tuple(FisherResult(prepared.p3d.registry, own) for own in own_data)
        if individual
        else None,
    )


def _validate_bins(bins):
    """Validate independent prepared bins in a shared global parameter basis.

    Parameters
    ----------
    bins : iterable of PreparedBin
        Nonempty sequence with unique bin IDs and disjoint redshift interiors.

    Returns
    -------
    bins : tuple of PreparedBin
        Input bin order retained after validation.

    Raises
    ------
    ValueError
        If bin types, IDs, registry identity, fiducials, h_fid, field
        identities, or redshift independence are inconsistent.
    """
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
    individual=False,
):
    """Run prepared independent bins in one registry; add prior exactly once.

    Parameters
    ----------
    bins : iterable of PreparedBin
        Independent prepared bins in the same exact registry and fiducial state.
    prior_fisher : array_like of shape (n_global, n_global), optional
        Shared prior information in inverse products of parameter units; default
        None applies no prior.
    batch_size : int, optional
        Positive maximum Fourier nodes per consecutive batch. Default None
        processes all nodes at once.
    steps : mapping of str to float, optional
        Absolute finite-difference steps by global ID in parameter units.
        Default None uses registry steps.
    step_scale : float, default=1.0
        Positive dimensionless multiplier of the numerical steps.
    numerical : bool, default=False
        Force numerical evaluation of otherwise analytic derivative columns.
    individual : bool, default=False
        Also accumulate each selected spectrum using its own one-spectrum
        covariance.

    Returns
    -------
    run : ForecastRun
        Ordered zero-prior bin results and their combined unmarginalized data
        information with the supplied prior added exactly once.

    Raises
    ------
    ValueError
        If bins, prior, derivative options, or any bin calculation fail.

    Notes
    -----
    No per-bin marginalization, implicit preparation/cache refresh or partial
    success return. Changing priors reuses the same prepared bins. The returned
    ForecastRun exposes each zero-prior bin result and ``combined``, where the
    caller's prior is applied once after summing independent-bin data Fisher
    matrices. ``individual=True`` also fills each BinRun's per-spectrum results.
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
            individual=individual,
        )
        for b in bins
    )
    return ForecastRun(
        tuple(b.id for b in bins),
        runs,
        combine_results([r.result for r in runs], prior_fisher=prior_fisher),
    )
