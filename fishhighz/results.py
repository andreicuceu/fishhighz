"""Named data/prior information, explicit fixing, and uncertainty diagnostics."""

from dataclasses import dataclass

import numpy as np

from ._arrays import readonly, scalar
from ._information import inspect_information
from .fisher import factor_covariance
from .kernels.fisher import _forward_substitute
from .parameters import ParameterRegistry


@dataclass(frozen=True, eq=False)
class FisherDiagnostics:
    """Rank in normalized coordinates x=scales*delta_theta (ordered by ids).

    null_directions has unit-Euclidean row vectors in x coordinates, including
    near-null directions at eigenvalue <= tolerance. Physical displacements are
    null_directions/scales[None,:]; they are not unit physical-coordinate vectors.
    scales=sqrt(F_total diagonal), using 1 in current units for exact zero rows.
    Eigenspace bases/signs are not unique. Arrays are owned/read-only float64.
    """

    ids: tuple[str, ...]
    rank: int
    condition: float
    tolerance: float
    scales: np.ndarray
    eigenvalues: np.ndarray
    null_directions: np.ndarray
    basis: str = "x = scales * delta_theta; unit Euclidean rows in x coordinates"


def _subset(registry, ids):
    if ids is None:
        return np.arange(len(registry.ids), dtype=np.int64)
    if isinstance(ids, str):
        raise ValueError("subset must be an ordered sequence of parameter IDs")
    ids = tuple(ids)
    if not ids or any(not isinstance(x, str) for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("subset must contain nonempty unique parameter IDs")
    lookup = {name: i for i, name in enumerate(registry.ids)}
    if any(name not in lookup for name in ids):
        raise ValueError("unknown parameter ID in subset")
    return np.array([lookup[name] for name in ids], dtype=np.int64)


def diagonal_prior(registry, sigmas):
    """Construct independent Gaussian prior precision in registry order.

    Parameters
    ----------
    registry : ParameterRegistry
        The exact global parameter basis.
    sigmas : mapping
        Explicit parameter ID to positive finite absolute width. Unlisted IDs
        receive zero information. Bounds and derivative steps are not priors.

    Returns
    -------
    prior : ndarray
        Owned C-contiguous float64 diagonal precision, not a covariance.
        Unrepresentable zero/infinite precision is rejected, not regularized.
    """
    if not isinstance(registry, ParameterRegistry):
        raise ValueError("registry must be a ParameterRegistry")
    prior = np.zeros((len(registry.ids), len(registry.ids)), dtype=np.float64)
    lookup = {name: i for i, name in enumerate(registry.ids)}
    for name, width in sigmas.items():
        if name not in lookup:
            raise ValueError(f"unknown prior parameter ID {name!r}")
        width = scalar(width, f"prior width for {name!r}")
        if width <= 0:
            raise ValueError(f"prior width for {name!r} must be positive")
        with np.errstate(over="ignore", under="ignore", divide="ignore"):
            precision = (np.float64(1) / width) ** 2
        if not np.isfinite(precision) or precision <= 0:
            raise ValueError(
                f"prior precision for {name!r} is not representable as positive finite float64"
            )
        prior[lookup[name], lookup[name]] = precision
    return prior


@dataclass(frozen=True, init=False, eq=False)
class FisherResult:
    """Information in one explicit registry, with no eager covariance inverse.

    Parameters
    ----------
    registry : ParameterRegistry
        Ordered IDs, fiducials, roles, bounds, and steps. Roles are metadata only.
    data_fisher : array_like, shape (n_global,n_global)
        Data information, possibly singular.
    prior_fisher : array_like, optional
        Separate finite symmetric PSD precision in the same units and at the
        same fiducial. Zero by default; bounds do not create priors.

    Notes
    -----
    Stored data_fisher, prior_fisher, total_fisher are owned read-only C-order
    float64. Negative diagonals and nonzero rows/columns with zero diagonal fail.
    R=F/scales/scales, scales=sqrt(diag(F)), or 1 for zero rows. Symmetry tolerance
    is elementwise 64*eps64*n*max(1,abs(Rij),abs(Rji)); accepted asymmetry is averaged
    only in owned copies. PSD allows lambda >= -t, rank counts lambda > t, with
    t=64*eps64*n*max(1,max(abs(eigenvalues(R)))). No eigenvalues are clipped.
    Diagnostics use total information; singular data or priors remain legal.
    """

    registry: ParameterRegistry
    data_fisher: np.ndarray
    prior_fisher: np.ndarray
    total_fisher: np.ndarray
    diagnostics: FisherDiagnostics

    def __init__(self, registry, data_fisher, *, prior_fisher=None):
        if not isinstance(registry, ParameterRegistry):
            raise ValueError("registry must be a ParameterRegistry")
        n = len(registry.ids)
        context = f"parameters {registry.ids!r}"
        data = inspect_information(data_fisher, f"data Fisher, {context}")[0]
        if data.shape != (n, n):
            raise ValueError("data Fisher shape must match registry")
        if prior_fisher is None:
            prior_fisher = np.zeros((n, n), dtype=np.float64)
        prior = inspect_information(prior_fisher, f"prior Fisher, {context}")[0]
        if prior.shape != (n, n):
            raise ValueError("prior Fisher shape must match registry")
        with np.errstate(over="ignore", invalid="ignore"):
            total = data + prior
        total, _, scales, values, vectors, tolerance = inspect_information(
            total, f"total Fisher, {context}"
        )
        rank = int(np.count_nonzero(values > tolerance))
        condition = float(values[-1] / values[0]) if rank == n else float("inf")
        diagnostics = FisherDiagnostics(
            ids=registry.ids,
            rank=rank,
            condition=condition,
            tolerance=float(tolerance),
            scales=readonly(scales, np.float64),
            eigenvalues=readonly(values, np.float64),
            null_directions=readonly(vectors[:, values <= tolerance].T, np.float64),
        )
        object.__setattr__(self, "registry", registry)
        for name, value in (
            ("data_fisher", data),
            ("prior_fisher", prior),
            ("total_fisher", total),
        ):
            object.__setattr__(self, name, readonly(value, np.float64))
        object.__setattr__(self, "diagnostics", diagnostics)

    def conditional_errors(self, ids=None):
        """Return ordered 1/sqrt(F_ii), fixing every other parameter.

        Parameters
        ----------
        ids : sequence of str, optional
            Requested IDs; all in registry order by default. Exactly zero
            information returns infinity, even for singular joint information.
        """
        index = _subset(self.registry, ids)
        diagonal = self.total_fisher.diagonal()[index]
        result = np.full(len(index), np.inf, dtype=np.float64)
        positive = diagonal > 0
        result[positive] = 1 / np.sqrt(diagonal[positive])
        return result

    def marginalized_covariance(self, ids=None):
        """Solve the full retained information, then select ordered covariance.

        Parameters
        ----------
        ids : sequence of str, optional
            Requested IDs; all by default. Other parameters remain free and
            marginalized over, not fixed. Singular full information raises with
            rank/null context, even if the requested sub-block is identifiable.

        Returns
        -------
        covariance : ndarray
            Owned C-contiguous float64 covariance of the requested parameters.
        """
        index = _subset(self.registry, ids)
        diag = self.diagnostics
        n = len(self.registry.ids)
        if diag.rank != n:
            raise ValueError(
                f"cannot marginalize singular Fisher: rank {diag.rank}/{n}, "
                f"IDs {diag.ids!r}, normalized threshold {diag.tolerance:.17g}; "
                f"null/near-null directions {diag.null_directions.tolist()} "
                "in x=scales*delta_theta (see diagnostics.scales); "
                "explicitly fix parameters or add a finite prior"
            )
        lower = factor_covariance(self.total_fisher[None])[0]
        solved = np.empty((n, n), dtype=np.float64)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            _forward_substitute(lower, np.eye(n), solved)
            covariance = solved.T @ solved
        if not np.all(np.isfinite(solved)) or not np.all(np.isfinite(covariance)):
            raise ValueError(
                f"nonfinite marginalized covariance for parameters {diag.ids!r}"
            )
        return np.array(
            covariance[np.ix_(index, index)], dtype=np.float64, order="C", copy=True
        )

    def marginalized_errors(self, ids=None):
        """Return sqrt of marginalized covariance diagonals, in requested order.

        Parameters
        ----------
        ids : sequence of str, optional
            Requested IDs; omitted parameters remain free, not fixed.
        """
        return np.sqrt(self.marginalized_covariance(ids).diagonal()).copy()

    def correlations(self, ids=None):
        """Return correlations derived from full marginalized covariance.

        Parameters
        ----------
        ids : sequence of str, optional
            Requested ordered IDs; defaults to all, without fixing others.
        """
        covariance = self.marginalized_covariance(ids)
        sigma = np.sqrt(covariance.diagonal())
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            correlation = covariance / sigma[:, None] / sigma[None, :]
        if not np.all(np.isfinite(correlation)):
            raise ValueError("nonfinite marginalized correlation")
        return correlation

    def fix_except(self, ids):
        """Keep ordered IDs free and explicitly fix their complement.

        Parameters
        ----------
        ids : nonempty sequence of str
            Retained free parameters. Principal submatrices preserve separate
            data/prior information and metadata. This is not marginalization.

        Returns
        -------
        result : FisherResult
            New result on the explicitly retained registry.
        """
        if ids is None:
            raise ValueError("fix_except requires an explicit nonempty ID subset")
        index = _subset(self.registry, ids)
        registry = ParameterRegistry([self.registry.parameters[i] for i in index])
        return FisherResult(
            registry,
            self.data_fisher[np.ix_(index, index)],
            prior_fisher=self.prior_fisher[np.ix_(index, index)],
        )


def combine_results(results, *, prior_fisher=None):
    """Sum independent unmarginalized data in an identical ordered registry.

    Parameters
    ----------
    results : nonempty iterable of FisherResult
        Data-only contributions. Nonzero existing priors are rejected. Every
        ordered parameter record (including fiducial, role, bounds, step) must
        match; singular individual contributions are legal.
    prior_fisher : array_like, optional
        One final precision matrix to apply once after summing all data.

    Returns
    -------
    result : FisherResult
        Combined data plus the one explicitly supplied prior. Combine before
        marginalizing shared nuisance parameters; no basis union is inferred.
    """
    results = tuple(results)
    if not results or any(not isinstance(r, FisherResult) for r in results):
        raise ValueError("combine_results requires nonempty FisherResult contributions")
    registry = results[0].registry
    total = np.zeros_like(results[0].data_fisher)
    for result in results:
        if result.registry.parameters != registry.parameters:
            raise ValueError(
                "contributions require identical ordered parameter metadata/fiducials"
            )
        if np.any(result.prior_fisher != 0):
            raise ValueError(
                "combine data-only results; supply a shared prior once at the end"
            )
        with np.errstate(over="ignore", invalid="ignore"):
            total += result.data_fisher
        if not np.all(np.isfinite(total)):
            raise ValueError("nonfinite combined data Fisher information")
    return FisherResult(registry, total, prior_fisher=prior_fisher)
