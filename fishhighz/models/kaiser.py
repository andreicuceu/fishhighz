"""One-bin intrinsic template Kaiser model for the plain P3D callable contract."""

from dataclasses import dataclass

import numpy as np

from .._arrays import integer, label, readonly, real_array, scalar
from ..fields import ObservedField
from ..kernels.kaiser import (
    _assemble,
    _coordinates,
    _damping,
    _field_factors,
    _resolve,
    _scales,
)
from .templates import PowerTemplate

_BASES = {
    "ap_at": ("ap", "at"),
    "alpha_phi": ("alpha", "phi"),
    "alpha_iso_epsilon": ("alpha_iso", "epsilon"),
    "alpha_iso_phi": ("alpha_iso", "phi"),
}


@dataclass(frozen=True, init=False)
class Scaling:
    """Explicit component basis and two fixed numbers or named local slots.

    Parameters
    ----------
    basis : {'ap_at', 'alpha_phi', 'alpha_iso_epsilon', 'alpha_iso_phi'}
        Coordinates respectively (ap,at), (alpha,phi), (alpha_iso,epsilon).
        These map to Vega ap_at, phi_alpha, aiso_epsilon conventions.
    **coordinates
        Exactly the two named coordinates in that basis. Numbers are fixed;
        strings name local free slots. No cross-basis ties are inferred.
    """

    basis: str
    values: tuple

    def __init__(self, basis, **coordinates):
        if basis not in _BASES or set(coordinates) != set(_BASES[basis]):
            raise ValueError(
                "scaling requires one supported basis and exactly its coordinate keys"
            )
        values = []
        for name in _BASES[basis]:
            value = coordinates[name]
            if isinstance(value, str):
                label(value, "scaling local slot")
            else:
                value = scalar(value, name)
                if value <= (-1 if name == "epsilon" else 0):
                    raise ValueError(f"{name}: invalid scaling coordinate")
            values.append(value)
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "values", tuple(values))


@dataclass(frozen=True, init=False, eq=False)
class KaiserModel:
    """Prepare a fixed-redshift intrinsic P3D callable.

    Parameters
    ----------
    template : PowerTemplate
        Accepted signed template in fixed h_fid units.
    fields : sequence of ObservedField
        Original field order, all supported by this model; no renumbering.
    biases : mapping
        Every field ID to a fixed finite number or a named local slot.
    betas : mapping
        Exactly the forest field IDs to fixed numbers or named local slots.
    widths : mapping
        Every field ID to explicit fixed (parallel,transverse) widths >=0,
        in Mpc/h_fid. Both values are required and may never be free slots.
    f : number or str, optional
        One shared galaxy growth rate, fixed or a named slot. Required with
        galaxies, forbidden without them. No galaxy beta parameter is accepted.
    local_names : sequence of str
        Explicit ordered free slots. Every slot must occur exactly once in
        settings; tie distinct slots via existing BoundParameters global IDs.
    smooth, wiggle : Scaling, optional
        Independent scaling blocks; omitted blocks are fixed identity ap_at.
    z, growth : float, optional
        One fixed bin and G power amplitude, following PowerTemplate's rules.
        Calls must use this exact z. G defaults to 1 only at template.z_ref.

    Notes
    -----
    Use as model(theta_local,z,k,mu,pairs) within P3DProvider. No production
    analytic Jacobian is supplied: Step 06 numerical derivatives use explicit
    steps and must remain in BOTH mapped template domains. Output excludes
    instrument response/noise; caller owns grid/template h_fid consistency.
    """

    template: PowerTemplate
    fields: tuple
    local_names: tuple
    fixed: np.ndarray
    destinations: np.ndarray
    sources: np.ndarray
    forest: np.ndarray
    widths: np.ndarray
    widths_squared: np.ndarray
    bases: tuple
    z: float
    growth: float

    def __init__(
        self,
        template,
        fields,
        *,
        biases,
        betas,
        widths,
        f=None,
        local_names=(),
        smooth=None,
        wiggle=None,
        z=None,
        growth=None,
    ):
        if not isinstance(template, PowerTemplate):
            raise ValueError("expected prepared PowerTemplate")
        fields = tuple(fields)
        if not fields or any(not isinstance(field, ObservedField) for field in fields):
            raise ValueError("require nonempty ObservedField sequence")
        ids = [field.id for field in fields]
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate field ID")
        forest = np.array([field.kind == "forest" for field in fields])
        if set(biases) != set(ids) or set(betas) != {
            field.id for field in fields if field.kind == "forest"
        }:
            raise ValueError(
                "biases must match all fields; betas must match only forest fields"
            )
        if (np.any(~forest) and f is None) or (np.all(forest) and f is not None):
            raise ValueError(
                "one shared f is required for galaxies and forbidden for forest-only models"
            )
        if set(widths) != set(ids):
            raise ValueError("fixed widths must include every model field")
        width_array = real_array([widths[name] for name in ids], "fixed widths")
        if width_array.shape != (len(ids), 2) or np.any(width_array < 0):
            raise ValueError(
                "require fixed nonnegative (parallel,transverse) widths per field"
            )
        with np.errstate(over="ignore"):
            width_squared = width_array**2
        if not np.all(np.isfinite(width_squared)) or np.any(
            (width_array != 0) & (width_squared == 0)
        ):
            raise ValueError("unrepresentable squared fixed widths")
        if isinstance(local_names, str):
            raise ValueError("local_names must be an ordered sequence")
        names = tuple(local_names)
        for name in names:
            label(name, "local slot")
        if len(set(names)) != len(names):
            raise ValueError("duplicate local slot name")
        smooth = Scaling("ap_at", ap=1, at=1) if smooth is None else smooth
        wiggle = Scaling("ap_at", ap=1, at=1) if wiggle is None else wiggle
        if not isinstance(smooth, Scaling) or not isinstance(wiggle, Scaling):
            raise ValueError("smooth and wiggle must be Scaling records")
        settings = (
            [biases[name] for name in ids]
            + [betas.get(name, 0) for name in ids]
            + [0 if f is None else f]
            + list(smooth.values)
            + list(wiggle.values)
        )
        fixed, destinations, sources, used = [], [], [], set()
        for destination, value in enumerate(settings):
            if isinstance(value, str):
                if value not in names or value in used:
                    raise ValueError(
                        f"unknown or duplicate local slot assignment {value!r}; use distinct slots with explicit global ties"
                    )
                used.add(value)
                destinations.append(destination)
                sources.append(names.index(value))
                fixed.append(0.0)
            else:
                fixed.append(scalar(value, f"setting {destination}"))
        if used != set(names):
            raise ValueError(f"unused local slots: {set(names) - used}")
        z = template.z_ref if z is None else scalar(z, "fixed redshift")
        # Reuse the accepted redshift/G validation, without preparing coefficients.
        template.evaluate([template.k[0]], z=z, growth=growth)
        growth = 1.0 if growth is None else scalar(growth, "fixed growth power G")
        for name, value in (
            ("template", template),
            ("fields", fields),
            ("local_names", names),
            ("fixed", readonly(fixed, np.float64)),
            ("destinations", readonly(destinations, np.int64)),
            ("sources", readonly(sources, np.int64)),
            ("forest", readonly(forest, np.bool_)),
            ("widths", readonly(width_array, np.float64)),
            ("widths_squared", readonly(width_squared, np.float64)),
            ("bases", (smooth.basis, wiggle.basis)),
            ("z", z),
            ("growth", growth),
        ):
            object.__setattr__(self, name, value)

    def __call__(self, theta_local, z, k, mu, pairs):
        """Return owned float64 intrinsic (node,pair) power in requested order.

        Real paired nodes, finite local state and original integer field indices
        are validated. Domain failures identify the component/scales. There is
        no mutable parameter cache and no cuts, noise or volume configuration.
        """
        theta = real_array(theta_local, "theta_local")
        if theta.shape != (len(self.local_names),):
            raise ValueError("theta_local length must match local_names")
        if scalar(z, "redshift") != self.z:
            raise ValueError(
                f"model prepared for z={self.z}; different redshift requested"
            )
        k, mu = real_array(k, "k"), real_array(mu, "mu")
        if (
            k.ndim != 1
            or not len(k)
            or mu.shape != k.shape
            or np.any(k <= 0)
            or np.any((mu < 0) | (mu > 1))
        ):
            raise ValueError("require nonempty paired k>0 and mu in [0,1]")
        raw = np.asarray(pairs, dtype=object)
        if raw.ndim != 2 or raw.shape[1] != 2 or not len(raw):
            raise ValueError("pairs must be nonempty (n_pair,2) integer field indices")
        pairs = np.array(
            [[integer(item, "field index") for item in pair] for pair in raw],
            dtype=np.int64,
        )
        if np.any(pairs >= len(self.fields)):
            raise ValueError("unknown field index")
        values = _resolve(self.fixed, self.destinations, self.sources, theta)
        n = len(self.fields)
        bias, beta, f = values[:n], values[n : 2 * n], values[2 * n]
        components, factors, prefactors = [], [], []
        for component, basis in enumerate(self.bases):
            name = ("smooth", "wiggle")[component]
            coordinates = values[2 * n + 1 + 2 * component : 2 * n + 3 + 2 * component]
            if coordinates[0] <= 0 or coordinates[1] <= (
                -1 if basis == "alpha_iso_epsilon" else 0
            ):
                raise ValueError(
                    f"{name}: invalid {basis} scaling values {coordinates.tolist()}"
                )
            with np.errstate(
                over="ignore", invalid="ignore", divide="ignore", under="ignore"
            ):
                ap, at, q = _scales(coordinates, tuple(_BASES).index(basis))
                if not np.all(np.isfinite([ap, at, q])) or min(ap, at, q) <= 0:
                    raise ValueError(
                        f"{name}: unrepresentable scaling factors/Q for {coordinates.tolist()}"
                    )
                mapped, angle, parallel, transverse = _coordinates(k, mu, ap, at)
            context = f"{name}, {basis}={coordinates.tolist()}, ap={ap}, at={at}"
            if not np.all(np.isfinite([mapped, angle, parallel, transverse])) or np.any(
                mapped <= 0
            ):
                raise ValueError(
                    f"{context}: nonfinite/unrepresentable mapped coordinates"
                )
            if mapped.min() < self.template.k[0] or mapped.max() > self.template.k[-1]:
                raise ValueError(
                    f"{context}: mapped range {(float(mapped.min()), float(mapped.max()))} outside template domain {self.template.domain}; pad template for all derivative stencils"
                )
            components.append(self.template.evaluate(mapped)[:, component])
            with np.errstate(over="ignore", invalid="ignore"):
                field_values = _field_factors(angle, bias, beta, f, self.forest)
            if not np.all(np.isfinite(field_values)):
                raise ValueError(f"{context}: nonfinite Kaiser factors")
            factors.append(field_values)
            prefactors.append(q)
        # parallel/transverse here are the WIGGLE component's coordinates.
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            damping, exponent = _damping(
                parallel, transverse, self.widths_squared, pairs
            )
            output = _assemble(
                *components, *factors, damping, pairs, *prefactors, self.growth
            )
        if not np.all(np.isfinite(exponent)) or not np.all(np.isfinite(output)):
            raise ValueError(
                f"{context}: nonfinite damping exponent or assembled power at z={self.z}"
            )
        return np.array(output, dtype=np.float64, order="C", copy=True)
