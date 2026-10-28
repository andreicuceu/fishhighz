"""Direct-cosmology EFT mean with one physical AP map and explicit gating."""

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

from .camb_snapshot import BASIS
from .physical_ap import physical_ap_coordinates


class CosmologicalEFTAdapter:
    """Evaluate shared physical tracers at paired fiducial-observed nodes.

    The snapshot provider owns transfer/growth/geometry. The caller supplies
    a pure tracer_parameters(values,growth_rate) mapping to eftkernels
    EFTTracerParameters, preserving its explicit bias convention. Repeated
    observed forest fields may point to the same physical tracer object.

    Intrinsic_additive(values,k,mu,pair,state,z_index), when provided, owns
    forecast-specific residual stochastic and higher-FoG terms in power units.
    Known observational noise and response products are always external.

    Features are currently permitted only in an explicitly undamped tree
    benchmark. Construction rejects a one-loop feature or BAO-only IR applied
    to a primordial feature: that production path awaits its scientific gate.
    Gaussian one-loop cosmology updates are supported independently.
    """

    def __init__(
        self,
        snapshot_provider,
        z,
        physical_tracers,
        parameter_names,
        defaults,
        tracer_parameters,
        *,
        order="loop",
        feature=None,
        no_ir_resummation=False,
        intrinsic_additive=None,
        nmax=256,
        master_integrals="bundled",
        num_mu=64,
        k_output=None,
    ):
        if order not in ("tree", "loop"):
            raise ValueError("order must be tree or loop")
        if feature is not None and (order != "tree" or not no_ir_resummation):
            raise ValueError(
                "primordial-feature one-loop/IR response is not validated; only explicit undamped tree benchmark is supported"
            )
        names = tuple(parameter_names)
        if len(set(names)) != len(names):
            raise ValueError("duplicate parameters")
        if any(n == "f" or n.startswith(("alpha_", "phi_", "f_sigma8")) for n in names):
            raise ValueError(
                "direct cosmology cannot include independent AP or growth coordinates"
            )
        extra = any(
            n.startswith(("Pshot_", "a0_", "a2_", "stoch_", "tilde_c_")) for n in names
        )
        if extra and intrinsic_additive is None:
            raise ValueError(
                "residual stochastic/FoG coordinates require an explicit intrinsic additive model"
            )
        if feature is None and any(n in names for n in ("a_sin", "a_cos")):
            raise ValueError("feature coordinates require an explicit feature template")
        self.provider, self.z, self.tracers = (
            snapshot_provider,
            float(z),
            tuple(physical_tracers),
        )
        self.names, self.defaults = names, {**snapshot_provider.defaults, **defaults}
        if not set(names) <= self.defaults.keys():
            raise ValueError("missing parameter fiducial")
        self.tracer_parameters = tracer_parameters
        self.order, self.feature, self.no_ir = order, feature, bool(no_ir_resummation)
        self.intrinsic_additive = intrinsic_additive
        self.fiducial = snapshot_provider.evaluate()
        matched = np.flatnonzero(
            np.isclose(self.fiducial.redshifts, self.z, rtol=0, atol=1e-12)
        )
        if len(matched) != 1:
            raise ValueError("adapter redshift absent from snapshot")
        self.z_index = int(matched[0])
        self.nmax, self.masters, self.num_mu = nmax, master_integrals, num_mu
        self.k_output = (
            np.geomspace(0.0001, 0.35, 600)
            if k_output is None
            else np.asarray(k_output, dtype=np.float64)
        )
        self.engine = None
        self._engine_state = None
        self._power = {}
        self._loop = {}
        self.update_count = 0

    def __call__(self, theta, z, k, mu, pairs):
        if abs(z - self.z) > 1e-12:
            raise ValueError("adapter redshift mismatch")
        theta = np.asarray(theta, dtype=np.float64)
        if theta.shape != (len(self.names),) or np.any(~np.isfinite(theta)):
            raise ValueError("invalid parameter vector")
        p = {**self.defaults, **dict(zip(self.names, theta))}
        state = self.provider.evaluate({n: p[n] for n in BASIS})
        i = self.z_index
        qp = state.DM[i] / self.fiducial.DM[i]
        ql = self.fiducial.H[i] / state.H[i]
        km, mm, j = physical_ap_coordinates(k, mu, qp, ql)
        if km.ndim != 1:
            raise ValueError("require paired one-dimensional k and mu")
        objects = self.tracer_parameters(p, float(state.growth_rate[i]))
        if not set(self.tracers) <= objects.keys():
            raise ValueError("missing physical tracer coefficients")
        pairs = np.asarray(pairs, dtype=int)
        if (
            pairs.ndim != 2
            or pairs.shape[1] != 2
            or np.any(pairs < 0)
            or np.any(pairs >= len(self.tracers))
        ):
            raise ValueError("invalid observed field pairs")
        if state.key not in self._power:
            self._power[state.key] = CubicSpline(
                np.log(state.k_reference), state.pk_linear[i]
            )
        if np.any(km < state.k_reference[0]) or np.any(km > state.k_reference[-1]):
            raise ValueError("AP query outside transfer snapshot")
        simple_tree = self.order == "tree" and self.no_ir
        if not simple_tree:
            from eftkernels import EFTModel

            if self.engine is None:
                self.engine = EFTModel(
                    {"name": "A", "type": "continuous"},
                    {"name": "B", "type": "continuous"},
                    tree_only=self.order == "tree",
                    no_ir_resummation=self.no_ir,
                    add_stochastic_terms=False,
                    nmax=self.nmax,
                    master_integrals=self.masters,
                    num_mu=self.num_mu,
                )
                self.engine.set_output_grid(self.k_output, self.num_mu)
            if self._engine_state != state.key:
                self.engine.update_fiducial(
                    state.k_reference,
                    state.pk_linear[i],
                    state.pk_smooth[i],
                    state.growth_rate[i],
                    state.r_drag_mpc * state.h_fid,
                )
                self._engine_state = state.key
                self.update_count += 1
        output = np.empty((len(km), len(pairs)))
        for col, (ia, ib) in enumerate(pairs):
            pair = self.tracers[ia], self.tracers[ib]
            a, b = (objects[t] for t in pair)
            if simple_tree:
                ka = a.bias_density - a.bias_eta * state.growth_rate[i] * mm**2
                kb = b.bias_density - b.bias_eta * state.growth_rate[i] * mm**2
                power = self._power[state.key](np.log(km))
                if self.feature is not None:
                    q = self.feature.quadratures(km, state.h_fid)
                    power = power * (
                        1 + p.get("a_sin", 0.0) * q[0] + p.get("a_cos", 0.0) * q[1]
                    )
                result = ka * kb * power
            else:
                from eftkernels.parameters import LEGACY_PARAMETER_NAMES

                key = (
                    state.key,
                    tuple(getattr(a, n) for n in LEGACY_PARAMETER_NAMES),
                    tuple(getattr(b, n) for n in LEGACY_PARAMETER_NAMES),
                )
                if key not in self._loop:
                    values = self.engine.compute((a, b, None))
                    angular = np.sqrt(self.engine.mu_powers[1])
                    spline = CubicSpline(angular, values, axis=0)
                    self._loop[key] = RectBivariateSpline(
                        np.r_[0.0, angular, 1.0],
                        np.log(self.k_output),
                        np.vstack((spline(0.0), values, spline(1.0))),
                    )
                if np.any(km < self.k_output[0]) or np.any(km > self.k_output[-1]):
                    raise ValueError("AP query outside padded EFT output grid")
                result = self._loop[key].ev(mm, np.log(km))
            if self.intrinsic_additive is not None:
                result = result + self.intrinsic_additive(p, km, mm, pair, state, i)
            output[:, col] = j * result
        return output
