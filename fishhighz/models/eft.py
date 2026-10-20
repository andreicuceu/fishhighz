"""Optional adapter for the prepared full ``eftkernels.EFTModel`` API.

The adapter owns the component AP maps and growth updates. Known observational
noise is supplied separately to FishHighz's fixed fiducial covariance. Model
stochasticity here has independent pair coefficients, avoiding the legacy
half-sum cross prescription in eftkernels. No covariance is differentiated.
"""

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline


def component_coordinates(k, mu, alpha, phi):
    """Map paired fiducial nodes using the FishHighz alpha/phi convention.

    k is in h/Mpc; mu, alpha and phi are dimensionless. Return mapped k, mu
    and the volume Jacobian. alpha=sqrt(q_parallel*q_perpendicular), rather
    than the volume-averaged dilation alpha_iso.
    """
    if alpha <= 0 or phi <= 0:
        raise ValueError("alpha and phi must be positive")
    q_parallel, q_perpendicular = alpha / np.sqrt(phi), alpha * np.sqrt(phi)
    parallel = k * mu / q_parallel
    transverse = k * np.sqrt(1 - mu**2) / q_perpendicular
    mapped_k = np.hypot(parallel, transverse)
    return mapped_k, parallel / mapped_k, 1 / (q_parallel * q_perpendicular**2)


class EFTModelAdapter:
    """Supply paired-node FF, FG and GG means to an external P3D provider.

    Parameters
    ----------
    fiducial : mapping
        Explicit k, pk_full, pk_nw, growth_rate and rd_mpc_h, as returned by
        eftkernels.read_template. Power is (Mpc/h)^3, k is h/Mpc.
    parameter_names : sequence of str
        Local vector order. ``f`` changes the prepared cosmology. ``alpha_s``,
        ``phi_s``, ``alpha_w`` and ``phi_w`` dilate smooth/wiggle separately.
        Other names are EFT dictionary keys or ``stoch_{FF,FG,GG}_{0,k2,mu2}``.
    defaults, scales : mapping
        Local fiducials and physical conversion factors. Dictionary coefficients
        equal local value times scale. This makes priors explicit in a chosen
        dimensionless basis. AP and f always retain their physical values.
    order : {'tree', 'loop'}
        IR-resummed tree or full one-loop calculation.
    k_output : array_like
        Padded evaluation grid, independent of forecast modes, in h/Mpc.
    num_mu, nmax : int
        EFT angular and FFTLog resolutions. Interpolation refinement must be
        checked for the actual forecast. Calls are serial, with mutable state.
    k_natural : float
        h/Mpc scale defining independent stochastic k^2 polynomials.

    Notes
    -----
    Fields are ordered F=0, G=1. The exact input template redshift is required.
    Stochastic terms belong to the smooth component. A pure polynomial is
    evaluated at mapped coordinates, with the same component Jacobian.
    Full wiggles are P_loop[full]-P_loop[smooth], retaining mixed loop terms.
    """

    def __init__(
        self,
        fiducial,
        parameter_names,
        defaults,
        scales=None,
        *,
        order="loop",
        k_output=None,
        num_mu=64,
        nmax=256,
        k_natural=1.0,
    ):
        from eftkernels import EFTModel

        if order not in ("tree", "loop"):
            raise ValueError("order must be tree or loop")
        self.fiducial = fiducial
        self.names = tuple(parameter_names)
        self.defaults, self.scales = dict(defaults), dict(scales or {})
        self.order, self.k_natural = order, float(k_natural)
        self.k = np.array(
            np.geomspace(0.005, 1.12, 420) if k_output is None else k_output,
            dtype=np.float64,
        )
        self.models = []
        forest, galaxy = (
            dict(name="F", type="continuous"),
            dict(name="G", type="discrete"),
        )
        self.pairs = ((0, 0), (0, 1), (1, 1))
        for a, b in ((forest, forest), (forest, galaxy), (galaxy, galaxy)):
            model = EFTModel(
                a,
                b,
                num_mu=num_mu,
                nmax=nmax,
                tree_only=order == "tree",
                add_stochastic_terms=False,
                master_integrals="analytic",
            )
            self.models.append(model)
        self._growth = None
        self.update_growth(fiducial["growth_rate"])
        for model in self.models:
            model.set_output_grid(self.k, num_mu)
        self.mu = self.models[0].muk_grid[:, 0].copy()
        self._tree_splines = {}
        if order == "tree":
            for component in ("smooth", "wiggle"):
                matter = self.models[0].compute(
                    {"bias_F": 1.0, "bias_eta_F": 0.0}, component=component
                )
                angular = CubicSpline(self.mu, matter, axis=0)
                self._tree_splines[component] = RectBivariateSpline(
                    np.r_[0.0, self.mu, 1.0],
                    np.log(self.k),
                    np.vstack((angular(0.0), matter, angular(1.0))),
                )

    def update_growth(self, growth):
        """Update loop kernels and IR damping through the explicit EFT API."""
        if growth == self._growth or (
            self.order == "tree" and self._growth is not None
        ):
            return
        d = self.fiducial
        for model in self.models:
            model.update_fiducial(
                d["k"], d["pk_full"], d["pk_nw"], growth, d["rd_mpc_h"]
            )
        self._growth = growth

    def __call__(self, theta_local, z, k, mu, pairs):
        """Return signed clustering means with shape (node, requested pair)."""
        if abs(z - self.fiducial["z_fiducial"]) > 1e-12:
            raise ValueError("EFT adapter requires the supplied template redshift")
        parameters = dict(self.defaults)
        parameters.update(zip(self.names, theta_local))
        self.update_growth(parameters.get("f", self.fiducial["growth_rate"]))
        physical = {
            name: value * self.scales.get(name, 1.0)
            for name, value in parameters.items()
        }
        query = []
        for suffix in ("s", "w"):
            query.append(
                component_coordinates(
                    k,
                    mu,
                    parameters.get("alpha_" + suffix, 1.0),
                    parameters.get("phi_" + suffix, 1.0),
                )
            )
        output = np.empty((len(k), len(pairs)), dtype=np.float64)
        for column, pair in enumerate(pairs):
            index = self.pairs.index(tuple(sorted(pair)))
            model = self.models[index]
            output[:, column] = 0.0
            for component, (mapped_k, mapped_mu, jacobian) in zip(
                ("smooth", "wiggle"), query
            ):
                if np.any(mapped_k < self.k[0]) or np.any(mapped_k > self.k[-1]):
                    raise ValueError("AP query outside padded EFT output grid")
                if self.order == "tree":
                    growth = parameters.get("f", self.fiducial["growth_rate"])
                    kernels = (
                        physical["bias_F"]
                        + physical["bias_eta_F"] * growth * mapped_mu**2,
                        physical["bias_G"] + growth * mapped_mu**2,
                    )
                    output[:, column] += (
                        jacobian
                        * kernels[pair[0]]
                        * kernels[pair[1]]
                        * self._tree_splines[component].ev(mapped_mu, np.log(mapped_k))
                    )
                    continue
                power = model.compute(physical, component=component)
                angular = CubicSpline(self.mu, power, axis=0)
                extended_mu = np.r_[0.0, self.mu, 1.0]
                extended_power = np.vstack((angular(0.0), power, angular(1.0)))
                spline = RectBivariateSpline(
                    extended_mu, np.log(self.k), extended_power, kx=3, ky=3
                )
                output[:, column] += jacobian * spline.ev(mapped_mu, np.log(mapped_k))
            mapped_k, mapped_mu, jacobian = query[0]
            label = ("FF", "FG", "GG")[index]
            output[:, column] += jacobian * (
                physical.get(f"stoch_{label}_0", 0.0)
                + (mapped_k / self.k_natural) ** 2
                * (
                    physical.get(f"stoch_{label}_k2", 0.0)
                    + mapped_mu**2 * physical.get(f"stoch_{label}_mu2", 0.0)
                )
            )
        return output
