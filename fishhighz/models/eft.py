"""Optional adapter for the prepared full ``eftkernels.EFTModel`` API.

The adapter owns the component AP maps and growth updates. Known observational
noise is supplied separately to FishHighz's fixed fiducial covariance. Model
stochasticity here has independent pair coefficients, avoiding the legacy
half-sum cross prescription in eftkernels. No covariance is differentiated.
"""

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

from fishhighz.kernels.kaiser import _coordinates, _scales


def component_coordinates(k, mu, alpha, phi):
    """Map paired fiducial nodes using the FishHighz alpha/phi convention.

    Parameters
    ----------
    k, mu : ndarray of shape (n_node,)
        Paired wavenumbers in h/Mpc and direction cosines on [0, 1].
    alpha, phi : float
        Positive dimensionless dilations. Here
        alpha=sqrt(q_parallel*q_perpendicular), rather than alpha_iso.

    Returns
    -------
    mapped_k, mapped_mu : ndarray of shape (n_node,)
        Dilated wavenumbers in h/Mpc and dimensionless direction cosines.
    jacobian : float
        Dimensionless inverse dilation volume.
    """
    if not np.isfinite([alpha, phi]).all() or alpha <= 0 or phi <= 0:
        raise ValueError("alpha and phi must be finite and positive")
    q_parallel, q_perpendicular, jacobian = _scales((alpha, phi), 1)
    mapped_k, mapped_mu, _, _ = _coordinates(k, mu, q_parallel, q_perpendicular)
    return mapped_k, mapped_mu, jacobian


class EFTModelAdapter:
    """Supply paired-node FF, FG and GG means to an external P3D provider.

    Parameters
    ----------
    fiducial : mapping
        Required keys: k, pk_full, pk_nw, growth_rate, rd_mpc_h and z_fiducial.
        Power is (Mpc/h)^3, k is h/Mpc. For eftkernels.read_template inputs,
        the FITS template must include PKNW; pk_smooth is not substituted.
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
        if order not in ("tree", "loop"):
            raise ValueError("order must be tree or loop")
        required = {"k", "pk_full", "pk_nw", "growth_rate", "rd_mpc_h", "z_fiducial"}
        missing = required.difference(fiducial)
        if missing:
            raise ValueError(
                "EFT fiducial missing required keys: "
                + ", ".join(sorted(missing))
                + ". A read_template FITS input must include PKNW for pk_nw."
            )
        self.fiducial = fiducial
        self.names = tuple(parameter_names)
        self.defaults, self.scales = dict(defaults), dict(scales or {})
        self.order, self.k_natural = order, float(k_natural)
        self.k = np.array(
            np.geomspace(0.005, 1.12, 420) if k_output is None else k_output,
            dtype=np.float64,
        )
        self._num_mu, self._nmax = num_mu, nmax
        self.pairs = ((0, 0), (0, 1), (1, 1))
        self._growth = None
        self._fiducial_models = None
        self._working_models = None
        self.update_growth(fiducial["growth_rate"])
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

    def _prepare_models(self, growth):
        """Prepare one serial model set through the public EFT API."""
        from eftkernels import EFTModel

        models = []
        forest, galaxy = (
            dict(name="F", type="continuous"),
            dict(name="G", type="discrete"),
        )
        for first, second in ((forest, forest), (forest, galaxy), (galaxy, galaxy)):
            model = EFTModel(
                first,
                second,
                num_mu=self._num_mu,
                nmax=self._nmax,
                tree_only=self.order == "tree",
                add_stochastic_terms=False,
            )
            self._update_model(model, growth)
            model.set_output_grid(self.k, self._num_mu)
            models.append(model)
        return models

    def _update_model(self, model, growth):
        """Install this adapter's fixed template at a new growth rate."""
        template = self.fiducial
        model.update_fiducial(
            template["k"],
            template["pk_full"],
            template["pk_nw"],
            growth,
            template["rd_mpc_h"],
        )

    def update_growth(self, growth):
        """Update loop kernels and IR damping through the public EFT API.

        Parameters
        ----------
        growth : float
            Finite dimensionless logarithmic growth rate, at fixed input power.

        Returns
        -------
        None

        Notes
        -----
        Keep the fiducial model set cached so returning from a growth stencil
        does no FFTLog work. A second set is reused for nonfiducial growth,
        bounding memory independently of the number of stencil points. This
        avoids assigning eftkernels' private cosmology-dependent state. Tree
        IR damping deliberately remains fixed at the fiducial growth rate.
        """
        if not np.isfinite(growth):
            raise ValueError("growth must be finite")
        if growth == self._growth or (
            self.order == "tree" and self._growth is not None
        ):
            return
        if self._fiducial_models is None:
            self._fiducial_models = self._prepare_models(growth)
            self.models = self._fiducial_models
        elif growth == self.fiducial["growth_rate"]:
            self.models = self._fiducial_models
        else:
            if self._working_models is None:
                self._working_models = self._prepare_models(growth)
            else:
                for model in self._working_models:
                    self._update_model(model, growth)
            self.models = self._working_models
        self._growth = growth

    def __call__(self, theta_local, z, k, mu, pairs):
        """Evaluate signed clustering means, excluding observational noise.

        Parameters
        ----------
        theta_local : array_like of shape (n_parameter,)
            Local parameter values in the order and units of parameter_names.
        z : float
            Dimensionless template redshift.
        k, mu : ndarray of shape (n_node,)
            Paired positive wavenumbers in h/Mpc and direction cosines on [0, 1].
        pairs : sequence of (int, int)
            Requested field pairs, with F=0, G=1. Reversed pairs are equivalent.

        Returns
        -------
        output : ndarray of float64, shape (n_node, n_pair)
            Signed mean power in (Mpc/h)^3, including pair stochasticity.

        Raises
        ------
        ValueError
            If redshift or nodes are invalid, or AP maps outside the output grid.
        """
        k, mu = np.asarray(k, dtype=np.float64), np.asarray(mu, dtype=np.float64)
        if (
            k.ndim != 1
            or mu.shape != k.shape
            or not np.isfinite(k).all()
            or not np.isfinite(mu).all()
            or np.any(k <= 0)
            or np.any((mu < 0) | (mu > 1))
        ):
            raise ValueError("k must be positive paired nodes and mu must be on [0, 1]")
        if abs(z - self.fiducial["z_fiducial"]) > 1e-12:
            raise ValueError("EFT adapter requires the supplied template redshift")
        parameters = dict(self.defaults)
        parameters.update(zip(self.names, theta_local))
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
        for mapped_k, _, _ in query:
            if np.any(mapped_k < self.k[0]) or np.any(mapped_k > self.k[-1]):
                raise ValueError("AP query outside padded EFT output grid")
        self.update_growth(parameters.get("f", self.fiducial["growth_rate"]))
        output = np.empty((len(k), len(pairs)), dtype=np.float64)
        pair_power = {}
        for column, pair in enumerate(pairs):
            canonical_pair = tuple(sorted(pair))
            if canonical_pair in pair_power:
                output[:, column] = pair_power[canonical_pair]
                continue
            index = self.pairs.index(canonical_pair)
            model = self.models[index]
            output[:, column] = 0.0
            for component, (mapped_k, mapped_mu, jacobian) in zip(
                ("smooth", "wiggle"), query
            ):
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
            pair_power[canonical_pair] = output[:, column]
        return output
