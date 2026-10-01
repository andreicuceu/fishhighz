"""Generated raw populations to two independent bins; entirely synthetic."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from fishhighz.adapters.legacy_inputs import (
    DensityReader,
    SNRReader,
    sample_forest_readers,
)
from fishhighz.fields import ObservedField, PairSelection
from fishhighz.forecast import prepare_bin, run_forecast
from fishhighz.geometry import LYA_REST_ANGSTROM, SPEED_LIGHT_KMS, prepare_geometry
from fishhighz.grids import gauss_legendre_grid
from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
from fishhighz.models.p1d import default_p1d
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import InstrumentResponse, pixel_width_angstrom_to_velocity
from fishhighz.results import diagonal_prior
from fishhighz.survey import BinSpec, ForestInput


def make_readers(root, population):
    """Write synthetic source counts and SNR tables and construct their readers.

    Parameters
    ----------
    root : str or pathlib.Path
        Existing directory receiving generated density and SNR text files.
    population : int
        Population index, either 0 or 1, controlling density and SNR amplitudes.

    Returns
    -------
    readers : tuple of DensityReader and SNRReader
        Readers for source counts per deg^2 per redshift/magnitude cell and
        dimensionless pixel SNR.

    Notes
    -----
    The source cells have explicit redshift spacing 0.5 and magnitude
    spacing 1. Files are generated locally without using bundled survey data.
    """
    root = Path(root)
    density = root / f"density-{population}.txt"
    redshift_grid = np.arange(2, 5, 0.5)
    magnitude_grid = np.arange(20, 24.0)
    np.savetxt(
        density,
        [
            [a, m, (1 + population) * (10 + a + (m - 20) ** 2) * 0.5]
            for a in redshift_grid
            for m in magnitude_grid
        ],
    )
    paths = []
    for m in magnitude_grid:
        path = root / f"snr-{population}-{m}.txt"
        wavelength_grid = np.arange(3500, 5501, 100.0)
        header = f"BAND= r MAG= {m} EXPTIME= 4000 NEXP= 4\nWave " + " ".join(
            f"SN(z={v})" for v in redshift_grid
        )
        np.savetxt(
            path,
            [
                [
                    w,
                    *[
                        (2 + population + 0.2 * a + 0.001 * w) / (1 + 0.1 * (m - 20))
                        for a in redshift_grid
                    ],
                ]
                for w in wavelength_grid
            ],
            header=header,
        )
        paths.append(path)
    return (
        DensityReader(
            density,
            semantics="cell_count_per_deg2",
            target_density=100 * (1 + population),
            z_norm_min=2.15,
            magnitude_bounds=(20, 23),
            label=f"population {population}",
        ),
        SNRReader(paths, smoothing="legacy", label=f"population {population}"),
    )


def run(root=None):
    """Forecast two synthetic bins with fixed weights and nuisance priors.

    Parameters
    ----------
    root : str or pathlib.Path, optional
        Existing directory for generated raw inputs. Default is None, which uses
        an automatically cleaned temporary directory.

    Returns
    -------
    report : dict
        Per-bin and combined Fisher matrices, dimensionless marginalized errors,
        source weights, volumes in (Mpc/h)^3, redshift bounds, and provider-call
        diagnostics.

    Raises
    ------
    AssertionError
        If derivative refinement changes the data Fisher beyond tolerance.

    Notes
    -----
    Writes synthetic raw inputs, prepares each redshift bin once, and compares
    three derivative step sizes with the same nuisance prior.
    """
    if root is None:
        with TemporaryDirectory(prefix="fishhighz-survey-") as directory:
            return run(directory)
    readers = [make_readers(root, i) for i in (0, 1)]
    registry = ParameterRegistry(
        [
            Parameter("A", 1, "target", step=0.001),
            Parameter("b_low", 1, "nuisance", step=0.001),
            Parameter("b_high", 1, "nuisance", step=0.001),
        ]
    )
    fields = [
        ObservedField("fq", "forest", "lya", background="qso"),
        ObservedField("fl", "forest", "lya", background="lbg"),
        ObservedField("g", "galaxy", "galaxy"),
    ]
    base = np.array([[4, 1, -0.5], [1, 3, -0.3], [-0.5, -0.3, 2.0]])

    def model(t, z, k, mu, pairs):
        """Evaluate the synthetic two-forest and galaxy spectra.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider's declared order.
        z : float
            Dimensionless evaluation redshift; unused by this synthetic model.
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines; unused where the model is
            isotropic.
        pairs : ndarray of int, shape (n_pairs, 2)
            Indices of the two observed fields in each requested spectrum.

        Returns
        -------
        power : ndarray of shape (n_nodes, n_pairs)
            Intrinsic three-dimensional power in (Mpc/h)^3; t contains amplitude and
            galaxy bias.
        """
        bias = np.array([1, 1, t[1]])
        biased_power = base * np.outer(bias, bias)
        return t[0] * (1 + k[:, None]) * biased_power[pairs[:, 0], pairs[:, 1]]

    bins = []
    for index, (lo, hi, z_eval) in enumerate([(2.2, 2.5, 2.3), (2.5, 2.7, 2.62)]):
        geometry = prepare_geometry(
            lo,
            hi,
            z_eval=z_eval,
            area_deg2=100,
            h_fid=0.7,
            z_order=8,
            hubble=lambda z: 70 * (1 + z) ** 1.5,
            transverse_distance=lambda z: (
                2 * SPEED_LIGHT_KMS / 70 * (1 - 1 / np.sqrt(1 + z))
            ),
        )
        grid = gauss_legendre_grid([0.02, 0.1, 0.2], k_order=2, mu_order=3, h_fid=0.7)
        selected = (
            [("fq", "fq"), ("fl", "g"), ("g", "g")]
            if index == 0
            else [("fl", "fl"), ("fq", "g"), ("g", "g")]
        )
        selection = PairSelection(fields, selected)
        binding = BoundParameters(
            registry,
            ("amplitude", "bias"),
            {"amplitude": "A", "bias": ("b_low", "b_high")[index]},
        )
        p3d = PreparedP3D(
            registry,
            selection,
            [P3DProvider("external", model, binding, selection.required_pairs)],
        )
        wavelength = LYA_REST_ANGSTROM * (1 + geometry.z_eval)
        responses = {"g": InstrumentResponse(0, 0)}
        forests = {}
        for i, name in enumerate(["fq", "fl"]):
            pixel = 0.8 + 0.2 * i
            responses[name] = InstrumentResponse(
                pixel_width_angstrom_to_velocity(pixel, lambda_obs_angstrom=wavelength),
                20 + 10 * i,
            )
            # Explicit illustrative rest limits, not an automatic source estimator.
            rest_min, rest_max = 1040.0, 1200.0
            z_source = wavelength / np.sqrt(rest_min * rest_max) - 1
            magnitude_grid = np.array([20.0, 21.0, 22.0, 23.0])
            sampled = sample_forest_readers(
                *readers[i],
                geometry,
                responses[name],
                z_source=z_source,
                magnitudes=magnitude_grid,
                pixel_width_angstrom=pixel,
                exposure_count=4,
            )
            forests[name] = ForestInput(
                dict(
                    z_source=z_source,
                    magnitudes=magnitude_grid,
                    quadrature=[0.5, 1, 1, 0.5],
                    rho=sampled["rho"],
                    variance=sampled["variance"],
                    length_velocity=SPEED_LIGHT_KMS * np.log(rest_max / rest_min),
                    method="legacy",
                    iterations=3,
                ),
                default_p1d,
                BoundParameters(registry, (), {}),
                registry.fiducials,
                auxiliary_coordinates=(2.4, 0.00035),
                provenance=sampled["provenance"],
            )
        galaxy_density = readers[0][0].local_galaxy_density(
            geometry, [20, 21, 22, 23], [0.5, 1, 1, 0.5]
        )
        bins.append(
            prepare_bin(
                BinSpec(
                    ("low", "high")[index],
                    geometry,
                    grid,
                    p3d,
                    responses,
                    forests=forests,
                    galaxies={"g": galaxy_density},
                    independent_sampling=True,
                )
            )
        )
    prior = diagonal_prior(registry, {"b_low": 0.5, "b_high": 0.5})
    results = [
        run_forecast(bins, batch_size=5, step_scale=s, prior_fisher=prior)
        for s in [1, 0.5, 0.25]
    ]
    unprior = run_forecast(bins)
    np.testing.assert_allclose(
        results[0].combined.data_fisher,
        results[-1].combined.data_fisher,
        rtol=5e-12,
        atol=0,
    )
    return dict(
        interpretation="synthetic illustration, not a DESI-2 forecast",
        assumptions="independent redshift interiors; common volume per bin; explicit local galaxy density; independent sampling",
        bin_ids=results[0].bin_ids,
        parameter_ids=registry.ids,
        per_bin_data=[b.result.data_fisher.tolist() for b in results[0].bins],
        combined_data=results[0].combined.data_fisher.tolist(),
        prior=prior.tolist(),
        total=results[0].combined.total_fisher.tolist(),
        marginalized_errors=results[0].combined.marginalized_errors().tolist(),
        data_null_directions=unprior.combined.diagnostics.null_directions.tolist(),
        preparation_calls=[dict(b.diagnostics["p3d_calls"]) for b in bins],
        derivative_calls=[
            [(c.provider, c.model, c.jacobian) for c in b.calls]
            for b in results[0].bins
        ],
        weights=[
            {name: w.weights.tolist() for name, w in b.weights.items()} for b in bins
        ],
        volume=[b.geometry.volume for b in bins],
        bin_bounds=[(b.geometry.z_min, b.geometry.z_max) for b in bins],
        evaluation_redshifts=[b.geometry.z_eval for b in bins],
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
