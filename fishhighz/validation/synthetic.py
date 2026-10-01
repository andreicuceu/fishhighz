"""Independent five-field covariance oracle and seven-case synthetic forecast."""

import numpy as np

from ..forecast import prepare_bin, run_bin
from ..geometry import prepare_geometry
from ..grids import gauss_legendre_grid
from ..models.external import BoundParameters, P3DProvider, PreparedP3D
from ..parameters import Parameter, ParameterRegistry
from ..response import InstrumentResponse
from ..survey import BinSpec
from .cases import CASE_IDS, selection

FIELD_ORDER = ("lya(qso)", "qso", "lbg", "lae", "lya(lbg)")


def run_case(case):
    """Compare correlated selected Fisher with independent full-field slicing.

    Parameters
    ----------
    case : str
        Identifier of one of the seven original DESI-2 validation
        configurations.

    Returns
    -------
    arrays : dict of str to ndarray
        Scalar amplitude Fisher matrix, shape (1, 1), and error, shape (1,).
    report : dict
        Original field/pair ordering and independent-oracle description.

    Notes
    -----
    Asserts covariance and information agreement against direct five-field reconstruction using tiny synthetic arrays.
    """
    select = selection(case)
    registry = ParameterRegistry([Parameter("A", 1, "target", step=0.001)])
    field_biases = np.array([-0.4, 2, 3, 1.5, -0.3])
    signal = np.outer(field_biases, field_biases) + np.diag([0.2, 0.3, 0.4, 0.5, 0.6])
    ids = np.array([FIELD_ORDER.index(f.id) for f in select.fields])
    local = signal[np.ix_(ids, ids)]

    def model(theta, z, k, mu, pairs):
        """Evaluate a correlated synthetic amplitude model on paired Fourier nodes.

        Parameters
        ----------
        theta : array_like, shape (n_parameter,)
            Model parameter values in the declared local parameter order.
        z : float
            Dimensionless evaluation redshift.
        k : array_like
            Comoving Fourier wavenumbers in h/Mpc; array shape follows the model or
            paired grid.
        mu : array_like
            Dimensionless line-of-sight direction cosines aligned with the Fourier
            grid.
        pairs : array_like, shape (n_pair, 2)
            Ordered pairs of integer field indices.

        Returns
        -------
        power : ndarray, shape (n_cell, n_pair)
            Synthetic pair power scaled by the dimensionless amplitude.
        """
        return theta[0] * (1 + k[:, None]) * local[pairs[:, 0], pairs[:, 1]]

    p3d = PreparedP3D(
        registry,
        select,
        [
            P3DProvider(
                "synthetic",
                model,
                BoundParameters(registry, ("A",), {"A": "A"}),
                select.required_pairs,
            )
        ],
    )
    geometry = prepare_geometry(
        2.47,
        2.705,
        z_eval=2.58,
        area_deg2=50,
        h_fid=0.7,
        hubble=lambda z: np.full_like(z, 250),
        transverse_distance=lambda z: np.full_like(z, 5000),
        z_order=4,
    )
    grid = gauss_legendre_grid([0.01, 0.1, 0.5], k_order=2, mu_order=3, h_fid=0.7)
    noise = np.tile(
        np.where(select.required_pairs[:, 0] == select.required_pairs[:, 1], 2, 0),
        (len(grid.k_flat), 1),
    )
    prepared = prepare_bin(
        BinSpec(
            case,
            geometry,
            grid,
            p3d,
            {f.id: InstrumentResponse(0, 0) for f in select.fields},
            full_noise=noise,
        )
    )
    result = run_bin(prepared)
    # Build the complete five-field matrix first, then slice named spectra.
    full_pairs = [(i, j) for i in range(5) for j in range(i, 5)]
    wanted = [(ids[i], ids[j]) for i, j in select.selected_pairs]
    indices = [full_pairs.index(tuple(sorted(pair))) for pair in wanted]
    fisher = 0.0
    covariances = []
    for k, modes in zip(grid.k_flat, prepared.modes):
        signal_power = (1 + k) * signal
        total = signal_power + 2 * np.eye(5)
        covariance = np.array(
            [
                [
                    (total[i, m] * total[j, n] + total[i, n] * total[j, m]) / modes
                    for m, n in full_pairs
                ]
                for i, j in full_pairs
            ]
        )
        selected_covariance = covariance[np.ix_(indices, indices)]
        amplitude_derivative = np.array([signal_power[i, j] for i, j in wanted])
        fisher += amplitude_derivative @ np.linalg.solve(
            selected_covariance, amplitude_derivative
        )
        covariances.append(selected_covariance)
    np.testing.assert_allclose(
        prepared.factors @ prepared.factors.swapaxes(-1, -2),
        covariances,
        rtol=5e-13,
        atol=0,
    )
    np.testing.assert_allclose(
        result.result.data_fisher, [[fisher]], rtol=5e-12, atol=0
    )
    return dict(
        fisher=np.array([[fisher]]), errors=np.array([1 / np.sqrt(fisher)])
    ), dict(
        case=case,
        fields=[f.id for f in select.fields],
        selected_pairs=select.selected_pairs.tolist(),
        required_pairs=select.required_pairs.tolist(),
        oracle="independent full-five-field covariance slicing",
        synthetic=True,
    )


def run():
    """All seven tiny synthetic forecasts, without local files or reference imports.

    Returns
    -------
    results : dict
        Seven case summaries and independent scalar amplitude information.
    """
    return {
        case: {**report, "F_AA": float(arrays["fisher"][0, 0])}
        for case in CASE_IDS
        for arrays, report in [run_case(case)]
    }


def evidence_payload(task, *, null=False):
    """Tiny semantically complete schema-2 oracle fixture, no external assets.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    null : bool
        Set the second parameter derivative to zero to produce a null direction.
        Default is ``False``.

    Returns
    -------
    arrays : dict of str to ndarray
        Tiny grid, power, derivative, covariance and Fisher evidence.
    report : dict
        Bound synthetic settings and validation metadata.
    """
    from .schema import assemble, grid_nodes

    settings = dict(
        profile=task["profile"],
        parameters=task["parameters"],
        bounds=task["bounds"],
        fields=[f["id"] for f in task["fields"]],
        grid=dict(
            kind="gauss_legendre",
            volume=1000.0,
            h_fid=0.7,
            k_intervals=2,
            k_order=2,
            mu_order=3,
        ),
    )
    k_grid, mu_grid, _ = grid_nodes(settings)
    n_fields = len(task["fields"])
    bias = np.arange(1, n_fields + 1, dtype=float)
    bias[0] *= -1
    matrix = np.outer(bias, bias) + np.eye(n_fields)
    required = np.array(task["required_pairs"])
    selected = np.array(task["selected_pairs"])
    total = np.tile(matrix[required[:, 0], required[:, 1]], (len(k_grid), 1))
    base = np.tile(matrix[selected[:, 0], selected[:, 1]], (len(k_grid), 1))
    targets = np.ones((len(k_grid), len(task["parameters"])))
    if len(task["parameters"]) == 2:
        targets[:, 1] = 0 if null else mu_grid**2
    return assemble(task, total, base[:, :, None] * targets[:, None, :], settings)


def convergence_payload(task, *, null=False):
    """Analytic constant-information control with distinct recorded refinements.

    Parameters
    ----------
    task : dict
        Declared case, bin, selected field pairs, parameter order and validation
        thresholds.
    null : bool
        Set the second derivative to zero to test constrained-coordinate
        handling. Default is ``False``.

    Returns
    -------
    arrays : dict of str to ndarray
        Constant-information payload with distinct recorded trial controls.
    report : dict
        Analytic refinement schedule and convergence outcomes.
    """
    from .schema import assemble, grid_nodes
    from .study import DEFAULT, study

    settings = dict(
        profile=task["profile"],
        parameters=task["parameters"],
        bounds=task["bounds"],
        fields=[f["id"] for f in task["fields"]],
        controls=dict(DEFAULT),
        grid=dict(
            kind="gauss_legendre",
            volume=1000.0,
            h_fid=0.7,
            k_intervals=128,
            k_order=4,
            mu_order=32,
        ),
    )
    k_grid, mu_grid, _ = grid_nodes(settings)
    n_fields = len(task["fields"])
    matrix = np.eye(n_fields) + np.ones((n_fields, n_fields))
    pairs = np.asarray(task["required_pairs"])
    total = np.tile(matrix[pairs[:, 0], pairs[:, 1]], (len(k_grid), 1))
    observed_jacobian = np.ones((len(k_grid), len(task["selected_pairs"]), 2))
    observed_jacobian[:, :, 1] = 0 if null else mu_grid[:, None] ** 2
    arrays, report = assemble(task, total, observed_jacobian, settings)

    class AnalyticStudy:
        selection = selection(task["case"])
        _prepared = {}

        def evaluate(self, task, controls):
            """Return the analytic constant-information state for a declared trial.

            Parameters
            ----------
            task : dict
                Declared case, bin, selected field pairs, parameter order and validation
                thresholds.
            controls : dict
                Quadrature orders, grid subdivisions, derivative step and applicable
                forest-weight convergence controls.

            Returns
            -------
            arrays : dict
                Shallow copy of the fixed synthetic numerical payload.
            report : dict
                Shallow copy of the fixed synthetic report.
            """
            return dict(arrays), dict(report)

    # This fixture asserts the analytic information limit, not a real survey.
    result, report = study(AnalyticStudy(), task)
    return result, report
