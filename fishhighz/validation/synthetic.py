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
    """Compare correlated selected Fisher with independent full-field slicing."""
    select = selection(case)
    registry = ParameterRegistry([Parameter("A", 1, "target", step=0.001)])
    b = np.array([-0.4, 2, 3, 1.5, -0.3])
    signal = np.outer(b, b) + np.diag([0.2, 0.3, 0.4, 0.5, 0.6])
    ids = np.array([FIELD_ORDER.index(f.id) for f in select.fields])
    local = signal[np.ix_(ids, ids)]

    def model(theta, z, k, mu, pairs):
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
        s = (1 + k) * signal
        total = s + 2 * np.eye(5)
        covariance = np.array(
            [
                [
                    (total[i, m] * total[j, n] + total[i, n] * total[j, m]) / modes
                    for m, n in full_pairs
                ]
                for i, j in full_pairs
            ]
        )
        sub = covariance[np.ix_(indices, indices)]
        jac = np.array([s[i, j] for i, j in wanted])
        fisher += jac @ np.linalg.solve(sub, jac)
        covariances.append(sub)
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
    """All seven tiny synthetic forecasts, without local files or reference imports."""
    return {
        case: {**report, "F_AA": float(arrays["fisher"][0, 0])}
        for case in CASE_IDS
        for arrays, report in [run_case(case)]
    }


def evidence_payload(task, *, null=False):
    """Tiny semantically complete schema-2 oracle fixture, no external assets."""
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
    k, mu, _ = grid_nodes(settings)
    n = len(task["fields"])
    bias = np.arange(1, n + 1, dtype=float)
    bias[0] *= -1
    matrix = np.outer(bias, bias) + np.eye(n)
    required = np.array(task["required_pairs"])
    selected = np.array(task["selected_pairs"])
    total = np.tile(matrix[required[:, 0], required[:, 1]], (len(k), 1))
    base = np.tile(matrix[selected[:, 0], selected[:, 1]], (len(k), 1))
    targets = np.ones((len(k), len(task["parameters"])))
    if len(task["parameters"]) == 2:
        targets[:, 1] = 0 if null else mu**2
    return assemble(task, total, base[:, :, None] * targets[:, None, :], settings)


def convergence_payload(task, *, null=False):
    """Analytic constant-information control with distinct recorded refinements."""
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
    k, mu, _ = grid_nodes(settings)
    n = len(task["fields"])
    matrix = np.eye(n) + np.ones((n, n))
    pairs = np.asarray(task["required_pairs"])
    total = np.tile(matrix[pairs[:, 0], pairs[:, 1]], (len(k), 1))
    j = np.ones((len(k), len(task["selected_pairs"]), 2))
    j[:, :, 1] = 0 if null else mu[:, None] ** 2
    arrays, report = assemble(task, total, j, settings)

    class AnalyticStudy:
        selection = selection(task["case"])
        _prepared = {}

        def evaluate(self, task, controls):
            return dict(arrays), dict(report)

    # This fixture asserts the analytic information limit, not a real survey.
    result, report = study(AnalyticStudy(), task)
    return result, report
