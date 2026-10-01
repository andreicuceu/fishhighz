"""Independent scalar limits and revised selected-forecast contracts."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from test_weights import prepare

from fishhighz.kernels.full_sum_weights import coefficients, seed, solve, update
from fishhighz.validation import trials
from fishhighz.validation.compatibility_weights import WeightInputs
from fishhighz.validation.profile_definitions import REVISION, forecast_selection
from fishhighz.validation.schema import request, validate_request
from fishhighz.validation.study import study


def inputs():
    """Construct a heterogeneous two-sample forest-weight fixture.

    Returns
    -------
    inputs : WeightInputs
        Source densities, quadrature, noise, and fixed auxiliary signal/P1D.
    """
    return WeightInputs(
        np.array([20.0, 21.0]),
        np.array([0.2, 0.7]),
        np.array([0.3, 0.8]),
        np.array([1.0, 7.0]),
        9.0,
        0.5,
        2.0,
        3.0,
    )


@pytest.mark.parametrize("variant", ["sum_historical", "sum_aliasing"])
def test_two_cell_scalar_update_and_amplitude(variant):
    """Check two cell scalar update and amplitude.

    Parameters
    ----------
    variant : str
        Forest-weight recurrence variant, supplied by pytest parametrization.
    """
    weight_inputs = inputs()
    source_weights = seed(weight_inputs)
    for _ in range(2):
        i1 = sum(
            float(d * q * a)
            for d, q, a in zip(
                weight_inputs.density, weight_inputs.quadrature, source_weights
            )
        )
        i2 = sum(
            float(d * q * a * a)
            for d, q, a in zip(
                weight_inputs.density, weight_inputs.quadrature, source_weights
            )
        )
        effective_signal = weight_inputs.signal + (
            weight_inputs.p1d / (weight_inputs.length * i1)
            if variant == "sum_historical"
            else weight_inputs.p1d * i2 / (weight_inputs.length * i1 * i1)
        )
        expected = [
            effective_signal
            / (
                effective_signal
                + weight_inputs.pixel * float(v) / (weight_inputs.length * i1)
            )
            for v in weight_inputs.variance
        ]
        source_weights = update(weight_inputs, source_weights, variant)
        np.testing.assert_allclose(source_weights, expected, rtol=5e-15)
    np.testing.assert_allclose(
        coefficients(weight_inputs, source_weights)[-2:],
        coefficients(weight_inputs, 2 * source_weights)[-2:],
        rtol=5e-15,
    )
    assert not np.allclose(
        update(weight_inputs, source_weights, variant),
        update(weight_inputs, 2 * source_weights, variant),
    )
    reverse = replace(
        weight_inputs,
        density=weight_inputs.density[::-1],
        quadrature=weight_inputs.quadrature[::-1],
        variance=weight_inputs.variance[::-1],
    )
    np.testing.assert_allclose(
        update(reverse, source_weights[::-1], variant)[::-1],
        update(weight_inputs, source_weights, variant),
        rtol=5e-15,
    )


def test_homogeneous_quadratic_fixed_point_and_zero_noise():
    """Check homogeneous quadratic fixed point and zero noise."""
    weight_inputs = inputs()
    weight_inputs = replace(weight_inputs, variance=np.full(2, 2.0))
    state = solve(weight_inputs, "sum_historical", rtol=1e-10)
    quadratic_coefficient = (
        weight_inputs.signal
        * weight_inputs.length
        * sum(weight_inputs.density * weight_inputs.quadrature)
    )
    # c w^2 + (B+lp*v-c) w - B = 0.
    linear = weight_inputs.p1d + weight_inputs.pixel * 2 - quadratic_coefficient
    root = (
        -linear + np.sqrt(linear**2 + 4 * quadratic_coefficient * weight_inputs.p1d)
    ) / (2 * quadratic_coefficient)
    assert state["status"] == "converged"
    np.testing.assert_allclose(state["weights"], root, rtol=1e-12)
    zero = solve(replace(weight_inputs, variance=np.zeros(2)), "sum_historical")
    assert zero["updates"] == 6
    np.testing.assert_array_equal(zero["weights"], np.ones(2))
    assert zero["coefficients"][-1] == 0
    for p in (0.0, 1e-12):
        small = replace(weight_inputs, signal=p)
        np.testing.assert_allclose(
            update(small, seed(small), "sum_historical"), seed(small), rtol=1e-11
        )
    assert (
        solve(replace(weight_inputs, signal=0.0), "sum_historical")["status"]
        == "ineligible"
    )
    assert solve(weight_inputs, "sum_historical", max_updates=3)["status"] == "capped"


def test_strict_public_methods_and_cap():
    """Check strict public methods and cap."""
    for method in ("early_lyaforecast", "mcdonald"):
        result = prepare(method=method, iterations=None)
        assert result.convergence["status"] == "converged"
        assert result.convergence["updates"] >= 6
        with pytest.raises(ValueError, match="capped"):
            prepare(method=method, iterations=None, max_updates=3)
        with pytest.raises(ValueError, match="nonnegative"):
            prepare(method=method, iterations=None, rho=[1, -1, 1])
        fixed = prepare(method=method, iterations=3)
        assert fixed.convergence["status"] == "fixed_count"


def test_selected_bin1_and_historical_identity():
    """Check selected bin1 and historical identity."""
    case = "lya_qso_lbg_lae_15x2pt"
    historical = request(case, 0, "accuracy")
    revised = request(case, 0, "accuracy", recipe_revision=REVISION)
    assert len(historical["selected_pairs"]) == 15
    assert len(revised["selected_pairs"]) == len(revised["required_pairs"]) == 3
    assert {
        revised["fields"][i]["id"].lower() for p in revised["selected_pairs"] for i in p
    } == {"lya(qso)", "qso"}
    assert sum(len(forecast_selection(case, i).selected_pairs) for i in range(6)) == 78
    validate_request(revised)


class AdaptiveStudy:
    weight_method = "early_lyaforecast"
    selection = SimpleNamespace(
        fields=[SimpleNamespace(kind="forest")], selected_pairs=np.array([[0, 0]])
    )

    def __init__(self):
        """Initialize the synthetic AdaptiveStudy fixture.

        Notes
        -----
        Sets the instance state used by the enclosing test; no scientific calculation is run.
        """
        self._prepared = {}

    def evaluate(self, task, controls):
        """Evaluate the synthetic Fisher trial for the supplied numerical controls.

        Parameters
        ----------
        task : dict
            Synthetic forecast request including case, redshift-bin index, profile,
            and pair selection.
        controls : dict
            Numerical quadrature and weighting controls for this synthetic trial.

        Returns
        -------
        payload : tuple
            Synthetic numerical arrays and validation report, including the
            requested test modification.

        Notes
        -----
        Uses small analytic matrices to exercise validation control flow; it does not run a survey forecast.
        """
        fisher = np.diag([2.0, 1.0])
        return dict(fisher=fisher, pair_fisher=fisher[None]), dict(
            settings=dict(
                controls=dict(controls),
                grid=dict(
                    volume=1.0,
                    k_intervals=controls["k_intervals"],
                    mu_order=controls["mu_order"],
                    k_order=4,
                ),
                forest_weighting=dict(
                    method=self.weight_method,
                    forests={"forest": dict(result=dict(status="converged"))},
                ),
            )
        )


def test_adaptive_trial_v3_replays_and_tightens_tolerance():
    """Check adaptive trial v3 replays and tightens tolerance."""
    arrays, report = study(AdaptiveStudy(), {})
    assert report["trial_contract"]["version"] == 3
    assert report["actual_levels"]["weights"] == [1e-4, 1e-5]
    assert report["final_controls"]["weight_rtol"] == 1e-5
    assert all("iterations" not in c for c in report["study_controls"])
    assert trials.validate(arrays, report)


def test_auxiliary_units_field_reuse_and_frozen_derivatives():
    """Check auxiliary units field reuse and frozen derivatives."""
    from test_weights import geometry

    from fishhighz.fields import ObservedField, PairSelection
    from fishhighz.forecast import prepare_bin, run_bin
    from fishhighz.grids import gauss_legendre_grid
    from fishhighz.models.external import BoundParameters, P3DProvider, PreparedP3D
    from fishhighz.parameters import Parameter, ParameterRegistry
    from fishhighz.response import InstrumentResponse, velocity_response
    from fishhighz.survey import BinSpec, ForestInput

    fields = [
        ObservedField("a", "forest", "lya", "qso"),
        ObservedField("b", "forest", "lya", "lbg"),
    ]
    selection = PairSelection(fields)
    registry = ParameterRegistry([Parameter("amplitude", 1.0, "target", step=0.001)])
    binding = BoundParameters(registry, ["A"], {"A": "amplitude"})
    p1d_binding = BoundParameters(registry, (), {})
    auxiliary_calls = []

    def model(t, z, k, mu, pairs):
        """Evaluate the synthetic spectrum used by the enclosing regression test.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        k : ndarray of shape (n_nodes,)
            Comoving wavenumbers in h/Mpc.
        mu : ndarray of shape (n_nodes,)
            Dimensionless line-of-sight direction cosines.
        pairs : ndarray of int, shape (n_pairs, 2)
            Observed-field indices defining the requested spectra.

        Returns
        -------
        power : ndarray of shape (n_nodes, n_pairs)
            Synthetic intrinsic power in (Mpc/h)^3 before response and noise.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        if len(k) == 1:
            auxiliary_calls.append((pairs.copy(), t.copy()))
        base = np.array([[4.0, 1.0], [1.0, 9.0]])
        return np.broadcast_to(
            t[0] * base[pairs[:, 0], pairs[:, 1]], (len(k), len(pairs))
        )

    p1d_calls = []

    def p1d(t, z, q):
        """Evaluate an independent synthetic one-dimensional forest spectrum.

        Parameters
        ----------
        t : ndarray of shape (n_parameters,)
            Local model parameters in the provider binding order.
        z : float or ndarray
            Dimensionless redshift.
        q : ndarray of shape (n_nodes,)
            Line-of-sight velocity wavenumbers in s/km.

        Returns
        -------
        power : ndarray of shape (n_nodes,)
            Intrinsic one-dimensional power in km/s.

        Notes
        -----
        Appends to the enclosing test call log so provider dispatch can be checked.
        """
        p1d_calls.append(q.copy())
        return np.full_like(q, 3.0)

    prepared_p3d = PreparedP3D(
        registry,
        selection,
        [P3DProvider("model", model, binding, selection.required_pairs)],
    )
    geom = geometry()
    responses = {"a": InstrumentResponse(0.5, 2.0), "b": InstrumentResponse(0.7, 5.0)}
    forests = {}
    for i, field in enumerate(fields):
        forests[field.id] = ForestInput(
            dict(
                method="early_lyaforecast",
                z_source=3.0,
                magnitudes=[20.0, 21.0],
                quadrature=[0.3, 0.8],
                rho=[0.2, 0.7],
                variance=[1.0 + i, 7.0 + i],
                length_velocity=9.0,
            ),
            p1d,
            p1d_binding,
            registry.fiducials,
            auxiliary_coordinates=(2.4, 0.00035),
        )
    grid = gauss_legendre_grid([0.02, 0.1], k_order=2, mu_order=2, h_fid=0.7)
    prepared = prepare_bin(
        BinSpec(
            "two",
            geom,
            grid,
            prepared_p3d,
            responses,
            forests=forests,
            galaxies={},
            independent_sampling=True,
        )
    )
    assert len(auxiliary_calls) == 2 and len(p1d_calls) == 4
    for i, field in enumerate(fields):
        weights = prepared.weights[field.id]
        response_amplitude = velocity_response(
            [0.00035],
            pixel_width_velocity=responses[field.id].pixel_width_velocity,
            gaussian_sigma_velocity=responses[field.id].gaussian_sigma_velocity,
        )[0]
        assert weights.signal == pytest.approx(
            [4.0, 9.0][i]
            * response_amplitude
            * response_amplitude
            * geom.a_v
            / geom.d_deg**2
        )
        assert weights.alias == pytest.approx(
            3.0 * response_amplitude * response_amplitude
        )
        noise = (
            (weights.A * 3.0 * prepared.response[:, i] ** 2 + weights.P_pixel)
            * geom.d_deg**2
            / geom.a_v
        )
        column = list(map(tuple, selection.required_pairs)).index((i, i))
        np.testing.assert_allclose(prepared.noise[:, column], noise, rtol=5e-15)
    snapshots = {
        name: source_weights.weights.copy()
        for name, source_weights in prepared.weights.items()
    }
    run_bin(prepared)
    assert len(auxiliary_calls) == 2 and len(p1d_calls) == 4
    for name, response_amplitude in prepared.weights.items():
        np.testing.assert_array_equal(response_amplitude.weights, snapshots[name])
    # Independent auto/cross Wick variances.
    total_power = prepared.total
    covariance = prepared.factors @ prepared.factors.swapaxes(-1, -2)
    np.testing.assert_allclose(
        covariance[:, 0, 0], 2 * total_power[:, 0] ** 2 / prepared.modes, rtol=5e-15
    )
    np.testing.assert_allclose(
        covariance[:, 1, 1],
        (total_power[:, 0] * total_power[:, 2] + total_power[:, 1] ** 2)
        / prepared.modes,
        rtol=5e-15,
    )


def test_accuracy_recipe_records_own_auxiliary_and_stopping(monkeypatch):
    """Check accuracy recipe records own auxiliary and stopping.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    from test_weighting_w12 import profile_fixture

    from fishhighz.forecast import prepare_bin
    from fishhighz.validation import accuracy

    recipe, _, calls = profile_fixture(monkeypatch)
    monkeypatch.setattr(accuracy, "prepare_bin", prepare_bin)
    recipe.weight_method = "early_lyaforecast"
    controls = dict(
        k_intervals=1,
        mu_order=2,
        z_order=2,
        magnitude_order=2,
        step=2.5e-4,
        weight_rtol=1e-5,
    )
    prepared, settings = recipe.prepare(0, controls)
    assert len(calls) == 3  # each auto auxiliary and the forecast grid
    rows = settings["forest_weighting"]["forests"]
    assert set(rows) == {"lya(qso)", "lya(lbg)"}
    for name, row in rows.items():
        assert row["result"]["status"] == "converged"
        assert row["iterations"]["stopping"]["rtol"] == 1e-5
        assert row["auxiliary"]["P"] == prepared.weights[name].signal
        assert row["auxiliary"]["B"] == prepared.weights[name].alias
        assert row["A"] == prepared.weights[name].A


def test_revised_runner_dispatches_exact_selected_profiles(tmp_path, monkeypatch):
    """Check revised runner dispatches exact selected profiles.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Isolated temporary directory supplied by pytest; generated test files
        are written here.
    monkeypatch : pytest.MonkeyPatch
        Fixture that restores patched callables, attributes, and environment
        variables after the test.
    """
    from fishhighz.validation import profiles, revised_compatibility
    from fishhighz.validation.evidence import modern_requests

    source = tmp_path / "source"
    source.write_text("fixture")
    ident = {"resources": {}, "reference": {"sources": {}}}
    monkeypatch.setattr(profiles, "verify_inventory", lambda *a: None)
    monkeypatch.setattr(profiles, "provenance", lambda *a: ident)
    monkeypatch.setattr(profiles, "background", lambda *a: (None, None))
    seen = []
    monkeypatch.setattr(
        revised_compatibility,
        "run",
        lambda task, *a: (seen.append(task) or {}, {"passed": True}),
    )

    class Recipe:
        def __init__(self, *a, **kw):
            """Initialize the synthetic Recipe fixture.

            Parameters
            ----------
            *a : tuple
                Positional arguments forwarded to the original callable or accepted by
                the test callback.
            **kw : dict
                Keyword options forwarded to the original callable or inspected by the
                test callback.

            Notes
            -----
            Sets the instance state used by the enclosing test; no scientific calculation is run.
            """
            assert kw == dict(
                weight_method="early_lyaforecast", recipe_revision=REVISION
            )
            self._samples = {}
            self._prepared = {}

        def study(self, task):
            """Record a synthetic task and return its fixed convergence result.

            Parameters
            ----------
            task : dict
                Synthetic forecast request including case, redshift-bin index, profile,
                and pair selection.

            Returns
            -------
            payload : tuple
                Synthetic numerical arrays and validation report, including the
                requested test modification.

            Notes
            -----
            Appends to the enclosing test call log so provider dispatch can be checked.
            """
            seen.append(task)
            return dict(
                fisher=np.eye(2), pair_fisher=np.tile(np.eye(2), (3, 1, 1))
            ), dict(passed=True, final_controls={})

    monkeypatch.setattr(profiles, "AccuracyRecipe", Recipe)

    def execute(out, **kw):
        """Dispatch only the requested synthetic profile tasks.

        Parameters
        ----------
        out : pathlib.Path
            Directory containing the synthetic evidence bundle.
        **kw : dict
            Keyword options forwarded to the original callable or inspected by the
            test callback.

        Returns
        -------
        tasks : list of dict
            Generated requests in dispatch order.
        """
        tasks = modern_requests(
            kw["suite"],
            kw["cases"],
            kw["bin_indices"],
            profiles=kw["profiles"],
            recipe_revision=kw["recipe_revision"],
        )
        for task in tasks:
            kw["worker"](task)
        return tasks

    monkeypatch.setattr(profiles, "execute", execute)
    tasks = profiles.run(
        tmp_path / "out",
        reference=tmp_path,
        template=source,
        reference_bundle=tmp_path,
        wheel=source,
        compatibility_bundle=tmp_path,
        bin_indices=[0],
        suite="full",
    )
    assert [t["profile"] for t in tasks] == [
        "full-compatibility",
        "fixed-compatibility",
        "accuracy",
    ]
    assert len(seen) == 3
    assert all(
        len(t["selected_pairs"]) == 3 and t["recipe_revision"] == REVISION for t in seen
    )


def test_angular_velocity_and_comoving_recurrence_equivalence():
    """Check angular velocity and comoving recurrence equivalence."""
    angular = inputs()
    a_v, d_deg = 73.0, 61.0
    comoving = replace(
        angular,
        density=angular.density * a_v / d_deg**2,
        length=angular.length / a_v,
        pixel=angular.pixel / a_v,
        signal=angular.signal * d_deg**2 / a_v,
        p1d=angular.p1d / a_v,
    )
    for variant in ("sum_historical", "sum_aliasing"):
        angular_solution = solve(angular, variant)
        comoving_solution = solve(comoving, variant)
        assert angular_solution["status"] == comoving_solution["status"] == "converged"
        assert angular_solution["updates"] == comoving_solution["updates"]
        np.testing.assert_allclose(
            angular_solution["weights"], comoving_solution["weights"], rtol=5e-15
        )
        np.testing.assert_allclose(
            comoving_solution["coefficients"][-2:],
            angular_solution["coefficients"][-2:] * [d_deg**2, d_deg**2 / a_v],
            rtol=5e-15,
        )
