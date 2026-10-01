"""Auto-only routes, independent P1D, immutable states and domain failures."""

import numpy as np
import pytest
from test_weights import FIELD, RESPONSE, geometry

from fishhighz.fields import ObservedField, PairSelection
from fishhighz.models.external import (
    BoundParameters,
    P3DProvider,
    PreparedP3D,
    evaluate_p3d,
)
from fishhighz.parameters import Parameter, ParameterRegistry
from fishhighz.response import InstrumentResponse, velocity_response
from fishhighz.weights import prepare_forest_weights, sample_auxiliary


def test_external_calls_and_independent_substitution():
    """Check external calls and independent substitution."""
    fields = [
        ObservedField("unused", "galaxy", "g"),
        FIELD,
        ObservedField("g", "galaxy", "g"),
    ]
    selection = PairSelection(fields, [("forest", "g")])
    registry = ParameterRegistry([Parameter("a", 2.0, "target", step=0.001)])
    binding = BoundParameters(registry, ["norm"], {"norm": "a"})
    calls, p1d_calls = [], []

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
        calls.append(pairs.copy())
        base = np.outer([1, -2, 3], [1, -2, 3])
        return np.broadcast_to(
            t[0] * base[pairs[:, 0], pairs[:, 1]], (len(k), len(pairs))
        )

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
        return np.full(q.shape, t[0])

    prepared = PreparedP3D(
        registry,
        selection,
        [P3DProvider("external", model, binding, selection.required_pairs)],
    )
    bin_geometry = geometry()
    original_auxiliary = sample_auxiliary(
        FIELD,
        bin_geometry,
        RESPONSE,
        prepared,
        [2],
        p1d_model=p1d,
        p1d_parameters=binding,
        theta_p1d=[3],
        k_t_deg=2.4,
        k_p_velocity=0.00035,
    )
    np.testing.assert_array_equal(calls[0], [[1, 1]])
    velocity_response_amplitude = velocity_response(
        [0.00035], pixel_width_velocity=0.5, gaussian_sigma_velocity=0
    )[0]
    assert original_auxiliary.signal == pytest.approx(
        8
        * velocity_response_amplitude
        * velocity_response_amplitude
        * bin_geometry.a_v
        / bin_geometry.d_deg**2,
        rel=5e-15,
    )
    assert original_auxiliary.alias == pytest.approx(
        3 * velocity_response_amplitude * velocity_response_amplitude, rel=5e-15
    )
    args = dict(
        z_source=3,
        magnitudes=[20, 21, 22],
        quadrature=[1] * 3,
        rho=[1, 2, 4],
        variance=[1, 4, 9],
        length_velocity=10,
    )
    for count in [0, 1, 3, 6]:
        prepare_forest_weights(
            FIELD,
            bin_geometry,
            RESPONSE,
            **args,
            method="legacy",
            iterations=count,
            auxiliary=original_auxiliary,
        )
    assert len(calls) == len(p1d_calls) == 1
    supplied = prepare_forest_weights(
        FIELD, bin_geometry, RESPONSE, **args, method="supplied", weights=[1, 2, 3]
    )
    assert len(calls) == len(p1d_calls) == 1
    changed_p1d_auxiliary = sample_auxiliary(
        FIELD,
        bin_geometry,
        RESPONSE,
        prepared,
        [2],
        p1d_model=p1d,
        p1d_parameters=binding,
        theta_p1d=[6],
        k_t_deg=2.4,
        k_p_velocity=0.00035,
    )
    assert changed_p1d_auxiliary.signal == original_auxiliary.signal
    assert changed_p1d_auxiliary.alias == 2 * original_auxiliary.alias
    old = prepare_forest_weights(
        FIELD,
        bin_geometry,
        RESPONSE,
        **args,
        method="legacy",
        iterations=0,
        auxiliary=original_auxiliary,
    )
    new = prepare_forest_weights(
        FIELD,
        bin_geometry,
        RESPONSE,
        **args,
        method="legacy",
        iterations=0,
        auxiliary=changed_p1d_auxiliary,
    )
    assert np.all(new.weights > old.weights)
    assert not np.allclose(new.weights, supplied.weights)
    intrinsic_power = evaluate_p3d(prepared, [2], 2.4, [0.1], [0.5])
    assert intrinsic_power[0, 1] < 0
    for state in (original_auxiliary.theta_p3d, original_auxiliary.theta_p1d):
        with pytest.raises(ValueError):
            state.flags.writeable = True
    with pytest.raises(ValueError, match="mismatch"):
        prepare_forest_weights(
            FIELD,
            geometry(h=0.8),
            RESPONSE,
            **args,
            method="legacy",
            iterations=3,
            auxiliary=original_auxiliary,
        )


def test_builtin_and_domain():
    """Check builtin and domain."""
    from test_kaiser import fields, model, template

    observed_fields = fields()
    registry = ParameterRegistry([Parameter("unused", 1, "target", step=0.001)])
    binding = BoundParameters(registry, (), {})
    selection = PairSelection(observed_fields)
    bin_geometry = geometry()
    for narrow in [False, True]:
        builtin = model(template(lo=0.1, hi=0.3)) if narrow else model()
        prepared = PreparedP3D(
            registry,
            selection,
            [P3DProvider("builtin", builtin, binding, selection.required_pairs)],
        )
        weight_options = dict(
            p1d_model=lambda t, z, q: np.full_like(q, 2),
            p1d_parameters=binding,
            theta_p1d=[1],
            k_t_deg=2.4,
            k_p_velocity=0.00035,
        )
        if narrow:
            with pytest.raises(ValueError, match="F.*auxiliary.*builtin"):
                sample_auxiliary(
                    observed_fields[0],
                    bin_geometry,
                    RESPONSE,
                    prepared,
                    [1],
                    **weight_options,
                )
        else:
            aux = sample_auxiliary(
                observed_fields[0],
                bin_geometry,
                RESPONSE,
                prepared,
                [1],
                **weight_options,
            )
            values = evaluate_p3d(prepared, [1], bin_geometry.z_eval, [aux.k], [aux.mu])
            velocity_response_amplitude = velocity_response(
                [aux.k_p_velocity], pixel_width_velocity=0.5, gaussian_sigma_velocity=0
            )[0]
            assert aux.signal == pytest.approx(
                values[0, 0]
                * velocity_response_amplitude
                * velocity_response_amplitude
                * bin_geometry.a_v
                / bin_geometry.d_deg**2,
                rel=5e-15,
            )
            assert values[0, 1] < 0


@pytest.mark.parametrize(
    "auto, response",
    [
        (0, RESPONSE),
        (-1, RESPONSE),
        (
            1,
            InstrumentResponse(0.5, 1e9),
        ),
    ],
)
def test_nonpositive_auxiliary(auto, response):
    """Check nonpositive auxiliary.

    Parameters
    ----------
    auto : int
        Auto-spectrum power or selection, supplied by pytest parametrization.
    response : InstrumentResponse
        Instrumental response input, supplied by pytest parametrization.
    """
    selection = PairSelection([FIELD])
    registry = ParameterRegistry([Parameter("a", 1, "target", step=0.001)])
    binding = BoundParameters(registry, (), {})
    prepared = PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "bad",
                lambda t, z, k, mu, pairs: np.full((len(k), len(pairs)), auto),
                binding,
                [(0, 0)],
            )
        ],
    )
    with pytest.raises(ValueError, match="forest.*auxiliary"):
        sample_auxiliary(
            FIELD,
            geometry(),
            response,
            prepared,
            [1],
            p1d_model=lambda t, z, q: np.ones_like(q),
            p1d_parameters=binding,
            theta_p1d=[1],
            k_t_deg=2.4,
            k_p_velocity=0.00035,
        )


def test_auxiliary_outside_cuts_and_substituted_noise():
    """Check auxiliary outside cuts and substituted noise."""
    from fishhighz.grids import gauss_legendre_grid
    from fishhighz.noise import forest_noise

    grid = gauss_legendre_grid([0.1, 0.2], k_order=2, mu_order=2, h_fid=0.7)
    snapshots = [x.tobytes() for x in [grid.k_flat, grid.mu_flat, grid.q_mode]]
    selection = PairSelection([FIELD])
    registry = ParameterRegistry([Parameter("a", 1, "target", step=0.001)])
    binding = BoundParameters(registry, (), {})
    prepared = PreparedP3D(
        registry,
        selection,
        [
            P3DProvider(
                "flat",
                lambda t, z, k, mu, pairs: np.full((len(k), len(pairs)), 3),
                binding,
                [(0, 0)],
            )
        ],
    )
    bin_geometry = geometry()
    args = dict(
        z_source=3,
        magnitudes=[20, 21],
        quadrature=[1, 1],
        rho=[1, 2],
        variance=[1, 9],
        length_velocity=10,
    )
    results = []
    for p in [2, 20]:
        aux = sample_auxiliary(
            FIELD,
            bin_geometry,
            RESPONSE,
            prepared,
            [1],
            p1d_model=lambda t, z, q: np.full_like(q, p),
            p1d_parameters=binding,
            theta_p1d=[1],
            k_t_deg=2.4,
            k_p_velocity=0.00035,
        )
        assert aux.k < 0.1
        weights = prepare_forest_weights(
            FIELD,
            bin_geometry,
            RESPONSE,
            **args,
            method="legacy",
            iterations=0,
            auxiliary=aux,
        )
        forest_noise_result = forest_noise(
            weights,
            FIELD,
            bin_geometry,
            RESPONSE,
            grid.k_flat,
            grid.mu_flat,
            np.full(4, p),
        )
        results.append((weights, forest_noise_result))
    fixed = prepare_forest_weights(
        FIELD,
        bin_geometry,
        RESPONSE,
        **args,
        method="supplied",
        weights=results[0][0].weights,
    )
    substituted = forest_noise(
        fixed, FIELD, bin_geometry, RESPONSE, grid.k_flat, grid.mu_flat, np.full(4, 20)
    )
    np.testing.assert_allclose(
        substituted.aliasing, 10 * results[0][1].aliasing, rtol=5e-15
    )
    np.testing.assert_array_equal(substituted.pixel, results[0][1].pixel)
    assert not np.allclose(substituted.total, results[1][1].total)
    assert snapshots == [x.tobytes() for x in [grid.k_flat, grid.mu_flat, grid.q_mode]]
