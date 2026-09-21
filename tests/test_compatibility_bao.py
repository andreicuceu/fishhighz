"""Availability, analytic information, and paired plotting cuts."""

import numpy as np

from fishhighz.validation.compatibility_bao import (
    dependencies,
    forecast,
    plot_errors,
    plot_ratio,
)


def test_plot_cut_masks_both_components_and_ratio_operands():
    x = np.array([[0.1, 0.21], [0.1, 0.2], [np.nan, 0.1], [0.1, 0.1]])
    y = np.array([[0.1, 0.1], [0.3, 0.1], [0.1, 0.1], [0.2, 0.2]])
    masked = plot_errors(x)
    assert np.isnan(masked[[0, 2]]).all()
    np.testing.assert_array_equal(masked[1], x[1])
    ratio = plot_ratio(x, y)
    assert np.isnan(ratio[:3]).all()
    np.testing.assert_array_equal(ratio[3], [-0.5, -0.5])
    np.testing.assert_array_equal(x[0], [0.1, 0.21])


def test_unavailable_forest_does_not_remove_galaxy_information():
    task = dict(
        fields=[dict(id="f", kind="forest"), dict(id="g", kind="galaxy")],
        required_pairs=[[0, 0], [0, 1], [1, 1]],
        selected_pairs=[[0, 0], [0, 1], [1, 1]],
    )
    assert dependencies(task, [0, 1]) == ["f_f"]
    arrays = dict(
        total=np.array([[10.0, 2.0, 3.0], [10.0, 2.0, 3.0]]),
        modes=np.array([2.0, 2.0]),
        observed_j=np.array([[[1.0, 0.0]] * 3, [[0.0, 1.0]] * 3]),
    )
    singles, joint, _, _, valid = forecast(arrays, task, {}, {})
    assert singles[:2] == [None, None] and joint is None
    np.testing.assert_array_equal(valid, [False, False, True])
    # Galaxy auto variance=2*3²/2=9; independent directions give sigma=3.
    np.testing.assert_array_equal(singles[2]["fisher"], np.eye(2) / 9)
    np.testing.assert_allclose(singles[2]["errors"], [3, 3], rtol=1e-15)


def test_joint_keeps_inter_spectrum_covariance():
    task = dict(
        fields=[dict(id="a", kind="galaxy"), dict(id="b", kind="galaxy")],
        required_pairs=[[0, 0], [0, 1], [1, 1]],
        selected_pairs=[[0, 0], [0, 1], [1, 1]],
    )
    arrays = dict(
        total=np.array([[2.0, 0.5, 3.0], [2.0, 0.5, 3.0]]),
        modes=np.array([2.0, 2.0]),
        observed_j=np.array(
            [[[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]], [[0.0, 1.0], [0.0, 2.0], [0.0, 3.0]]]
        ),
    )
    singles, joint, _, covariance, _ = forecast(arrays, task, {}, {})
    expected = np.array([[4, 1, 0.25], [1, 3.125, 1.5], [0.25, 1.5, 9]])
    np.testing.assert_array_equal(covariance[0], expected)
    value = np.array([1, 2, 3]) @ np.linalg.solve(expected, [1, 2, 3])
    np.testing.assert_allclose(joint["fisher"], np.eye(2) * value, rtol=1e-14)
    assert not np.allclose(joint["fisher"], sum(x["fisher"] for x in singles))
