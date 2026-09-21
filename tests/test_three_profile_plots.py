"""Synthetic selection gaps, exact joint extraction and paired plot cuts."""

import numpy as np

from fishhighz.validation.three_profile_plots import PROFILES, REVISION, build_tables


def synthetic_records():
    records = []
    for profile in PROFILES:
        for index in range(6):
            pairs = (
                [(i, j) for i in range(5) for j in range(i, 5)]
                if index
                else [(0, 0), (0, 1), (1, 1)]
            )
            count = len(pairs)
            arrays = dict(
                fisher=np.eye(2) * 100,
                covariance=np.eye(2) / 100,
                errors=np.array([0.1, 0.1]),
                correlation=np.eye(2),
                rank=np.array([2]),
                pair_fisher=np.tile(np.eye(2) * 25, (count, 1, 1)),
                pair_covariance=np.tile(np.eye(2) / 25, (count, 1, 1)),
                pair_errors=np.full((count, 2), 0.2),
                pair_correlation=np.tile(np.eye(2), (count, 1, 1)),
                pair_rank=np.full(count, 2),
            )
            task = dict(
                profile=profile,
                bin=index,
                recipe_revision=REVISION,
                fields=[dict(id=f"field{i}") for i in range(5)],
                selected_pairs=pairs,
                bounds=[index + 2, index + 3],
            )
            records.append((dict(task=task), arrays))
    return records


def test_selection_and_joint_preserved():
    table = build_tables(synthetic_records())
    for profile in PROFILES:
        rows = [r for r in table["rows"] if r["profile"] == profile]
        assert (
            sum(
                r["availability"] == "available" and r["spectrum"] != "joint"
                for r in rows
            )
            == 78
        )
        assert sum(r["availability"] == "excluded_by_selection" for r in rows) == 12
        joint = [r for r in rows if r["spectrum"] == "joint"]
        assert len(joint) == 6
        assert all(r["fisher"] == [[100.0, 0.0], [0.0, 100.0]] for r in joint)
        assert all(r["numerical_status"] == "unqualified" for r in rows)


def test_hidden_values_and_both_ratio_operands():
    records = synthetic_records()
    records[0][1]["pair_errors"][0] = [0.3, 0.1]
    table = build_tables(records)
    hidden = table["rows"][0]
    assert hidden["sigma_parallel"] == 0.3
    assert hidden["sigma_transverse"] == 0.1
    assert hidden["plot_omission"] == "sigma_above_0.2"
    ratios = [
        r
        for r in table["ratios"]
        if r["bin"] == 0 and r["spectrum"] == hidden["spectrum"]
    ]
    assert ratios[0]["plot_omission"] == ratios[1]["plot_omission"] == "sigma_above_0.2"
    assert ratios[2]["plot_omission"] is None
    np.testing.assert_allclose(ratios[0]["sigma_parallel_fraction"], 0.2 / 0.3 - 1)


def test_failure_distinct_from_exclusion():
    records = synthetic_records()
    records[0] = (records[0][0], None)
    rows = build_tables(records)["rows"][:16]
    assert sum(r["availability"] == "failed_or_unconstrained" for r in rows) == 4
    assert sum(r["availability"] == "excluded_by_selection" for r in rows) == 12


def test_render_saved_weights_and_coefficients(tmp_path):
    from fishhighz.validation.three_profile_plots import render

    records = synthetic_records()
    records[6][0]["report"] = dict(
        settings=dict(
            magnitude_nodes=[20.0, 21.0],
            forest_weighting=dict(
                forests=dict(
                    field0=dict(
                        weights=[1.0, 2.0],
                        coefficients=[1.0, 2.0, 3.0, 4.0, 5.0],
                        status="converged",
                    )
                )
            ),
        )
    )
    table = build_tables(records)
    render(table, tmp_path)
    assert len(list(tmp_path.glob("*.png"))) == 10
    assert (tmp_path / "three-profile-tables.json").is_file()
    assert table["weighting"][6]["fields"][0]["A"] == 4.0
    assert table["weighting"][6]["fields"][0]["P_pixel"] == 5.0
