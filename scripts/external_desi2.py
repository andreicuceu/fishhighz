"""One actual intrinsic-P3D amplitude demonstration with schema-2 evidence."""

import argparse
from pathlib import Path

import numpy as np

from fishhighz.adapters.lyaforecast import IntrinsicP3D
from fishhighz.validation.evidence import check, execute
from fishhighz.validation.profiles import provenance
from fishhighz.validation.real import RealRecipe
from fishhighz.validation.schema import assemble


def run(reference, template, wheel, output):
    identity = provenance(reference, template, wheel)
    recipe = RealRecipe(
        reference,
        template,
        "lya_qso_lbg_lae_15x2pt",
        negative_policy="floor_negative",
        bin_indices=[2],
    )
    oracle = recipe.amplitude(2)
    prepared, _ = recipe.prepare(2, k_intervals=2, mu_order=3, z_order=8)
    sel = recipe.selection
    routes = {
        (i, j): f"{sel.fields[i].id}_{sel.fields[j].id}" for i, j in sel.required_pairs
    }
    bridge = IntrinsicP3D(
        recipe.external,
        routes=routes,
        n_fields=len(sel.fields),
        h_source=recipe.h,
        h_fid=recipe.h,
        k_domain=(0.01, 0.5),
        z_domain=(1.8, 4.5),
    )
    intrinsic = bridge(
        [], prepared.geometry.z_eval, prepared.k, prepared.mu, sel.required_pairs
    )
    signal = prepared.products * intrinsic

    def worker(task):
        settings = dict(
            profile="accuracy",
            parameters=["A"],
            bounds=task["bounds"],
            fields=[f.id for f in sel.fields],
            grid=dict(
                kind="gauss_legendre",
                volume=prepared.geometry.volume,
                h_fid=recipe.h,
                k_intervals=2,
                k_order=4,
                mu_order=3,
            ),
            response_ownership=task["response_ownership"],
            estimator="External amplitude only, separate from the primary BAO accuracy profile; labeled legacy input/response preparation",
        )
        arrays, report = assemble(
            task,
            signal + prepared.noise,
            signal[:, sel.selected_to_required, None],
            settings,
        )
        np.testing.assert_allclose(
            arrays["fisher"][0, 0], oracle["oracle"], rtol=5e-12, atol=0
        )
        report.update(
            provenance=identity, external_amplitude=oracle, analytic_amplitude=True
        )
        return arrays, report

    m = execute(
        output,
        suite="quick",
        kind="external_amplitude",
        worker=worker,
        inputs=[
            wheel,
            template,
            *identity["resources"],
            *identity["reference"]["sources"],
        ],
    )
    check(output)
    return m


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("reference", "template", "wheel", "output"):
        p.add_argument("--" + name, required=True, type=Path)
    a = p.parse_args()
    run(a.reference, a.template, a.wheel, a.output)
    print(
        "Actual external amplitude, trace oracle and schema-2 independent reconstruction passed"
    )
