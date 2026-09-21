"""Bounded refinement controller, reusable offline without model imports."""

from collections import OrderedDict

import numpy as np

from .numerics import change
from .schema import LIMITS, add_metrics

DEFAULT = dict(
    k_intervals=128,
    mu_order=32,
    z_order=32,
    magnitude_order=16,
    iterations=12,
    step=2.5e-4,
)


def study(self, task, *, payload_factory=OrderedDict):
    """Bounded isolated/combined convergence; unresolved controls never pass."""
    method = getattr(self, "weight_method", "legacy")
    if method not in ("legacy", "inverse_variance", "early_lyaforecast", "mcdonald"):
        raise ValueError("unknown accuracy forest-weight method")
    controls = dict(DEFAULT)
    if method != "legacy":
        controls.pop("iterations")
    adaptive = method in ("early_lyaforecast", "mcdonald")
    if adaptive:
        controls["weight_rtol"] = 1e-5
    trials = {}
    outcomes = {}
    failures = []
    # Scope is this single study: task, recipe, options and providers cannot
    # change between calls. Keys contain every refinement control, including
    # the derivative step. Keep at most three payloads; failed calls enter
    # neither cache. No external callable is cached across public calls.
    payloads = payload_factory()
    retained_key = None

    def evaluate(c):
        key = tuple(sorted(c.items()))
        if key in payloads:
            payloads.move_to_end(key)
            return payloads[key]
        from .trials import trial_id

        identity = trial_id(c)
        try:
            a, r = self.evaluate(task, c)
        except (ValueError, FloatingPointError) as error:
            outcomes.setdefault(
                identity,
                dict(id=identity, controls=dict(c), outcome="failed", error=str(error)),
            )
            raise
        outcomes.setdefault(
            identity,
            dict(
                id=identity,
                controls=dict(c),
                outcome="success",
                array_index=len(trials),
            ),
        )
        if adaptive:
            outcomes[identity]["forest_weighting"] = r.get("settings", {}).get(
                "forest_weighting", {}
            )
        payloads[key] = (a, r)
        while len(payloads) > 3:
            # Preserve the current selected control candidate, so later
            # lower-level comparisons cannot evict its final payload.
            if next(iter(payloads)) == retained_key:
                payloads.move_to_end(retained_key)
            payloads.popitem(last=False)
        if key not in trials:
            trials[key] = (
                dict(
                    fisher=a["fisher"],
                    pair_fisher=a["pair_fisher"],
                    volume=r["settings"]["grid"]["volume"],
                ),
                dict(c),
            )
        return a, r

    def summary(c):
        key = tuple(sorted(c.items()))
        if key not in trials:
            evaluate(c)
        return trials[key][0]

    def metric(a, b):
        return change(
            a["fisher"],
            b["fisher"],
            a["pair_fisher"],
            b["pair_fisher"],
            a["volume"],
            b["volume"],
        )

    def passes(m):
        return all(np.isfinite(m[k]) and m[k] <= LIMITS[k] for k in LIMITS)

    # Finite iteration policy belongs only to the historical cumulative method.
    successful = []
    weight_levels = []
    if method == "legacy":
        weight_levels = [3, 6, 12]
        forest = any(
            f.kind == "forest" and i in np.unique(self.selection.selected_pairs)
            for i, f in enumerate(self.selection.fields)
        )
        for count in weight_levels:
            c = {**controls, "iterations": count}
            try:
                successful.append((count, summary(c)))
            except (ValueError, FloatingPointError) as error:
                failures.append(dict(control="weights", level=count, error=str(error)))
        if not successful:
            raise ValueError(f"all prescribed weight counts failed: {failures}")
        controls["iterations"] = successful[-1][0]
        if len(successful) >= 2 and not passes(
            metric(successful[-2][1], successful[-1][1])
        ):
            try:
                c = {**controls, "iterations": 24}
                successful.append((24, summary(c)))
                controls["iterations"] = 24
                weight_levels.append(24)
            except (ValueError, FloatingPointError) as error:
                failures.append(dict(control="weights", level=24, error=str(error)))
        if not forest:
            controls["iterations"] = 12
    retained_key = tuple(sorted(controls.items()))
    levels = dict(
        k=[32, 64, 128],
        mu=[8, 16, 32],
        magnitude=[4, 8, 16],
        volume=[8, 16, 32],
        step=[1e-3, 5e-4, 2.5e-4],
    )
    keys = dict(
        k="k_intervals",
        mu="mu_order",
        magnitude="magnitude_order",
        volume="z_order",
        step="step",
    )
    if adaptive:
        levels["weights"] = [1e-4, 1e-5]
        keys["weights"] = "weight_rtol"
    for name, values in levels.items():
        key = keys[name]
        for extra in range(3):
            c1 = {**controls, key: values[-2]}
            c2 = {**controls, key: values[-1]}
            try:
                m = metric(summary(c1), summary(c2))
            except (ValueError, FloatingPointError) as error:
                failures.append(dict(control=name, level=values[-1], error=str(error)))
                break
            controls[key] = values[-1]
            retained_key = tuple(sorted(controls.items()))
            if passes(m) or extra == 2 or name == "weights":
                break
            values.append(values[-1] * (0.5 if name == "step" else 2))
    # Re-evaluate every isolated level at the same final other controls.
    base_a, base_r = evaluate(controls)
    base = summary(controls)
    names = []
    mf = []
    mp = []
    mv = []
    actual = {}
    for name, values in levels.items():
        key = keys[name]
        rows = []
        for level in values:
            c = {**controls, key: level}
            try:
                rows.append((level, summary(c)))
            except (ValueError, FloatingPointError) as error:
                failures.append(dict(control=name, level=level, error=str(error)))
        if len(rows) < 2:
            raise ValueError(f"not enough valid levels for {name}: {failures}")
        a, b = rows[-2][1], rows[-1][1]
        names.append(name)
        mf.append([a["fisher"], b["fisher"]])
        mp.append([a["pair_fisher"], b["pair_fisher"]])
        mv.append([a["volume"], b["volume"]])
        actual[name] = [x[0] for x in rows]
    weight_rows = []
    if method == "legacy":
        for count, _ in successful:
            try:
                weight_rows.append((count, summary({**controls, "iterations": count})))
            except (ValueError, FloatingPointError) as error:
                failures.append(dict(control="weights", level=count, error=str(error)))
        if len(weight_rows) >= 2:
            a, b = weight_rows[-2][1], weight_rows[-1][1]
        else:
            a = b = base
            failures.append(
                dict(
                    control="weights",
                    error="no valid second count; repeated placeholder is not convergence evidence",
                )
            )
        names.append("weights")
        mf.append([a["fisher"], b["fisher"]])
        mp.append([a["pair_fisher"], b["pair_fisher"]])
        mv.append([a["volume"], b["volume"]])
        actual["weights"] = [r[0] for r in weight_rows]
    lower = dict(controls)
    for name, values in actual.items():
        if name != "weights" or adaptive:
            lower[keys[name]] = values[-2]
    if method == "legacy" and len(weight_rows) >= 2:
        lower["iterations"] = weight_rows[-2][0]
    try:
        combined = summary(lower)
    except (ValueError, FloatingPointError) as error:
        combined = base
        failures.append(dict(control="combined", error=str(error)))
    names.append("combined")
    mf.append([combined["fisher"], base["fisher"]])
    mp.append([combined["pair_fisher"], base["pair_fisher"]])
    mv.append([combined["volume"], base["volume"]])
    add_metrics(base_a, base_r, names, mf, mp, mv)
    base_r.update(
        final_controls=controls,
        actual_levels=actual,
        unresolved_controls=failures,
        study_controls=[v[1] for v in trials.values()],
        combined_lower_controls=lower,
        interpretation="Numerical accuracy within the fixed existing model; retained input policies are not calibrated physics",
    )
    base_a["study_fisher"] = np.stack([v[0]["fisher"] for v in trials.values()])
    base_a["study_pair_fisher"] = np.stack(
        [v[0]["pair_fisher"] for v in trials.values()]
    )
    base_a["study_volume"] = np.array([v[0]["volume"] for v in trials.values()])
    from .trials import bind

    bind(base_a, base_r, list(outcomes.values()), method=method)
    if failures:
        base_r["passed"] = False
    self._prepared.clear()
    return base_a, base_r
