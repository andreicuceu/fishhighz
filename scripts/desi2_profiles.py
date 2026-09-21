"""Explicit revision-2 profile validation; ordinary default is one-bin quick."""

import argparse
import json

from fishhighz.validation.evidence import check

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["run", "check"])
    p.add_argument("--output", required=True)
    p.add_argument("--suite", choices=["quick", "full"], default="quick")
    p.add_argument(
        "--profiles",
        nargs="+",
        choices=[
            "compatibility",
            "full-compatibility",
            "fixed-compatibility",
            "accuracy",
        ],
        default=["full-compatibility", "fixed-compatibility", "accuracy"],
    )
    p.add_argument("--cases", nargs="+")
    p.add_argument("--bins", nargs="+", type=int)
    p.add_argument("--reference")
    p.add_argument("--template")
    p.add_argument("--reference-bundle")
    p.add_argument("--compatibility-bundle")
    p.add_argument(
        "--accuracy-method",
        choices=["early_lyaforecast", "mcdonald", "inverse_variance", "legacy"],
        default="early_lyaforecast",
    )
    p.add_argument("--wheel")
    p.add_argument("--reuse-completed")
    a = p.parse_args()
    if a.mode == "check":
        print(json.dumps(dict(complete=check(a.output)["complete"])))
    else:
        if not all([a.reference, a.template, a.reference_bundle, a.wheel]):
            p.error("run requires --reference --template --reference-bundle --wheel")
        from fishhighz.validation.profiles import run

        m = run(
            a.output,
            reference=a.reference,
            template=a.template,
            reference_bundle=a.reference_bundle,
            compatibility_bundle=a.compatibility_bundle,
            accuracy_method=a.accuracy_method,
            wheel=a.wheel,
            suite=a.suite,
            cases=a.cases,
            bin_indices=a.bins,
            profiles=a.profiles,
            reuse_completed=a.reuse_completed,
        )
        print(
            json.dumps(
                dict(
                    execution_finished=m["execution_finished"], complete=m["complete"]
                ),
                indent=2,
            )
        )
        if not m["complete"]:
            raise SystemExit(1)
