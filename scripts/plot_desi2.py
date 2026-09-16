"""Generate failure-preserving DESI-2 tables and plots offline."""

import argparse
import json
from pathlib import Path

from fishhighz.validation import plots

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--bundle", required=True)
p.add_argument("--output", required=True)
a = p.parse_args()
table = plots.tables(a.bundle, a.output)
# This command renders only the cases in its checked input inventory.
plots.CASE_IDS = tuple(dict.fromkeys(row["case"] for row in table["rows"]))
plots.plot(table, a.output)

# The inherited grouped panels predate failure-preserving sensitivity plots.
# Retain them separately; canonical diagnostics use plot_desi2_sensitivities.py.
out = Path(a.output)
manifest_path = out / "plots-manifest.json"
manifest = json.loads(manifest_path.read_text())
legacy = out / "superseded-diagnostics"
legacy.mkdir()
for name in list(manifest["files"]):
    if name.endswith(("-diagnostics.png", "-diagnostics.pdf")):
        (out / name).rename(legacy / name)
        del manifest["files"][name]
manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
