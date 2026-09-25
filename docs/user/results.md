# Results and saved output

`Forecast.run()` returns a `SurveyResult`. Keep the distinction between
individual-spectrum and joint constraints when interpreting forecast precision:

| Attribute | Content |
| --- | --- |
| `individual` / `individual_spectra` | One constraint per selected spectrum, calculated with its own covariance. |
| `joint` / `joint_constraints` | One constraint per bin from the full selected-spectrum covariance. |
| `excluded` | Explicit records for pairs omitted by the bin selection. |
| `constraints` | Individual records followed by joint records; excludes omitted pairs. |
| `fisher_results` | Fisher objects in `constraints` order. |
| `combined` | Combined joint Fisher result across independent bins, in the full named parameter basis. |
| `convergence` | Recorded forest-weight iteration diagnostics by bin and field. |
| `resolved_settings` | Expanded configuration, input identities, parameter order and execution settings. |
| `config`, `prepared` | Parsed configuration and retained prepared state. |

Each `SpectrumConstraint` gives `bin_index` (zero-based), `bin_id`, redshift
`bounds`, `z_eval`, `kind`, `pair`, `parameter_ids`, `status`, `fisher` (also
`result`), `sigma_ap`, `sigma_at` and `correlation`. Joint records have
`pair=None`. INI bin sections are one-based; parameter IDs and result indices
are zero-based (`ap_0`, `at_0` for the first bin).

In BAO mode, `status="available"` means the two-parameter Fisher information has full rank;
`available` is the corresponding Boolean property. Rank-deficient constraints
have `status="unavailable"`, retain their Fisher object, and report errors and
correlation as `None`. Excluded records have `status="excluded"` and no Fisher
object or uncertainties. Do not interpret either status as a zero uncertainty.

The two BAO errors are marginalized within each record's AP pair. Individual
constraints are not blocks of the inverse joint covariance. Correlations between
spectra remain in the joint calculation even though the native sampling-noise
model is diagonal. Across bins, `combined` retains separate AP parameters; it
does not impose a shared dilation parameter.

```python
for item in result.joint:
    print(item.bin_index, item.status, item.sigma_ap, item.sigma_at)

output = result.save("accuracy-desi2")
```

`save()` creates a new directory, including missing parents, and returns its
resolved path. Existing destinations raise an error. Two files are written:

- `settings.json`: result schema `fishhighz-forecast-result`, version 1;
  expanded settings and input identities/checksums captured during construction
  and preparation; convergence records; constraint identities, statuses,
  uncertainties and counts. Non-finite JSON numbers are rejected.
- `results.npz`: `individual_fisher`, `joint_fisher` and `combined_fisher` arrays
  of data Fisher information. Individual and joint arrays follow their record
  order; excluded records have no matrix. Pair matrices use their recorded
  `parameter_ids`; the combined matrix uses saved `parameter_order`.

```python
import json
import numpy as np

with open("accuracy-desi2/settings.json") as stream:
    metadata = json.load(stream)
with np.load("accuracy-desi2/results.npz", allow_pickle=False) as arrays:
    combined = arrays["combined_fisher"]
```

This is an inspectable saved result, not a serialization of the full live
`SurveyResult`: no public inverse `load()` reconstructs prepared objects. Input
hashes identify the inputs used; the output directory does not contain copies
of every source table. Retain external inputs separately for reproducibility.

For direct Python analyses, `FisherResult` keeps data and prior information
separate and exposes named fixing, marginalization and rank diagnostics. Bounds
and derivative steps are not priors. A converged forest iteration does not
establish quadrature convergence or physical accuracy; consult the
[research prescription](../research/RESEARCH_BASELINE.md) for the qualified
calculation and remaining limitations.

## BAO constraints with marginalized nuisances

With `[model] mode=bao_marginalized`, each joint or individual record contains
`target_ids=(alpha_iso_N,phi_N)` in that order, where `N` is the original
zero-based bin index. `target_fiducials=(1,1)`, and `target_covariance`,
`target_errors` and `target_correlations` describe the two targets after
marginalizing the active tracer biases and forest beta. The galaxy growth rate
is fixed to the background's fiducial `f(z)` in each bin and is recorded in
`resolved_settings.fixed_growth_by_bin`.

Saves use schema `fishhighz-bao-marginalized-result`, version 1, with
`settings.json` and `results.npz` following the named-array layout below.
Each record includes its active `parameter_ids`, Fisher matrix, registry
fiducials, rank and null directions. A singular active Fisher matrix yields
status `unavailable` with no target covariance; an empty selected bin has
status `excluded`. Individual spectra use their own covariance closure.

## Full-shape constraints and output

With `[model] mode=full_shape`, each record retains the complete active internal
Fisher matrix in `fisher`, with `parameter_ids` specifying its order. Its registry
contains fiducials, target/nuisance roles and derivative steps. Only parameters
from other bins and nuisances of absent tracers are removed before rank testing.
The default retains all five targets, including a zero-sensitivity growth
target where the selected spectra cannot constrain it. The forest-only
`target_set=dilation_only` option retains the four dilation targets and fixes
growth to its fiducial value.

| Attribute | Full-shape content |
| --- | --- |
| `target_ids` | Ordered `alpha_w`, `phi_w`, `alpha_s`, `phi_s` for `alpha_phi`, or `alpha_iso_w`, `phi_w`, `alpha_iso_s`, `phi_s` for `alpha_iso_phi`; the five-target fit adds `fsigma8`. Each has the original zero-based bin suffix. |
| `target_fiducials` | Fiducials in the reported target basis. |
| `target_errors` | Marginalized standard deviations in `target_ids` order. |
| `target_covariance`, `target_correlations` | Marginalized target covariance and correlation matrices. |
| `sigma8_fid` | Fixed fiducial sigma8 used to report growth. |

Active nuisances are marginalized before extracting the target covariance.
Internally growth is `f` with fixed template normalization. Reported covariance
is `J @ C_internal @ J.T`, with `J=diag(1,1,1,1,sigma8_fid)`; the transformation
includes growth cross-covariances when growth is a target. The four-target
forest-only fit has the identity reporting Jacobian. The BAO-specific `sigma_ap`, `sigma_at` and
`correlation` fields are `None` for full-shape records.

`resolved_settings.grid_by_bin` stores actual radial and angular quadrature
nodes and weights for each nonempty bin. Its category intervals give the
observed-cut bounds, radial-node counts and selected spectra used in each
joint covariance block. These arrays make category boundaries and refinement
grids inspectable from the saved result.

If the active Fisher is singular, the entire constraint is `unavailable` and
all target uncertainty arrays are `None`. The Fisher matrix and its rank,
threshold, eigenvalues, scales and null directions remain available; no
regularization or pseudoinverse is applied. Individual-spectrum constraints use
their own covariance and marginalize their own active nuisances. `combined`
retains the independent-bin internal Fisher basis with explicit bin identities.

```python
for item in result.joint:
    print(item.bin_index, item.status, item.target_ids, item.target_errors)
```

Full-shape saves use schema `fishhighz-full-shape-result`, version 1.
`settings.json` contains active parameter metadata, reported target identities
and fiducials, statuses, rank diagnostics, requested Fourier cuts and effective
quadrature. Each record's `arrays` mapping names its numeric arrays in
`results.npz`, such as `record_0000_fisher` and
`record_0000_target_covariance`. Different active matrix sizes have separate keys,
so `numpy.load(..., allow_pickle=False)` remains sufficient. Excluded records
have no matrices; unavailable records omit target uncertainty arrays.
`combined_parameters` gives the exact order of `combined_fisher`.
`resolved_settings` records effective category cutoffs, selected pairs by
original bin index, the target basis, and whether individual calculations ran.
An empty bin has one joint record with status `excluded`, no Fisher matrix,
and separate excluded pair records. Its parameter IDs are empty. The
`individuals=False` run keeps these exclusion records.
The existing BAO output schema and accessors are unchanged.
