"""Transfer-aware immutable CAMB snapshots in fixed fiducial h units.

This optional provider preserves a caller-supplied fixed-species CAMB model.
It intentionally supports only five flat-LambdaCDM coordinates; neutrino mass,
radiation and dark-energy extensions need a different validated model.
"""

import hashlib
import inspect
import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

BASIS = ("omega_b", "omega_c", "h", "n_s", "ln_1e10_A_s")


def evaluate_linear_power(interpolator, redshift, k, *, low_k_policy="raise"):
    """Evaluate inside transfer support, or explicitly continue its low-k slope.

    CAMB's multi-redshift spline can silently clamp outside its k support.
    ``matched_power_law`` instead matches the value and the logarithmic slope
    measured over 0.1% above the actual transfer minimum. It is a numerical
    soft-integral continuation, not extra computed transfer modes. The default
    rejects unsupported nodes. No high-k extrapolation is permitted.
    """
    if low_k_policy not in ("raise", "matched_power_law"):
        raise ValueError("unsupported low-k continuation policy")
    k = np.asarray(k, dtype=np.float64)
    lower, upper = float(interpolator.kmin), float(interpolator.kmax)
    if np.any(~np.isfinite(k)) or np.any(k <= 0) or np.any(k > upper):
        raise ValueError("positive k within the upper transfer support required")
    below = k < lower
    if np.any(below) and low_k_policy == "raise":
        raise ValueError("k below CAMB transfer support; choose an explicit policy")
    values = np.asarray(interpolator.P(float(redshift), np.maximum(k, lower)))
    if np.any(below):
        anchor = float(interpolator.P(float(redshift), lower))
        sample = float(interpolator.P(float(redshift), lower * 1.001))
        slope = np.log(sample / anchor) / np.log(1.001)
        if not np.isfinite(slope) or not 0 < slope < 2:
            raise ValueError("unexpected low-k transfer slope")
        values = np.where(below, anchor * (k / lower) ** slope, values)
    return values


def _immutable(values):
    result = np.array(values, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class CosmologySnapshot:
    """Linear/smooth power (nz,nk), growth and geometry, caller redshift order.

    k_reference: h_fid/Mpc; powers: (Mpc/h_fid)^3; H: km/s/Mpc;
    DM and r_drag: physical Mpc. f and sigma8 are CAMB's fixed-species
    total-matter scalar-growth approximation, not exact f(k,z).
    """

    key: str
    parameters: tuple
    redshifts: np.ndarray
    k_reference: np.ndarray
    pk_linear: np.ndarray
    pk_smooth: np.ndarray
    growth_rate: np.ndarray
    sigma8: np.ndarray
    H: np.ndarray
    DM: np.ndarray
    r_drag_mpc: float
    h_fid: float
    metadata_json: str


class CambSnapshotProvider:
    """Cache one actual transfer solve per cosmology for all requested bins.

    base_parameters is an explicit CAMBparams instance. smooth_spectrum is a
    pure (k_reference,pk_linear)->pk_smooth callable. No optional CAMB import
    occurs until construction. Covariance, survey preparation and responses
    are deliberately outside this provider and are never differentiated.
    """

    def __init__(
        self,
        base_parameters,
        k_reference,
        redshifts,
        h_fid,
        smooth_spectrum,
        *,
        cache_dir=None,
        accuracy_boost=1.0,
        transfer_variable="delta_tot",
        low_k_policy="raise",
    ):
        import camb

        self._camb = camb
        self._base = base_parameters.copy()
        if (
            self._base.omk != 0
            or self._base.DarkEnergy.w != -1
            or self._base.DarkEnergy.wa != 0
        ):
            raise ValueError(
                "this provider supports flat LambdaCDM with fixed species only"
            )
        k = _immutable(k_reference)
        z = _immutable(redshifts)
        if (
            k.ndim != 1
            or len(k) < 4
            or np.any(~np.isfinite(k))
            or np.any(k <= 0)
            or np.any(np.diff(k) <= 0)
        ):
            raise ValueError("require a finite increasing positive k grid")
        if (
            z.ndim != 1
            or not len(z)
            or np.any(~np.isfinite(z))
            or np.any(z < 0)
            or len(np.unique(z)) != len(z)
        ):
            raise ValueError("require explicit unique nonnegative redshifts")
        if (
            not np.isfinite(h_fid)
            or h_fid <= 0
            or not np.isfinite(accuracy_boost)
            or accuracy_boost <= 0
        ):
            raise ValueError("reference h and accuracy must be positive")
        if transfer_variable != "delta_tot":
            raise ValueError(
                "only the declared total-matter approximation is supported"
            )
        if low_k_policy not in ("raise", "matched_power_law"):
            raise ValueError("unsupported low-k continuation policy")
        self.k_reference, self.redshifts, self.h_fid = k, z, float(h_fid)
        self.smooth_spectrum = smooth_spectrum
        self.accuracy_boost = float(accuracy_boost)
        self.transfer_variable = transfer_variable
        self.low_k_policy = low_k_policy
        self.cache_dir = None if cache_dir is None else Path(cache_dir)
        if self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.defaults = dict(
            zip(
                BASIS,
                (
                    self._base.ombh2,
                    self._base.omch2,
                    self._base.H0 / 100.0,
                    self._base.InitPower.ns,
                    np.log(self._base.InitPower.As * 1e10),
                ),
            )
        )
        smooth_source = inspect.getsource(smooth_spectrum)
        self._signature = {
            "provider_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "smoothing_sha256": hashlib.sha256(smooth_source.encode()).hexdigest(),
            "base_parameters": str(self._base),
            "camb_version": camb.__version__,
            "accuracy_boost": self.accuracy_boost,
            "transfer_variable": transfer_variable,
            "low_k_policy": low_k_policy,
            "h_fid": self.h_fid,
            "k_reference": k.tolist(),
            "redshifts": z.tolist(),
            "growth": "get_fsigma8/get_sigma8, scalar total-matter approximation, fixed neutrino/radiation content",
        }
        self._memory = {}
        self.solve_count = self.memory_hits = self.disk_hits = 0
        self.last_results = None

    def _resolve(self, values):
        p = dict(self.defaults)
        if values is not None:
            if not set(values) <= set(BASIS):
                raise ValueError("unsupported cosmology coordinate")
            p.update(values)
        if any(not np.isfinite(v) for v in p.values()) or any(
            p[n] <= 0 for n in ("omega_b", "omega_c", "h", "n_s")
        ):
            raise ValueError("invalid cosmological parameter values")
        return p

    def evaluate(self, values=None):
        """Return a snapshot; nuisance-only calls reuse the same immutable state."""
        p = self._resolve(values)
        meta = {**self._signature, "parameters": p}
        metadata = json.dumps(meta, sort_keys=True)
        key = hashlib.sha256(metadata.encode()).hexdigest()
        if key in self._memory:
            self.memory_hits += 1
            return self._memory[key]
        path = None if self.cache_dir is None else self.cache_dir / (key + ".npz")
        if path is not None and path.exists():
            with np.load(path, allow_pickle=False) as a:
                if str(a["metadata_json"].item()) != metadata:
                    raise ValueError("cosmology cache provenance mismatch")
                arrays = {
                    n: _immutable(a[n])
                    for n in (
                        "redshifts",
                        "k_reference",
                        "pk_linear",
                        "pk_smooth",
                        "growth_rate",
                        "sigma8",
                        "H",
                        "DM",
                    )
                }
                rdrag = float(a["r_drag_mpc"])
            self.disk_hits += 1
        else:
            camb = self._camb
            pars = self._base.copy()
            # Preserve every fixed neutrino, radiation, helium and DE field.
            pars.ombh2, pars.omch2, pars.H0 = p["omega_b"], p["omega_c"], 100 * p["h"]
            pars.InitPower.ns = p["n_s"]
            pars.InitPower.As = float(np.exp(p["ln_1e10_A_s"]) * 1e-10)
            pars.WantCls = pars.Want_CMB = pars.Want_CMB_lensing = False
            pars.NonLinear = camb.model.NonLinear_none
            pars.set_accuracy(AccuracyBoost=self.accuracy_boost)
            pars.set_matter_power(
                redshifts=sorted(self.redshifts, reverse=True),
                kmax=1.05 * self.k_reference[-1] * self.h_fid,
            )
            results = camb.get_results(pars)
            self.last_results = results
            self.solve_count += 1
            interpolator = results.get_matter_power_interpolator(
                nonlinear=False,
                var1=self.transfer_variable,
                var2=self.transfer_variable,
                hubble_units=False,
                k_hunit=False,
            )
            full = np.array(
                [
                    evaluate_linear_power(
                        interpolator,
                        z,
                        self.k_reference * self.h_fid,
                        low_k_policy=self.low_k_policy,
                    )
                    * self.h_fid**3
                    for z in self.redshifts
                ]
            )
            smooth = np.array(
                [self.smooth_spectrum(self.k_reference, power) for power in full]
            )
            actual_z = np.array(results.Params.Transfer.PK_redshifts)
            indices = []
            for z in self.redshifts:
                where = np.flatnonzero(np.isclose(actual_z, z, atol=1e-12, rtol=0))
                if len(where) != 1:
                    raise ValueError("CAMB redshift ordering mismatch")
                indices.append(int(where[0]))
            sigma = np.array(results.get_sigma8())[indices]
            growth = np.array(results.get_fsigma8())[indices] / sigma
            arrays = {
                "redshifts": self.redshifts,
                "k_reference": self.k_reference,
                "pk_linear": full,
                "pk_smooth": smooth,
                "sigma8": sigma,
                "growth_rate": growth,
                "H": results.hubble_parameter(self.redshifts),
                "DM": (1 + self.redshifts)
                * results.angular_diameter_distance(self.redshifts),
            }
            rdrag = float(results.get_derived_params()["rdrag"])
            arrays = {n: _immutable(v) for n, v in arrays.items()}
            if path is not None:
                temporary = path.with_name(path.stem + f".{os.getpid()}.partial.npz")
                np.savez_compressed(
                    temporary, **arrays, r_drag_mpc=rdrag, metadata_json=metadata
                )
                temporary.replace(path)
        if arrays["pk_linear"].shape != (len(self.redshifts), len(self.k_reference)):
            raise ValueError("snapshot power shape mismatch")
        if any(
            np.any(~np.isfinite(v)) or np.any(v <= 0)
            for n, v in arrays.items()
            if n not in ("redshifts", "k_reference", "DM")
        ):
            raise ValueError("nonpositive or nonfinite snapshot data")
        if np.any(~np.isfinite(arrays["DM"])) or np.any(arrays["DM"] < 0):
            raise ValueError("invalid comoving angular diameter distance")
        result = CosmologySnapshot(
            key,
            tuple(p.items()),
            **arrays,
            r_drag_mpc=rdrag,
            h_fid=self.h_fid,
            metadata_json=metadata,
        )
        self._memory[key] = result
        return result
