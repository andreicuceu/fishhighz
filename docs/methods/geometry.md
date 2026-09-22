# Survey geometry, response, and intrinsic P1D

These composable primitives use supplied backgrounds and synthetic noise; they
are not a survey loader or realistic survey forecast. Run the small three-field
example with `python examples/survey_primitives.py`. The retained external and
built-in examples remain available.

`fishhighz.geometry.prepare_geometry(z_min, z_max, *, z_eval, area_deg2, h_fid,
hubble, transverse_distance, z_order)` prepares one common bin. Supply explicit
`0 <= z_min < z_max`, positive `z_eval` inside the bin, area in square degrees
(up to the full sky), and positive fixed `h_fid`. Ordinary callables receive
immutable 1D redshift arrays and return exactly matching finite positive arrays:
H in km/s/Mpc and transverse comoving D_M in Mpc. D_M is neither angular-diameter
distance nor an already h-scaled distance. Both backgrounds are sampled in
batches at interior quadrature nodes and at `z_eval`, during preparation only.

The immutable `BinGeometry` retains bin/evaluation redshifts, area/solid angle,
`h_fid`, `z_order`, `z_nodes`, dz weights `w_z`, sampled `hubble_nodes` and
`transverse_distance_nodes`, evaluation values, and conversion scalars. It
integrates `volume = Omega*h_fid**3*sum(w_z*c*D_M**2/H)` in `(Mpc/h_fid)^3`.
The explicit Gauss–Legendre order controls volume accuracy; check its convergence
separately from Fourier-grid convergence. No bin-centre approximation, flatness,
background consistency, or interpolation is inferred. `speed_light_kms` records
299792.458. The record contains no cosmology object or background callable.

`prepare_astropy_geometry(cosmology, z_min, z_max, *, z_eval, area_deg2, h_fid,
z_order)` accepts a caller-created Astropy FLRW and delegates to the same
integrator, using unit-converted H and **transverse** comoving distance even for
curved cosmologies. It never substitutes `cosmology.h` for `h_fid` or selects a
default cosmology. Install `pip install 'fishhighz[cosmology]'` for this optional
Astropy/SciPy preparation path. Development acceptance uses
`pip install -e '.[dev,templates,cosmology]'`. Both extras remain optional;
NumPy alone supports supplied backgrounds, response, P1D, and external models.

At `z_eval`, `a_v = H/((1+z_eval)*h_fid)` is in `(km/s)/(Mpc/h_fid)` and
`d_deg = h_fid*D_M*pi/180` is in `(Mpc/h_fid)/degree`. `a_v` is not the scale
factor. Explicit geometry helpers take keyword `a_v`:

- `wavenumber_comoving_to_velocity(k_parallel_comoving, *, a_v)` divides by a_v;
  `wavenumber_velocity_to_comoving` is its inverse.
- `p1d_velocity_to_comoving(power_velocity, *, a_v)` divides km/s power by a_v;
  `p1d_comoving_to_velocity` is its inverse.
- `width_velocity_to_comoving(width_velocity, *, a_v)` divides km/s widths by a_v;
  `width_comoving_to_velocity` is its inverse. Neither changes full-width/sigma
  meaning. These conversion helpers accept finite real arrays or scalars.
- `mode_counts(geometry, grid)` returns positive `volume*grid.q_mode`, shape
  `(node,)`, and rejects unequal geometry/grid h_fid. It does not change the
  grid, cuts, weights, node order, or conjugate-mode normalization.

In `fishhighz.response`, `InstrumentResponse(pixel_width_velocity,
gaussian_sigma_velocity)` requires nonnegative full top-hat pixel width and
Gaussian **one-sigma** in km/s. Use explicit `InstrumentResponse(0, 0)` for
identity, including galaxies. `prepare_response(fields, k, mu, *, a_v, settings)`
requires a mapping from every field ID to its own settings; shared physical
labels imply no sharing. It returns immutable `W(node,field)` in field order.
Inputs are arbitrary paired 1D observed k (h_fid/Mpc) and mu in [0,1], with
slicing and zero k supported. No AP transformation is applied. With q=k*mu/a_v,
`W = np.sinc(q*pixel_width_velocity/(2*pi))*exp(-(q*sigma_velocity)**2/2)`.
Zero q gives exactly one; negative sinc lobes are preserved. Gaussian
attenuation may underflow to zero; unrepresentable arguments fail.
`velocity_response(q, *, pixel_width_velocity, gaussian_sigma_velocity)` returns
the corresponding `(node,)` transfer directly for nonnegative s/km q.

`pair_response(W, selection)` returns `W_i*W_j` in original required-pair order.
Multiply intrinsic P3D by these products, and intrinsic Jacobians by
`products[:, :, None]`, then gather `selection.selected_to_required` for means
and mean derivatives. Supplied noise receives **no automatic response factor**.
The explicit P1D consumer applies W² and the 1/a_v conversion once. Keep geometry,
W, modes, noise, and covariance factors fixed during mean perturbations.

Instrument conversion helpers have deliberately distinct meanings:

| Helper in `fishhighz.response` | Conversion to km/s |
| --- | --- |
| `pixel_width_angstrom_to_velocity(width, *, lambda_obs_angstrom)` | c*full_pixel_width/lambda_obs |
| `gaussian_sigma_angstrom_to_velocity(sigma, *, lambda_obs_angstrom)` | c*sigma_lambda/lambda_obs |
| `gaussian_fwhm_velocity_to_sigma(fwhm_velocity)` | FWHM/(2*sqrt(2*log(2))) |
| `resolving_power_fwhm_to_sigma(resolving_power_fwhm)` | c/(R_FWHM*2*sqrt(2*log(2))) |
| `legacy_resolving_power_to_sigma(resolving_power_legacy)` | c/R_legacy, explicitly interpreted as sigma |

Wavelength inputs are Angstrom and use the local narrow-width approximation.
Supply observed wavelength explicitly; the Lyα convention is
`lambda_obs = 1215.67*(1+z_eval)` Angstrom. FWHM converts once. The labeled legacy
helper uses the new c=299792.458; lyaforecast used 299800, so their values differ
by the known ratio. Direct velocity sigma has no light-speed conversion.

`fishhighz.models.p1d.default_p1d(theta_local, z, k_parallel_velocity)` is a
plain-compatible intrinsic PD2013 callable. Its local vector must be empty;
`BoundParameters(registry, (), {})` connects it to `evaluate_p1d`. The registry
remains nonempty. Direct calls accept finite scalar z>-1; the existing evaluator
retains its nonnegative-redshift contract. k is nonempty finite nonnegative
1D in s/km, including zero. Output is positive `(node,)` km/s power, with no
response, noise, G, BAO, or template normalization. `p1d_floor(z)` exposes the
original redshift-dependent stationary floor. Below it the spectrum is flat;
its first k derivative joins continuously at zero, while the second derivative
need not. Mathematical validity does not establish empirical accuracy outside
the fit's calibration range. Constants and floor are adapted from lyaforecast's
`analytic_p1d_PD2013.py` under its GPLv3 terms; see LICENSE. External P1D uses its
own explicit callable/binding and the same downstream response/conversion.
Changing P3D never selects, integrates, or modifies P1D.

