# Bundled DESI Run-2 reference inputs

This directory contains byte-for-byte copies of the DESI-2 inputs used by the
FishHighz input/adaptor work.  The files are data, not regenerated numerical
products.  Their inventory is tested from the installed package; the source
distribution includes this directory through `MANIFEST.in`.

## Origins and provenance

The QSO density table `dn_dzdr_qso_desi_2.dat` and the DESI-2 SNR tables were
copied from the read-only `lyaforecast` checkout at commit
`5abe8bcf8d12cc31d1f5ecd89c87e739c6a14b81` (density table) and
`46f5515283542c0ffb698fbb7e458aefe7bbf6a3` (DESI-2 SNR tables).  The upstream
data notes are preserved verbatim in `licenses/LYAFORECAST-DATA-README.md`.
Those notes identify the density-table DESI-2 SRD origins, but do not provide
generation settings beyond the dataset provenance.  The SNR headers identify the
files as Lyman-alpha forest S/N per Angstrom tables and record the band,
magnitude, exposure time, number of exposures, source-redshift nodes, and
wavelength nodes; no independent copyright or redistribution statement was
supplied.

**LBG and LAE density tables (replaced 2026-10-01).**  The bundled
`lbg_matched_dndzdr.txt` and `lae_matched_dndzdr.txt` were originally copied from
the same `lyaforecast` commit (`5abe8bcf...`).  Those tables had rounded redshift
centres (2.38, 2.60, 2.83, 3.07, 3.29), which do not tile under the 0.22
first-spacing width and are incompatible with piecewise-constant density cells.
They were replaced on 2026-10-01 by the regular-grid DESI-2 SRD v2 tables, copied
byte for byte from the read-only `desi-2-srd` checkout
(`/global/cfs/cdirs/desicollab/users/acuceu/vega_dev/lib/desi-2-srd/data/output_fishhighz_v2/{lbg,lae}_matched_dndzdr.txt`).
The replaced files had SHA-256
`995e1b06e09735e4babdcdfb6f0d253371c9d108b06c5826eb92a9774defcdef` (LBG) and
`f56f7d3ba3144ae0fdfaaeba93c557d337e74f8de993110133aaaba3d7f1302e` (LAE), 19125
bytes each.  The replacements have redshift centres 2.375 + 0.23 k (k = 0..4),
i.e. cell edges 2.26 to 3.41 with zero tiling residual, and the same magnitude
grid and file format.

Provenance of the replacement files: the `desi-2-srd` checkout was at commit
`04f5a716267306f8dbd76787c7497ffadaa283a4`; the last commit touching these two
files is `0470a1f` ("Fix z-bin edges").  The copied files are the **uncommitted
working-tree versions** (modification time 2026-10-01 20:14), which differ from
the committed blobs at that HEAD (committed SHA-256 `d80665f7...` for LBG and
`a5b11d5f...` for LAE).  The copied files are 19125 bytes each with SHA-256
`5c955c62ad2f57c8728bb6f578bb89d6e99ef1db1ad8efb9781f5e61c5a4d420` (LBG) and
`45a691768ec56ed2740385145a80440e4330f419965cb72732f767791201445d` (LAE).  Copied
on 2026-10-01.  Because no commit identifies them, they should be re-pinned to a
`desi-2-srd` commit once the SRD v2 inputs are committed there.  The QSO table
is byte-identical to `desi-2-srd/data/qso_target_selection/dn_dzdr_qso_desi_2.dat`.
The bundled `desi2_accuracy.ini` `[survey] z_edges` was changed in the same
update to the SRD v2 edges `2.0, 2.26, 2.52, 2.73, 2.93, 3.15, 3.41`.  The
historical validation code reads the LBG/LAE tables from a `lyaforecast`
checkout, not from this directory, so saved S2--S4 evidence is unaffected.

`camb_configs/Planck18.ini` is copied from `lyaforecast` commit
`40f31bb101e80c7f0354584a2d0cccea3f5b325f` (the file was introduced in the
earlier `ecd92e9fe548ac1b2f1fb47d3e3038000a20e94b` history).  It is a CAMB
parameter file in CAMB's native conventions; individual parameter units are
those defined by CAMB, not reinterpreted here.

`templates/Planck18_z_2.406.fits` is a byte-for-byte copy from the read-only
`vega` checkout at commit `ecd92e9fe548ac1b2f1fb47d3e3038000a20e94b` (the
commit that added this template).  The FITS binary table is named
`PK` and has columns `K`, `PK`, and `PKSB`; its header records `ZREF=2.406`,
`H0=67.36`, and the fiducial cosmological metadata.  The `K` and power units
follow the Vega template convention (`h/Mpc` and `(Mpc/h)^3` respectively).
No separate FITS-data copyright or redistribution statement was supplied by
the source checkout.

The source repository checkout used for this copy was
`01f1a9784767e43e68c89ee4464c2bdcb15420c0`; the file-history commit and the
immutable SHA-256 are both recorded here so that provenance can be checked
without relying on the current checkout state.

Both source checkouts carry `GPL-3.0-or-later` repository licenses.  Exact
copies are included as `licenses/LYAFORECAST-LICENSE.txt` and
`licenses/VEGA-LICENSE.txt`.  A repository-level GPL notice does not by itself
establish separate rights for every scientific data product; maintainers must
confirm redistribution of the SNR and FITS assets before publishing a release.
This is the only unresolved packaging/legal gate; no provenance or licensing
fact has been inferred.

## Units and semantics

* `dn_dzdr_qso_desi_2.dat`, `lbg_matched_dndzdr.txt`, and
  `lae_matched_dndzdr.txt` are three-column `(z, r-band magnitude, count)`
  tables.  The count is per redshift/magnitude cell per square degree, as
  stated by the upstream notes; the density adapter applies the explicit cell
  widths and any caller-requested normalization.  The LBG/LAE tables have
  contiguous redshift cells of width 0.23 (centres 2.375 to 3.295); piecewise-
  constant interpolation rejects tables whose cells do not tile the axis.
* `DESI-2-QSO/*.dat` and `DESI-2-LBG/*.dat` are header-driven S/N tables.  The
  row axis is wavelength in Angstrom, the column axis is source redshift, and
  the values are S/N per Angstrom for the stated r-band magnitude, exposure
  time, and number of exposures.  Headers are retained unchanged.
* `camb_configs/Planck18.ini` is a CAMB configuration in CAMB-native units and
  syntax.
* `templates/Planck18_z_2.406.fits` is a Vega-format `PK` binary table.  `PK`
  is the full power and `PKSB` the smooth/no-wiggle component; the signed
  decomposition is retained exactly.

## Exact inventory

The sizes and SHA-256 values below refer to the files in this directory.

| Relative path | Bytes | SHA-256 |
| --- | ---: | --- |
| `dn_dzdr_qso_desi_2.dat` | 162000 | `78914d3a9ad0307ad94072d1b39a54b41d4043ebca5a0e4d1d4d202a3d489c8d` |
| `lbg_matched_dndzdr.txt` | 19125 | `5c955c62ad2f57c8728bb6f578bb89d6e99ef1db1ad8efb9781f5e61c5a4d420` |
| `lae_matched_dndzdr.txt` | 19125 | `45a691768ec56ed2740385145a80440e4330f419965cb72732f767791201445d` |
| `camb_configs/Planck18.ini` | 2422 | `45a04472fb946a2306b0c9668081922e6b0bddd11dfe28f7634aac34d6db9199` |
| `templates/Planck18_z_2.406.fits` | 28800 | `b4a73103e1105b7f9cdb59bbe8133ff0f752b05bc27580dd526f80b9c2fdc26a` |

All QSO files are 809829 bytes and all LBG files are 809811 bytes.  Their
individual checksums are listed here to make accidental population or
magnitude substitution detectable:

| Population | Magnitude | SHA-256 |
| --- | ---: | --- |
| LBG | 19.25 | `130e262131de9497e57fa0ed7700199e6512ef4e4e27b54fac2cb03fff2cf9f2` |
| LBG | 19.75 | `957b92007b050318cd91b58ca969d3162e27124d814520469de7f4f2f691fb64` |
| LBG | 20.25 | `943567fca00b2f0e8794c020a711cf0b0a778e05fdeae54aa80a951e21e5a4e0` |
| LBG | 20.75 | `2e302b76c186adbaedb6e74aa75de8ff1990967eb23629a12201e22dcc0ee52d` |
| LBG | 21.25 | `6e4b2e58526028e6c7755c3fbfe43b94fe34228786f8d89fb4261dc86c021e05` |
| LBG | 21.75 | `04117d795a79d4d5572bdafa41775cbe4a6f6d02f93285832d65930772a8e932` |
| LBG | 22.25 | `e427fe3bbbc3d3ab812195672dedfd562ab94c696cc91c9cf1af2604b706cab7` |
| LBG | 22.75 | `e8557e5a0832f5833a2b4b2052092e5d646d92280be573c74368b34a0cac4ee1` |
| LBG | 23.25 | `a3e412023e4b08869cfa0148a00d59e7653c82477d4cb19a3792b932b9569004` |
| LBG | 23.75 | `e71a0c927f78f69b1a0a719eb192c28ccd3a33ac5783f2ac92b56c8014396c25` |
| LBG | 24.25 | `0073ce17431a8d24d8a662c2167c0f1afbecff691fcc2d1b8cf2d4e284fc620c` |
| LBG | 24.75 | `d04fc1974d1cb8531d41dd06b984a628118bc1ef26fb35357ccc64fc44c1efc9` |
| QSO | 19.25 | `480e22143980f5e35ccb9999b3ea77934e64c60c934898b45d285af8c4fe744c` |
| QSO | 19.75 | `36695a2e74ff92cdf40dc58c37b5b05bc4017847f2a0ad79f98376ae60526794` |
| QSO | 20.25 | `d771ac42e1592abbfba6ce9a6145770b5d471d98ae6219e35d0916b0096339e6` |
| QSO | 20.75 | `4efb34233109e6f348a4b8c287768677b0363364f064a5db3436d0ef2df8f323` |
| QSO | 21.25 | `e39c0ab8412dd430985a93a79441dccb34ca4df976347b9cd5225f6503e30024` |
| QSO | 21.75 | `53b118ce9e3d0b096ff5cceb4e2e049e237d20ec742010f47e2e33c28ee2e75c` |
| QSO | 22.25 | `73f162ff17706fa9b470c332cc581eaf04ad2be4a8b799889ff8a0c0c75110fb` |
| QSO | 22.75 | `4b1ab8a7999a6673c9690b45c42d62ef45c624ed04cca6d23d8f1512dd78d612` |
| QSO | 23.25 | `1d806d637a7640b4e2e3b9aae06c5e4525c8101700c2a9c91bf94713939c02a0` |
| QSO | 23.75 | `a3a729c0e27f1e26e116e63ef4160fc23139d0b8bd7b1781b50beeb55f368473` |
| QSO | 24.25 | `a86acce2ce1ed62bca55f6150205e31f16080ff44b1407844437e81d52d82db4` |
| QSO | 24.75 | `779fd1296f53b764450fce82d6e497eb90165436a4806dab209cd315ed50197b` |
