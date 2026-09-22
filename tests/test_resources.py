"""Deterministic checks for the bundled DESI-2 input contract."""

import hashlib
import zipfile
from importlib import resources
from pathlib import Path

import fishhighz.resources as resource_api
from fishhighz.resources import (
    bundled_path,
    bundled_paths,
    bundled_resource,
    resolve_input_path,
)

EXPECTED = {
    "desi2_accuracy.ini": (
        2054,
        "c2053fdfe3b9c89fb62bdf527fb9f053b24ee17416a6350629f80db8c6eb94e6",
    ),
    "README.md": (
        6586,
        "d97bb80fa2574c975e97243ffa2dce2f0dae49eda1ad5b345615e9f8c41fa8bc",
    ),
    "DESI-2-LBG/snr-r19.25-t4000-nexp4.dat": (
        809811,
        "130e262131de9497e57fa0ed7700199e6512ef4e4e27b54fac2cb03fff2cf9f2",
    ),
    "DESI-2-LBG/snr-r19.75-t4000-nexp4.dat": (
        809811,
        "957b92007b050318cd91b58ca969d3162e27124d814520469de7f4f2f691fb64",
    ),
    "DESI-2-LBG/snr-r20.25-t4000-nexp4.dat": (
        809811,
        "943567fca00b2f0e8794c020a711cf0b0a778e05fdeae54aa80a951e21e5a4e0",
    ),
    "DESI-2-LBG/snr-r20.75-t4000-nexp4.dat": (
        809811,
        "2e302b76c186adbaedb6e74aa75de8ff1990967eb23629a12201e22dcc0ee52d",
    ),
    "DESI-2-LBG/snr-r21.25-t4000-nexp4.dat": (
        809811,
        "6e4b2e58526028e6c7755c3fbfe43b94fe34228786f8d89fb4261dc86c021e05",
    ),
    "DESI-2-LBG/snr-r21.75-t4000-nexp4.dat": (
        809811,
        "04117d795a79d4d5572bdafa41775cbe4a6f6d02f93285832d65930772a8e932",
    ),
    "DESI-2-LBG/snr-r22.25-t4000-nexp4.dat": (
        809811,
        "e427fe3bbbc3d3ab812195672dedfd562ab94c696cc91c9cf1af2604b706cab7",
    ),
    "DESI-2-LBG/snr-r22.75-t4000-nexp4.dat": (
        809811,
        "e8557e5a0832f5833a2b4b2052092e5d646d92280be573c74368b34a0cac4ee1",
    ),
    "DESI-2-LBG/snr-r23.25-t4000-nexp4.dat": (
        809811,
        "a3e412023e4b08869cfa0148a00d59e7653c82477d4cb19a3792b932b9569004",
    ),
    "DESI-2-LBG/snr-r23.75-t4000-nexp4.dat": (
        809811,
        "e71a0c927f78f69b1a0a719eb192c28ccd3a33ac5783f2ac92b56c8014396c25",
    ),
    "DESI-2-LBG/snr-r24.25-t4000-nexp4.dat": (
        809811,
        "0073ce17431a8d24d8a662c2167c0f1afbecff691fcc2d1b8cf2d4e284fc620c",
    ),
    "DESI-2-LBG/snr-r24.75-t4000-nexp4.dat": (
        809811,
        "d04fc1974d1cb8531d41dd06b984a628118bc1ef26fb35357ccc64fc44c1efc9",
    ),
    "DESI-2-QSO/snr-r19.25-t4000-nexp4.dat": (
        809829,
        "480e22143980f5e35ccb9999b3ea77934e64c60c934898b45d285af8c4fe744c",
    ),
    "DESI-2-QSO/snr-r19.75-t4000-nexp4.dat": (
        809829,
        "36695a2e74ff92cdf40dc58c37b5b05bc4017847f2a0ad79f98376ae60526794",
    ),
    "DESI-2-QSO/snr-r20.25-t4000-nexp4.dat": (
        809829,
        "d771ac42e1592abbfba6ce9a6145770b5d471d98ae6219e35d0916b0096339e6",
    ),
    "DESI-2-QSO/snr-r20.75-t4000-nexp4.dat": (
        809829,
        "4efb34233109e6f348a4b8c287768677b0363364f064a5db3436d0ef2df8f323",
    ),
    "DESI-2-QSO/snr-r21.25-t4000-nexp4.dat": (
        809829,
        "e39c0ab8412dd430985a93a79441dccb34ca4df976347b9cd5225f6503e30024",
    ),
    "DESI-2-QSO/snr-r21.75-t4000-nexp4.dat": (
        809829,
        "53b118ce9e3d0b096ff5cceb4e2e049e237d20ec742010f47e2e33c28ee2e75c",
    ),
    "DESI-2-QSO/snr-r22.25-t4000-nexp4.dat": (
        809829,
        "73f162ff17706fa9b470c332cc581eaf04ad2be4a8b799889ff8a0c0c75110fb",
    ),
    "DESI-2-QSO/snr-r22.75-t4000-nexp4.dat": (
        809829,
        "4b1ab8a7999a6673c9690b45c42d62ef45c624ed04cca6d23d8f1512dd78d612",
    ),
    "DESI-2-QSO/snr-r23.25-t4000-nexp4.dat": (
        809829,
        "1d806d637a7640b4e2e3b9aae06c5e4525c8101700c2a9c91bf94713939c02a0",
    ),
    "DESI-2-QSO/snr-r23.75-t4000-nexp4.dat": (
        809829,
        "a3a729c0e27f1e26e116e63ef4160fc23139d0b8bd7b1781b50beeb55f368473",
    ),
    "DESI-2-QSO/snr-r24.25-t4000-nexp4.dat": (
        809829,
        "a86acce2ce1ed62bca55f6150205e31f16080ff44b1407844437e81d52d82db4",
    ),
    "DESI-2-QSO/snr-r24.75-t4000-nexp4.dat": (
        809829,
        "779fd1296f53b764450fce82d6e497eb90165436a4806dab209cd315ed50197b",
    ),
    "camb_configs/Planck18.ini": (
        2422,
        "45a04472fb946a2306b0c9668081922e6b0bddd11dfe28f7634aac34d6db9199",
    ),
    "dn_dzdr_qso_desi_2.dat": (
        162000,
        "78914d3a9ad0307ad94072d1b39a54b41d4043ebca5a0e4d1d4d202a3d489c8d",
    ),
    "lae_matched_dndzdr.txt": (
        19125,
        "f56f7d3ba3144ae0fdfaaeba93c557d337e74f8de993110133aaaba3d7f1302e",
    ),
    "lbg_matched_dndzdr.txt": (
        19125,
        "995e1b06e09735e4babdcdfb6f0d253371c9d108b06c5826eb92a9774defcdef",
    ),
    "templates/Planck18_z_2.406.fits": (
        28800,
        "b4a73103e1105b7f9cdb59bbe8133ff0f752b05bc27580dd526f80b9c2fdc26a",
    ),
    "licenses/LYAFORECAST-DATA-README.md": (
        1429,
        "2ddef050a78d69badc0ce31214ceb3087f8a466fb0d9bce16efa62d28bb2a93e",
    ),
    "licenses/LYAFORECAST-LICENSE.txt": (
        35141,
        "589ed823e9a84c56feb95ac58e7cf384626b9cbf4fda2a907bc36e103de1bad2",
    ),
    "licenses/VEGA-LICENSE.txt": (
        35141,
        "589ed823e9a84c56feb95ac58e7cf384626b9cbf4fda2a907bc36e103de1bad2",
    ),
}


def _inventory(resource):
    result = []
    for child in resource.iterdir():
        relative = child.name
        if child.is_dir():
            result.extend(f"{relative}/{name}" for name in _inventory(child))
        else:
            result.append(relative)
    return result


def test_exact_bundled_inventory_and_checksums():
    root = resources.files("fishhighz").joinpath("data")
    assert set(_inventory(root)) == set(EXPECTED)
    for name, (size, digest) in EXPECTED.items():
        content = bundled_resource(name).read_bytes()
        assert len(content) == size
        assert hashlib.sha256(content).hexdigest() == digest


def test_bundled_resource_path_works_in_installed_style_context():
    with bundled_path("camb_configs/Planck18.ini") as path:
        assert path.is_file()
        assert (
            path.read_bytes()
            == bundled_resource("camb_configs/Planck18.ini").read_bytes()
        )


def test_bundled_snr_directory_materializes_each_file_with_live_context():
    with bundled_paths("DESI-2-QSO") as paths:
        assert len(paths) == 12
        assert all(path.is_file() and path.stat().st_size > 0 for path in paths)
        assert [path.name for path in paths] == sorted(path.name for path in paths)


def test_bundled_snr_directory_materializes_zip_traversable_files(
    tmp_path, monkeypatch
):
    names = (
        "DESI-2-QSO/snr-r19.25-t4000-nexp4.dat",
        "DESI-2-QSO/snr-r19.75-t4000-nexp4.dat",
    )
    expected = tuple(bundled_resource(name).read_bytes() for name in names)
    archive_path = tmp_path / "fishhighz-data.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        for name in names:
            archive.writestr(
                f"fishhighz/data/{name}", bundled_resource(name).read_bytes()
            )
    with zipfile.ZipFile(archive_path) as archive:
        monkeypatch.setattr(
            resource_api,
            "_DATA_ROOT",
            zipfile.Path(archive, at="fishhighz/data"),
        )
        with resource_api.bundled_paths("DESI-2-QSO") as paths:
            assert len(paths) == len(names)
            assert tuple(path.read_bytes() for path in paths) == expected
            assert all(path.is_file() and path.stat().st_size > 0 for path in paths)


def test_ini_relative_paths_do_not_depend_on_cwd(tmp_path, monkeypatch):
    ini = tmp_path / "survey" / "run.ini"
    ini.parent.mkdir()
    other = tmp_path / "other"
    other.mkdir()
    monkeypatch.chdir(other)
    expected = (ini.parent / "resources/density.dat").resolve()
    assert resolve_input_path("resources/density.dat", ini) == expected
    absolute = tmp_path / "absolute.dat"
    assert resolve_input_path(absolute, ini) == absolute


def test_distribution_declarations_include_package_data():
    root = Path(__file__).parents[1]
    assert 'fishhighz = ["data/**/*"]' in (root / "pyproject.toml").read_text()
    assert "recursive-include fishhighz/data *" in (root / "MANIFEST.in").read_text()
