import re
import shutil
import subprocess
import sys

import pytest

pytest.importorskip("ase")


@pytest.fixture(scope="session")
def tmp_dir_fixture(tmp_path_factory):
    """Fixture that provides a temporary directory for tests."""
    tmp_path = tmp_path_factory.mktemp("cli_tests")
    return tmp_path


@pytest.fixture(scope="session")
def caffeine_xyz_file_fixture(tmp_dir_fixture, caffeine_molecule_xyz):
    """Fixture that creates a caffeine.xyz file in the temporary directory."""
    file_path = tmp_dir_fixture / "caffeine.xyz"
    with open(file_path, "w") as f:
        f.write(caffeine_molecule_xyz)
    return file_path


@pytest.fixture(scope="session")
def caffeine_xyz_rotated_file_fixture(tmp_dir_fixture, caffeine_molecule_xyz_rotated):
    """Fixture that creates a caffeine_rotated.xyz file in the temporary
    directory."""
    file_path = tmp_dir_fixture / "caffeine_rotated.xyz"
    with open(file_path, "w") as f:
        f.write(caffeine_molecule_xyz_rotated)
    return file_path


@pytest.fixture(scope="session")
def caffeine_xyz_obabel_file_fixture(tmp_dir_fixture, caffeine_molecule_xyz_obabel):
    """Fixture that creates a caffeine_obabel.xyz file in the temporary
    directory."""
    file_path = tmp_dir_fixture / "caffeine_obabel.xyz"
    with open(file_path, "w") as f:
        f.write(caffeine_molecule_xyz_obabel)
    return file_path


# FIXME: better share data between cli and api tests
# TODO: test file creation this would also test rotated and shifted molecule
@pytest.mark.parametrize(
    "file1_fixture,file2_fixture,additional_cli_args,expected_rmsd",
    [
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_file_fixture",
            [],
            0.0,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_file_fixture",
            ["--heavy"],
            0.0,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_rotated_file_fixture",
            [],
            0.0,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_rotated_file_fixture",
            ["--heavy"],
            0.0,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_obabel_file_fixture",
            [],
            0.1393943463,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_obabel_file_fixture",
            ["--heavy"],
            0.17258155,
        ),
    ],
)
def test_cli_rmsd(
    file1_fixture,
    file2_fixture,
    additional_cli_args,
    expected_rmsd,
    request,
):
    file1 = request.getfixturevalue(file1_fixture)
    file2 = request.getfixturevalue(file2_fixture)

    exe = shutil.which("irmsd")
    assert exe, "CLI not found on PATH (editable install?)"
    out = subprocess.check_output(
        [
            exe,
            "--rmsd",
            *additional_cli_args,
            str(file1),
            str(file2),
        ],
        text=True,
    )
    match = re.search(r"Cartesian RMSD:\s*([-+]?\d*\.\d+|\d+) Å", out)
    assert match, f"Cartesian RMSD not found in output: {out}"

    rmsd_value = float(match.group(1))
    assert pytest.approx(rmsd_value, abs=1e-6) == expected_rmsd


# def test_cli_runs():
#    exe = shutil.which("irmsd")
#    assert exe, "CLI not found on PATH (editable install?)"
#    out = subprocess.check_output([exe, "--a","0.5","--x","1","2","--y","10","20"], text=True)
#    assert out.strip() == "10.5 21.0"
