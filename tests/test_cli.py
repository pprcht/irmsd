import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from irmsd.cli import main

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
def test_cli_compare_quarternion(
    file1_fixture,
    file2_fixture,
    additional_cli_args,
    expected_rmsd,
    request,
    capsys,
):
    file1 = request.getfixturevalue(file1_fixture)
    file2 = request.getfixturevalue(file2_fixture)

    mains_args = [
        "compare",
        "--quaternion",
        *additional_cli_args,
        str(file1),
        str(file2),
    ]
    retcode = main(mains_args)
    out, err = capsys.readouterr()

    match = re.search(r"Cartesian RMSD:\s*([-+]?\d*\.\d+|\d+) Å", out)
    assert match, f"Cartesian RMSD not found in output: {out}"

    rmsd_value = float(match.group(1))
    assert pytest.approx(rmsd_value, abs=1e-6) == expected_rmsd


@pytest.mark.parametrize(
    "file1_fixture,file2_fixture,additional_cli_args,expected_irmsd",
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
            "caffeine_xyz_obabel_file_fixture",
            [],
            2.2149021075,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_obabel_file_fixture",
            ["--heavy"],
            2.2149021075,
        ),
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_obabel_file_fixture",
            # only check that files are created
            ["--heavy", "--output"],
            2.2149021075,
        ),
    ],
)
def test_cli_compare_irmsd(
    file1_fixture,
    file2_fixture,
    additional_cli_args,
    expected_irmsd,
    request,
    capsys,
):
    file1 = request.getfixturevalue(file1_fixture)
    file2 = request.getfixturevalue(file2_fixture)

    if "--output" in additional_cli_args:
        # insert filename after --output
        output_filename = str(file1.parent / "output.xyz")
        additional_cli_args.insert(
            additional_cli_args.index("--output") + 1, output_filename
        )

    mains_args = [
        "compare",
        *additional_cli_args,
        str(file1),
        str(file2),
    ]
    retcode = main(mains_args)
    out, err = capsys.readouterr()

    match = re.search(r"iRMSD:\s*([-+]?\d*\.\d+|\d+) Å", out)
    assert match, f"iRMSD not found in output: {out}"

    rmsd_value = float(match.group(1))
    assert pytest.approx(rmsd_value, abs=1e-6) == expected_irmsd

    if "--output" in additional_cli_args:
        outfile_path = Path(output_filename)
        outfile_path_ref = outfile_path.with_stem(outfile_path.stem + "_ref")
        outfile_path_aligned = outfile_path.with_stem(outfile_path.stem + "_aligned")

        assert outfile_path_ref.exists(), "Reference output file not created"
        assert outfile_path_aligned.exists(), "Aligned output file not created"


@pytest.mark.parametrize(
    "file_fixtures,additional_cli_args,expected_delta_rmsds",
    [
        (
            [
                "caffeine_xyz_file_fixture",
                "caffeine_xyz_file_fixture",
                "caffeine_xyz_obabel_file_fixture",
            ],
            [],
            [0.0, 0.1394],
        ),
        (
            [
                "caffeine_xyz_file_fixture",
                "caffeine_xyz_file_fixture",
                "caffeine_xyz_obabel_file_fixture",
            ],
            ["--align"],
            [0.0, 0.0, 2.2149],
        ),
    ],
)
def test_cli_sorter(
    file_fixtures, additional_cli_args, expected_delta_rmsds, request, capfd
):
    files = [request.getfixturevalue(fixture) for fixture in file_fixtures]

    mains_args = [
        "sort",
        *additional_cli_args,
        *(str(f) for f in files),
    ]
    retcode = main(mains_args)
    out, err = capfd.readouterr()

    match_no_structures = re.search(r"number of structures\s*:\s*(\d+)", out)
    no_structures = int(match_no_structures.group(1))
    assert no_structures == len(
        files
    ), f"Number of structures mismatch: expected {len(files)}, got {no_structures}"

    delta_rmsd_values = re.findall(
        r"(?<=\s)\d+\.\d{4}(?=\s*$)", out, flags=re.MULTILINE
    )

    assert len(delta_rmsd_values) == len(
        expected_delta_rmsds
    ), f"Number of delta RMSD values mismatch: expected {len(expected_delta_rmsds)}, got {len(delta_rmsd_values)}"
    for i, delta_rmsd_str in enumerate(delta_rmsd_values):
        delta_rmsd_value = float(delta_rmsd_str)
        assert (
            pytest.approx(delta_rmsd_value, abs=1e-4) == expected_delta_rmsds[i]
        ), f"Delta RMSD value mismatch at index {i}: expected {expected_delta_rmsds[i]}, got {delta_rmsd_value}"


# def test_cli_runs():
#    exe = shutil.which("irmsd")
#    assert exe, "CLI not found on PATH (editable install?)"
#    out = subprocess.check_output([exe, "--a","0.5","--x","1","2","--y","10","20"], text=True)
#    assert out.strip() == "10.5 21.0"
