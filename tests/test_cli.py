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


@pytest.fixture(scope="session")
def alanine_glycine_xyz_file_fixture(tmp_dir_fixture, alanine_glycine_conformers_xyz):
    """Fixture that creates an alanine.xyz file in the temporary directory."""
    file_path = tmp_dir_fixture / "alanine_glycine.xyz"
    with open(file_path, "w") as f:
        f.write(alanine_glycine_conformers_xyz)
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
        (
            "caffeine_xyz_file_fixture",
            "caffeine_xyz_obabel_file_fixture",
            ["--heavy", "--output"],
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

    if "--output" in additional_cli_args:
        # insert filename after --output
        output_filename = str(file1.parent / "output_rmsd.xyz")
        additional_cli_args.insert(
            additional_cli_args.index("--output") + 1, output_filename
        )

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

    if "--output" in additional_cli_args:
        assert Path(output_filename).exists(), "Output file not created"


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
        output_filename = str(file1.parent / "output_irmsd.xyz")
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


@pytest.mark.parametrize(
    "data,additional_cli_args",
    [
        (
            [
                # file fixture, no of structures, expected delta RMSDs
                [["caffeine_xyz_file_fixture"], 1, [0]],
                # file fixture, no of structures, expected delta RMSDs
                [
                    ["alanine_glycine_xyz_file_fixture"],
                    4,
                    [0.0000, 1.0423, 0.7799, 1.9073],
                ],
            ],
            [],
        ),
        (
            [
                # file fixture, no of structures, expected delta RMSDs
                [["caffeine_xyz_file_fixture", "caffeine_xyz_file_fixture"], 2, [0]],
                # file fixture, no of structures, expected delta RMSDs
                [
                    ["alanine_glycine_xyz_file_fixture"],
                    4,
                    [0.0000, 1.0423, 0.7799, 1.9073],
                ],
            ],
            ["--output"],
        ),
        (
            [
                # file fixture, no of structures, expected delta RMSDs
                [["caffeine_xyz_file_fixture", "caffeine_xyz_file_fixture"], 2, [0, 0]],
                # file fixture, no of structures, expected delta RMSDs
                [
                    ["alanine_glycine_xyz_file_fixture"],
                    4,
                    [0.0000, 1.0423, 0.7799, 1.9073],
                ],
            ],
            ["--align"],
        ),
        (
            [
                # file fixture, no of structures, expected delta RMSDs
                [["caffeine_xyz_file_fixture", "caffeine_xyz_file_fixture"], 2, [0, 0]],
                # file fixture, no of structures, expected delta RMSDs
                [
                    ["alanine_glycine_xyz_file_fixture"],
                    4,
                    [0.0000, 1.0423, 0.7799, 1.9073],
                ],
            ],
            ["--align", "--output"],
        ),
    ],
)
def test_cli_sorter_several_molecules(data, additional_cli_args, request, capfd):
    files = []
    no_structures = []
    expected_delta_rmsds = []
    for file, n_structures, delta_rmsds in data:
        files += [request.getfixturevalue(fname) for fname in file]
        no_structures.append(n_structures)
        expected_delta_rmsds.append(delta_rmsds)
    if "--output" in additional_cli_args:
        # insert filename after --output
        output_filename = str(files[0].parent / "output_several.xyz")
        additional_cli_args.insert(
            additional_cli_args.index("--output") + 1, output_filename
        )

    mains_args = [
        "sort",
        *additional_cli_args,
        *(str(f) for f in files),
    ]
    retcode = main(mains_args)
    out, err = capfd.readouterr()

    match_no_structures = re.findall(r"number of structures\s*:\s*(\d+)", out)
    for i, match in enumerate(match_no_structures):
        n_structures = int(match)
        assert n_structures == no_structures[i], (
            f"Number of structures mismatch for file {files[i]}: "
            f"expected {no_structures[i]}, got {n_structures}"
        )

    delta_rmsd_values = re.findall(
        r"(?<=\s)\d+\.\d{4}(?=\s*$)", out, flags=re.MULTILINE
    )
    current_index = 0
    for i, expected_deltas in enumerate(expected_delta_rmsds):
        for j, expected_delta in enumerate(expected_deltas):
            delta_rmsd_str = delta_rmsd_values[current_index]
            delta_rmsd_value = float(delta_rmsd_str)
            assert (
                pytest.approx(delta_rmsd_value, abs=1e-4) == expected_delta
            ), f"Delta RMSD value mismatch for file {files[i]} at index {j}: expected {expected_delta}, got {delta_rmsd_value}"
            current_index += 1

    if "--output" in additional_cli_args:
        outfile_path = Path(output_filename).parent
        assert len(list(outfile_path.glob("*output_several_*.xyz"))) > 0


@pytest.mark.parametrize(
    "file_fixtures,additional_cli_args,property_printout",
    [
        (
            [
                "caffeine_xyz_file_fixture",
            ],
            ["--cn"],
            """Atom Symbol             CN
---- ------ --------------
   1      C       3.687259
   2      N       2.867262
   3      C       3.150725
   4      N       1.893570
   5      C       3.203212
   6      C       3.080100
   7      C       2.766851
   8      O       0.858278
   9      N       2.742942
  10      C       2.714023
  11      O       0.857489
  12      N       2.737923
  13      C       3.693785
  14      C       3.700082
  15      H       0.924157
  16      H       0.924170
  17      H       0.924168
  18      H       0.925322
  19      H       0.924123
  20      H       0.923952
  21      H       0.923952
  22      H       0.924191
  23      H       0.923896
  24      H       0.923906
""",
        ),
        (
            [
                "caffeine_xyz_file_fixture",
            ],
            ["--canonical"],
            """Atom Symbol Canonical Rank
---- ------ --------------
   1      C              1
   2      N              8
   3      C              5
   4      N              7
   5      C             13
   6      C             14
   7      C             12
   8      O              6
   9      N             10
  10      C              9
  11      O              2
  12      N             11
  13      C              4
  14      C              3
  15      H             15
  16      H             15
  17      H             15
  18      H             18
  19      H             17
  20      H             17
  21      H             17
  22      H             16
  23      H             16
  24      H             16
""",
        ),
        (
            [
                "caffeine_xyz_file_fixture",
            ],
            ["--rot"],
            """Rotational constants (MHz):
1068.0731    710.7118    430.2612

Average momentum a.u. (10⁻⁴⁷kg m²): 1.305645e-44

Rotation matrix:
  0.6712     -0.7413      0.0012
  0.7413      0.6712     -0.0011
  0.0000      0.0017      1.0000
""",
        ),
    ],
)
def test_cli_prop(
    file_fixtures, additional_cli_args, property_printout, request, capfd
):
    files = [request.getfixturevalue(fixture) for fixture in file_fixtures]

    mains_args = [
        "prop",
        *additional_cli_args,
        *(str(f) for f in files),
    ]
    retcode = main(mains_args)
    out, err = capfd.readouterr()

    assert (
        property_printout in out
    ), f"Expected property printout not found in output:\n{out}"


def test_cli_help(caffeine_xyz_file_fixture, capfd):
    mains_args = ["prop", str(caffeine_xyz_file_fixture)]
    retcode = main(mains_args)
    out, err = capfd.readouterr()

    assert (
        retcode == 1
    ), f"Expected return code 1 when no property selected, got {retcode}"
    assert (
        "CLI to read an arbitrary number of structures with ASE and run selected" in out
    ), f"Help message not found in output:\n{out}"
