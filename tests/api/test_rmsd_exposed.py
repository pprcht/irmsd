import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.rmsd_exposed import get_quaternion_rmsd_fortran


# TODO: test mask
def test_get_quaternion_rmsd_fortran(caffeine_rmsd_test_data):
    (
        conformer1,
        conformer2,
        expected_rmsd,
        expected_aligned_conformer,
        expected_Umat,
    ) = caffeine_rmsd_test_data
    mask = None
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)
    atom_numbers2, positions2 = get_atom_num_and_pos_from_xyz(conformer2)
    expected_positions = get_atom_num_and_pos_from_xyz(expected_aligned_conformer)[1]

    rmsd, new_positions, umat = get_quaternion_rmsd_fortran(
        atom_numbers1,
        positions1,
        atom_numbers2,
        positions2,
        mask,
    )

    assert pytest.approx(rmsd, abs=1e-6) == expected_rmsd
    assert pytest.approx(new_positions, abs=1e-6) == expected_positions
    assert pytest.approx(umat, abs=1e-6) == expected_Umat
