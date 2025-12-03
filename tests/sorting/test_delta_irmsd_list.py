import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd import delta_irmsd_list


def test_delta_irmsd_list_caffeine(caffeine_irmsd_test_data):
    # Same fixture as for get_irmsd
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )

    # Convert XYZs to arrays
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)
    atom_numbers2, positions2 = get_atom_num_and_pos_from_xyz(conformer2)

    # Build the list inputs for delta_irmsd_list
    atom_numbers_list = [atom_numbers1, atom_numbers2]
    positions_list = [positions1, positions2]
    nat = atom_numbers1.shape[0]

    delta, xyz_structs, Z_structs = delta_irmsd_list(
        atom_numbers_list,
        positions_list,
        nat=nat,
        iinversion=1,
        allcanon=True,
        printlvl=0,
    )

    # Basic sanity checks on shapes
    assert delta.shape[0] == 2
    assert len(xyz_structs) == 2
    assert len(Z_structs) == 2

    # The iRMSD between structure 1 and 1 should be zero
    assert pytest.approx(0.0, abs=1e-6) == delta[0]
    # The iRMSD between structure 1 and 2 should be as the ref data
    assert pytest.approx(expected_irmsd, abs=1e-6) == delta[1]
