import numpy as np
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd import sorter_irmsd  # adjust if needed


def test_sorter_irmsd(caffeine_irmsd_test_data):
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )

    # Extract coordinates
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)
    atom_numbers2, positions2 = get_atom_num_and_pos_from_xyz(conformer2)

    atom_numbers_list = [atom_numbers1, atom_numbers2]
    positions_list = [positions1, positions2]

    nat = atom_numbers1.shape[0]
    rthr = 0.125

    groups, xyz_structs, Z_structs = sorter_irmsd(
        atom_numbers_list,
        positions_list,
        nat=nat,
        rthr=rthr,
        iinversion=1,
        allcanon=True,
        printlvl=0,
    )

    # Basic structure checks
    assert isinstance(groups, np.ndarray)
    assert groups.dtype == np.int32

    assert len(xyz_structs) == 2
    assert len(Z_structs) == 2

    # Groups should be non-negative integers
    # and both conformers should belong to the same group
    assert np.all(groups >= 0)
    assert groups.shape == (2,)
    np.testing.assert_array_equal(groups, [1, 1])
