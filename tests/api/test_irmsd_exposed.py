import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.irmsd_exposed import get_irmsd


def test_get_irmsd(caffeine_irmsd_test_data):
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)
    atom_numbers2, positions2 = get_atom_num_and_pos_from_xyz(conformer2)
    expected_atom_numbers, expected_positions = get_atom_num_and_pos_from_xyz(
        expected_aligned_conformer
    )

    (
        rmsd,
        atom_numbers_aligned1,
        positions_aligned1,
        atom_numbers_aligned2,
        positions_aligned2,
    ) = get_irmsd(
        atom_numbers1, positions1, atom_numbers2, positions2, iinversion=1
    )

    assert pytest.approx(expected_irmsd, abs=1e-6) == rmsd
