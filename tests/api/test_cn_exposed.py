import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.cn_exposed import get_cn_fortran


def test_get_cn_fortran(caffeine_cn_test_data):
    conformer1, expected_cn = caffeine_cn_test_data
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)

    cn = get_cn_fortran(
        atom_numbers1,
        positions1,
    )

    assert pytest.approx(expected_cn, abs=1e-6) == cn
