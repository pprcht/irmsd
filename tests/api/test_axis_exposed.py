import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.axis_exposed import get_axis_fortran


def test_get_axis_fortran(caffeine_axis_test_data):
    conformer1, expected_rot, expected_avmom, expected_evec = caffeine_axis_test_data
    (
        atom_numbers1,
        positions1,
    ) = get_atom_num_and_pos_from_xyz(conformer1)

    rot, avmom, evec = get_axis_fortran(
        atom_numbers1,
        positions1,
    )

    assert pytest.approx(expected_rot, abs=1e-6) == rot
    assert pytest.approx(expected_avmom, abs=1e-6) == avmom
    assert pytest.approx(expected_evec, abs=1e-6) == evec
