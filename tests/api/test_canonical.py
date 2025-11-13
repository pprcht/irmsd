import numpy as np
import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.canonical_exposed import (
    get_canonical_fortran,
    get_canonical_from_connect_fortran,
)


@pytest.fixture(
    scope="module",
    params=[
        # conformer1, connectivity_matrix, heavy, canonical_atom_id
        # fmt: off
        (
        "caffeine_molecule_xyz",
        np.asarray([
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0],
        [1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
        [0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
            ],
        # fmt: on
            dtype=np.int32,
        ),
        False,
        np.asarray(
            # fmt: off
            [1,12,6,7,10,11,8,4,13,9,5,14,3,2,15,15,15,18,17,17,17,16,16,16,],
            # fmt: on
            dtype=int,
        ),
        ),
    ],
)
def caffeine_canonical_from_connect_test_data(request):
    return (
        request.getfixturevalue(request.param[0]),
        request.param[1],
        request.param[2],
        request.param[3],
    )


def test_get_canonical_fortran(caffeine_canonical_test_data):
    conformer1, heavy, expected_canonical_atom_id = caffeine_canonical_test_data
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)

    canonical_atom_id = get_canonical_fortran(
        atom_numbers1,
        positions1,
        heavy=heavy,
    )
    assert all(canonical_atom_id == expected_canonical_atom_id)


def test_get_canonical_from_connect_fortran(caffeine_canonical_from_connect_test_data):
    conformer1, connectivity, heavy, expected_canonical_atom_id = (
        caffeine_canonical_from_connect_test_data
    )

    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)

    canonical_atom_id = get_canonical_from_connect_fortran(
        atom_numbers1,
        connectivity,
        heavy=heavy,
    )
    assert all(canonical_atom_id == expected_canonical_atom_id)
