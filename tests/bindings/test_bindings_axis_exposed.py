import numpy as np
import pytest

from irmsd.bindings.axis_exposed import get_axis_fortran_raw


def make_valid_inputs(n=4):
    types = np.arange(n, dtype=np.int32)
    coords = np.arange(3 * n, dtype=np.float64)
    rot = np.zeros(3, dtype=np.float64)
    avmom = np.zeros(1, dtype=np.float64)
    evec = np.zeros((3, 3), dtype=np.float64)
    return types, coords, rot, avmom, evec


# ------------------------------------------------
# WRONG TYPES / DTYPE / CONTIGUITY
# ------------------------------------------------


def test_types_wrong_dtype():
    types, coords, rot, avmom, evec = make_valid_inputs()
    types = types.astype(np.int64)
    with pytest.raises(TypeError, match="types must be int32"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_types_not_c_contiguous():
    types, coords, rot, avmom, evec = make_valid_inputs()
    types = np.arange(8, dtype=np.int32).reshape(2, 4)[
        :, 1
    ]  # slicing breaks contiguity
    assert not types.flags.c_contiguous
    with pytest.raises(TypeError, match="types must be int32 and C-contiguous"):
        get_axis_fortran_raw(types.size, types, coords, rot, avmom, evec)


def test_coords_wrong_dtype():
    types, coords, rot, avmom, evec = make_valid_inputs()
    coords = coords.astype(np.float32)
    with pytest.raises(TypeError, match="coords_flat must be float64"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_coords_not_c_contiguous():
    types, coords, rot, avmom, evec = make_valid_inputs()
    # slicing like coords[::2] breaks contiguity (size is wrong too; we fix size afterward)
    c = np.arange(12, dtype=np.float64).reshape(3, 4)[:, ::2]  # non-contiguous 1D
    with pytest.raises(TypeError, match="coords_flat must be float64 and C-contiguous"):
        get_axis_fortran_raw(types.size, types, c, rot, avmom, evec)


def test_rot_wrong_dtype():
    types, coords, rot, avmom, evec = make_valid_inputs()
    rot = rot.astype(np.float32)
    with pytest.raises(TypeError, match="rot must be float64"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_rot_wrong_size():
    types, coords, rot, avmom, evec = make_valid_inputs()
    rot = np.zeros(2, dtype=np.float64)
    with pytest.raises(
        TypeError, match="rot must be float64, C-contiguous, shape \\(3\\)"
    ):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_rot_not_c_contiguous():
    types, coords, rot, avmom, evec = make_valid_inputs()
    rot = np.arange(9, dtype=np.float64).reshape(3, 3)[
        :, 0
    ]  # a slice → non-contiguous 1D
    assert not rot.flags.c_contiguous
    with pytest.raises(TypeError, match="rot must be float64, C-contiguous"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_evec_wrong_dtype():
    types, coords, rot, avmom, evec = make_valid_inputs()
    evec = evec.astype(np.float32)
    with pytest.raises(TypeError, match="evec must be float64"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_evec_wrong_shape():
    types, coords, rot, avmom, evec = make_valid_inputs()
    evec = np.zeros((3, 2), dtype=np.float64)
    with pytest.raises(
        TypeError, match="evec must be float64, C-contiguous, shape \\(3, 3\\)"
    ):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_evec_not_c_contiguous():
    types, coords, rot, avmom, evec = make_valid_inputs()
    # transpose breaks C-contiguity for 2D arrays
    evec = evec.T
    assert not evec.flags.c_contiguous
    with pytest.raises(TypeError, match="evec must be float64, C-contiguous"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_avmom_wrong_dtype():
    types, coords, rot, avmom, evec = make_valid_inputs()
    avmom = avmom.astype(np.float32)
    with pytest.raises(TypeError, match="avmom must be float64"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_avmom_wrong_size():
    types, coords, rot, avmom, evec = make_valid_inputs()
    avmom = np.zeros(2, dtype=np.float64)
    with pytest.raises(TypeError, match="avmom must be float64, C-contiguous, size 1"):
        get_axis_fortran_raw(len(types), types, coords, rot, avmom, evec)


def test_avmom_not_c_contiguous():
    types, coords, rot, avmom, evec = make_valid_inputs()
    # slice a 2-element vector to break contiguity
    x = np.arange(4, dtype=np.float64)
    avmom = x[::2]  # size 1, but non-contiguous
    with pytest.raises(TypeError, match="avmom must be float64, C-contiguous"):
        get_axis_fortran_raw(len(x), types[:1], coords[:3], rot, avmom, evec)


# ------------------------------------------------
# VALUE ERROR BRANCHES
# ------------------------------------------------


def test_coords_wrong_length():
    types, coords, rot, avmom, evec = make_valid_inputs(n=4)
    coords = coords[:-1]  # length 11 instead of 12
    with pytest.raises(ValueError, match="coords_flat length must be 3\\*natoms"):
        get_axis_fortran_raw(4, types, coords, rot, avmom, evec)


def test_types_wrong_length():
    types, coords, rot, avmom, evec = make_valid_inputs(n=4)
    types = types[:-1]
    with pytest.raises(ValueError, match="types length must be natoms"):
        get_axis_fortran_raw(4, types, coords, rot, avmom, evec)
