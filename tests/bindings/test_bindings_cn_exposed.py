# test_get_cn_fortran_raw_input_validation.py
import numpy as np
import pytest

from irmsd.bindings.cn_exposed import get_cn_fortran_raw


def make_valid_inputs(n=4):
    types = np.arange(n, dtype=np.int32)
    coords = np.arange(3 * n, dtype=np.float64)
    cn = np.zeros(n, dtype=np.float64)
    return types, coords, cn


# ------------------------------
# TYPE CHECKS
# ------------------------------


def test_types_wrong_dtype():
    types, coords, cn = make_valid_inputs()
    types = types.astype(np.int64)  # wrong dtype
    with pytest.raises(TypeError, match="types must be int32"):
        get_cn_fortran_raw(len(types), types, coords, cn)


def test_types_not_c_contiguous():
    types, coords, cn = make_valid_inputs()
    # actually break C-contiguous requirement by slicing with a step of 2 and then taking every other element
    types = types[::2]
    with pytest.raises(TypeError, match="types must be int32 and C-contiguous"):
        get_cn_fortran_raw(len(types), types, coords, cn)


def test_coords_wrong_dtype():
    types, coords, cn = make_valid_inputs()
    coords = coords.astype(np.float32)
    with pytest.raises(TypeError, match="coords_flat must be float64"):
        get_cn_fortran_raw(len(types), types, coords, cn)


def test_coords_not_c_contiguous():
    types, coords, cn = make_valid_inputs()
    # actually break C-contiguous requirement by slicing with a step of 2 and then taking every other element
    coords = coords[::2]
    with pytest.raises(TypeError, match="coords_flat must be float64 and C-contiguous"):
        get_cn_fortran_raw(len(types), types, coords, cn)


def test_cn_wrong_dtype():
    types, coords, cn = make_valid_inputs()
    cn = cn.astype(np.float32)
    with pytest.raises(TypeError, match="cn_flat must be float64"):
        get_cn_fortran_raw(len(types), types, coords, cn)


def test_cn_not_c_contiguous():
    types, coords, cn = make_valid_inputs()
    # actually break C-contiguous requirement by slicing with a step of 2 and then taking every other element
    cn = cn[::2]
    with pytest.raises(TypeError, match="cn_flat must be float64, C-contiguous"):
        get_cn_fortran_raw(len(types), types, coords, cn)


def test_cn_wrong_size():
    types, coords, cn = make_valid_inputs(n=4)
    cn = np.zeros(3, dtype=np.float64)  # wrong size
    with pytest.raises(TypeError, match="shape \\(natoms\\)"):
        get_cn_fortran_raw(4, types, coords, cn)


# ------------------------------
# VALUE CHECKS
# ------------------------------


def test_coords_wrong_length():
    types, coords, cn = make_valid_inputs(n=4)
    coords = coords[:-1]  # length = 11 != 12
    with pytest.raises(ValueError, match="coords_flat length must be 3\\*natoms"):
        get_cn_fortran_raw(4, types, coords, cn)


def test_types_wrong_length():
    types, coords, cn = make_valid_inputs(n=4)
    types = types[:-1]  # length = 3 != 4
    with pytest.raises(ValueError, match="types length must be natoms"):
        get_cn_fortran_raw(4, types, coords, cn)
