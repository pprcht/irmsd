import numpy as np
import pytest

from irmsd.bindings.rmsd_exposed import get_quaternion_rmsd_fortran_raw

# ------------------------------------------------------
# Utility: correct inputs
# ------------------------------------------------------


def make_valid(n1=3, n2=4):
    types1 = np.arange(n1, dtype=np.int32)
    types2 = np.arange(n2, dtype=np.int32)

    coords1 = np.arange(3 * n1, dtype=np.float64)
    coords2 = np.arange(3 * n2, dtype=np.float64)

    Umat_F = np.zeros((3, 3), dtype=np.float64, order="F")

    mask = np.ones(n1, dtype=bool)

    return types1, coords1, types2, coords2, Umat_F, mask


# ======================================================
# TYPE CHECKS — dtype & contiguity
# ======================================================


def test_types1_wrong_dtype():
    t1, c1, t2, c2, U, m = make_valid()
    with pytest.raises(TypeError, match="types1 must be int32"):
        get_quaternion_rmsd_fortran_raw(3, t1.astype(np.int64), c1, 4, t2, c2, U)


def test_types1_not_c_contiguous():
    t1, c1, t2, c2, U, m = make_valid()
    bad = np.arange(6, dtype=np.int32).reshape(2, 3)[:, 1]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="types1 must be int32 and C-contiguous"):
        get_quaternion_rmsd_fortran_raw(bad.size, bad, c1[: 3 * bad.size], 4, t2, c2, U)


def test_types2_wrong_dtype():
    t1, c1, t2, c2, U, m = make_valid()
    with pytest.raises(TypeError, match="types2 must be int32"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2.astype(np.int64), c2, U)


def test_types2_not_c_contiguous():
    t1, c1, _, _, U, m = make_valid(n1=4, n2=3)
    bad = np.arange(6, dtype=np.int32).reshape(2, 3)[:, 1]
    assert not bad.flags.c_contiguous
    # adjust arrays to match sizes
    n2 = bad.size
    c2 = np.arange(3 * n2, dtype=np.float64)
    with pytest.raises(TypeError, match="types2 must be int32 and C-contiguous"):
        get_quaternion_rmsd_fortran_raw(4, t1, c1, n2, bad, c2, U)


def test_coords1_wrong_dtype():
    t1, c1, t2, c2, U, m = make_valid()
    with pytest.raises(TypeError, match="coords1_flat must be float64"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1.astype(np.float32), 4, t2, c2, U)


def test_coords1_not_c_contiguous():
    t1, c1, t2, c2, U, m = make_valid()
    bad = np.arange(c1.size, dtype=np.float64).reshape(3, -1)[:, 0]
    assert not bad.flags.c_contiguous
    with pytest.raises(
        TypeError, match="coords1_flat must be float64 and C-contiguous"
    ):
        get_quaternion_rmsd_fortran_raw(t1.size, t1, bad, 4, t2, c2, U)


def test_coords2_wrong_dtype():
    t1, c1, t2, c2, U, m = make_valid()
    with pytest.raises(TypeError, match="coords2_flat must be float64"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2.astype(np.float32), U)


def test_coords2_not_c_contiguous():
    t1, c1, t2, c2, U, m = make_valid()
    bad = np.arange(c2.size, dtype=np.float64).reshape(3, -1)[:, 0]
    assert not bad.flags.c_contiguous
    with pytest.raises(
        TypeError, match="coords2_flat must be float64 and C-contiguous"
    ):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, t2.size, t2, bad, U)


# ======================================================
# Umat_F: must be float64, shape (3,3), F-contiguous
# ======================================================


def test_Umat_wrong_dtype():
    t1, c1, t2, c2, U, m = make_valid()
    U2 = U.astype(np.float32)
    with pytest.raises(TypeError, match="Umat_F must be float64"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2, U2)


def test_Umat_wrong_shape():
    t1, c1, t2, c2, U, m = make_valid()
    bad = np.zeros((3, 2), dtype=np.float64, order="F")
    with pytest.raises(TypeError, match="Umat_F must be float64, shape \\(3,3\\)"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2, bad)


def test_Umat_not_F_contiguous():
    t1, c1, t2, c2, U, m = make_valid()
    bad = U.copy(order="C")  # 3x3 but C-contiguous
    assert not bad.flags.f_contiguous
    with pytest.raises(TypeError, match="Fortran-contiguous"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2, bad)


# ======================================================
# VALUE CHECKS — length mismatches
# ======================================================


def test_types1_wrong_length():
    t1, c1, t2, c2, U, m = make_valid(n1=3)
    with pytest.raises(ValueError, match="types1 length must be n1"):
        get_quaternion_rmsd_fortran_raw(3, t1[:-1], c1, 4, t2, c2, U)


def test_coords1_wrong_length():
    t1, c1, t2, c2, U, m = make_valid(n1=3)
    with pytest.raises(ValueError, match="coords1_flat length must be 3\\*n1"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1[:-1], 4, t2, c2, U)


def test_types2_wrong_length():
    t1, c1, t2, c2, U, m = make_valid(n2=4)
    with pytest.raises(ValueError, match="types2 length must be n2"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2[:-1], c2, U)


def test_coords2_wrong_length():
    t1, c1, t2, c2, U, m = make_valid(n2=4)
    with pytest.raises(ValueError, match="coords2_flat length must be 3\\*n2"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2[:-1], U)


# ======================================================
# MASK — optional, must be bool, C-contiguous, size n1
# ======================================================


def test_mask_wrong_dtype():
    t1, c1, t2, c2, U, m = make_valid()
    bad = m.astype(np.int32)
    with pytest.raises(TypeError, match="mask must be bool"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2, U, mask=bad)


def test_mask_not_c_contiguous():
    t1, c1, t2, c2, U, m = make_valid()
    bad = np.asarray([True, False, True, True, True, False]).reshape(2, 3)[
        :, 1
    ]  # slice → non-C
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="mask must be bool, C-contiguous"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2, U, mask=bad)


def test_mask_wrong_size():
    t1, c1, t2, c2, U, m = make_valid(n1=3)
    bad = np.ones(2, dtype=bool)
    with pytest.raises(TypeError, match="mask must be bool.*size natoms"):
        get_quaternion_rmsd_fortran_raw(3, t1, c1, 4, t2, c2, U, mask=bad)
