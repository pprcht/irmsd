import numpy as np
import pytest

from irmsd.bindings.irmsd_exposed import get_irmsd_fortran_raw

# -------------------------------------------------------------
# Utility: produce valid baseline test arrays
# -------------------------------------------------------------


def make_valid(n1=3, n2=4):
    types1 = np.arange(n1, dtype=np.int32)
    types2 = np.arange(n2, dtype=np.int32)

    coords1 = np.arange(3 * n1, dtype=np.float64)
    coords2 = np.arange(3 * n2, dtype=np.float64)

    types_out1 = np.zeros(n2, dtype=np.int32)
    types_out2 = np.zeros(n2, dtype=np.int32)

    coords_out1 = np.zeros(3 * n2, dtype=np.float64)
    coords_out2 = np.zeros(3 * n2, dtype=np.float64)

    return (
        types1,
        coords1,
        types2,
        coords2,
        types_out1,
        coords_out1,
        types_out2,
        coords_out2,
    )


# -------------------------------------------------------------
# TYPE CHECKS — wrong dtype or wrong contiguity
# -------------------------------------------------------------


def test_types1_wrong_dtype():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = t1.astype(np.int64)
    with pytest.raises(TypeError, match="types1 must be int32"):
        get_irmsd_fortran_raw(3, bad, c1, 4, t2, c2, 0, to1, co1, to2, co2)


def test_types1_not_c_contiguous():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = np.arange(6, dtype=np.int32).reshape(2, 3)[:, 1]  # slice → non-contiguous
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="types1 must be int32 and C-contiguous"):
        get_irmsd_fortran_raw(
            bad.size, bad, c1[: 3 * bad.size], 4, t2, c2, 0, to1, co1, to2, co2
        )


def test_types2_wrong_dtype():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = t2.astype(np.int64)
    with pytest.raises(TypeError, match="types2 must be int32"):
        get_irmsd_fortran_raw(3, t1, c1, 4, bad, c2, 0, to1, co1, to2, co2)


def test_types2_not_c_contiguous():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=4, n2=3)
    bad = np.arange(6, dtype=np.int32).reshape(2, 3)[:, 1]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="types2 must be int32 and C-contiguous"):
        get_irmsd_fortran_raw(
            4,
            t1,
            c1,
            bad.size,
            bad,
            c2[: 3 * bad.size],
            0,
            to1[: bad.size],
            co1[: 3 * bad.size],
            to2[: bad.size],
            co2[: 3 * bad.size],
        )


def test_types_out1_wrong_dtype():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = to1.astype(np.int64)
    with pytest.raises(TypeError, match="types_out must be int32"):
        get_irmsd_fortran_raw(3, t1, c1, 4, t2, c2, 0, bad, co1, to2, co2)


def test_types_out1_not_c_contiguous():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = np.arange(6, dtype=np.int32).reshape(2, 3)[:, 1]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="types_out must be int32 and C-contiguous"):
        get_irmsd_fortran_raw(
            3,
            t1,
            c1,
            bad.size,
            t2[: bad.size],
            c2[: 3 * bad.size],
            0,
            bad,
            co1[: 3 * bad.size],
            to2[: bad.size],
            co2[: 3 * bad.size],
        )


def test_coords1_wrong_dtype():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    with pytest.raises(TypeError, match="coords1_flat must be float64"):
        get_irmsd_fortran_raw(
            3, t1, c1.astype(np.float32), 4, t2, c2, 0, to1, co1, to2, co2
        )


def test_coords1_not_c_contiguous():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = np.arange(c1.size, dtype=np.float64).reshape(3, -1)[:, 0]
    assert not bad.flags.c_contiguous
    with pytest.raises(
        TypeError, match="coords1_flat must be float64 and C-contiguous"
    ):
        get_irmsd_fortran_raw(t1.size, t1, bad, 4, t2, c2, 0, to1, co1, to2, co2)


def test_coords2_wrong_dtype():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    with pytest.raises(TypeError, match="coords2_flat must be float64"):
        get_irmsd_fortran_raw(
            3, t1, c1, 4, t2, c2.astype(np.float32), 0, to1, co1, to2, co2
        )


def test_coords2_not_c_contiguous():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = np.arange(c2.size, dtype=np.float64).reshape(3, -1)[:, 0]
    assert not bad.flags.c_contiguous
    with pytest.raises(
        TypeError, match="coords2_flat must be float64 and C-contiguous"
    ):
        get_irmsd_fortran_raw(3, t1, c1, t2.size, t2, bad, 0, to1, co1, to2, co2)


def test_coords_out1_wrong_dtype():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    with pytest.raises(TypeError, match="coords_out_flat must be float64"):
        get_irmsd_fortran_raw(
            3, t1, c1, 4, t2, c2, 0, to1, co1.astype(np.float32), to2, co2
        )


def test_coords_out1_not_c_contiguous():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid()
    bad = np.arange(co1.size, dtype=np.float64).reshape(3, -1)[:, 0]
    assert not bad.flags.c_contiguous
    with pytest.raises(
        TypeError, match="coords_out_flat must be float64 and C-contiguous"
    ):
        get_irmsd_fortran_raw(3, t1, c1, 4, t2, c2, 0, to1, bad, to2, co2)


# -------------------------------------------------------------
# VALUE CHECKS — wrong lengths or mismatched sizes
# -------------------------------------------------------------


def test_types1_wrong_length():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=3, n2=4)
    bad = t1[:-1]
    with pytest.raises(ValueError, match="types1 length must be n1"):
        get_irmsd_fortran_raw(3, bad, c1, 4, t2, c2, 0, to1, co1, to2, co2)


def test_coords1_wrong_length():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=3, n2=4)
    bad = c1[:-1]
    with pytest.raises(ValueError, match="coords1_flat length must be 3\\*n1"):
        get_irmsd_fortran_raw(3, t1, bad, 4, t2, c2, 0, to1, co1, to2, co2)


def test_types2_wrong_length():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=3, n2=4)
    bad = t2[:-1]
    with pytest.raises(ValueError, match="types2 length must be n2"):
        get_irmsd_fortran_raw(3, t1, c1, 4, bad, c2, 0, to1, co1, to2, co2)


def test_coords2_wrong_length():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=3, n2=4)
    bad = c2[:-1]
    with pytest.raises(ValueError, match="coords2_flat length must be 3\\*n2"):
        get_irmsd_fortran_raw(3, t1, c1, 4, t2, bad, 0, to1, co1, to2, co2)


def test_types_out1_wrong_length():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=3, n2=4)
    bad = to1[:-1]
    with pytest.raises(
        ValueError, match="types_out length must be equal to types2 length"
    ):
        get_irmsd_fortran_raw(3, t1, c1, 4, t2, c2, 0, bad, co1, to2, co2)


def test_coords_out1_wrong_length():
    t1, c1, t2, c2, to1, co1, to2, co2 = make_valid(n1=3, n2=4)
    bad = co1[:-1]
    with pytest.raises(
        ValueError, match="coords_out_flat length must be equal to coords2_flat length"
    ):
        get_irmsd_fortran_raw(3, t1, c1, 4, t2, c2, 0, to1, bad, to2, co2)
