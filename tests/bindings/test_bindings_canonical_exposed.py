import numpy as np
import pytest

from irmsd.bindings.canonical_exposed import (
    get_canonical_sorter_fortran_raw,
    get_ids_from_connect_fortran_raw,
)

# ===============================================================
# Utility: valid baseline inputs
# ===============================================================


def make_valid_inputs(n=4):
    types = np.arange(n, dtype=np.int32)
    coords = np.arange(3 * n, dtype=np.float64)
    rank = np.zeros(n, dtype=np.int32)
    wbo = np.zeros((n, n), dtype=np.float64)
    return types, coords, rank, wbo


# ===============================================================
# Tests for get_canonical_sorter_fortran_raw
# ===============================================================


def test_types_wrong_dtype():
    types, coords, rank, wbo = make_valid_inputs()
    types = types.astype(np.int64)
    with pytest.raises(TypeError, match="types must be int32"):
        get_canonical_sorter_fortran_raw(len(types), types, coords, rank)


def test_types_not_c_contiguous():
    types, coords, rank, wbo = make_valid_inputs()
    bad = np.arange(8, dtype=np.int32).reshape(2, 4)[:, 1]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="types must be int32 and C-contiguous"):
        get_canonical_sorter_fortran_raw(bad.size, bad, coords, rank)


def test_coords_wrong_dtype():
    types, coords, rank, wbo = make_valid_inputs()
    with pytest.raises(TypeError, match="coords_flat must be float64"):
        get_canonical_sorter_fortran_raw(
            len(types), types, coords.astype(np.float32), rank
        )


def test_coords_not_c_contiguous():
    types, coords, rank, wbo = make_valid_inputs()
    # slicing from a reshaped array breaks contiguity
    bad = np.arange(coords.size, dtype=np.float64).reshape(3, -1)[:, 0]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="coords_flat must be float64 and C-contiguous"):
        get_canonical_sorter_fortran_raw(types.size, types, bad, rank)


def test_coords_wrong_length():
    types, coords, rank, wbo = make_valid_inputs(n=4)
    coords = coords[:-1]
    with pytest.raises(ValueError, match="coords_flat length must be 3\\*natoms"):
        get_canonical_sorter_fortran_raw(4, types, coords, rank)


def test_types_wrong_length():
    types, coords, rank, wbo = make_valid_inputs(n=4)
    types = types[:-1]
    with pytest.raises(ValueError, match="types length must be natoms"):
        get_canonical_sorter_fortran_raw(4, types, coords, rank)


# ---------------- Rank array checks ----------------


def test_rank_wrong_dtype():
    types, coords, rank, wbo = make_valid_inputs()
    bad = rank.astype(np.int64)
    with pytest.raises(TypeError, match="rank must be int32"):
        get_canonical_sorter_fortran_raw(len(types), types, coords, bad)


def test_rank_not_c_contiguous():
    types, coords, rank, wbo = make_valid_inputs()
    bad = np.arange(8, dtype=np.int32).reshape(2, 4)[:, 1]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="rank must be int32, C-contiguous"):
        get_canonical_sorter_fortran_raw(
            bad.size, types[: bad.size], coords[: 3 * bad.size], bad
        )


def test_rank_wrong_size():
    types, coords, rank, wbo = make_valid_inputs(n=4)
    bad = np.zeros(3, dtype=np.int32)
    with pytest.raises(TypeError, match="size natoms"):
        get_canonical_sorter_fortran_raw(4, types, coords, bad)


# ---------------- wbo optional matrix ----------------


def test_wbo_wrong_dtype():
    types, coords, rank, wbo = make_valid_inputs()
    bad = wbo.astype(np.float32)
    with pytest.raises(TypeError, match="wbo must be float64"):
        get_canonical_sorter_fortran_raw(len(types), types, coords, rank, wbo=bad)


def test_wbo_wrong_shape():
    types, coords, rank, wbo = make_valid_inputs(n=4)
    bad = np.zeros((3, 4), dtype=np.float64)
    with pytest.raises(TypeError, match="wbo must be float64, C-contiguous, shape"):
        get_canonical_sorter_fortran_raw(4, types, coords, rank, wbo=bad)


def test_wbo_not_c_contiguous():
    types, coords, rank, wbo = make_valid_inputs()
    bad = wbo.T
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="wbo must be float64, C-contiguous"):
        get_canonical_sorter_fortran_raw(len(types), types, coords, rank, wbo=bad)


# ===============================================================
# Tests for get_ids_from_connect_fortran_raw
# ===============================================================


def make_connectivity(n=4):
    return np.zeros((n, n), dtype=np.int32, order="F")


def test_connect_types_wrong_dtype():
    n = 4
    types = np.arange(n, dtype=np.int64)
    conn = make_connectivity(n)
    rank = np.zeros(n, dtype=np.int32)
    with pytest.raises(TypeError, match="types must be int32"):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)


def test_connect_types_not_c_contiguous():
    n = 4
    t = np.arange(8, dtype=np.int32).reshape(2, 4)[:, 1]
    assert not t.flags.c_contiguous
    conn = make_connectivity(t.size)
    rank = np.zeros(t.size, dtype=np.int32)
    with pytest.raises(TypeError, match="types must be int32 and C-contiguous"):
        get_ids_from_connect_fortran_raw(t.size, t, conn, rank)


def test_connect_flat_not_f_contiguous():
    n = 4
    types = np.arange(n, dtype=np.int32)
    conn = np.zeros((n, n), dtype=np.int32)  # C-order, not F-order
    rank = np.zeros(n, dtype=np.int32)
    assert not conn.flags.f_contiguous
    with pytest.raises(
        TypeError, match="connectivity_flat must be int32 and F-contiguous"
    ):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)


def test_connectivity_wrong_dtype():
    n = 4
    types = np.arange(n, dtype=np.int32)
    conn = np.zeros((n, n), dtype=np.float64, order="F")
    rank = np.zeros(n, dtype=np.int32)
    with pytest.raises(TypeError, match="connectivity_flat must be int32"):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)


def test_connectivity_wrong_length():
    n = 4
    types = np.arange(n, dtype=np.int32)
    conn = np.zeros((n * n - 1,), dtype=np.int32, order="F")
    rank = np.zeros(n, dtype=np.int32)
    with pytest.raises(
        ValueError, match="connectivity_flat length must be natoms\\*natoms"
    ):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)


def test_connect_rank_wrong_dtype():
    n = 4
    types = np.arange(n, dtype=np.int32)
    conn = make_connectivity(n)
    rank = np.zeros(n, dtype=np.int64)
    with pytest.raises(TypeError, match="rank must be int32"):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)


def test_connect_rank_wrong_size():
    n = 4
    types = np.arange(n, dtype=np.int32)
    conn = make_connectivity(n)
    rank = np.zeros(n - 1, dtype=np.int32)
    with pytest.raises(TypeError, match="size natoms"):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)


def test_connect_rank_not_c_contiguous():
    n = 4
    types = np.arange(n, dtype=np.int32)
    conn = make_connectivity(n)
    bad = np.arange(8, dtype=np.int32).reshape(2, 4)[:, 1]
    assert not bad.flags.c_contiguous
    with pytest.raises(TypeError, match="rank must be int32, C-contiguous"):
        get_ids_from_connect_fortran_raw(n, types, conn, bad)


def test_connect_types_wrong_length():
    n = 4
    types = np.arange(n - 1, dtype=np.int32)
    conn = make_connectivity(n)
    rank = np.zeros(n, dtype=np.int32)
    with pytest.raises(ValueError, match="types length must be natoms"):
        get_ids_from_connect_fortran_raw(n, types, conn, rank)
