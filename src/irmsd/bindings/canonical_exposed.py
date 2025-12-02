from __future__ import annotations

import ctypes as ct

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

# We expose the canonical function as:
#   get_canonical_sorter_fortran_allargs(natoms: c_int,
#                  types: int32[C_CONTIGUOUS](natoms),
#                  coords: float64[C_CONTIGUOUS](3*natoms),
#                  wbo: float64[F_CONTIGUOUS](natoms,natoms) or None,
#                  invtype:
#                  heavy: boolean
#                  rank: int32[C_CONTIGUOUS](natoms),
#                  invariants: int32[C_CONTIGUOUS](natoms),
LIB.get_canonical_sorter_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.POINTER(ct.c_double),  # wbo can be None
    ct.c_char_p,
    ct.c_bool,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
]
LIB.get_canonical_sorter_fortran.restype = None


def get_canonical_sorter_fortran_raw(
    natoms: int,
    types: np.ndarray,
    coords_flat: np.ndarray,
    rank: np.ndarray,
    heavy: bool = False,
    wbo: np.ndarray | None = None,
    invtype: str = "apsp+",
) -> None:
    """Low-level call that matches the Fortran signature exactly. Operates IN-
    PLACE on rank.

    Parameters
    ----------
    natoms : int
        Number of atoms (must be consistent with array lengths).
    types : (natoms,) int32, C-contiguous
        Atomic numbers (or type IDs).
    coords_flat : (3*natoms,) float64, C-contiguous
        Flat coordinates [x1,y1,z1,x2,y2,z2,...].
    rank : (natoms,) int32, C-contiguous
        Output array for rank.
    heavy : bool, optional
        Whether to consider only heavy atoms (default: False).
    wbo: (natoms, natoms) float64, C-contiguous, optional
        Optional Wiberg bond order matrix, required if invtype is 'cangen', ignored in case of 'apsp+'.
    invtype : str, optional
        alogrithm type for invariants calculation (default: apsp+), alternativly 'cangen'.

    Raises
    ------
    TypeError
        If any of the arrays do not have the expected dtype or memory layout.
    ValueError
        If array sizes do not match natoms.
    """
    # light validation to catch mismatches early
    if types.dtype != np.int32 or not types.flags.c_contiguous:
        raise TypeError("types must be int32 and C-contiguous")
    if coords_flat.dtype != np.float64 or not coords_flat.flags.c_contiguous:
        raise TypeError("coords_flat must be float64 and C-contiguous")
    if coords_flat.size != 3 * natoms:
        raise ValueError("coords_flat length must be 3*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")
    if rank.dtype != np.int32 or not rank.flags.c_contiguous or rank.size != natoms:
        raise TypeError("rank must be int32, C-contiguous, size natoms")

    if wbo is not None:
        if (
            wbo.dtype != np.float64
            or not wbo.flags.c_contiguous
            or wbo.shape != (natoms, natoms)
        ):
            raise TypeError("wbo must be float64, C-contiguous, shape (natoms, natoms)")
        wbo_ptr = wbo.ctypes.data_as(ct.POINTER(ct.c_double))
    else:
        wbo_ptr = None

    invtype_bytes = invtype.encode("utf-8")
    LIB.get_canonical_sorter_fortran(
        int(natoms),
        types,
        coords_flat,
        wbo_ptr,
        invtype_bytes,
        bool(heavy),
        rank,
    )


# We expose the canonical function as:
#   get_ids_from_connect_fortran(natoms: c_int,
#                  types: int32[C_CONTIGUOUS](natoms),
#                  connectivity_flat: int32[F_CONTIGUOUS](natoms,natoms) or None,
#                  heavy: boolean
#                  rank: int32[C_CONTIGUOUS](natoms)
LIB.get_ids_from_connect_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),
    ct.c_bool,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
]
LIB.get_ids_from_connect_fortran.restype = None


def get_ids_from_connect_fortran_raw(
    natoms: int,
    types: np.ndarray,
    connectivity_flat: np.ndarray,
    rank: np.ndarray,
    heavy: bool = False,
) -> None:
    """Low-level call that matches the Fortran signature exactly.

    Parameters
    ----------
    natoms : int
        Number of atoms (must be consistent with array lengths).
    types : (natoms,) int32, C-contiguous
        Atomic numbers (or type IDs).
    connectivity_flat : (natoms*natoms,) int32, F-contiguous
        indicates whether a bond between two atoms exists
    rank : (natoms,) int32, C-contiguous
        Output array for rank.
    heavy : bool, optional
        Whether to consider only heavy atoms (default: False).

    Raises
    ------
    TypeError
        If any of the arrays do not have the expected dtype or memory layout.
    ValueError
        If array sizes do not match natoms.
    """
    # light validation to catch mismatches early
    if types.dtype != np.int32 or not types.flags.c_contiguous:
        raise TypeError("types must be int32 and C-contiguous")
    if connectivity_flat.dtype != np.int32 or not connectivity_flat.flags.f_contiguous:
        raise TypeError("connectivity_flat must be int32 and F-contiguous")
    if connectivity_flat.size != natoms * natoms:
        raise ValueError("connectivity_flat length must be natoms*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")
    if rank.dtype != np.int32 or not rank.flags.c_contiguous or rank.size != natoms:
        raise TypeError("rank must be int32, C-contiguous, size natoms")

    LIB.get_ids_from_connect_fortran(
        int(natoms),
        types,
        connectivity_flat,
        bool(heavy),
        rank,
    )
