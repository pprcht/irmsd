from __future__ import annotations

import ctypes as ct

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

# We expose the CN as:
#   get_cn_fortran(natoms: c_int,
#                  types: int32[C_CONTIGUOUS](natoms),
#                  coords: float64[C_CONTIGUOUS](3*natoms),
#                  cn: float64[F_CONTIGUOUS](natoms))
LIB.get_cn_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
LIB.get_cn_fortran.restype = None


def get_cn_fortran_raw(
    natoms: int,
    types: np.ndarray,
    coords_flat: np.ndarray,
    cn_flat: np.ndarray,
) -> None:
    """Low-level call that matches the Fortran signature exactly. Operates IN-
    PLACE on cn_flat.

    Parameters
    ----------
    natoms : int
        Number of atoms (must be consistent with array lengths).
    types : (natoms,) int32, C-contiguous
        Atomic numbers (or type IDs).
    coords_flat : (3*natoms,) float64, C-contiguous
        Flat coordinates [x1,y1,z1,x2,y2,z2,...].
    cn_flat : (natoms) float64, C-contiguous
        Output CN array written by Fortran.

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
    if (
        cn_flat.dtype != np.float64
        or not cn_flat.flags.c_contiguous
        or cn_flat.size != natoms
    ):
        raise TypeError("cn_flat must be float64, C-contiguous, shape (natoms)")
    if coords_flat.size != 3 * natoms:
        raise ValueError("coords_flat length must be 3*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")

    LIB.get_cn_fortran(int(natoms), types, coords_flat, cn_flat)
