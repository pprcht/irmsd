# python/bindings/xyz_bridge.py
from __future__ import annotations
from numpy.ctypeslib import ndpointer
import ctypes as ct
import numpy as np
from .._lib import LIB

# Fortran: subroutine xyz_to_fortran(natoms, types_ptr, coords_ptr, mat_ptr) bind(C)
# We expose it as:
#   xyz_to_fortran(natoms: c_int,
#                  types: int32[C_CONTIGUOUS](natoms),
#                  coords: float64[C_CONTIGUOUS](3*natoms),
#                  mat: float64[F_CONTIGUOUS](3,3))
LIB.xyz_to_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32,   flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
]
LIB.xyz_to_fortran.restype = None


def xyz_to_fortran_raw(
    natoms: int,
    types: np.ndarray,
    coords_flat: np.ndarray,
    mat_3x3_F: np.ndarray,
) -> None:
    """
    Low-level call that matches the Fortran signature exactly.
    Operates IN-PLACE on coords_flat and mat_3x3_F.

    Parameters
    ----------
    natoms : int
        Number of atoms (must be consistent with array lengths).
    types : (natoms,) int32, C-contiguous
        Atomic numbers (or type IDs).
    coords_flat : (3*natoms,) float64, C-contiguous
        Flat coordinates [x1,y1,z1,x2,y2,z2,...].
    mat_3x3_F : (3,3) float64, Fortran-contiguous
        Output matrix written by Fortran (column-major).
    """
    # light validation to catch mismatches early
    if types.dtype != np.int32 or not types.flags.c_contiguous:
        raise TypeError("types must be int32 and C-contiguous")
    if coords_flat.dtype != np.float64 or not coords_flat.flags.c_contiguous:
        raise TypeError("coords_flat must be float64 and C-contiguous")
    if mat_3x3_F.dtype != np.float64 or not mat_3x3_F.flags.f_contiguous or mat_3x3_F.shape != (3, 3):
        raise TypeError("mat_3x3_F must be float64, Fortran-contiguous, shape (3,3)")
    if coords_flat.size != 3 * natoms:
        raise ValueError("coords_flat length must be 3*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")

    LIB.xyz_to_fortran(int(natoms), types, coords_flat, mat_3x3_F)

