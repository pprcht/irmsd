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


# void xyz_to_fortran_pair(int n1, int* types1, double* coords1,
#                          int n2, int* types2, double* coords2,
#                          double* mat1(3x3 F), double* mat2(3x3 F))
LIB.xyz_to_fortran_pair.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32,   flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.c_int,
    ndpointer(dtype=np.int32,   flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
]
LIB.xyz_to_fortran_pair.restype = None


def xyz_to_fortran_pair_raw(
    n1: int,
    types1: np.ndarray,         # (n1,) int32 C
    coords1_flat: np.ndarray,   # (3*n1,) float64 C
    n2: int,
    types2: np.ndarray,         # (n2,) int32 C
    coords2_flat: np.ndarray,   # (3*n2,) float64 C
    mat1_F: np.ndarray,         # (3,3) float64 F
    mat2_F: np.ndarray,         # (3,3) float64 F
) -> None:
    # Validate buffers to catch ABI mismatches early
    if types1.dtype != np.int32 or not types1.flags.c_contiguous:
        raise TypeError("types1 must be int32 and C-contiguous")
    if types2.dtype != np.int32 or not types2.flags.c_contiguous:
        raise TypeError("types2 must be int32 and C-contiguous")

    if coords1_flat.dtype != np.float64 or not coords1_flat.flags.c_contiguous:
        raise TypeError("coords1_flat must be float64 and C-contiguous")
    if coords2_flat.dtype != np.float64 or not coords2_flat.flags.c_contiguous:
        raise TypeError("coords2_flat must be float64 and C-contiguous")

    if mat1_F.dtype != np.float64 or mat1_F.shape != (3,3) or not mat1_F.flags.f_contiguous:
        raise TypeError("mat1_F must be float64, shape (3,3), Fortran-contiguous")
    if mat2_F.dtype != np.float64 or mat2_F.shape != (3,3) or not mat2_F.flags.f_contiguous:
        raise TypeError("mat2_F must be float64, shape (3,3), Fortran-contiguous")

    if types1.size != n1:        raise ValueError("types1 length must be n1")
    if coords1_flat.size != 3*n1: raise ValueError("coords1_flat length must be 3*n1")
    if types2.size != n2:        raise ValueError("types2 length must be n2")
    if coords2_flat.size != 3*n2: raise ValueError("coords2_flat length must be 3*n2")

    LIB.xyz_to_fortran_pair(
        int(n1), types1, coords1_flat,
        int(n2), types2, coords2_flat,
        mat1_F, mat2_F
    )    
