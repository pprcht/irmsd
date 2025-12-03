from __future__ import annotations

import ctypes as ct

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

# We expose the axis as:
#   get_axis_fortran(natoms: c_int,
#                  types: int32[C_CONTIGUOUS](natoms),
#                  coords: float64[C_CONTIGUOUS](3*natoms),
#                  rot: float64[F_CONTIGUOUS](3))
#                  avmom: c_double)
#                  evec: float64[F_CONTIGUOUS](3,3))
LIB.get_axis0_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
LIB.get_axis0_fortran.restype = None


def get_axis_fortran_raw(
    natoms: int,
    types: np.ndarray,
    coords_flat: np.ndarray,
    rot: np.ndarray,
    avmom: np.ndarray,
    evec: np.ndarray,
) -> None:
    """Low-level call that matches the Fortran signature exactly. Operates IN-
    PLACE on rot, avmom, evec.

    Parameters
    ----------
    natoms : int
        Number of atoms (must be consistent with array lengths).
    types : (natoms,) int32, C-contiguous
        Atomic numbers (or type IDs).
    coords_flat : (3*natoms,) float64, C-contiguous
        Flat coordinates [x1,y1,z1,x2,y2,z2,...].
    rot : (3,) float64, C-contiguous
        Output rotation axis.
    avmom : (1,) float64, C-contiguous
        Output average moment.
    evec : (3,3) float64, C-contiguous
        Output eigenvectors.

    Raises
    ------
    TypeError
        If any of the arrays do not have the expected dtype or memory layout.
    """
    # light validation to catch mismatches early
    if types.dtype != np.int32 or not types.flags.c_contiguous:
        raise TypeError("types must be int32 and C-contiguous")
    if coords_flat.dtype != np.float64 or not coords_flat.flags.c_contiguous:
        raise TypeError("coords_flat must be float64 and C-contiguous")
    if rot.dtype != np.float64 or not rot.flags.c_contiguous or rot.size != 3:
        raise TypeError("rot must be float64, C-contiguous, shape (3)")
    if evec.dtype != np.float64 or not evec.flags.c_contiguous or evec.shape != (3, 3):
        raise TypeError("evec must be float64, C-contiguous, shape (3, 3)")
    if avmom.dtype != np.float64 or not evec.flags.c_contiguous or avmom.size != 1:
        raise TypeError("avmom must be float64, C-contiguous, size 1")
    if coords_flat.size != 3 * natoms:
        raise ValueError("coords_flat length must be 3*natoms")
    if types.size != natoms:
        raise ValueError("types length must be natoms")

    LIB.get_axis0_fortran(int(natoms), types, coords_flat, rot, avmom, evec)
