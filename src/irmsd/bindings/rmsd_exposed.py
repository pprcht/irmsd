from __future__ import annotations

import ctypes as ct

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

# void get_quaternion_rmsd_fortran(int n1, int* types1, double* coords1,
#                          int n2, int* types2, double* coords2,
#                          double rmsd, double* Umat(3x3 F), bool* mask)
LIB.get_quaternion_rmsd_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),
    ct.POINTER(ct.c_bool),  # wbo can be None
]
LIB.get_quaternion_rmsd_fortran.restype = None


def get_quaternion_rmsd_fortran_raw(
    n1: int,
    types1: np.ndarray,  # (n1,) int32 C
    coords1_flat: np.ndarray,  # (3*n1,) float64 C
    n2: int,
    types2: np.ndarray,  # (n2,) int32 C
    coords2_flat: np.ndarray,  # (3*n2,) float64 C
    Umat_F: np.ndarray,  # (3,3) float64 F
    mask: np.ndarray | None = None,  # (bool,
) -> float:
    """Low-level call that matches the Fortran signature exactly.

    Parameters
    ----------
    n1 : int
        Number of atoms in structure 1 (must be consistent with array lengths).
    types1 : (n1,) int32, C-contiguous
        Atomic numbers for structure 1.
    coords1_flat : (3*n1,) float64, C-contiguous
        Flat coordinates [x1,y1,z1,x2,y2,z2,...] for structure 1.
    n2 : int
        Number of atoms in structure 2 (must be consistent with array lengths).
    types2 : (n2,) int32, C-contiguous
        Atomic numbers for structure 2.
    coords2_flat : (3*n2,) float64, C-contiguous
        Flat coordinates [x1,y1,z1,x2,y2,z2,...] for structure 2.
    Umat_F : (3,3) float64, Fortran-contiguous
        Output rotation matrix.
    mask : (n1,) bool, C-contiguous, optional
        Optional mask to include only a subset of atoms from structure 1.

    Returns
    -------
    float
        The computed RMSD.

    Raises
    ------
    TypeError
        If any of the arrays do not have the expected dtype or memory layout.
    ValueError
        If array sizes do not match n1 or n2.
    """
    # Validate buffers to catch ABI mismatches early
    if types1.dtype != np.int32 or not types1.flags.c_contiguous:
        raise TypeError("types1 must be int32 and C-contiguous")
    if types2.dtype != np.int32 or not types2.flags.c_contiguous:
        raise TypeError("types2 must be int32 and C-contiguous")

    if coords1_flat.dtype != np.float64 or not coords1_flat.flags.c_contiguous:
        raise TypeError("coords1_flat must be float64 and C-contiguous")
    if coords2_flat.dtype != np.float64 or not coords2_flat.flags.c_contiguous:
        raise TypeError("coords2_flat must be float64 and C-contiguous")

    if (
        Umat_F.dtype != np.float64
        or Umat_F.shape != (3, 3)
        or not Umat_F.flags.f_contiguous
    ):
        raise TypeError("Umat_F must be float64, shape (3,3), Fortran-contiguous")

    if types1.size != n1:
        raise ValueError("types1 length must be n1")
    if coords1_flat.size != 3 * n1:
        raise ValueError("coords1_flat length must be 3*n1")
    if types2.size != n2:
        raise ValueError("types2 length must be n2")
    if coords2_flat.size != 3 * n2:
        raise ValueError("coords2_flat length must be 3*n2")
    rmsd = np.zeros(1, dtype=np.float64)

    if mask is not None:
        if mask.dtype != np.bool or not mask.flags.c_contiguous or mask.size != n1:
            raise TypeError("mask must be bool, C-contiguous, size natoms")
        mask_ptr = mask.ctypes.data_as(ct.POINTER(ct.c_bool))
    else:
        mask_ptr = None

    LIB.get_quaternion_rmsd_fortran(
        int(n1),
        types1,
        coords1_flat,
        int(n2),
        types2,
        coords2_flat,
        rmsd,
        Umat_F,
        mask_ptr,
    )
    return rmsd[0]
