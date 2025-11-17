from __future__ import annotations

import ctypes as ct

import numpy as np
from numpy.ctypeslib import ndpointer

from .._lib import LIB

# void get_irmsd_fortran(int n1, int* types1, double* coords1,
#                          int n2, int* types2, double* coords2,
#                          int iinversion, double rmsd, int* types_out, double* coords3)
LIB.get_irmsd_fortran.argtypes = [
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.c_int,
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ct.c_int,
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
LIB.get_irmsd_fortran.restype = None


def get_irmsd_fortran_raw(
    n1: int,
    types1: np.ndarray,  # (n1,) int32 C
    coords1_flat: np.ndarray,  # (3*n1,) float64 C
    n2: int,
    types2: np.ndarray,  # (n2,) int32 C
    coords2_flat: np.ndarray,  # (3*n2,) float64 C
    iinversion: int,
    types_out1: np.ndarray,  # (n2,) int32 C
    coords_out1_flat: np.ndarray,  # (3*n2,) float64 C
    types_out2: np.ndarray,  # (n2,) int32 C
    coords_out2_flat: np.ndarray,  # (3*n2,) float64 C
) -> float:
    # Validate buffers to catch ABI mismatches early
    if types1.dtype != np.int32 or not types1.flags.c_contiguous:
        raise TypeError("types1 must be int32 and C-contiguous")
    if types2.dtype != np.int32 or not types2.flags.c_contiguous:
        raise TypeError("types2 must be int32 and C-contiguous")
    if types_out1.dtype != np.int32 or not types_out1.flags.c_contiguous:
        raise TypeError("types_out must be int32 and C-contiguous")

    if coords1_flat.dtype != np.float64 or not coords1_flat.flags.c_contiguous:
        raise TypeError("coords1_flat must be float64 and C-contiguous")
    if coords2_flat.dtype != np.float64 or not coords2_flat.flags.c_contiguous:
        raise TypeError("coords2_flat must be float64 and C-contiguous")
    if coords_out1_flat.dtype != np.float64 or not coords_out1_flat.flags.c_contiguous:
        raise TypeError("coords_out_flat must be float64 and C-contiguous")

    if types1.size != n1:
        raise ValueError("types1 length must be n1")
    if coords1_flat.size != 3 * n1:
        raise ValueError("coords1_flat length must be 3*n1")
    if types2.size != n2:
        raise ValueError("types2 length must be n2")
    if types_out1.size != types2.size:
        raise ValueError("types_out length must be equal to types2 length")
    if coords2_flat.size != 3 * n2:
        raise ValueError("coords2_flat length must be 3*n2")
    if coords_out1_flat.size != coords2_flat.size:
        raise ValueError("coords_out_flat length must be equal to coords2_flat length")
    rmsd = np.zeros(1, dtype=np.float64)

    LIB.get_irmsd_fortran(
        int(n1),
        types1,
        coords1_flat,
        int(n2),
        types2,
        coords2_flat,
        int(iinversion),
        rmsd,
        types_out1,
        coords_out1_flat,
        types_out2,
        coords_out2_flat,
    )
    return rmsd[0]
