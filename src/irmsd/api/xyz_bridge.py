from __future__ import annotations
from typing import Tuple
import numpy as np
from ..bindings import xyz_bridge as _F
from ..bindings.xyz_bridge import xyz_to_fortran_pair_raw 

def xyz_to_fortran(atom_numbers: np.ndarray, positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Core API: call the Fortran routine using plain arrays (no ASE dependency).

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.

    Returns
    -------
    new_positions : (N, 3) float64 ndarray
        Positions after the Fortran subroutine overwrote them.
    mat : (3, 3) float64 ndarray (Fortran-ordered)
        3×3 matrix produced by the Fortran routine.
    """
    atom_numbers = np.ascontiguousarray(atom_numbers, dtype=np.int32)
    pos = np.ascontiguousarray(positions, dtype=np.float64)

    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = int(pos.shape[0])

    coords_flat = pos.reshape(-1).copy(order="C")
    mat = np.zeros((3, 3), dtype=np.float64, order="F")

    _F.xyz_to_fortran_raw(n, atom_numbers, coords_flat, mat)

    new_positions = coords_flat.reshape(n, 3)
    return new_positions, mat

###############################################################################

def xyz_to_fortran_pair(
    atom_numbers1: np.ndarray, positions1: np.ndarray,
    atom_numbers2: np.ndarray, positions2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Pair API: call the Fortran routine on TWO structures.

    Parameters
    ----------
    atom_numbers1 : (N1,) int32-like
    positions1    : (N1, 3) float64-like
    atom_numbers2 : (N2,) int32-like
    positions2    : (N2, 3) float64-like

    Returns
    -------
    new_positions1 : (N1, 3) float64
    new_positions2 : (N2, 3) float64
    mat1           : (3, 3) float64 (Fortran-ordered)
    mat2           : (3, 3) float64 (Fortran-ordered)
    """
    Z1 = np.ascontiguousarray(atom_numbers1, dtype=np.int32)
    Z2 = np.ascontiguousarray(atom_numbers2, dtype=np.int32)

    P1 = np.ascontiguousarray(positions1, dtype=np.float64)
    P2 = np.ascontiguousarray(positions2, dtype=np.float64)

    if P1.ndim != 2 or P1.shape[1] != 3:
        raise ValueError("positions1 must have shape (N1, 3)")
    if P2.ndim != 2 or P2.shape[1] != 3:
        raise ValueError("positions2 must have shape (N2, 3)")

    n1 = int(P1.shape[0]); n2 = int(P2.shape[0])

    c1 = P1.reshape(-1).copy(order="C")
    c2 = P2.reshape(-1).copy(order="C")
    M1 = np.zeros((3, 3), dtype=np.float64, order="F")
    M2 = np.zeros((3, 3), dtype=np.float64, order="F")

    xyz_to_fortran_pair_raw(n1, Z1, c1, n2, Z2, c2, M1, M2)

    new_P1 = c1.reshape(n1, 3)
    new_P2 = c2.reshape(n2, 3)
    return new_P1, new_P2, M1, M2

