from __future__ import annotations
from typing import Tuple
import numpy as np
from ..bindings import xyz_bridge as _F

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

