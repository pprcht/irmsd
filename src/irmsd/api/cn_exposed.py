from __future__ import annotations

import numpy as np

from ..bindings import cn_exposed as _F


def get_cn_fortran(atom_numbers: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """
    Core API: call the Fortran routine to calculate CN

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.

    Returns
    -------
    cn : (N) float64 ndarray
        array with coordination numbers

    Raises
    ------
    ValueError
        If positions does not have shape (N, 3).
    """
    atom_numbers = np.ascontiguousarray(atom_numbers, dtype=np.int32)
    pos = np.ascontiguousarray(positions, dtype=np.float64)

    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = int(pos.shape[0])

    coords_flat = pos.reshape(-1).copy(order="C")
    cn = np.zeros(n, dtype=np.float64, order="C")

    _F.get_cn_fortran_raw(n, atom_numbers, coords_flat, cn)

    return cn
