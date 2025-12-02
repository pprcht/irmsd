from __future__ import annotations

from typing import Tuple

import numpy as np

from ..bindings import axis_exposed as _F


def get_axis(
    atom_numbers: np.ndarray, positions: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Core API: call the Fortran routine to calculate the rotation axis, average moment, and eigenvectors

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.

    Returns
    -------
    rot : (3,) float64 ndarray
        Rotation axis.
    avmom : (1,) float64 ndarray
        Average moment.
    evec : (3, 3) float64 ndarray
        Eigenvectors.

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

    rot = np.ascontiguousarray(np.zeros(3), dtype=np.float64)
    avmom = np.ascontiguousarray(np.zeros(1), dtype=np.float64)
    evec = np.ascontiguousarray(np.zeros((3, 3)), dtype=np.float64)

    _F.get_axis_fortran_raw(n, atom_numbers, coords_flat, rot, avmom, evec)

    return rot, avmom, evec
