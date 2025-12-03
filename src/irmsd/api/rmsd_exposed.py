from __future__ import annotations

from typing import Tuple

import numpy as np

from ..bindings import rmsd_exposed as _F


def get_quaternion_rmsd_fortran(
    atom_numbers1: np.ndarray,
    positions1: np.ndarray,
    atom_numbers2: np.ndarray,
    positions2: np.ndarray,
    mask: np.ndarray | None = None,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Pair API: call the Fortran routine on TWO structures.

    Parameters
    ----------
    atom_numbers1 : (N1,) int32-like
    positions1    : (N1, 3) float64-like
    atom_numbers2 : (N2,) int32-like
    positions2    : (N2, 3) float64-like
    mask          : (N1,) bool-like or None

    Returns
    -------
    rmsdval        : float64
    new_positions2 : (N2, 3) float64, (positions2 @ Umat.T)
    Umat           : (3, 3) float64 (Fortran-ordered)

    Notes
    -----
    The returned new_positions2 is aligned onto positions1.
    1. If mask is provided, only the atoms where mask==True in structure 1 are used to compute the RMSD.
    2. The rotation matrix Umat is Fortran-ordered, i.e., to rotate positions2, do: new_positions2 = positions2 @ Umat.T
    3. The returned new_positions2 is also translated to have the same barycenter as positions1.

    Raises
    ------
    ValueError
        If the input arrays do not have the correct shapes or types.
    """
    Z1 = np.ascontiguousarray(atom_numbers1, dtype=np.int32)
    Z2 = np.ascontiguousarray(atom_numbers2, dtype=np.int32)

    P1 = np.ascontiguousarray(positions1, dtype=np.float64)
    P2 = np.ascontiguousarray(positions2, dtype=np.float64)

    if P1.ndim != 2 or P1.shape[1] != 3:
        raise ValueError("positions1 must have shape (N1, 3)")
    if P2.ndim != 2 or P2.shape[1] != 3:
        raise ValueError("positions2 must have shape (N2, 3)")

    n1 = int(P1.shape[0])
    n2 = int(P2.shape[0])

    if mask is not None:
        mask = np.ascontiguousarray(mask, dtype=bool)
        if mask.ndim != 1 or mask.shape[0] != n1:
            raise ValueError("mask must have shape (N1,)")

    c1 = P1.reshape(-1).copy(order="C")
    c2 = P2.reshape(-1).copy(order="C")
    M2 = np.zeros((3, 3), dtype=np.float64, order="F")

    rmsdval = _F.get_quaternion_rmsd_fortran_raw(n1, Z1, c1, n2, Z2, c2, M2, mask=mask)

    new_P2 = c2.reshape(n2, 3) @ M2.T

    bc1 = P1.mean(axis=0)  # barycenter of reference structure
    bc2 = new_P2.mean(axis=0)  # barycenter of rotated structure
    new_P2 = new_P2 + (bc1 - bc2)

    return rmsdval, new_P2, M2
