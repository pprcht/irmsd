from __future__ import annotations

from typing import Tuple

import numpy as np

from ..bindings import canonical_exposed as _F


def get_canonical_fortran(
    atom_numbers: np.ndarray,
    positions: np.ndarray,
    wbo: np.ndarray | None = None,
    invtype: str = "apsp+",
    heavy: bool = False,
) -> np.ndarray:
    """
    Core API: call the Fortran routine to calculate the canonical ranking of atoms

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.
    heavy : bool, optional
        Whether to consider only heavy atoms (default: False).
    wbo: (natoms, natoms) float64, C-contiguous, optional
        Optional Wiberg bond order matrix, required if invtype is 'cangen', ignored in case of 'apsp+'.
    invtype : str, optional
        alogrithm type for invariants calculation (default: apsp+), alternativly 'cangen'.

    Returns
    -------
    rank : (N,) int32
        Rank array.

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

    rank = np.ascontiguousarray(np.zeros(n), dtype=np.int32)

    _F.get_canonical_sorter_fortran_raw(
        n,
        atom_numbers,
        coords_flat,
        rank,
        heavy=heavy,
        wbo=wbo,
        invtype=invtype,
    )

    return rank


def get_canonical_from_connect_fortran(
    atom_numbers: np.ndarray,
    connectivity: np.ndarray,
    heavy: bool = False,
) -> np.ndarray:
    """
    Core API: call the Fortran routine to calculate the canonical ranking of atoms from connectivity matrix

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    connectivity : (N, N) int32-like
        Connectivity matrix
    heavy : bool, optional
        Whether to consider only heavy atoms (default: False).

    Returns
    -------
    rank : (N,) int32
        Rank array.

    Raises
    ------
    ValueError
        If connectivity does not have shape (N, N).
    """
    atom_numbers = np.ascontiguousarray(atom_numbers, dtype=np.int32)
    connect = np.ascontiguousarray(connectivity, dtype=np.int32)

    n = int(connect.shape[0])
    if connect.ndim != 2 or connect.shape[1] != n:
        raise ValueError("Connectivity must have shape (N, N)")

    connect_flat = connect.reshape(-1).copy(order="F")

    rank = np.ascontiguousarray(np.zeros(n), dtype=np.int32)

    _F.get_ids_from_connect_fortran_raw(
        n,
        atom_numbers,
        connect_flat,
        rank,
        heavy=heavy,
    )

    return rank
