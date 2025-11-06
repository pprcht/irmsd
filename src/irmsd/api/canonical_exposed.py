from __future__ import annotations

from typing import Tuple

import numpy as np

from ..bindings import canonical_exposed as _F


# TODO: disucss with phillip how to deal with the invtype string
#       the fortran implementation will not throw an error if a wrong
#       string is passed, it will just ignore it and use the default
#       should we check it here?
def get_canonical_fortran(
    atom_numbers: np.ndarray,
    positions: np.ndarray,
    wbo: np.ndarray | None = None,
    invtype: str = "apsp+",
    heavy: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Core API: call the Fortran routine to calculate CN

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.
    heavy : bool, optional
        Whether to consider only heavy atoms (default: False).
    wbo: (natoms, natoms) float64, C-contiguous, optional
        Optional Wiberg bond order matrix.
    invtype : str, optional
        alogrithm type for invariants calculation (default: None <=> apsp+), alternativly 'cangen'.

    Returns
    -------
    rank : (N,) int32
        Rank array.
    invariants : (N,) int32
        Invariants array.
    """
    atom_numbers = np.ascontiguousarray(atom_numbers, dtype=np.int32)
    pos = np.ascontiguousarray(positions, dtype=np.float64)

    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = int(pos.shape[0])

    coords_flat = pos.reshape(-1).copy(order="C")

    rank = np.ascontiguousarray(np.zeros(n), dtype=np.int32)
    invariants = np.ascontiguousarray(np.zeros(n), dtype=np.int32)

    _F.get_canonical_sorter_fortran_raw(
        n,
        atom_numbers,
        coords_flat,
        rank,
        invariants,
        heavy=heavy,
        wbo=wbo,
        invtype=invtype,
    )

    return rank, invariants
