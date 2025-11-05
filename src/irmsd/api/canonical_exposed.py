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
    invtype: str = "cangen",
    heavy: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
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
    TODO
    """
    atom_numbers = np.ascontiguousarray(atom_numbers, dtype=np.int32)
    pos = np.ascontiguousarray(positions, dtype=np.float64)

    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = int(pos.shape[0])

    coords_flat = pos.reshape(-1).copy(order="C")

    # FIXME: get coorect wbo somehow not just identity
    wbo = np.ascontiguousarray(np.identity(n), dtype=np.float64)
    rank = np.ascontiguousarray(np.zeros(n), dtype=np.int32)
    invariants = np.ascontiguousarray(np.zeros(n), dtype=np.int32)

    _F.get_canonical_sorter_fortran_raw(
        n, atom_numbers, coords_flat, wbo, invtype, heavy, rank, invariants
    )

    return rank, invariants
