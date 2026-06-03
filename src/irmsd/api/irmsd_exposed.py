from __future__ import annotations

from typing import Tuple

import numpy as np

from ..bindings import irmsd_exposed as _F


def get_irmsd(
    atom_numbers1: np.ndarray,
    positions1: np.ndarray,
    atom_numbers2: np.ndarray,
    positions2: np.ndarray,
    iinversion: int = 0,
    ranks1: np.ndarray | None = None,
    ranks2: np.ndarray | None = None,
) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Core API: call the Fortran routine to calculate the iRMSD between two structures

    Parameters
    ----------
    atom_numbers1 : (N1,) int32-like
        Atomic numbers (or types) of structure 1.
    positions1 : (N1, 3) float64-like
        Cartesian coordinates in Å of structure 1.
    atom_numbers2 : (N2,) int32-like
        Atomic numbers (or types) of structure 2.
    positions2 : (N2, 3) float64-like
        Cartesian coordinates in Å of structure 2.
    iinversion : int, optional
        Whether to consider inversion symmetry. Default is 0 (auto). Set to 1 to use inversion, set
        to 2 to disable inversion.
    ranks1 : (N1,) int-like, optional
        Externally supplied per-atom canonical ranks for structure 1 (e.g. read
        from a file). Used directly only if both ``ranks1`` and ``ranks2`` are
        given, contain no zeros, and pass the internal consistency check;
        otherwise the canonical ranks are recomputed. A zero entry is the
        "not provided" sentinel.
    ranks2 : (N2,) int-like, optional
        Externally supplied per-atom canonical ranks for structure 2. See
        ``ranks1``.

    Returns
    -------
    rmsdval : float
        The calculated iRMSD value.
    Z3 : (N1,) int32 ndarray
        Atomic numbers of the aligned structure 1.
    P3 : (N1, 3) float64 ndarray
        Aligned and centered coordinates of structure 1.
    Z4 : (N2,) int32 ndarray
        Atomic numbers of the aligned structure 2.
    P4 : (N2, 3) float64 ndarray
        Aligned and centered coordinates of structure 2.

    Notes
    -----
    The returned coordinates P3 and P4 are centered at the origin.

    Raises
    ------
    ValueError
        If positions1 or positions2 do not have shape (Ni, 3).
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

    # Externally supplied ranks: an all-zero array is the "not provided"
    # sentinel the Fortran side falls back on.
    if ranks1 is None:
        R1 = np.zeros(n1, dtype=np.int32)
    else:
        R1 = np.ascontiguousarray(ranks1, dtype=np.int32)
        if R1.shape != (n1,):
            raise ValueError("ranks1 must have shape (N1,)")
    if ranks2 is None:
        R2 = np.zeros(n2, dtype=np.int32)
    else:
        R2 = np.ascontiguousarray(ranks2, dtype=np.int32)
        if R2.shape != (n2,):
            raise ValueError("ranks2 must have shape (N2,)")

    c1 = P1.reshape(-1).copy(order="C")
    c2 = P2.reshape(-1).copy(order="C")
    Z3 = np.zeros_like(Z2)
    c3 = np.zeros_like(c2)

    Z4 = np.zeros_like(Z2)
    c4 = np.zeros_like(c2)

    rmsdval = _F.get_irmsd_fortran_raw(
        n1,
        Z1,
        c1,
        n2,
        Z2,
        c2,
        iinversion,
        Z3,
        c3,
        Z4,
        c4,
        ranks1=R1,
        ranks2=R2,
    )

    P3 = c3.reshape(n1, 3)
    P4 = c4.reshape(n2, 3)

    center3 = P3.mean(axis=0)
    center4 = P4.mean(axis=0)

    P3 -= center3
    P4 -= center4

    return rmsdval, Z3, P3, Z4, P4
