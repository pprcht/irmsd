from __future__ import annotations

from typing import Sequence, Tuple, List
import numpy as np

from ..bindings import sorter_exposed as _F


def sorter_irmsd(
    atom_numbers_list: Sequence[np.ndarray],
    positions_list: Sequence[np.ndarray],
    nat: int,
    rthresh: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    """
    High-level API: call the sorter_exposed_xyz_fortran Fortran routine.

    Parameters
    ----------
    atom_numbers_list : sequence of (N,) int32 arrays
        Per-structure atom numbers.
    positions_list : sequence of (N,3) float64 arrays
        Per-structure coordinates.
    nat : int
        Number of atoms for which the groups array is defined.
        Must satisfy 1 <= nat <= N.
    rthresh : float
        Distance threshold for the Fortran sorter.
    iinversion : int
        Inversion symmetry flag.
    allcanon : bool
        Canonicalization flag.
    printlvl : int
        Verbosity level.

    Returns
    -------
    groups : (nat,) int32
        Group index for each of the first `nat` atoms.
    xyz_structs : list of (N,3) float64 arrays
        Updated coordinates for each structure.
    Z_structs : list of (N,) int32 arrays
        Updated atom numbers for each structure.
    """
    # Check basic structure count
    if len(atom_numbers_list) == 0:
        raise ValueError("atom_numbers_list must contain at least one structure")
    if len(atom_numbers_list) != len(positions_list):
        raise ValueError("atom_numbers_list and positions_list must have same length")

    nall = len(atom_numbers_list)

    # Validate shapes using the first structure
    Z0 = np.ascontiguousarray(atom_numbers_list[0], dtype=np.int32)
    P0 = np.ascontiguousarray(positions_list[0], dtype=np.float64)

    if P0.ndim != 2 or P0.shape[1] != 3:
        raise ValueError("Each positions array must have shape (N,3)")
    if Z0.ndim != 1 or Z0.shape[0] != P0.shape[0]:
        raise ValueError("Each atom_numbers array must be shape (N,) matching positions")

    N = int(Z0.shape[0])
    if not (1 <= nat <= N):
        raise ValueError(f"nat must satisfy 1 ≤ nat ≤ N; got nat={nat}, N={N}")

    # Normalize all structures and enforce equal lengths
    at_list = []
    xyz_list = []
    for i, (Zi, Pi) in enumerate(zip(atom_numbers_list, positions_list)):
        Zi = np.ascontiguousarray(Zi, dtype=np.int32)
        Pi = np.ascontiguousarray(Pi, dtype=np.float64)

        if Pi.shape != (N, 3):
            raise ValueError(f"positions_list[{i}] must have shape (N,3)")
        if Zi.shape != (N,):
            raise ValueError(f"atom_numbers_list[{i}] must have shape (N,)")

        at_list.append(Zi)
        xyz_list.append(Pi)

    # Pack into contiguous meta-arrays
    atall = np.stack(at_list, axis=0)     # (nall, N)
    xyzall = np.stack(xyz_list, axis=0)   # (nall, N, 3)

    # Allocate groups
    groups = np.empty(nat, dtype=np.int32)

    # ---- Raw Fortran call ----
    _F.sorter_exposed_xyz_fortran_raw(
        int(nat),
        int(nall),
        xyzall,      # flattened buffer (nall*N*3)
        atall,       # flattened buffer (nall*N)
        groups,      # length nat
        float(rthresh),
        int(iinversion),
        bool(allcanon),
        int(printlvl),
    )

    # ---- Extract back into per-structure arrays ----
    xyz_structs = [xyzall[i, :, :].copy() for i in range(nall)]
    Z_structs   = [atall[i, :].copy() for i in range(nall)]

    return groups, xyz_structs, Z_structs

