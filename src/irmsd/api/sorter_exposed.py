from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

from ..bindings import sorter_exposed as _F


def sorter_irmsd(
    atom_numbers_list: Sequence[np.ndarray],
    positions_list: Sequence[np.ndarray],
    nat: int,
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    energies_list: Sequence[np.ndarray] | None = None,
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
    rthr : float
        Distance threshold for the Fortran sorter. In Angström.
    iinversion : int
        Inversion symmetry flag.
    allcanon : bool
        Canonicalization flag.
    printlvl : int
        Verbosity level.
    ethr: float | None
        Inter-conformer energy threshold (optional). In Hartree.
    energies_list: sequence of (Nall,) floats | None
        List of energies for the passed structures (optional). In Hartree.

    Returns
    -------
    groups : (nat,) int32
        Group index for each of the first `nat` atoms.
    xyz_structs : list of (N,3) float64 arrays
        Updated coordinates for each structure.
    Z_structs : list of (N,) int32 arrays
        Updated atom numbers for each structure.

    Raises
    ------
    ValueError
        If input arrays have inconsistent shapes or invalid parameters.
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
        raise ValueError(
            "Each atom_numbers array must be shape (N,) matching positions"
        )

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
    atall = np.stack(at_list, axis=0)  # (nall, N)
    xyzall = np.stack(xyz_list, axis=0)  # (nall, N, 3)

    # Handling of the optional energy threshold (can't be None for Fortran/C-bindings)
    if ethr is None:
        ethr = -99999

    # Similar for energies, we can't pass None to the Fortran/C-bindings
    if energies_list is None:
        energies = np.zeros(nall, dtype=np.float64)
    else:
        # user provided a list/sequence of scalars → make a contiguous 1D array
        energies = np.asarray(energies_list, dtype=np.float64)
        if energies.ndim != 1:
            raise ValueError("energies_list must be 1D (one energy per structure)")

    # Allocate groups
    groups = np.empty(nall, dtype=np.int32)

    # ---- Raw Fortran call ----
    _F.sorter_exposed_xyz_fortran_raw(
        int(nat),
        int(nall),
        xyzall,  # flattened buffer (nall*N*3)
        atall,  # flattened buffer (nall*N)
        groups,  # length nat
        float(rthr),
        int(iinversion),
        bool(allcanon),
        int(printlvl),
        float(ethr),
        energies,
    )

    # ---- Extract back into per-structure arrays ----
    xyz_structs = [xyzall[i, :, :].copy() for i in range(nall)]
    Z_structs = [atall[i, :].copy() for i in range(nall)]

    return groups, xyz_structs, Z_structs


def delta_irmsd_list(
    atom_numbers_list: Sequence[np.ndarray],
    positions_list: Sequence[np.ndarray],
    nat: int,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    """
    High-level API: call the delta_irmsd_list_fortran Fortran routine.

    Parameters
    ----------
    atom_numbers_list : sequence of (N,) int32 arrays
        Per-structure atom numbers.
    positions_list : sequence of (N,3) float64 arrays
        Per-structure coordinates.
    nat : int
        Number of atoms for which the groups array is defined.
        Must satisfy 1 <= nat <= N.
    iinversion : int
        Inversion symmetry flag.
    allcanon : bool
        Canonicalization flag.
    printlvl : int
        Verbosity level.

    Returns
    -------
    delta : (nall,) float64
        iRMSD value between structure x_i and x_i-1
    xyz_structs : list of (N,3) float64 arrays
        Updated coordinates for each structure.
    Z_structs : list of (N,) int32 arrays
        Updated atom numbers for each structure.

    Raises
    ------
    ValueError
        If input arrays have inconsistent shapes or invalid parameters.
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
        raise ValueError(
            "Each atom_numbers array must be shape (N,) matching positions"
        )

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
    atall = np.stack(at_list, axis=0)  # (nall, N)
    xyzall = np.stack(xyz_list, axis=0)  # (nall, N, 3)

    # Allocate delta iRMSD storage
    delta = np.empty(nall, dtype=np.float64)

    # ---- Raw Fortran call ----
    _F.delta_irmsd_list_fortran_raw(
        int(nat),
        int(nall),
        xyzall,  # flattened buffer (nall*N*3)
        atall,  # flattened buffer (nall*N)
        int(iinversion),
        delta,  # length nall
        bool(allcanon),
        int(printlvl),
    )

    # ---- Extract back into per-structure arrays ----
    xyz_structs = [xyzall[i, :, :].copy() for i in range(nall)]
    Z_structs = [atall[i, :].copy() for i in range(nall)]

    return delta, xyz_structs, Z_structs


def cregen_raw(
    atom_numbers_list: Sequence[np.ndarray],
    positions_list: Sequence[np.ndarray],
    energies_list: Sequence[np.ndarray], 
    nat: int,
    rthr: float = 0.125,
    ethr: float = 7.96800686e-5, # == 0.05 kcal/mol
    bthr: float = 0.01,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    """
    High-level API: call the cregen_exposed_fortran Fortran routine.
    Note, this routine assumes (and checks) if all structures have
    the same sequence of atomic numbers

    Parameters
    ----------
    atom_numbers_list : sequence of (N,) int32 arrays
        Per-structure atom numbers.
    positions_list : sequence of (N,3) float64 arrays
        Per-structure coordinates.
    energies_list: sequence of (Nall,) float64
        List of energies for the passed structures (optional). In Hartree. 
    nat : int
        Number of atoms for which the groups array is defined.
        Must satisfy 1 <= nat <= N.
    rthr : float
        Distance threshold for the Fortran sorter. In Angström.
    ethr: float
        Inter-conformer energy threshold. In Hartree. 
    bthr: float
        Inter-conformer rotational constant threshold (relative)
    printlvl : int
        Verbosity level.

    Returns
    -------
    groups : (nat,) int32
        Group index for each of the first `nat` atoms.
    xyz_structs : list of (N,3) float64 arrays
        Updated coordinates for each structure.
    energies_list:
        Updated energies_list

    Raises
    ------
    ValueError
        If input arrays have inconsistent shapes or invalid parameters.
    """
    # Check basic structure count
    if len(atom_numbers_list) == 0:
        raise ValueError("atom_numbers_list must contain at least one structure")
    if len(atom_numbers_list) != len(positions_list):
        raise ValueError("atom_numbers_list and positions_list must have same length")
    if len(positions_list) != len(energies_list):
        raise ValueError("positions_list and energies_list must have same length") 

    ref = atom_numbers_list[0]
    for i, arr in enumerate(atom_numbers_list[1:], start=1):
        if not np.array_equal(arr, ref):
            raise ValueError(f"atom_numbers_list[{i}] differs from the first array")

    nall = len(atom_numbers_list)

    # Validate shapes using the first structure
    Z0 = np.ascontiguousarray(atom_numbers_list[0], dtype=np.int32)
    P0 = np.ascontiguousarray(positions_list[0], dtype=np.float64)

    if P0.ndim != 2 or P0.shape[1] != 3:
        raise ValueError("Each positions array must have shape (N,3)")
    if Z0.ndim != 1 or Z0.shape[0] != P0.shape[0]:
        raise ValueError(
            "Each atom_numbers array must be shape (N,) matching positions"
        )

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
    atall = np.stack(at_list, axis=0)  # (nall, N)
    xyzall = np.stack(xyz_list, axis=0)  # (nall, N, 3)
    energies = np.ascontiguousarray(energies_list, dtype=np.float64)

    # Allocate groups
    groups = np.empty(nall, dtype=np.int32)

    # ---- Raw Fortran call ----
    _F.cregen_fortran_raw(
        int(nat),
        int(nall),
        xyzall,  # flattened buffer (nall*N*3)
        atall,  # flattened buffer (nall*N)
        groups,  # length nat
        float(rthr),
        float(ethr),
        float(bthr),
        int(printlvl),
        energies,
    )

    # ---- Extract back into per-structure arrays ----
    xyz_structs = [xyzall[i, :, :].copy() for i in range(nall)]
    energies_list = [energies[i].copy() for i in range(nall)]

    return groups, xyz_structs, energies_list
