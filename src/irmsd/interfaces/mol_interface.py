from __future__ import annotations

from typing import Sequence, Tuple, List

import numpy as np

from ..core.molecule import Molecule
from ..sorting import first_by_assignment


##########################################################################################
# Energy utilities
##########################################################################################


def get_energies_from_molecule_list(
    molecule_list: Sequence[Molecule],
) -> np.ndarray:
    """
    Collect potential energies from a sequence of Molecule objects.

    For each Molecule, this function calls ``get_potential_energy()`` and
    stores the result in a 1D NumPy array of dtype float. If the energy
    is not available (for example, if ``get_potential_energy()`` raises
    ``AttributeError`` or returns ``None``), the corresponding entry is
    set to 0.0.

    Parameters
    ----------
    molecule_list : Sequence[Molecule]
        Sequence of Molecule instances.

    Returns
    -------
    energies : np.ndarray
        Array of shape (n_structures,) containing one energy per Molecule.
    """
    energies: list[float] = []

    for mol in molecule_list:
        if not isinstance(mol, Molecule):
            raise TypeError(
                "get_energies_from_molecule_list expects a sequence of Molecule objects"
            )

        # Our Molecule.get_potential_energy() may raise if energy is unset;
        # treat missing energies as 0.0.
        try:
            e = mol.get_potential_energy()
        except AttributeError:
            e = None

        energies.append(0.0 if e is None else float(e))

    return np.asarray(energies, dtype=float)


def prune_by_energy_window(
    molecule_list: Sequence[Molecule],
    energy_window: float | None,
) -> List[Molecule]:
    """
    Prune a list of Molecule objects based on a relative energy window.

    This function identifies the lowest-energy structure in ``molecule_list``
    and retains only those molecules whose energies lie within
    ``energy_window`` Hartree above this minimum. Energies are obtained via
    ``get_energies_from_molecule_list`` and are assumed to be aligned with the
    order of ``molecule_list``.

    If ``energy_window`` is ``None``, the input list is returned unchanged.

    Parameters
    ----------
    molecule_list : list of Molecule
        The structures to be filtered.
    energy_window : float or None
        Allowed energy range above the minimum structure, in Hartree. If
        ``None``, no pruning is performed.

    Returns
    -------
    pruned : list of Molecule
        A new list containing only the structures within the specified
        energy window.

    Notes
    -----
    The reference energy is the minimum energy found in the provided list.
    No sorting of the returned molecules is performed; the original order
    (minus removed entries) is preserved.
    """
    if energy_window is None:
        return molecule_list

    energies = get_energies_from_molecule_list(molecule_list)
    e_min = float(np.min(energies))
    threshold = e_min + energy_window

    mask = energies <= threshold
    pruned = [mol for mol, keep in zip(molecule_list, mask) if keep]

    return pruned


##########################################################################################
# RMSD utilities
##########################################################################################


def get_rmsd_molecule(
    molecule1: Molecule,
    molecule2: Molecule,
    mask=None,
) -> Tuple[float, Molecule, np.ndarray]:
    """
    Compute the RMSD between two Molecule objects using the quaternion-based
    RMSD backend, and return an aligned copy of the second Molecule.

    Parameters
    ----------
    molecule1 : Molecule
        Reference structure.
    molecule2 : Molecule
        Structure to be rotated/translated onto ``molecule1``.
    mask : array-like of bool, optional
        Optional mask selecting which atoms of ``molecule1`` participate
        in the RMSD. Must be broadcastable / compatible with the Fortran
        backend's mask semantics.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation in Ångström.
    new_molecule2 : Molecule
        Copy of ``molecule2`` with its positions replaced by the aligned
        coordinates returned by the backend.
    rotation_matrix : np.ndarray
        3×3 rotation matrix applied to align ``molecule2`` onto ``molecule1``.

    Raises
    ------
    TypeError
        If either input is not a Molecule.
    """
    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(molecule1, Molecule) or not isinstance(molecule2, Molecule):
        raise TypeError("get_rmsd_molecule expects two irmsd.Molecule objects")

    Z1 = molecule1.get_atomic_numbers()  # (N1,)
    P1 = molecule1.get_positions()  # (N1, 3)
    Z2 = molecule2.get_atomic_numbers()  # (N2,)
    P2 = molecule2.get_positions()  # (N2, 3)

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    # Work on a copy; do not modify the original second Molecule in-place.
    new_molecule2 = molecule2.copy()
    new_molecule2.set_positions(new_P2)

    return rmsdval, new_molecule2, umat


##########################################################################################
# iRMSD utilities
##########################################################################################


def get_irmsd_molecule(
    molecule1: Molecule,
    molecule2: Molecule,
    iinversion: int = 0,
) -> Tuple[float, Molecule, Molecule]:
    """
    Compute the iRMSD between two Molecule objects using the iRMSD backend.

    The backend may reorder atoms and/or change atomic numbers according to
    its canonicalization / matching logic. This wrapper returns *copies* of
    both input Molecules with the updated atomic numbers and positions.

    Parameters
    ----------
    molecule1 : Molecule
        First input structure.
    molecule2 : Molecule
        Second input structure.
    iinversion : int, optional
        Inversion flag passed directly to the Fortran backend. See the
        backend documentation for allowed values and meanings.

    Returns
    -------
    irmsd : float
        iRMSD value in Ångström.
    new_molecule1 : Molecule
        Copy of ``molecule1`` with updated atomic numbers and positions.
    new_molecule2 : Molecule
        Copy of ``molecule2`` with updated atomic numbers and positions.

    Raises
    ------
    TypeError
        If either input is not a Molecule.
    """
    from ..api.irmsd_exposed import get_irmsd

    if not isinstance(molecule1, Molecule) or not isinstance(molecule2, Molecule):
        raise TypeError("get_irmsd_molecule expects two irmsd.Molecule objects")

    Z1 = molecule1.get_atomic_numbers()  # (N1,)
    P1 = molecule1.get_positions()  # (N1, 3)
    Z2 = molecule2.get_atomic_numbers()  # (N2,)
    P2 = molecule2.get_positions()  # (N2, 3)

    irmsdval, new_Z1, new_P1, new_Z2, new_P2 = get_irmsd(
        Z1, P1, Z2, P2, iinversion=iinversion
    )

    new_molecule1 = molecule1.copy()
    new_molecule1.set_atomic_numbers(new_Z1)
    new_molecule1.set_positions(new_P1)

    new_molecule2 = molecule2.copy()
    new_molecule2.set_atomic_numbers(new_Z2)
    new_molecule2.set_positions(new_P2)

    return irmsdval, new_molecule1, new_molecule2


def sorter_irmsd_molecule(
    molecule_list: Sequence[Molecule],
    rthr: float = 0.125,  # aligned with '--rthr' in src/irmsd/cli.py
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
) -> Tuple[np.ndarray, List[Molecule]]:
    """
    High-level wrapper around the Fortran-backed ``sorter_irmsd`` that
    operates directly on Molecule objects.

    Parameters
    ----------
    molecule_list : Sequence[Molecule]
        Sequence of Molecule objects. All molecules must have the same
        number of atoms.
    rthr : float
        Distance threshold for the sorter (passed through to the backend).
    iinversion : int, optional
        Inversion symmetry flag, passed through to the backend.
    allcanon : bool, optional
        Canonicalization flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.
    ethr : float | None
        Optional energy threshold to accelerate by pre-sorting
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    groups : np.ndarray
        Integer array of shape (nat,) with group indices for the first
        ``nat`` atoms (as defined by the backend).
    new_molecule_list : list[Molecule]
        New Molecule objects reconstructed from the sorted atomic numbers
        and positions returned by the backend. The list has the same length
        and ordering as ``molecule_list``.

    Raises
    ------
    TypeError
        If ``molecule_list`` does not contain Molecule instances.
    ValueError
        If ``molecule_list`` is empty or if the Molecules do not all have
        the same number of atoms.
    """
    from ..api.sorter_exposed import sorter_irmsd

    # --- Basic checks on molecule_list ---
    if not isinstance(molecule_list, (list, tuple)):
        raise TypeError(
            "sorter_irmsd_molecule expects a sequence (list/tuple) of Molecule objects"
        )

    if len(molecule_list) == 0:
        raise ValueError("molecule_list must contain at least one Molecule object")

    for i, mol in enumerate(molecule_list):
        if not isinstance(mol, Molecule):
            raise TypeError(
                "sorter_irmsd_molecule expects a sequence of Molecule objects; "
                f"item {i} has type {type(mol)}"
            )

    # --- Check that all Molecules have the same number of atoms and define nat ---
    nat = len(molecule_list[0])
    for i, mol in enumerate(molecule_list):
        if len(mol) != nat:
            raise ValueError(
                "All Molecule objects must have the same number of atoms; "
                f"item 0 has {nat} atoms, item {i} has {len(mol)} atoms"
            )

    # --- Remove structures too high in energy, if needed
    if printlvl > 0 and ewin is not None:
        print(f"EWIN energy window : {ewin} Ha ({ewin*627.5095:.2f} kcal/mol)")
    n_orig = len(molecule_list)
    molecule_list = prune_by_energy_window(molecule_list, ewin)
    n_new = len(molecule_list)
    n_diff = n_orig - n_new
    if printlvl > 0 and ewin is not None:
        print(f" --> removed {n_diff} of {n_orig} structures.\n")

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for mol in molecule_list:
        Z = np.asarray(mol.get_atomic_numbers(), dtype=np.int32)  # (nat,)
        P = np.asarray(mol.get_positions(), dtype=np.float64)  # (nat, 3)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Molecule positions array must have shape (nat, 3); "
                f"got {P.shape}"
            )

        atom_numbers_list.append(Z)
        positions_list.append(P)

    energies_list = get_energies_from_molecule_list(molecule_list)

    # --- Call the Fortran-backed sorter_irmsd ---
    groups, xyz_structs, Z_structs = sorter_irmsd(
        atom_numbers_list=atom_numbers_list,
        positions_list=positions_list,
        nat=nat,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
        ethr=ethr,
        energies_list=energies_list,
    )

    # --- Reconstruct new Molecule objects ---
    new_molecule_list: List[Molecule] = []
    for mol_orig, Z_new, P_new in zip(molecule_list, Z_structs, xyz_structs):
        # Start from a copy to preserve metadata (cell, pbc, info, energy, etc.)
        new_mol = mol_orig.copy()
        # Update atomic numbers and positions according to sorter output
        new_mol.set_atomic_numbers(Z_new)
        new_mol.set_positions(P_new)
        new_molecule_list.append(new_mol)

    return groups, new_molecule_list


def delta_irmsd_list_molecule(
    molecule_list: Sequence[Molecule],
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List[Molecule]]:
    """
    High-level wrapper around the Fortran-backed ``delta_irmsd_list`` that
    operates directly on Molecule objects.

    Parameters
    ----------
    molecule_list : Sequence[Molecule]
        Sequence of Molecule objects. All molecules must have the same
        number of atoms.
    iinversion : int, optional
        Inversion symmetry flag, passed through to the backend.
    allcanon : bool, optional
        Canonicalization flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.

    Returns
    -------
    delta : np.ndarray
        Float array returned by the backend (see ``delta_irmsd_list`` for
        detailed semantics).
    new_molecule_list : list[Molecule]
        New Molecule objects reconstructed from the atomic numbers and
        positions returned by the backend. The list has the same length
        and ordering as ``molecule_list``.

    Raises
    ------
    TypeError
        If ``molecule_list`` does not contain Molecule instances.
    ValueError
        If ``molecule_list`` is empty or if the Molecules do not all have
        the same number of atoms.
    """
    from ..api.sorter_exposed import delta_irmsd_list

    # --- Basic checks on molecule_list ---
    if not isinstance(molecule_list, (list, tuple)):
        raise TypeError(
            "delta_irmsd_list_molecule expects a sequence (list/tuple) of Molecule objects"
        )

    if len(molecule_list) == 0:
        raise ValueError("molecule_list must contain at least one Molecule object")

    for i, mol in enumerate(molecule_list):
        if not isinstance(mol, Molecule):
            raise TypeError(
                "delta_irmsd_list_molecule expects a sequence of Molecule objects; "
                f"item {i} has type {type(mol)}"
            )

    # --- Check that all Molecules have the same number of atoms and define nat ---
    nat = len(molecule_list[0])
    for i, mol in enumerate(molecule_list):
        if len(mol) != nat:
            raise ValueError(
                "All Molecule objects must have the same number of atoms; "
                f"item 0 has {nat} atoms, item {i} has {len(mol)} atoms"
            )

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for mol in molecule_list:
        Z = np.asarray(mol.get_atomic_numbers(), dtype=np.int32)
        P = np.asarray(mol.get_positions(), dtype=np.float64)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Molecule positions array must have shape (nat, 3); "
                f"got {P.shape}"
            )

        atom_numbers_list.append(Z)
        positions_list.append(P)

    # --- Call the Fortran-backed delta_irmsd_list ---
    delta, xyz_structs, Z_structs = delta_irmsd_list(
        atom_numbers_list=atom_numbers_list,
        positions_list=positions_list,
        nat=nat,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    # --- Reconstruct new Molecule objects ---
    new_molecule_list: List[Molecule] = []
    for mol_orig, Z_new, P_new in zip(molecule_list, Z_structs, xyz_structs):
        new_mol = mol_orig.copy()
        new_mol.set_atomic_numbers(Z_new)
        new_mol.set_positions(P_new)
        new_molecule_list.append(new_mol)

    return delta, new_molecule_list


def cregen(
    molecule_list: Sequence[Molecule],
    rthr: float = 0.125,
    ethr: float = 7.96800686e-5,  # == 0.05 kcal/mol
    bthr: float = 0.01,
    printlvl: int = 0,
    ewin: float | None = None,
) -> List[Molecule]:
    """
    High-level wrapper around the Fortran-backed ``cregen_raw`` that
    operates directly on Molecule objects.
    Returns a pruned & energy-sorted list of structures.

    Parameters
    ----------
    molecule_list : Sequence[Molecule]
        Sequence of Molecule objects. All molecules must have the same
        number of atoms.
    rthr : float
        Distance threshold for the sorter (passed through to the backend).
    ethr : float
        Inter-conformer energy threshold (in Hartree)
    bthr : float
        Inter-conformer rotational constant threshold (fractional)
    iinversion : int, optional
        Inversion symmetry flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.


    Returns
    -------
    new_molecule_list : list[Molecule]
        New Molecule objects reconstructed from the sorted atomic numbers
        and positions returned by the backend. The list contains only n_structures
        defined as ``unique`` according to the selected thresholds.

    Raises
    ------
    TypeError
        If ``molecule_list`` does not contain Molecule instances.
    ValueError
        If ``molecule_list`` is empty or if the Molecules do not all have
        the same number of atoms.
    """
    from ..api.sorter_exposed import cregen_raw

    # --- Basic checks on molecule_list ---
    if not isinstance(molecule_list, (list, tuple)):
        raise TypeError(
            "sorter_irmsd_molecule expects a sequence (list/tuple) of Molecule objects"
        )

    if len(molecule_list) == 0:
        raise ValueError("molecule_list must contain at least one Molecule object")

    for i, mol in enumerate(molecule_list):
        if not isinstance(mol, Molecule):
            raise TypeError(
                "sorter_irmsd_molecule expects a sequence of Molecule objects; "
                f"item {i} has type {type(mol)}"
            )

    # --- Check that all Molecules have the same number of atoms and define nat ---
    nat = len(molecule_list[0])
    for i, mol in enumerate(molecule_list):
        if len(mol) != nat:
            raise ValueError(
                "All Molecule objects must have the same number of atoms; "
                f"item 0 has {nat} atoms, item {i} has {len(mol)} atoms"
            )

    # --- Check same order of atomic numbers
    ref = molecule_list[0].get_atomic_numbers()
    for i, mol in enumerate(molecule_list[1:], start=1):
        arr = mol.get_atomic_numbers()
        if not np.array_equal(arr, ref):
            raise ValueError(
                f"Molecule {i} has different atomic numbers than molecule 0"
            )

    # --- Remove structures too high in energy, if needed
    if printlvl > 0 and ewin is not None:
        print(f"EWIN energy window : {ewin} Ha ({ewin*627.5095:.2f} kcal/mol)")
    n_orig = len(molecule_list)
    molecule_list = prune_by_energy_window(molecule_list, ewin)
    n_new = len(molecule_list)
    n_diff = n_orig - n_new
    if printlvl > 0 and ewin is not None:
        print(f" --> removed {n_diff} of {n_orig} structures.\n")

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for mol in molecule_list:
        Z = np.asarray(mol.get_atomic_numbers(), dtype=np.int32)  # (nat,)
        P = np.asarray(mol.get_positions(), dtype=np.float64)  # (nat, 3)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Molecule positions array must have shape (nat, 3); "
                f"got {P.shape}"
            )

        atom_numbers_list.append(Z)
        positions_list.append(P)

    energies_list = get_energies_from_molecule_list(molecule_list)

    # --- Call the Fortran-backed sorter_irmsd ---
    groups, xyz_structs, energies_list = cregen_raw(
        atom_numbers_list=atom_numbers_list,
        positions_list=positions_list,
        energies_list=energies_list,
        nat=nat,
        rthr=rthr,
        ethr=ethr,
        bthr=bthr,
        printlvl=printlvl,
    )

    # --- Reconstruct new Molecule objects ---
    new_molecule_list: List[Molecule] = []
    for mol_orig, P_new, E_new in zip(molecule_list, xyz_structs, energies_list):
        # Start from a copy to preserve metadata (cell, pbc, info, energy, etc.)
        new_mol = mol_orig.copy()
        # Update atomic numbers and positions according to sorter output
        # new_mol.set_atomic_numbers(Z_new)
        new_mol.set_positions(P_new)
        new_mol.set_potential_energy(E_new)
        new_molecule_list.append(new_mol)

    new_molecule_list = first_by_assignment(new_molecule_list, groups)

    return new_molecule_list


def prune(
    molecule_list: Sequence[Molecule],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
) -> List[Molecule]:
    """
    High-level wrapper around the Fortran-backed ``sorter_irmsd`` that
    operates directly on Molecule objects. Returns a pruned list of structures.

    Parameters
    ----------
    molecule_list : Sequence[Molecule]
        Sequence of Molecule objects. All molecules must have the same
        number of atoms.
    rthr : float
        Distance threshold for the sorter (passed through to the backend).
    iinversion : int, optional
        Inversion symmetry flag, passed through to the backend.
    allcanon : bool, optional
        Canonicalization flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.
    ethr : float | None
        Optional energy threshold to accelerate by pre-sorting
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    new_molecule_list : list[Molecule]
        New Molecule objects reconstructed from the sorted atomic numbers
        and positions returned by the backend. The list contains only n_structures
        defined as ``unique`` according to the selected thresholds.

    Raises
    ------
    TypeError
        If ``molecule_list`` does not contain Molecule instances.
    ValueError
        If ``molecule_list`` is empty or if the Molecules do not all have
        the same number of atoms.
    """

    # --- Basic checks on molecule_list ---
    if not isinstance(molecule_list, (list, tuple)):
        raise TypeError(
            "sorter_irmsd_molecule expects a sequence (list/tuple) of Molecule objects"
        )

    if len(molecule_list) == 0:
        raise ValueError("molecule_list must contain at least one Molecule object")

    for i, mol in enumerate(molecule_list):
        if not isinstance(mol, Molecule):
            raise TypeError(
                "sorter_irmsd_molecule expects a sequence of Molecule objects; "
                f"item {i} has type {type(mol)}"
            )

    # --- Check that all Molecules have the same number of atoms and define nat ---
    nat = len(molecule_list[0])
    for i, mol in enumerate(molecule_list):
        if len(mol) != nat:
            raise ValueError(
                "All Molecule objects must have the same number of atoms; "
                f"item 0 has {nat} atoms, item {i} has {len(mol)} atoms"
            )

    # --- Remove structures too high in energy, if needed
    if printlvl > 0 and ewin is not None:
        print(f"EWIN energy window : {ewin} Ha ({ewin*627.5095:.2f} kcal/mol)")
    n_orig = len(molecule_list)
    molecule_list = prune_by_energy_window(molecule_list, ewin)
    n_new = len(molecule_list)
    n_diff = n_orig - n_new
    if printlvl > 0 and ewin is not None:
        print(f" --> removed {n_diff} of {n_orig} structures.\n")

    groups, new_molecule_list = sorter_irmsd_molecule(
        molecule_list, rthr, iinversion, allcanon, printlvl, ethr
    )

    new_molecule_list = first_by_assignment(new_molecule_list, groups)

    return new_molecule_list
