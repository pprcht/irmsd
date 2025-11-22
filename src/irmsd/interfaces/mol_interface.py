from __future__ import annotations

from typing import Tuple

import numpy as np

from ..core import Molecule

#####################################################################################


def get_energies_from_molecule_list(molecule_list):
    """
    Given a list of irmsd Molecule objects, call get_potential_energy() for each,
    collect the energies into a float NumPy array, and replace any `None`
    returned by the energy function with 0.0.
    """
    energies = []
    for molecule in molecule_list:
        e = molecule.get_potential_energy()  # get energy, if stored
        energies.append(0.0 if e is None else float(e))
    return np.array(energies, dtype=float)


#########################################################################################


def get_rmsd_molecule(
    molecule1, molecule2, mask=None
) -> Tuple[float, "Molecule", np.ndarray]:
    """
    Optional irmsd utility: operate on TWO irmsd.Molecule. Returns the RMSD in Angström,
    the modified second Molecule plus the 3×3 rotation matrix produced by
    the Fortran routine.
    """
    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(molecule1, Molecule) or not isinstance(molecule2, Molecule):
        raise TypeError("get_rmsd_molecule expects two irmsd.Molecule objects")

    Z1 = molecule1.get_atomic_numbers()  # (N1,)
    P1 = molecule1.get_positions()  # (N1, 3)
    Z2 = molecule2.get_atomic_numbers()  # (N2,)
    P2 = molecule2.get_positions()  # (N2, 3)

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    # Return ONLY the modified second structure, as requested
    new_molecule2 = molecule2.copy()
    new_molecule2.set_positions(new_P2, apply_constraint=False)

    return rmsdval, new_molecule2, umat


#########################################################################################


def get_irmsd_molecule(
    molecule1, molecule2, iinversion=0
) -> Tuple[float, "Molecule", "Molecule"]:
    """
    Optional irmsd utility: operate on TWO irmsd Molecule. Returns the iRMSD in Angström,
    the modified second Molecule plus the 3×3 rotation matrix produced by
    the Fortran routine.
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
    new_molecule1.set_positions(new_P1, apply_constraint=False)
    # Return ONLY the modified second structure, as requested
    new_molecule2 = molecule2.copy()
    new_molecule2.set_atomic_numbers(new_Z2)
    new_molecule2.set_positions(new_P2, apply_constraint=False)

    return irmsdval, new_molecule1, new_molecule2


####################################################################################


def sorter_irmsd_molecule(
    molecule_list: Sequence["Molecule"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["Molecule"]]:
    """
    irmsd wrapper around sorter_irmsd.

    Parameters
    ----------
    molecule_list : sequence of irmsd.Molecule
        List/sequence of irmsd Molecule objects. All must have the same number of molecule.
    rthr : float
        Distance threshold for sorter_irmsd.
    iinversion : int, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.

    Returns
    -------
    groups : (nat,) int32
        Group indices for the first `nat` molecule.
    new_molecule_list : list of irmsd.Molecule
        New Molecule objects reconstructed from the sorted atom types and positions.
    """
    from ..api.sorter_exposed import sorter_irmsd

    # --- Basic checks on molecule_list ---
    if not isinstance(molecule_list, (list, tuple)):
        raise TypeError(
            "sorter_irmsd_molecule expects a sequence (list/tuple) of irmsd.Molecule objects"
        )

    if len(molecule_list) == 0:
        raise ValueError(
            "molecule_list must contain at least one irmsd.Molecule object"
        )

    for i, at in enumerate(molecule_list):
        if not isinstance(at, Molecule):
            raise TypeError(
                f"sorter_irmsd_molecule expects a sequence of irmsd.Molecule objects; "
                f"item {i} has type {type(at)}"
            )

    # --- Check that all Molecule have the same number of molecule and define nat ---
    nat = len(molecule_list[0])
    for i, at in enumerate(molecule_list):
        if len(at) != nat:
            raise ValueError(
                "All Molecule objects must have the same number of molecule; "
                f"item 0 has {nat} molecule, item {i} has {len(at)} molecule"
            )

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for at in molecule_list:
        Z = np.asarray(at.numbers, dtype=np.int32)
        P = np.asarray(at.get_positions(), dtype=np.float64)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Molecule positions array must have shape (nat, 3); "
                f"got {P.shape}"
            )

        atom_numbers_list.append(Z)
        positions_list.append(P)

    # --- Call the Fortran-backed sorter_irmsd ---
    groups, xyz_structs, Z_structs = sorter_irmsd(
        atom_numbers_list=atom_numbers_list,
        positions_list=positions_list,
        nat=nat,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    # --- Reconstruct new irmsd Molecule objects ---
    new_molecule_list: List[Molecule] = []
    for at_orig, Z_new, P_new in zip(molecule_list, Z_structs, xyz_structs):
        # Preserve cell and PBC; other metadata can be copied as needed
        new_at = Molecule(
            numbers=Z_new,
            positions=P_new,
            cell=at_orig.cell,
            pbc=at_orig.pbc,
        )

        # Optionally preserve info and constraints
        new_at.info = dict(at_orig.info)
        new_at.calc = at_orig.calc
        if getattr(at_orig, "constraints", None):
            new_at.set_constraint(at_orig.constraints)

        new_molecule_list.append(new_at)

    return groups, new_molecule_list


####################################################################################


def delta_irmsd_list_molecule(
    molecule_list: Sequence["Molecule"],
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["Molecule"]]:
    """
    irmsd wrapper around delta_irmsd_list.

    Parameters
    ----------
    molecule_list : sequence of irmsd.Molecule
        List/sequence of irmsd Molecule objects. All must have the same number of molecule.
    iinversion : int, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.

    Returns
    -------
    delta : (nat,) float64
        Group indices for the first `nat` molecule.
    new_molecule_list : list of irmsd.Molecule
        New Molecule objects reconstructed from the sorted atom types and positions.
    """
    from ..api.sorter_exposed import delta_irmsd_list

    # --- Basic checks on molecule_list ---
    if not isinstance(molecule_list, (list, tuple)):
        raise TypeError(
            "delta_irmsd_list expects a sequence (list/tuple) of irmsd.Molecule objects"
        )

    if len(molecule_list) == 0:
        raise ValueError(
            "molecule_list must contain at least one irmsd.Molecule object"
        )

    for i, at in enumerate(molecule_list):
        if not isinstance(at, Molecule):
            raise TypeError(
                f"delta_irmsd_list expects a sequence of irmsd.Molecule objects; "
                f"item {i} has type {type(at)}"
            )

    # --- Check that all Molecule have the same number of molecule and define nat ---
    nat = len(molecule_list[0])
    for i, at in enumerate(molecule_list):
        if len(at) != nat:
            raise ValueError(
                "All Molecule objects must have the same number of molecule; "
                f"item 0 has {nat} molecule, item {i} has {len(at)} molecule"
            )

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for at in molecule_list:
        Z = np.asarray(at.numbers, dtype=np.int32)
        P = np.asarray(at.get_positions(), dtype=np.float64)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Molecule positions array must have shape (nat, 3); "
                f"got {P.shape}"
            )

        atom_numbers_list.append(Z)
        positions_list.append(P)

    # --- Call the Fortran-backed sorter_irmsd ---
    delta, xyz_structs, Z_structs = delta_irmsd_list(
        atom_numbers_list=atom_numbers_list,
        positions_list=positions_list,
        nat=nat,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    # --- Reconstruct new irmsd Molecule objects ---
    new_molecule_list: List[Molecule] = []
    for at_orig, Z_new, P_new in zip(molecule_list, Z_structs, xyz_structs):
        # Preserve cell and PBC; other metadata can be copied as needed
        new_at = Molecule(
            numbers=Z_new,
            positions=P_new,
            cell=at_orig.cell,
            pbc=at_orig.pbc,
        )

        # Optionally preserve info and constraints
        new_at.info = dict(at_orig.info)
        new_at.calc = at_orig.calc
        if getattr(at_orig, "constraints", None):
            new_at.set_constraint(at_orig.constraints)

        new_molecule_list.append(new_at)

    return delta, new_molecule_list
