from __future__ import annotations

from typing import Tuple

import numpy as np

from ..core import Mol

#####################################################################################


def get_energy_ase(mol):
    """
    Return the energy stored in an ASE Mol object.

    Checks, in order:
        1. mol.info["energy"]
        2. mol.calc.get_potential_energy()
        3. mol.calc.results["energy"]

    Returns None if nothing is found.
    """
    # 1. Info dict (extended XYZ metadata)
    E = mol.info.get("energy")
    if isinstance(E, (int, float)):
        return float(E)

    # 2. Calculator energy (the ASE way)
    calc = mol.calc
    if calc is not None:
        try:
            return float(mol.get_potential_energy())
        except Exception:
            pass

        # 3. Direct access to stored results, if present
        try:
            E = calc.results.get("energy")
            if isinstance(E, (int, float)):
                return float(E)
        except Exception:
            pass

    # 4. Nothing found
    return None


def get_energies_from_mol_list(mol_list):
    """
    Given a list of ASE Mol objects, call `get_energy_ase(mol)` for each,
    collect the energies into a float NumPy array, and replace any `None`
    returned by the energy function with 0.0.
    """
    energies = []
    for mol in mol_list:
        e = get_energy_ase(mol)  # user-defined energy calculator
        energies.append(0.0 if e is None else float(e))
    return np.array(energies, dtype=float)


def get_cn_ase(mol) -> Tuple["Mol", np.ndarray]:
    """
    Optional utility: accepts ASE Mol, adapts to core get_cn_fortran, and
    returns a numpy array with the coordination numbers per atom
    """

    require_ase()

    from ase import Mol

    from ..api.cn_exposed import get_cn_fortran

    if not isinstance(mol, Mol):
        raise TypeError("get_cn_ase expects an ase.Mol object")

    Z = mol.get_atomic_numbers()  # (N,)
    pos = mol.get_positions()  # (N, 3) float64

    new_cn = get_cn_fortran(Z, pos)

    return new_cn


def get_axis_mol(mol) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optional utility: accepts MOL object, adapts to core get_axis, and
    returns a numpy array with the rotation constants in Mhz, the average  momentum in a.u.
    and the rot. matrix.
    """
    require_ase()

    from ase import Mol

    from ..api.axis_exposed import get_axis

    if not isinstance(mol, Mol):
        raise TypeError("get_cn_ase expects an ase.Mol object")

    Z = mol.get_atomic_numbers()  # (N,)
    pos = mol.get_positions()  # (N, 3) float64

    rot, avmom, evec = get_axis(Z, pos)

    return rot, avmom, evec


def get_canonical_ase(
    mol, wbo=None, invtype="apsp+", heavy: bool = False
) -> np.ndarray:
    """
    Optional utility: accepts ASE Mol, adapts to core get_canonical_fortran, and
    returns rank and invariants arrays.
    """

    require_ase()

    from ase import Mol

    from ..api.canonical_exposed import get_canonical_fortran

    if not isinstance(mol, Mol):
        raise TypeError("get_canonical_ase expects an ase.Mol object")

    Z = mol.get_atomic_numbers()  # (N,)
    pos = mol.get_positions()  # (N, 3) float64

    rank = get_canonical_fortran(Z, pos, wbo=wbo, invtype=invtype, heavy=heavy)

    return rank


#########################################################################################


def get_rmsd_ase(mol1, mol2, mask=None) -> Tuple[float, "Mol", np.ndarray]:
    """
    Optional ASE utility: operate on TWO ASE Mol. Returns the RMSD in Angström,
    the modified second Mol plus the 3×3 rotation matrix produced by
    the Fortran routine.
    """

    require_ase()
    from ase import Mol

    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(mol1, Mol) or not isinstance(mol2, Mol):
        raise TypeError("ase_quaternion_rmsd expects two ase.Mol objects")

    Z1 = mol1.get_atomic_numbers()  # (N1,)
    P1 = mol1.get_positions()  # (N1, 3)
    Z2 = mol2.get_atomic_numbers()  # (N2,)
    P2 = mol2.get_positions()  # (N2, 3)

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    # Return ONLY the modified second structure, as requested
    new_mol2 = mol2.copy()
    new_mol2.set_positions(new_P2, apply_constraint=False)

    return rmsdval, new_mol2, umat


def get_irmsd_ase(mol1, mol2, iinversion=0) -> Tuple[float, "Mol", "Mol"]:
    """
    Optional ASE utility: operate on TWO ASE Mol. Returns the iRMSD in Angström,
    the modified second Mol plus the 3×3 rotation matrix produced by
    the Fortran routine.
    """

    require_ase()
    from ase import Mol

    from ..api.irmsd_exposed import get_irmsd

    if not isinstance(mol1, Mol) or not isinstance(mol2, Mol):
        raise TypeError("ase_quaternion_rmsd expects two ase.Mol objects")

    Z1 = mol1.get_atomic_numbers()  # (N1,)
    P1 = mol1.get_positions()  # (N1, 3)
    Z2 = mol2.get_atomic_numbers()  # (N2,)
    P2 = mol2.get_positions()  # (N2, 3)

    irmsdval, new_Z1, new_P1, new_Z2, new_P2 = get_irmsd(
        Z1, P1, Z2, P2, iinversion=iinversion
    )
    new_mol1 = mol1.copy()
    new_mol1.set_atomic_numbers(new_Z1)
    new_mol1.set_positions(new_P1, apply_constraint=False)
    # Return ONLY the modified second structure, as requested
    new_mol2 = mol2.copy()
    new_mol2.set_atomic_numbers(new_Z2)
    new_mol2.set_positions(new_P2, apply_constraint=False)

    return irmsdval, new_mol1, new_mol2


####################################################################################
def sorter_irmsd_ase(
    mol_list: Sequence["Mol"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["Mol"]]:
    """
    ASE wrapper around sorter_irmsd.

    Parameters
    ----------
    mol_list : sequence of ase.Mol
        List/sequence of ASE Mol objects. All must have the same number of mol.
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
        Group indices for the first `nat` mol.
    new_mol_list : list of ase.Mol
        New Mol objects reconstructed from the sorted atom types and positions.
    """
    require_ase()
    from ase import Mol
    from ..api.sorter_exposed import sorter_irmsd

    # --- Basic checks on mol_list ---
    if not isinstance(mol_list, (list, tuple)):
        raise TypeError(
            "ase_sorter_irmsd expects a sequence (list/tuple) of ase.Mol objects"
        )

    if len(mol_list) == 0:
        raise ValueError("mol_list must contain at least one ase.Mol object")

    for i, at in enumerate(mol_list):
        if not isinstance(at, Mol):
            raise TypeError(
                f"ase_sorter_irmsd expects a sequence of ase.Mol objects; "
                f"item {i} has type {type(at)}"
            )

    # --- Check that all Mol have the same number of mol and define nat ---
    nat = len(mol_list[0])
    for i, at in enumerate(mol_list):
        if len(at) != nat:
            raise ValueError(
                "All Mol objects must have the same number of mol; "
                f"item 0 has {nat} mol, item {i} has {len(at)} mol"
            )

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for at in mol_list:
        Z = np.asarray(at.numbers, dtype=np.int32)
        P = np.asarray(at.get_positions(), dtype=np.float64)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Mol positions array must have shape (nat, 3); " f"got {P.shape}"
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

    # --- Reconstruct new ASE Mol objects ---
    new_mol_list: List[Mol] = []
    for at_orig, Z_new, P_new in zip(mol_list, Z_structs, xyz_structs):
        # Preserve cell and PBC; other metadata can be copied as needed
        new_at = Mol(
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

        new_mol_list.append(new_at)

    return groups, new_mol_list


####################################################################################
def delta_irmsd_list_ase(
    mol_list: Sequence["Mol"],
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["Mol"]]:
    """
    ASE wrapper around delta_irmsd_list.

    Parameters
    ----------
    mol_list : sequence of ase.Mol
        List/sequence of ASE Mol objects. All must have the same number of mol.
    iinversion : int, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.

    Returns
    -------
    delta : (nat,) float64
        Group indices for the first `nat` mol.
    new_mol_list : list of ase.Mol
        New Mol objects reconstructed from the sorted atom types and positions.
    """
    require_ase()
    from ase import Mol
    from ..api.sorter_exposed import delta_irmsd_list

    # --- Basic checks on mol_list ---
    if not isinstance(mol_list, (list, tuple)):
        raise TypeError(
            "delta_irmsd_list expects a sequence (list/tuple) of ase.Mol objects"
        )

    if len(mol_list) == 0:
        raise ValueError("mol_list must contain at least one ase.Mol object")

    for i, at in enumerate(mol_list):
        if not isinstance(at, Mol):
            raise TypeError(
                f"delta_irmsd_list expects a sequence of ase.Mol objects; "
                f"item {i} has type {type(at)}"
            )

    # --- Check that all Mol have the same number of mol and define nat ---
    nat = len(mol_list[0])
    for i, at in enumerate(mol_list):
        if len(at) != nat:
            raise ValueError(
                "All Mol objects must have the same number of mol; "
                f"item 0 has {nat} mol, item {i} has {len(at)} mol"
            )

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for at in mol_list:
        Z = np.asarray(at.numbers, dtype=np.int32)
        P = np.asarray(at.get_positions(), dtype=np.float64)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Mol positions array must have shape (nat, 3); " f"got {P.shape}"
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

    # --- Reconstruct new ASE Mol objects ---
    new_mol_list: List[Mol] = []
    for at_orig, Z_new, P_new in zip(mol_list, Z_structs, xyz_structs):
        # Preserve cell and PBC; other metadata can be copied as needed
        new_at = Mol(
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

        new_mol_list.append(new_at)

    return delta, new_mol_list
