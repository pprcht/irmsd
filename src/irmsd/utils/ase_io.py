from __future__ import annotations

from typing import Tuple

import numpy as np

from .utils import require_ase


def ase_to_fortran(atoms) -> Tuple["Atoms", np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core xyz_to_fortran, and returns:
      - a NEW Atoms with updated positions
      - the 3×3 matrix.
    """

    require_ase()

    from ..api.xyz_bridge import xyz_to_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("ase_to_fortran expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    new_pos, mat = xyz_to_fortran(Z, pos)

    new_atoms = atoms.copy()
    new_atoms.set_positions(new_pos, apply_constraint=False)
    return new_atoms, mat


###############################################################################


def ase_to_fortran_pair(atoms1, atoms2) -> Tuple["Atoms", np.ndarray, np.ndarray]:
    """
    Optional ASE utility: operate on TWO ASE Atoms. Returns ONLY the modified second Atoms
    plus both 3×3 matrices produced by the Fortran routine.

    Parameters
    ----------
    atoms1 : ase.Atoms
    atoms2 : ase.Atoms

    Returns
    -------
    new_atoms2 : ase.Atoms
        Copy of `atoms2` with updated positions (atoms1 is NOT returned).
    mat1 : (3, 3) float64 ndarray (Fortran-ordered)
        Matrix corresponding to the first structure.
    mat2 : (3, 3) float64 ndarray (Fortran-ordered)
        Matrix corresponding to the second structure.

    Notes
    -----
    - The corresponding core API is `irmsd.api.xyz_bridge.xyz_to_fortran_pair`, which takes
      plain NumPy arrays and has no ASE dependency.
    """

    require_ase()

    from ..api.xyz_bridge import xyz_to_fortran_pair

    if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
        raise TypeError("ase_to_fortran_pair expects two ase.Atoms objects")

    Z1 = atoms1.get_atomic_numbers()  # (N1,)
    P1 = atoms1.get_positions()  # (N1, 3)
    Z2 = atoms2.get_atomic_numbers()  # (N2,)
    P2 = atoms2.get_positions()  # (N2, 3)

    new_P1, new_P2, M1, M2 = xyz_to_fortran_pair(Z1, P1, Z2, P2)

    # Return ONLY the modified second structure, as requested
    new_atoms2 = atoms2.copy()
    new_atoms2.set_positions(new_P2, apply_constraint=False)

    return new_atoms2, M1, M2


#####################################################################################

def get_energy_ase(atoms):
    """
    Return the energy stored in an ASE Atoms object.

    Checks, in order:
        1. atoms.info["energy"]
        2. atoms.calc.get_potential_energy()
        3. atoms.calc.results["energy"]

    Returns None if nothing is found.
    """
    # 1. Info dict (extended XYZ metadata)
    E = atoms.info.get("energy")
    if isinstance(E, (int, float)):
        return float(E)

    # 2. Calculator energy (the ASE way)
    calc = atoms.calc
    if calc is not None:
        try:
            return float(atoms.get_potential_energy())
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


def get_cn_ase(atoms) -> Tuple["Atoms", np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core get_cn_fortran, and
    returns a numpy array with the coordination numbers per atom
    """

    require_ase()

    from ase import Atoms

    from ..api.cn_exposed import get_cn_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("get_cn_ase expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    new_cn = get_cn_fortran(Z, pos)

    return new_cn


def get_axis_ase(atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core get_axis_fortran, and
    returns a numpy array with the rotation constants in Mhz, the average  momentum in a.u.
    and the rot. matrix.
    """
    require_ase()

    from ase import Atoms

    from ..api.axis_exposed import get_axis_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("get_cn_ase expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    rot, avmom, evec = get_axis_fortran(Z, pos)

    return rot, avmom, evec


def get_canonical_ase(
    atoms, wbo=None, invtype="apsp+", heavy: bool = False
) -> np.ndarray:
    """
    Optional utility: accepts ASE Atoms, adapts to core get_canonical_fortran, and
    returns rank and invariants arrays.
    """

    require_ase()

    from ase import Atoms

    from ..api.canonical_exposed import get_canonical_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("get_canonical_ase expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    rank = get_canonical_fortran(Z, pos, wbo=wbo, invtype=invtype, heavy=heavy)

    return rank


#########################################################################################


def get_rmsd_ase(atoms1, atoms2, mask=None) -> Tuple[float, "Atoms", np.ndarray]:
    """
    Optional ASE utility: operate on TWO ASE Atoms. Returns the RMSD in Angström,
    the modified second Atoms plus the 3×3 rotation matrix produced by
    the Fortran routine.
    """

    require_ase()
    from ase import Atoms

    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
        raise TypeError("ase_quaternion_rmsd expects two ase.Atoms objects")

    Z1 = atoms1.get_atomic_numbers()  # (N1,)
    P1 = atoms1.get_positions()  # (N1, 3)
    Z2 = atoms2.get_atomic_numbers()  # (N2,)
    P2 = atoms2.get_positions()  # (N2, 3)

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    # Return ONLY the modified second structure, as requested
    new_atoms2 = atoms2.copy()
    new_atoms2.set_positions(new_P2, apply_constraint=False)

    return rmsdval, new_atoms2, umat


def get_irmsd_ase(atoms1, atoms2, iinversion=0) -> Tuple[float, "Atoms", "Atoms"]:
    """
    Optional ASE utility: operate on TWO ASE Atoms. Returns the iRMSD in Angström,
    the modified second Atoms plus the 3×3 rotation matrix produced by
    the Fortran routine.
    """

    require_ase()
    from ase import Atoms

    from ..api.irmsd_exposed import get_irmsd

    if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
        raise TypeError("ase_quaternion_rmsd expects two ase.Atoms objects")

    Z1 = atoms1.get_atomic_numbers()  # (N1,)
    P1 = atoms1.get_positions()  # (N1, 3)
    Z2 = atoms2.get_atomic_numbers()  # (N2,)
    P2 = atoms2.get_positions()  # (N2, 3)

    irmsdval, new_Z1, new_P1, new_Z2, new_P2 = get_irmsd(
        Z1, P1, Z2, P2, iinversion=iinversion
    )
    new_atoms1 = atoms1.copy()
    new_atoms1.set_atomic_numbers(new_Z1)
    new_atoms1.set_positions(new_P1, apply_constraint=False)
    # Return ONLY the modified second structure, as requested
    new_atoms2 = atoms2.copy()
    new_atoms2.set_atomic_numbers(new_Z2)
    new_atoms2.set_positions(new_P2, apply_constraint=False)

    return irmsdval, new_atoms1, new_atoms2


####################################################################################
def sorter_irmsd_ase(
    atoms_list: Sequence["Atoms"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["Atoms"]]:
    """
    ASE wrapper around sorter_irmsd.

    Parameters
    ----------
    atoms_list : sequence of ase.Atoms
        List/sequence of ASE Atoms objects. All must have the same number of atoms.
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
        Group indices for the first `nat` atoms.
    new_atoms_list : list of ase.Atoms
        New Atoms objects reconstructed from the sorted atom types and positions.
    """
    require_ase()
    from ase import Atoms
    from ..api.sorter_exposed import sorter_irmsd

    # --- Basic checks on atoms_list ---
    if not isinstance(atoms_list, (list, tuple)):
        raise TypeError(
            "ase_sorter_irmsd expects a sequence (list/tuple) of ase.Atoms objects"
        )

    if len(atoms_list) == 0:
        raise ValueError("atoms_list must contain at least one ase.Atoms object")

    for i, at in enumerate(atoms_list):
        if not isinstance(at, Atoms):
            raise TypeError(
                f"ase_sorter_irmsd expects a sequence of ase.Atoms objects; "
                f"item {i} has type {type(at)}"
            )

    # --- Check that all Atoms have the same number of atoms and define nat ---
    nat = len(atoms_list[0])
    for i, at in enumerate(atoms_list):
        if len(at) != nat:
            raise ValueError(
                "All Atoms objects must have the same number of atoms; "
                f"item 0 has {nat} atoms, item {i} has {len(at)} atoms"
            )

    # --- Build atom_numbers_list and positions_list ---
    atom_numbers_list: List[np.ndarray] = []
    positions_list: List[np.ndarray] = []

    for at in atoms_list:
        Z = np.asarray(at.numbers, dtype=np.int32)
        P = np.asarray(at.get_positions(), dtype=np.float64)

        if P.shape != (nat, 3):
            raise ValueError(
                "Each Atoms positions array must have shape (nat, 3); " f"got {P.shape}"
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

    # --- Reconstruct new ASE Atoms objects ---
    new_atoms_list: List[Atoms] = []
    for at_orig, Z_new, P_new in zip(atoms_list, Z_structs, xyz_structs):
        # Preserve cell and PBC; other metadata can be copied as needed
        new_at = Atoms(
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

        new_atoms_list.append(new_at)

    return groups, new_atoms_list
