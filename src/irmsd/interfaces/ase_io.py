from __future__ import annotations

from typing import overload, Sequence, Tuple, List
import numpy as np

from ..utils.utils import require_ase
from ..core.molecule import Molecule
from .mol_interface import (
    get_rmsd_molecule,
    get_irmsd_molecule,
    sorter_irmsd_molecule,
    delta_irmsd_list_molecule,
)

# -------------------------------------------------------------------
# Some I/O
# -------------------------------------------------------------------


def get_energy_ase(atoms):
    """
    Return the energy stored in an ASE Atoms object.

    Checks, in order:
        1. atoms.info["energy"]
        2. calc.results["energy"] / ["free_energy"] / ["enthalpy"]
        3. atoms.get_potential_energy() *only if it will NOT trigger a calculation*

    Returns None if nothing is found.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_energy_ase expects an ase.Atoms object")

    # 1. info["energy"]
    E = atoms.info.get("energy")
    if isinstance(E, (int, float)):
        return float(E)

    # 2. calculator results
    calc = getattr(atoms, "calc", None)
    if calc is None:
        return None

    results = getattr(calc, "results", None)
    if isinstance(results, dict):
        for key in ("energy", "free_energy", "enthalpy"):
            val = results.get(key)
            if isinstance(val, (int, float)):
                return float(val)

    # 3. get_potential_energy only if it won't trigger a calculation
    try:
        if hasattr(calc, "calculation_required"):
            if calc.calculation_required(atoms):
                return None
        return float(atoms.get_potential_energy())
    except Exception:
        pass

    return None


@overload
def ase_to_molecule(atoms: "ase.Atoms") -> Molecule: ...
@overload
def ase_to_molecule(atoms: Sequence["ase.Atoms"]) -> list[Molecule]: ...


def ase_to_molecule(atoms):
    """
    Convert an ASE `Atoms` object (or a sequence of them) into the internal
    `irmsd.core.Molecule` type.

    This function is intentionally non-invasive: it does not trigger any new
    ASE calculator evaluations. It merely extracts whatever structural and
    metadata information is already present in the ASE object.

    Parameters
    ----------
    atoms : ase.Atoms or Sequence[ase.Atoms]
        A single ASE Atoms instance or a sequence of them.

    Returns
    -------
    Molecule or list[Molecule]
        - If `atoms` is a single Atoms object, a single Molecule is returned.
        - If `atoms` is a sequence of Atoms objects, a list of Molecules is
          returned in the same order.

    Notes
    -----
    - This routine requires ASE to be installed. If ASE is missing, a clear
      and controlled error message is raised via `require_ase()`.
    - This routine does not modify either the input Atoms object or its
      attached calculator.
    - The returned Molecule is guaranteed to be fully self-contained and
      ASE-independent.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If the input is neither an ASE Atoms instance nor a sequence of them.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    def _one(a):
        # check type
        if not isinstance(a, ASEAtoms):
            raise TypeError("ase_to_molecule expects ase.Atoms or a sequence thereof")

        symbols = a.get_chemical_symbols()
        positions = a.get_positions()

        cell_array = None
        try:
            cell = a.get_cell()
            cell_array = np.asarray(cell, float)
            if cell_array.shape != (3, 3):
                cell_array = None
        except Exception:
            pass

        pbc = tuple(bool(x) for x in getattr(a, "pbc", (False, False, False)))

        info = dict(getattr(a, "info", {}))
        energy = get_energy_ase(a)

        return Molecule(
            symbols=symbols,
            positions=positions,
            energy=energy,
            info=info,
            cell=cell_array,
            pbc=pbc,
        )

    # sequence vs single
    if isinstance(atoms, ASEAtoms):
        return _one(atoms)
    return [_one(a) for a in atoms]


@overload
def molecule_to_ase(molecules: Molecule) -> "ase.Atoms": ...
@overload
def molecule_to_ase(molecules: Sequence[Molecule]) -> list["ase.Atoms"]: ...


def molecule_to_ase(
    molecules: Molecule | Sequence[Molecule],
):
    """
    Convert an internal `irmsd.core.Molecule` instance (or a sequence of them)
    into ASE `Atoms` objects.

    This routine performs a purely structural and metadata-level conversion:
    it does not create or attach any calculator, nor does it trigger any new
    ASE calculations.

    Parameters
    ----------
    molecules : Molecule or Sequence[Molecule]
        A single Molecule or a sequence of Molecule objects.

    Returns
    -------
    ase.Atoms or list[ase.Atoms]
        - If `molecules` is a single Molecule, a single ASE Atoms object is
          returned.
        - If `molecules` is a sequence of Molecule objects, a list of ASE
          Atoms objects is returned in the same order.

    Notes
    -----
    - This routine requires ASE to be installed. If ASE is missing, a clear
      RuntimeError is raised via `require_ase()`.
    - The returned Atoms objects are structurally independent copies; further
      modifications to the original Molecule will not affect them.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If the input is neither a Molecule instance nor a sequence of Molecule
        instances.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    def _one(mol: Molecule) -> "ase.Atoms":  # type: ignore[name-defined]
        if not isinstance(mol, Molecule):
            raise TypeError(
                "molecule_to_ase expects Molecule or a sequence of Molecule"
            )

        symbols = mol.get_chemical_symbols()
        positions = mol.get_positions(copy=True)

        # Cell: either a proper (3,3) array or None
        cell = None
        if mol.cell is not None:
            cell_arr = np.asarray(mol.cell, dtype=float)
            if cell_arr.shape == (3, 3):
                cell = cell_arr

        # PBC: pass through if set, otherwise False
        pbc = mol.pbc if mol.pbc is not None else False

        # Info: shallow copy to avoid mutating the original
        info = dict(mol.info)

        # Energy: only set info["energy"] if it is not already present
        if mol.energy is not None and "energy" not in info:
            info["energy"] = float(mol.energy)

        atoms = ASEAtoms(
            symbols=symbols,
            positions=positions,
            cell=cell,
            pbc=pbc,
            info=info,
        )
        return atoms

    if isinstance(molecules, Molecule):
        return _one(molecules)

    # Treat as sequence
    try:
        return [_one(m) for m in molecules]
    except TypeError as exc:
        raise TypeError(
            "molecule_to_ase expects either a single Molecule or a sequence of Molecule objects"
        ) from exc


def get_energies_from_atoms_list(atoms_list):
    """
    Given a list of ASE Atoms objects, call `get_energy_ase(atoms)` for each,
    collect the energies into a float NumPy array, and replace any `None`
    returned by the energy function with 0.0.
    """
    energies = []
    for atoms in atoms_list:
        e = get_energy_ase(atoms)  # user-defined energy calculator
        energies.append(0.0 if e is None else float(e))
    return np.array(energies, dtype=float)


# -----------------------------------------------------------------------------
# Callable functions
# -----------------------------------------------------------------------------


def get_cn_ase(atoms) -> np.ndarray:
    """
    High-level utility: accepts an ASE Atoms object, converts it into an
    internal Molecule instance, and returns the coordination-number array
    as computed by `Molecule.get_cn()`.

    This routine does *not* trigger any new ASE calculator evaluation.

    Parameters
    ----------
    atoms : ase.Atoms
        A single ASE Atoms object.

    Returns
    -------
    np.ndarray
        Array of coordination numbers with shape (N,).
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_cn_ase expects a single ASE Atoms object")

    # Convert ASE → Molecule using conversion routine
    mol: Molecule = ase_to_molecule(atoms)

    # Delegate to Molecule API
    return mol.get_cn()


def get_axis_ase(atoms) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    High-level utility: accepts an ASE Atoms object, converts it into an
    internal Molecule instance, and returns rotation constants, average
    angular momentum, and eigenvectors via `Molecule.get_axis()`.

    This routine never triggers a new ASE calculator evaluation.

    Parameters
    ----------
    atoms : ase.Atoms
        A single ASE Atoms object.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        rot_constants_MHz, avg_momentum_au, rotation_matrix
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_axis_ase expects a single ASE Atoms object")

    # Convert ASE → Molecule
    mol: Molecule = ase_to_molecule(atoms)

    # Delegate to Molecule's high-level API
    return mol.get_axis()


def get_canonical_ase(
    atoms,
    wbo: np.ndarray | None = None,
    invtype: str = "apsp+",
    heavy: bool = False,
) -> np.ndarray:
    """
    High-level utility: accepts an ASE Atoms object, converts it into an
    internal Molecule instance, and returns the canonicalization rank /
    invariants as computed by `Molecule.get_canonical()`.

    This routine does not trigger any new ASE calculator evaluation.

    Parameters
    ----------
    atoms : ase.Atoms
        A single ASE Atoms object.
    wbo : np.ndarray or None, optional
        Optional Wiberg bond order matrix or similar, passed through to
        `Molecule.get_canonical()` and ultimately to the Fortran backend.
    invtype : str, optional
        Invariant type selector, e.g. "apsp+" (default). Forwarded directly
        to the canonicalization backend.
    heavy : bool, optional
        If True, restricts invariants to heavy atoms only, as defined by the
        underlying backend. Defaults to False.

    Returns
    -------
    np.ndarray
        Canonicalization rank / invariants array as returned by
        `Molecule.get_canonical()`.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If `atoms` is not an ASE Atoms instance.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_canonical_ase expects a single ASE Atoms object")

    # Convert ASE → Molecule
    mol: Molecule = ase_to_molecule(atoms)

    # Delegate to the Molecule API
    return mol.get_canonical(wbo=wbo, invtype=invtype, heavy=heavy)


# -----------------------------------------------------------------------------
# Comparison functions
# -----------------------------------------------------------------------------
#
def get_rmsd_ase(atoms1, atoms2, mask=None) -> Tuple[float, "ase.Atoms", np.ndarray]:
    """
    ASE wrapper for ``get_rmsd_molecule``.

    Converts two ASE ``Atoms`` objects to internal Molecule objects, calls
    ``get_rmsd_molecule``, and converts the aligned second structure back to
    an ASE ``Atoms`` object.

    Parameters
    ----------
    atoms1 : ase.Atoms
        Reference structure.
    atoms2 : ase.Atoms
        Structure to be rotated/translated onto ``atoms1``.
    mask : array-like of bool, optional
        Optional mask selecting which atoms in the first structure participate
        in the RMSD (forwarded to the backend via ``get_rmsd_molecule``).

    Returns
    -------
    rmsd : float
        RMSD value in Ångström.
    new_atoms2 : ase.Atoms
        New ASE Atoms object with coordinates aligned to ``atoms1``.
    rotation_matrix : np.ndarray
        3×3 rotation matrix used for the alignment.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If inputs are not ASE Atoms.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms1, ASEAtoms) or not isinstance(atoms2, ASEAtoms):
        raise TypeError("get_rmsd_ase expects two ASE Atoms objects")

    mol1, mol2 = ase_to_molecule([atoms1, atoms2])  # sequence form

    rmsd, new_mol2, umat = get_rmsd_molecule(mol1, mol2, mask=mask)
    new_atoms2 = molecule_to_ase(new_mol2)

    return rmsd, new_atoms2, umat


def get_irmsd_ase(
    atoms1,
    atoms2,
    iinversion: int = 0,
) -> Tuple[float, "ase.Atoms", "ase.Atoms"]:
    """
    ASE wrapper for ``get_irmsd_molecule``.

    Converts two ASE ``Atoms`` objects to Molecules, calls
    ``get_irmsd_molecule``, and converts both resulting Molecules back to
    ASE ``Atoms`` objects.

    Parameters
    ----------
    atoms1 : ase.Atoms
        First structure.
    atoms2 : ase.Atoms
        Second structure.
    iinversion : int, optional
        Inversion flag passed through to the backend.

    Returns
    -------
    irmsd : float
        iRMSD value in Ångström.
    new_atoms1 : ase.Atoms
        New ASE Atoms object corresponding to the transformed first Molecule.
    new_atoms2 : ase.Atoms
        New ASE Atoms object corresponding to the transformed second Molecule.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms1, ASEAtoms) or not isinstance(atoms2, ASEAtoms):
        raise TypeError("get_irmsd_ase expects two ASE Atoms objects")

    mol1, mol2 = ase_to_molecule([atoms1, atoms2])

    irmsd, new_mol1, new_mol2 = get_irmsd_molecule(mol1, mol2, iinversion=iinversion)

    new_atoms1 = molecule_to_ase(new_mol1)
    new_atoms2 = molecule_to_ase(new_mol2)

    return irmsd, new_atoms1, new_atoms2


def sorter_irmsd_ase(
    atoms_list: Sequence["ase.Atoms"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["ase.Atoms"]]:
    """
    ASE wrapper for ``sorter_irmsd_molecule``.

    Converts a sequence of ASE ``Atoms`` objects to Molecules, calls
    ``sorter_irmsd_molecule``, and converts the resulting Molecules back
    to ASE ``Atoms`` objects.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        Sequence of ASE Atoms objects. All must have the same number of atoms.
    rthr : float
        Distance threshold for the sorter.
    iinversion : int, optional
        Inversion symmetry flag.
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.

    Returns
    -------
    groups : np.ndarray
        Integer array of shape (nat,) with group indices as returned by
        ``sorter_irmsd_molecule`` / backend.
    new_atoms_list : list[ase.Atoms]
        New ASE Atoms objects reconstructed from the sorted Molecules.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms_list, (list, tuple)):
        raise TypeError("sorter_irmsd_ase expects a sequence (list/tuple) of ASE Atoms")

    for i, at in enumerate(atoms_list):
        if not isinstance(at, ASEAtoms):
            raise TypeError(
                "sorter_irmsd_ase expects a sequence of ASE Atoms; "
                f"item {i} has type {type(at)}"
            )

    mols = ase_to_molecule(atoms_list)  # returns list[Molecule]

    groups, new_mols = sorter_irmsd_molecule(
        molecule_list=mols,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    new_atoms_list = [molecule_to_ase(m) for m in new_mols]

    return groups, new_atoms_list


def delta_irmsd_list_ase(
    atoms_list: Sequence["ase.Atoms"],
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["ase.Atoms"]]:
    """
    ASE wrapper for ``delta_irmsd_list_molecule``.

    Converts a sequence of ASE ``Atoms`` objects to Molecules, calls
    ``delta_irmsd_list_molecule``, and converts the resulting Molecules
    back to ASE ``Atoms`` objects.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        Sequence of ASE Atoms objects. All must have the same number of atoms.
    iinversion : int, optional
        Inversion symmetry flag.
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.

    Returns
    -------
    delta : np.ndarray
        Float array returned by the backend (see ``delta_irmsd_list`` for
        detailed semantics).
    new_atoms_list : list[ase.Atoms]
        New ASE Atoms objects reconstructed from the transformed Molecules.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms_list, (list, tuple)):
        raise TypeError(
            "delta_irmsd_list_ase expects a sequence (list/tuple) of ASE Atoms"
        )

    for i, at in enumerate(atoms_list):
        if not isinstance(at, ASEAtoms):
            raise TypeError(
                "delta_irmsd_list_ase expects a sequence of ASE Atoms; "
                f"item {i} has type {type(at)}"
            )

    mols = ase_to_molecule(atoms_list)

    delta, new_mols = delta_irmsd_list_molecule(
        molecule_list=mols,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    new_atoms_list = [molecule_to_ase(m) for m in new_mols]

    return delta, new_atoms_list


#########################################################################################


# def get_rmsd_ase(atoms1, atoms2, mask=None) -> Tuple[float, "Atoms", np.ndarray]:
#     """
#     Optional ASE utility: operate on TWO ASE Atoms. Returns the RMSD in Angström,
#     the modified second Atoms plus the 3×3 rotation matrix produced by
#     the Fortran routine.
#     """
#
#     require_ase()
#     from ase import Atoms
#
#     from ..api.rmsd_exposed import get_quaternion_rmsd_fortran
#
#     if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
#         raise TypeError("ase_quaternion_rmsd expects two ase.Atoms objects")
#
#     Z1 = atoms1.get_atomic_numbers()  # (N1,)
#     P1 = atoms1.get_positions()  # (N1, 3)
#     Z2 = atoms2.get_atomic_numbers()  # (N2,)
#     P2 = atoms2.get_positions()  # (N2, 3)
#
#     rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)
#
#     # Return ONLY the modified second structure, as requested
#     new_atoms2 = atoms2.copy()
#     new_atoms2.set_positions(new_P2, apply_constraint=False)
#
#     return rmsdval, new_atoms2, umat
#
#
# def get_irmsd_ase(atoms1, atoms2, iinversion=0) -> Tuple[float, "Atoms", "Atoms"]:
#    """
#    Optional ASE utility: operate on TWO ASE Atoms. Returns the iRMSD in Angström,
#    the modified second Atoms plus the 3×3 rotation matrix produced by
#    the Fortran routine.
#    """
#
#    require_ase()
#    from ase import Atoms
#
#    from ..api.irmsd_exposed import get_irmsd
#
#    if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
#        raise TypeError("ase_quaternion_rmsd expects two ase.Atoms objects")
#
#    Z1 = atoms1.get_atomic_numbers()  # (N1,)
#    P1 = atoms1.get_positions()  # (N1, 3)
#    Z2 = atoms2.get_atomic_numbers()  # (N2,)
#    P2 = atoms2.get_positions()  # (N2, 3)
#
#    irmsdval, new_Z1, new_P1, new_Z2, new_P2 = get_irmsd(
#        Z1, P1, Z2, P2, iinversion=iinversion
#    )
#    new_atoms1 = atoms1.copy()
#    new_atoms1.set_atomic_numbers(new_Z1)
#    new_atoms1.set_positions(new_P1, apply_constraint=False)
#    # Return ONLY the modified second structure, as requested
#    new_atoms2 = atoms2.copy()
#    new_atoms2.set_atomic_numbers(new_Z2)
#    new_atoms2.set_positions(new_P2, apply_constraint=False)
#
#    return irmsdval, new_atoms1, new_atoms2
#
#
#####################################################################################
# def sorter_irmsd_ase(
#    atoms_list: Sequence["Atoms"],
#    rthr: float,
#    iinversion: int = 0,
#    allcanon: bool = True,
#    printlvl: int = 0,
# ) -> Tuple[np.ndarray, List["Atoms"]]:
#    """
#    ASE wrapper around sorter_irmsd.
#
#    Parameters
#    ----------
#    atoms_list : sequence of ase.Atoms
#        List/sequence of ASE Atoms objects. All must have the same number of atoms.
#    rthr : float
#        Distance threshold for sorter_irmsd.
#    iinversion : int, optional
#        Inversion symmetry flag, passed through.
#    allcanon : bool, optional
#        Canonicalization flag, passed through.
#    printlvl : int, optional
#        Verbosity level, passed through.
#
#    Returns
#    -------
#    groups : (nat,) int32
#        Group indices for the first `nat` atoms.
#    new_atoms_list : list of ase.Atoms
#        New Atoms objects reconstructed from the sorted atom types and positions.
#    """
#    require_ase()
#    from ase import Atoms
#    from ..api.sorter_exposed import sorter_irmsd
#
#    # --- Basic checks on atoms_list ---
#    if not isinstance(atoms_list, (list, tuple)):
#        raise TypeError(
#            "ase_sorter_irmsd expects a sequence (list/tuple) of ase.Atoms objects"
#        )
#
#    if len(atoms_list) == 0:
#        raise ValueError("atoms_list must contain at least one ase.Atoms object")
#
#    for i, at in enumerate(atoms_list):
#        if not isinstance(at, Atoms):
#            raise TypeError(
#                f"ase_sorter_irmsd expects a sequence of ase.Atoms objects; "
#                f"item {i} has type {type(at)}"
#            )
#
#    # --- Check that all Atoms have the same number of atoms and define nat ---
#    nat = len(atoms_list[0])
#    for i, at in enumerate(atoms_list):
#        if len(at) != nat:
#            raise ValueError(
#                "All Atoms objects must have the same number of atoms; "
#                f"item 0 has {nat} atoms, item {i} has {len(at)} atoms"
#            )
#
#    # --- Build atom_numbers_list and positions_list ---
#    atom_numbers_list: List[np.ndarray] = []
#    positions_list: List[np.ndarray] = []
#
#    for at in atoms_list:
#        Z = np.asarray(at.numbers, dtype=np.int32)
#        P = np.asarray(at.get_positions(), dtype=np.float64)
#
#        if P.shape != (nat, 3):
#            raise ValueError(
#                "Each Atoms positions array must have shape (nat, 3); " f"got {P.shape}"
#            )
#
#        atom_numbers_list.append(Z)
#        positions_list.append(P)
#
#    # --- Call the Fortran-backed sorter_irmsd ---
#    groups, xyz_structs, Z_structs = sorter_irmsd(
#        atom_numbers_list=atom_numbers_list,
#        positions_list=positions_list,
#        nat=nat,
#        rthr=rthr,
#        iinversion=iinversion,
#        allcanon=allcanon,
#        printlvl=printlvl,
#    )
#
#    # --- Reconstruct new ASE Atoms objects ---
#    new_atoms_list: List[Atoms] = []
#    for at_orig, Z_new, P_new in zip(atoms_list, Z_structs, xyz_structs):
#        # Preserve cell and PBC; other metadata can be copied as needed
#        new_at = Atoms(
#            numbers=Z_new,
#            positions=P_new,
#            cell=at_orig.cell,
#            pbc=at_orig.pbc,
#        )
#
#        # Optionally preserve info and constraints
#        new_at.info = dict(at_orig.info)
#        new_at.calc = at_orig.calc
#        if getattr(at_orig, "constraints", None):
#            new_at.set_constraint(at_orig.constraints)
#
#        new_atoms_list.append(new_at)
#
#    return groups, new_atoms_list
#
#
#####################################################################################
# def delta_irmsd_list_ase(
#    atoms_list: Sequence["Atoms"],
#    iinversion: int = 0,
#    allcanon: bool = True,
#    printlvl: int = 0,
# ) -> Tuple[np.ndarray, List["Atoms"]]:
#    """
#    ASE wrapper around delta_irmsd_list.
#
#    Parameters
#    ----------
#    atoms_list : sequence of ase.Atoms
#        List/sequence of ASE Atoms objects. All must have the same number of atoms.
#    iinversion : int, optional
#        Inversion symmetry flag, passed through.
#    allcanon : bool, optional
#        Canonicalization flag, passed through.
#    printlvl : int, optional
#        Verbosity level, passed through.
#
#    Returns
#    -------
#    delta : (nat,) float64
#        Group indices for the first `nat` atoms.
#    new_atoms_list : list of ase.Atoms
#        New Atoms objects reconstructed from the sorted atom types and positions.
#    """
#    require_ase()
#    from ase import Atoms
#    from ..api.sorter_exposed import delta_irmsd_list
#
#    # --- Basic checks on atoms_list ---
#    if not isinstance(atoms_list, (list, tuple)):
#        raise TypeError(
#            "delta_irmsd_list expects a sequence (list/tuple) of ase.Atoms objects"
#        )
#
#    if len(atoms_list) == 0:
#        raise ValueError("atoms_list must contain at least one ase.Atoms object")
#
#    for i, at in enumerate(atoms_list):
#        if not isinstance(at, Atoms):
#            raise TypeError(
#                f"delta_irmsd_list expects a sequence of ase.Atoms objects; "
#                f"item {i} has type {type(at)}"
#            )
#
#    # --- Check that all Atoms have the same number of atoms and define nat ---
#    nat = len(atoms_list[0])
#    for i, at in enumerate(atoms_list):
#        if len(at) != nat:
#            raise ValueError(
#                "All Atoms objects must have the same number of atoms; "
#                f"item 0 has {nat} atoms, item {i} has {len(at)} atoms"
#            )
#
#    # --- Build atom_numbers_list and positions_list ---
#    atom_numbers_list: List[np.ndarray] = []
#    positions_list: List[np.ndarray] = []
#
#    for at in atoms_list:
#        Z = np.asarray(at.numbers, dtype=np.int32)
#        P = np.asarray(at.get_positions(), dtype=np.float64)
#
#        if P.shape != (nat, 3):
#            raise ValueError(
#                "Each Atoms positions array must have shape (nat, 3); " f"got {P.shape}"
#            )
#
#        atom_numbers_list.append(Z)
#        positions_list.append(P)
#
#    # --- Call the Fortran-backed sorter_irmsd ---
#    delta, xyz_structs, Z_structs = delta_irmsd_list(
#        atom_numbers_list=atom_numbers_list,
#        positions_list=positions_list,
#        nat=nat,
#        iinversion=iinversion,
#        allcanon=allcanon,
#        printlvl=printlvl,
#    )
#
#    # --- Reconstruct new ASE Atoms objects ---
#    new_atoms_list: List[Atoms] = []
#    for at_orig, Z_new, P_new in zip(atoms_list, Z_structs, xyz_structs):
#        # Preserve cell and PBC; other metadata can be copied as needed
#        new_at = Atoms(
#            numbers=Z_new,
#            positions=P_new,
#            cell=at_orig.cell,
#            pbc=at_orig.pbc,
#        )
#
#        # Optionally preserve info and constraints
#        new_at.info = dict(at_orig.info)
#        new_at.calc = at_orig.calc
#        if getattr(at_orig, "constraints", None):
#            new_at.set_constraint(at_orig.constraints)
#
#        new_atoms_list.append(new_at)
#
#    return delta, new_atoms_list
