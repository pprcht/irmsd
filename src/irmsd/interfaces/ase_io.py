from __future__ import annotations

from typing import List, Sequence, Tuple, overload

import numpy as np

from ..core.molecule import Molecule
from ..utils.utils import require_ase
from .mol_interface import (
    delta_irmsd_list_molecule,
    get_irmsd_molecule,
    get_rmsd_molecule,
    sorter_irmsd_molecule,
    cregen,
    prune,
)

# -------------------------------------------------------------------
# Some I/O
# -------------------------------------------------------------------


def get_energy_ase(atoms):
    """Retrieve the energy associated with an ASE Atoms object.

    This function attempts to extract the energy of the given ASE Atoms object
    through several common avenues, in the following order:
    1. Check if the energy is stored in `atoms.info["energy"]`.
    2. If a calculator is attached, check its `results` dictionary for keys
         "energy", "free_energy", or "enthalpy".
    3. Call `atoms.get_potential_energy()` only if no calculation is needed.

    Parameters
    ----------
    atoms : ase.Atoms
        An ASE Atoms object from which to retrieve the energy.

    Returns
    -------
    float or None
        The energy value if found, otherwise `None`.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If the input is not an ASE Atoms object.
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
    """Convert an ASE `Atoms` object (or a sequence of them) into the internal
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
    """Convert an internal `irmsd.core.Molecule` instance (or a sequence of
    them) into ASE `Atoms` objects.

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


def get_energies_from_atoms_list(atoms_list: Sequence["ase.Atoms"]) -> np.ndarray:
    """Given a list of ASE Atoms objects, call `get_energy_ase(atoms)` for
    each, collect the energies into a float NumPy array, and replace any `None`
    returned by the energy function with 0.0.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        Sequence of ASE Atoms objects.

    Returns
    -------
    np.ndarray
        Float array of energies with shape (N,).
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
    """ASE wrapper for ``get_rmsd_molecule``.

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
    """ASE wrapper for ``get_irmsd_molecule``.

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
        Inversion flag passed through to the backend. (0 = 'auto', 1 = 'on', 2 = 'off')

    Returns
    -------
    irmsd : float
        iRMSD value in Ångström.
    new_atoms1 : ase.Atoms
        New ASE Atoms object corresponding to the transformed first Molecule.
    new_atoms2 : ase.Atoms
        New ASE Atoms object corresponding to the transformed second Molecule.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
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
    rthr: float = 0.125,  # aligned with '--rthr' in src/irmsd/cli.py
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
) -> Tuple[np.ndarray, List["ase.Atoms"]]:
    """ASE wrapper for ``sorter_irmsd_molecule``.

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
        Inversion symmetry flag. (0 = 'auto', 1 = 'on', 2 = 'off')
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.
    ethr : float | None
        Optional energy threshold to accelerate by pre-sorting. In Hartree.
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

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
        ethr=ethr,
        ewin=ewin,
    )

    new_atoms_list = molecule_to_ase(new_mols)

    return groups, new_atoms_list


def delta_irmsd_list_ase(
    atoms_list: Sequence["ase.Atoms"],
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["ase.Atoms"]]:
    """ASE wrapper for ``delta_irmsd_list_molecule``.

    Converts a sequence of ASE ``Atoms`` objects to Molecules, calls
    ``delta_irmsd_list_molecule``, and converts the resulting Molecules
    back to ASE ``Atoms`` objects.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        Sequence of ASE Atoms objects. All must have the same number of atoms.
    iinversion : int, optional
        Inversion symmetry flag. (0 = 'auto', 1 = 'on', 2 = 'off')
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

    new_atoms_list = molecule_to_ase(new_mols)

    return delta, new_atoms_list


def cregen_ase(
    atoms_list: Sequence["ase.Atoms"],
    rthr: float = 0.125,
    ethr: float = 8.0e-5,
    bthr: float = 0.01,
    printlvl: int = 0,
    ewin: float | None = None,
) -> List["ase.Atoms"]:
    """ASE wrapper for ``cregen()`` from mol_interface.

    Converts a sequence of ASE ``Atoms`` objects to Molecules, calls
    ``cregen()``, and converts the resulting Molecules back
    to ASE ``Atoms`` objects.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        Sequence of ASE Atoms objects. All must have the same number of atoms.
    rthr : float
        Distance threshold for the sorter. In Angström.
    ethr : float                                                   
        Energy threshold to accelerate by pre-sorting. In Hartree. 
    bthr : float
        Rotational constant comparison threshold. Relative value (default: 0.01)
    printlvl : int, optional
        Verbosity level.
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    new_atoms_list : list[ase.Atoms]
        New ASE Atoms objects reconstructed from the sorted Molecules.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms_list, (list, tuple)):
        raise TypeError("prune_ase expects a sequence (list/tuple) of ASE Atoms")

    for i, at in enumerate(atoms_list):
        if not isinstance(at, ASEAtoms):
            raise TypeError(
                "prune_ase expects a sequence of ASE Atoms; "
                f"item {i} has type {type(at)}"
            )

    mols = ase_to_molecule(atoms_list)  # returns list[Molecule]

    new_mols = cregen(
        molecule_list=mols,
        rthr=rthr,
        printlvl=printlvl,
        ethr=ethr,
        bthr=bthr,
        ewin=ewin,
    )

    new_atoms_list = molecule_to_ase(new_mols)
    return new_atoms_list


def prune_ase(
    atoms_list: Sequence["ase.Atoms"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
) -> List["ase.Atoms"]:
    """ASE wrapper for ``prune()`` from mol_interface.

    Converts a sequence of ASE ``Atoms`` objects to Molecules, calls
    ``prune()``, and converts the resulting Molecules back
    to ASE ``Atoms`` objects.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        Sequence of ASE Atoms objects. All must have the same number of atoms.
    rthr : float
        Distance threshold for the sorter.
    iinversion : int, optional
        Inversion symmetry flag. (0 = 'auto', 1 = 'on', 2 = 'off')
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.
    ethr : float | None
        Optional energy threshold to accelerate by pre-sorting. In Hartree.
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    new_atoms_list : list[ase.Atoms]
        New ASE Atoms objects reconstructed from the sorted Molecules.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms_list, (list, tuple)):
        raise TypeError("prune_ase expects a sequence (list/tuple) of ASE Atoms")

    for i, at in enumerate(atoms_list):
        if not isinstance(at, ASEAtoms):
            raise TypeError(
                "prune_ase expects a sequence of ASE Atoms; "
                f"item {i} has type {type(at)}"
            )

    mols = ase_to_molecule(atoms_list)  # returns list[Molecule]

    new_mols = prune(
        molecule_list=mols,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
        ethr=ethr,
        ewin=ewin,
    )

    new_atoms_list = molecule_to_ase(new_mols)

    return new_atoms_list
