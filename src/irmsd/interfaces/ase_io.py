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

def get_energy_ase(atoms):
    """Return the energy already stored on an ASE Atoms object, or None.

    Looks in ``atoms.info["energy"]``, then the calculator's ``results``
    ("energy", "free_energy", "enthalpy"), then ``get_potential_energy()``
    only if that needs no new calculation. The value is in ASE units (eV).

    Parameters
    ----------
    atoms : ase.Atoms

    Returns
    -------
    float or None

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If `atoms` is not an ASE Atoms object.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_energy_ase expects an ase.Atoms object")

    E = atoms.info.get("energy")
    if isinstance(E, (int, float)):
        return float(E)

    calc = getattr(atoms, "calc", None)
    if calc is None:
        return None

    results = getattr(calc, "results", None)
    if isinstance(results, dict):
        for key in ("energy", "free_energy", "enthalpy"):
            val = results.get(key)
            if isinstance(val, (int, float)):
                return float(val)

    try:
        if hasattr(calc, "calculation_required"):
            if calc.calculation_required(atoms):
                return None
        return float(atoms.get_potential_energy())
    except Exception:
        pass

    return None


# Energy units: ASE uses eV, Molecule uses Hartree.

#: Info key (case-insensitive) carrying an explicit energy-unit declaration.
_ENERGY_UNITS_KEY = "energy_units"
#: Values (case-insensitive) marking an energy as already in Hartree.
_HARTREE_UNIT_ALIASES = frozenset(
    {"hartree", "hartrees", "ha", "au", "a.u.", "eh", "e_h", "atomic"}
)


def _declared_hartree(info) -> bool:
    """True if ``info`` carries ``energy_units=Hartree`` (as our extxyz writers emit)."""
    if not isinstance(info, dict):
        return False
    for key, val in info.items():
        if str(key).lower() == _ENERGY_UNITS_KEY and isinstance(val, str):
            return val.strip().lower() in _HARTREE_UNIT_ALIASES
    return False


def _strip_energy_units(info: dict) -> dict:
    """Copy of ``info`` without ``energy_units``; avoids double conversion on round-trip."""
    return {k: v for k, v in info.items() if str(k).lower() != _ENERGY_UNITS_KEY}


def _ase_energy_to_hartree(atoms) -> float | None:
    """Energy of `atoms` in Hartree; eV unless an ``energy_units=Hartree`` marker says otherwise."""
    e = get_energy_ase(atoms)
    if e is None:
        return None
    if _declared_hartree(getattr(atoms, "info", {})):
        return float(e)
    ase = require_ase()
    return float(e) / ase.units.Hartree  # type: ignore[attr-defined]


def _hartree_to_ev(energy: float | None) -> float | None:
    """Hartree to eV, passing None through."""
    if energy is None:
        return None
    ase = require_ase()
    return float(energy) * ase.units.Hartree  # type: ignore[attr-defined]


@overload
def ase_to_molecule(atoms: "ase.Atoms") -> Molecule: ...
@overload
def ase_to_molecule(atoms: Sequence["ase.Atoms"]) -> list[Molecule]: ...


def ase_to_molecule(atoms):
    """Convert ASE Atoms (single or sequence) to `irmsd.core.Molecule`.

    Triggers no calculator evaluation and does not modify the input. Energies
    are converted from eV to Hartree unless ``info`` declares
    ``energy_units=Hartree``; the marker is dropped from the result's ``info``.
    A per-atom ``canonical_id`` array is carried over as ``Molecule.ids``.

    Parameters
    ----------
    atoms : ase.Atoms or Sequence[ase.Atoms]

    Returns
    -------
    Molecule or list[Molecule]
        A list, in input order, if `atoms` is a sequence.

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

        energy = _ase_energy_to_hartree(a)
        info = _strip_energy_units(dict(getattr(a, "info", {})))

        # e.g. from an extxyz canonical_id:I:1 column
        ids = None
        try:
            if a.has("canonical_id"):
                ids = np.asarray(a.get_array("canonical_id"), dtype=np.int32)
        except Exception:
            ids = None

        return Molecule(
            symbols=symbols,
            positions=positions,
            energy=energy,
            info=info,
            cell=cell_array,
            pbc=pbc,
            ids=ids,
        )

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
    """Convert Molecule(s) to ASE Atoms without attaching a calculator.

    The Hartree energy is written to ``info["energy"]`` in eV unless
    ``info`` already has an ``energy`` entry; any ``energy_units`` marker is
    dropped. ``Molecule.ids`` becomes a per-atom ``canonical_id`` array.
    The returned objects do not share state with the input.

    Parameters
    ----------
    molecules : Molecule or Sequence[Molecule]

    Returns
    -------
    ase.Atoms or list[ase.Atoms]
        A list, in input order, if `molecules` is a sequence.

    Raises
    ------
    RuntimeError
        If ASE is not installed.
    TypeError
        If the input is neither a Molecule nor a sequence of Molecules.
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

        cell = None
        if mol.cell is not None:
            cell_arr = np.asarray(mol.cell, dtype=float)
            if cell_arr.shape == (3, 3):
                cell = cell_arr

        pbc = mol.pbc if mol.pbc is not None else False

        # stale marker would mislabel the eV energy written below
        info = _strip_energy_units(dict(mol.info))

        if mol.energy is not None and "energy" not in info:
            info["energy"] = _hartree_to_ev(mol.energy)

        atoms = ASEAtoms(
            symbols=symbols,
            positions=positions,
            cell=cell,
            pbc=pbc,
            info=info,
        )

        # ase.io.write then emits a canonical_id column in extxyz Properties
        if mol.ids is not None:
            atoms.set_array("canonical_id", np.asarray(mol.ids, dtype=int))

        return atoms

    if isinstance(molecules, Molecule):
        return _one(molecules)

    try:
        return [_one(m) for m in molecules]
    except TypeError as exc:
        raise TypeError(
            "molecule_to_ase expects either a single Molecule or a sequence of Molecule objects"
        ) from exc


def get_energies_from_atoms_list(atoms_list: Sequence["ase.Atoms"]) -> np.ndarray:
    """Energies of `atoms_list` in Hartree, shape (N,); 0.0 where none is available.

    Units follow :func:`_ase_energy_to_hartree`.
    """
    energies = []
    for atoms in atoms_list:
        e = _ase_energy_to_hartree(atoms)
        energies.append(0.0 if e is None else float(e))
    return np.array(energies, dtype=float)


def get_cn_ase(atoms) -> np.ndarray:
    """Coordination numbers via `Molecule.get_cn()`, shape (N,).

    Parameters
    ----------
    atoms : ase.Atoms
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_cn_ase expects a single ASE Atoms object")

    mol: Molecule = ase_to_molecule(atoms)
    return mol.get_cn()


def get_axis_ase(atoms) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Principal-axis data via `Molecule.get_axis()`.

    Parameters
    ----------
    atoms : ase.Atoms

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        rot_constants_MHz, avg_momentum_au, rotation_matrix
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_axis_ase expects a single ASE Atoms object")

    mol: Molecule = ase_to_molecule(atoms)
    return mol.get_axis()


def get_point_group_ase(atoms, **settings) -> str | None:
    """Schoenflies point group via `Molecule.get_point_group()`.

    Parameters
    ----------
    atoms : ase.Atoms
    **settings
        Analyzer settings (threshold, primary_threshold, max_axis_order,
        max_opt_cycles, max_atoms), see :func:`irmsd.get_point_group`.

    Returns
    -------
    str or None
        None if the analysis was skipped.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms, ASEAtoms):
        raise TypeError("get_point_group_ase expects a single ASE Atoms object")

    mol: Molecule = ase_to_molecule(atoms)
    return mol.get_point_group(**settings)


def get_canonical_ase(
    atoms,
    wbo: np.ndarray | None = None,
    invtype: str = "apsp+",
    heavy: bool = False,
) -> np.ndarray:
    """Canonical ranks / invariants via `Molecule.get_canonical()`.

    Parameters
    ----------
    atoms : ase.Atoms
    wbo : np.ndarray or None, optional
        Wiberg bond order matrix or similar, forwarded to the Fortran backend.
    invtype : str, optional
        Invariant type selector forwarded to the backend.
    heavy : bool, optional
        Restrict invariants to heavy atoms.

    Returns
    -------
    np.ndarray

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

    mol: Molecule = ase_to_molecule(atoms)
    return mol.get_canonical(wbo=wbo, invtype=invtype, heavy=heavy)


def get_rmsd_ase(atoms1, atoms2, mask=None) -> Tuple[float, "ase.Atoms", np.ndarray]:
    """ASE wrapper for ``get_rmsd_molecule``.

    Parameters
    ----------
    atoms1 : ase.Atoms
        Reference structure.
    atoms2 : ase.Atoms
        Structure aligned onto `atoms1`.
    mask : array-like of bool, optional
        Atoms of the first structure that enter the RMSD.

    Returns
    -------
    rmsd : float
        In Angstrom.
    new_atoms2 : ase.Atoms
        New object with coordinates aligned to `atoms1`.
    rotation_matrix : np.ndarray
        3x3 rotation used for the alignment.

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

    mol1, mol2 = ase_to_molecule([atoms1, atoms2])

    rmsd, new_mol2, umat = get_rmsd_molecule(mol1, mol2, mask=mask)
    new_atoms2 = molecule_to_ase(new_mol2)

    return rmsd, new_atoms2, umat


def get_irmsd_ase(
    atoms1,
    atoms2,
    iinversion: int = 0,
) -> Tuple[float, "ase.Atoms", "ase.Atoms"]:
    """ASE wrapper for ``get_irmsd_molecule``.

    Parameters
    ----------
    atoms1, atoms2 : ase.Atoms
    iinversion : int, optional
        0 = 'auto', 1 = 'on', 2 = 'off'.

    Returns
    -------
    irmsd : float
        In Angstrom.
    new_atoms1, new_atoms2 : ase.Atoms
        New objects for the transformed structures.

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

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        List or tuple; all structures must have the same atom count.
    rthr : float
        Distance threshold for the sorter.
    iinversion : int, optional
        0 = 'auto', 1 = 'on', 2 = 'off'.
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.
    ethr : float or None
        Energy threshold for pre-sorting, in Hartree.
    ewin : float or None
        Energy window above the lowest structure, in Hartree.

    Returns
    -------
    groups : np.ndarray
        Integer group index per structure.
    new_atoms_list : list[ase.Atoms]
        Rebuilt from the sorted Molecules.
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

    mols = ase_to_molecule(atoms_list)

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

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        List or tuple; all structures must have the same atom count.
    iinversion : int, optional
        0 = 'auto', 1 = 'on', 2 = 'off'.
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.

    Returns
    -------
    delta : np.ndarray
        Float array from the backend; see ``delta_irmsd_list``.
    new_atoms_list : list[ase.Atoms]
        Rebuilt from the transformed Molecules.
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
    """ASE wrapper for ``cregen()``.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        List or tuple; all structures must have the same atom count.
    rthr : float
        Distance threshold in Angstrom.
    ethr : float
        Energy threshold for pre-sorting, in Hartree.
    bthr : float
        Relative threshold for comparing rotational constants.
    printlvl : int, optional
        Verbosity level.
    ewin : float or None
        Energy window above the lowest structure, in Hartree.

    Returns
    -------
    list[ase.Atoms]
        Rebuilt from the sorted Molecules.
    """
    ase = require_ase()
    ASEAtoms = ase.Atoms  # type: ignore[attr-defined]

    if not isinstance(atoms_list, (list, tuple)):
        raise TypeError("cregen_ase expects a sequence (list/tuple) of ASE Atoms")

    for i, at in enumerate(atoms_list):
        if not isinstance(at, ASEAtoms):
            raise TypeError(
                "cregen_ase expects a sequence of ASE Atoms; "
                f"item {i} has type {type(at)}"
            )

    mols = ase_to_molecule(atoms_list)

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
    """ASE wrapper for ``prune()``.

    Parameters
    ----------
    atoms_list : Sequence[ase.Atoms]
        List or tuple; all structures must have the same atom count.
    rthr : float
        Distance threshold for the sorter.
    iinversion : int, optional
        0 = 'auto', 1 = 'on', 2 = 'off'.
    allcanon : bool, optional
        Canonicalization flag.
    printlvl : int, optional
        Verbosity level.
    ethr : float or None
        Energy threshold for pre-sorting, in Hartree.
    ewin : float or None
        Energy window above the lowest structure, in Hartree.

    Returns
    -------
    list[ase.Atoms]
        Rebuilt from the sorted Molecules.
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

    mols = ase_to_molecule(atoms_list)

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
