from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

import os
import numpy as np

try:
    from .ase_io import (
        get_energy_ase,
        get_energies_from_atoms_list,
        get_axis_ase,
        get_canonical_ase,
        get_cn_ase,
        get_irmsd_ase,
        get_rmsd_ase,
        sorter_irmsd_ase,
    )
except Exception:  # pragma: no cover
    get_cn_ase = None  # type: ignore
    get_axis_ase = None  # type: ignore
    get_canonical_ase = None  # type: ignore

from .utils import print_array, print_structur, require_ase
from .printouts import print_structure_summary
from ..sorting import first_by_assignment, group_by, sort_by_value

if TYPE_CHECKING:
    from ase import Atoms  # type: ignore

__all__ = ["compute_cn_and_print", "compute_axis_and_print"]


def compute_cn_and_print(atoms_list: Sequence["Atoms"]) -> List[np.ndarray]:
    """Compute coordination numbers for each structure and print them.

    Parameters
    ----------
    atoms_list : list[ase.Atoms]
        Structures to analyze.

    Returns
    -------
    list[np.ndarray]
        One integer array per structure, same order as ``atoms_list``.
    """
    # Ensure ASE is present only when this command is actually invoked
    require_ase()

    results: List[np.ndarray] = []
    for i, atoms in enumerate(atoms_list, start=1):
        if get_cn_ase is not None:
            cn_vec = get_cn_ase(atoms)
        else:
            cn_vec = None
        results.append(cn_vec)
        print_array(f"CN[ structure {i} ] (n={len(atoms)})", cn_vec)
    return results


def compute_axis_and_print(
    atoms_list: Sequence["Atoms"],
) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Compute rotational constants, averge momentum and rotation matrix for
    each structure and prints them.

    Parameters
    ----------
    atoms_list : list[ase.Atoms]
        Structures to analyze.

    Returns
    -------
    list[np.ndarray, np.ndarray, np.ndarray]
        One float array with the 3 rotational constants, one float with the average momentum
        and one float array with the rotation matrix (3, 3) per structure, same order as ``atoms_list``.
    """
    # Ensure ASE is present only when this command is actually invoked
    require_ase()

    results: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for i, atoms in enumerate(atoms_list, start=1):
        if get_axis_ase is not None:
            rot, avmom, evec = get_axis_ase(atoms)
        else:
            rot, avmom, evec = None, None, None
        results.append((rot, avmom, evec))
        print_array(f"Rotational constants (MHz) for structure {i}", rot)
        print(f"Average momentum a.u. (10⁻⁴⁷kg m²) for structure {i}: {avmom}")
        print_array(f"Rotation matrix for structure {i}", evec)
    return results


def compute_canonical_and_print(
    atoms_list: Sequence["Atoms"], heavy=False
) -> List[np.ndarray]:
    """Computes the canonical atom identifiers for each structure and prints
    them.

    Parameters
    ----------
    atoms_list : list[ase.Atoms]
        Structures to analyze.

    Returns
    -------
    list[np.ndarray]
        One integer array with the canonical ranks per structure, same order as ``atoms_list``.
    """
    # Ensure ASE is present only when this command is actually invoked
    require_ase()

    results: List[np.ndarray] = []
    for i, atoms in enumerate(atoms_list, start=1):
        if get_canonical_ase is not None:
            rank = get_canonical_ase(atoms, heavy=heavy)
        else:
            rank = None
        results.append(rank)
        print_array(f"Canonical rank for structure {i}", rank)
    return results


def compute_quaternion_rmsd_and_print(
    atoms_list: Sequence["Atoms"], heavy=False, outfile=None
) -> None:
    """Computes the canonical atom identifiers for a SINGLE PAIR of molecules
    and print the RMSD in Angström between them.

    Parameters
    ----------
    atoms_list : list[ase.Atoms]
        Structures to analyze. Must contain exactly two strucutres

    Returns
    -------
    list[np.ndarray]
        One integer array with the canonical ranks per structure, same order as ``atoms_list``.
    """
    # Ensure ASE is present only when this command is actually invoked
    require_ase()
    from ase.io import write as asewrite

    print("Reference structure:")
    print_structur(atoms_list[0])
    print("Structure to align:")
    print_structur(atoms_list[1])
    if heavy:
        mask0 = atoms_list[0].get_atomic_numbers() > 1
    else:
        mask0 = None
    rmsd, new_atoms, umat = get_rmsd_ase(atoms_list[0], atoms_list[1], mask=mask0)

    if outfile is not None:
        print(f"\nAligned structure written to {outfile}")
        asewrite(outfile, new_atoms)
    else:
        print("Aligned structure:")
        print_structur(new_atoms)

    print_array("\nU matrix (Fortran order)", umat)
    print(f"Cartesian RMSD: {rmsd:.10f} Å")


def compute_irmsd_and_print(
    atoms_list: Sequence["Atoms"], inversion=None, outfile=None
) -> None:
    """Computes the iRMSD between a SINGLE PAIR of molecules and print the
    iRMSD value.

    Parameters
    ----------
    atoms_list : list[ase.Atoms]
        Structures to analyze. Must contain exactly two strucutres
    inversion :
        parameter to instruct inversion in iRMSD routine

    Returns
    -------
    None
    """
    # Ensure ASE is present only when this command is actually invoked
    require_ase()
    from ase.io import write as asewrite

    if inversion is not None:
        print(f"Inversion check: {inversion}\n")

    print("Reference structure:")
    print_structur(atoms_list[0])
    print("Structure to align:")
    print_structur(atoms_list[1])

    if inversion is not None:
        iinversion = {"auto": 0, "on": 1, "off": 2}[inversion]

    irmsd_value, new_atoms_ref, new_atoms_aligned = get_irmsd_ase(
        atoms_list[0], atoms_list[1], iinversion=iinversion
    )

    if outfile is not None:
        print(f"\nAligned reference structure written to {outfile}")
        outfile_ref = outfile
        outfile_ref.stem += "_ref"
        asewrite(outfile_ref, new_atoms_ref)
        print(f"\nAligned probe structure written to {outfile}")
        outfile_aligned = outfile
        outfile_aligned.stem += "_ref"
        asewrite(outfile_aligned, new_atoms_aligned)
    else:
        print("Aligned reference structure:")
        print_structur(new_atoms_ref)
        print()
        print("Aligned probe structure:")
        print_structur(new_atoms_aligned)

    print(f"\niRMSD: {irmsd_value:.10f} Å")


def sort_structures_and_print(
    atoms_list: Sequence["Atoms"],
    rthr: float,
    inversion: str = None,
    allcanon: bool = True,
    printlvl: int = 0,
    outfile: str | None = None,
) -> None:
    """
    Convenience wrapper around presorted_sort_structures_and_print:

    - Analyzes the atoms_list to separate them by composition
    - Sorts by energy if applicable.
    - Calls presorted_sort_structures_and_print for each group

    Parameters
    ----------
    atoms_list : sequence of ase.Atoms
        Input structures.
    rthresh : float
        Distance threshold for sorter_irmsd_ase.
    iinversion : int, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.
    outfile : str or None, optional
        If not None, write all resulting structures to this file
        (e.g. 'sorted.xyz') using ASE's write function.
        Gets automatic name appendage if there are more than one
        type of molecule in the atoms_list
    """
    require_ase()
    from ase.io import write  # local import after require_ase

    # sort the atoms_list by chemical sum formula
    mol_dict = group_by(atoms_list, key=lambda a: a.get_chemical_formula(mode="hill"))

    if len(mol_dict) == 1:
        # Exactly one molecule type
        key, atoms_list = next(iter(mol_dict.items()))
        # Sort by energy (if possible)
        energies = get_energies_from_atoms_list(atoms_list)
        atoms_list, energies = sort_by_value(atoms_list,energies)
        mol_dict[key] = Presorted_sort_structures_and_print(
            atoms_list, rthr, inversion, allcanon, printlvl, outfile
        )
        irmsdvals = np.zeros(len(atoms_list))
        print_structure_summary(key,energies,irmsdvals,max_rows=25)

    else:
        # Multiple molecule types
        for key, atoms_list in mol_dict.items():
            if outfile is not None:
                root, ext = os.path.splitext(outfile)
                outfile_key = f"{root}_{key}{ext}"
            else:
                outfile_key = None
            # Sort by energy (if possible)
            energies = get_energies_from_atoms_list(atoms_list)
            atoms_list, energies = sort_by_value(atoms_list,energies) 
            mol_dict[key] = Presorted_sort_structures_and_print(
                atoms_list, rthr, inversion, allcanon, printlvl, outfile_key
            )


def Presorted_sort_structures_and_print(
    atoms_list: Sequence["Atoms"],
    rthr: float,
    inversion: str = None,
    allcanon: bool = True,
    printlvl: int = 0,
    outfile: str | None = None,
) -> None:
    """
    Convenience wrapper around sorter_irmsd_ase:

    - Calls sorter_irmsd_ase on the given list of ASE Atoms.
    - Prints the resulting groups array.
    - Optionally writes all resulting structures to `outfile` via ASE.

    Parameters
    ----------
    atoms_list : sequence of ase.Atoms
        Input structures.
    rthresh : float
        Distance threshold for sorter_irmsd_ase.
    iinversion : int, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.
    outfile : str or None, optional
        If not None, write all resulting structures to this file
        (e.g. 'sorted.xyz') using ASE's write function.
    """
    require_ase()
    from ase.io import write  # local import after require_ase

    if inversion is not None:
        iinversion = {"auto": 0, "on": 1, "off": 2}[inversion]

    # Call the ASE-level sorter
    groups, new_atoms_list = sorter_irmsd_ase(
        atoms_list=atoms_list,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    # Print groups to screen
    repr = np.max(groups)
    print(f"List of structures was processed: {repr} group{'s' if repr != 1 else ''}.")

    new_atoms_list = first_by_assignment(new_atoms_list, groups)

    # Optionally write all resulting structures to file (e.g. multi-structure XYZ)
    if outfile is not None:
        write(outfile, new_atoms_list)
        print(
            f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile}"
        )

    return new_atoms_list
