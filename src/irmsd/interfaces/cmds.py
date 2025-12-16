from __future__ import annotations

import os
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

from ..core import Molecule
from ..sorting import first_by_assignment, group_by, sort_by_value
from ..utils.io import write_structures
from ..utils.printouts import (
    print_conformer_structures,
    print_molecule_summary,
    print_pretty_array,
    print_structure_summary,
)
from .mol_interface import (
    cregen,
    delta_irmsd_list_molecule,
    get_energies_from_molecule_list,
    get_irmsd_molecule,
    get_rmsd_molecule,
    sorter_irmsd_molecule,
)

# ------------------------------------------------------
# CMDs for "prop" runtypes
# ------------------------------------------------------


def compute_cn_and_print(
    molecule_list: Sequence["Molecule"],
    run_multiple: bool = False,
) -> List[np.ndarray]:
    """Compute coordination numbers for each structure and print them.

    Parameters
    ----------
    molecule_list : list[irmsd.Molecule]
        Structures to analyze.

    Returns
    -------
    list[np.ndarray]
        One integer array per structure, same order as ``molecule_list``.
    """

    results: List[np.ndarray] = []
    for i, mol in enumerate(molecule_list, start=1):
        cn_vec = mol.get_cn()
        results.append(cn_vec)
    if not run_multiple:
        print_molecule_summary(molecule_list, **{"CN": results})
    return results


def compute_axis_and_print(
    molecule_list: Sequence["Molecule"],
    run_multiple: bool = False,
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Compute rotational constants, averge momentum and rotation matrix for
    each structure and prints them.

    Parameters
    ----------
    molecule_list : list[irmsd.Molecule]
        Structures to analyze.

    Returns
    -------
    list[np.ndarray, np.ndarray, np.ndarray]
        One float array with the 3 rotational constants, one float with the average momentum
        and one float array with the rotation matrix (3, 3) per structure, same order as ``molecule_list``.
    """

    results: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for i, mol in enumerate(molecule_list, start=1):
        axd = dict()
        rot, avmom, evec = mol.get_axis()
        axd["Rotational constants (MHz)"] = rot
        axd["Rotation matrix"] = evec
        results.append(axd)
    if not run_multiple:
        print_molecule_summary(molecule_list, axis=results)
    return results


def compute_canonical_and_print(
    molecule_list: Sequence["Molecule"],
    heavy: bool = False,
    run_multiple: bool = False,
) -> List[np.ndarray]:
    """Computes the canonical atom identifiers for each structure and prints
    them.

    Parameters
    ----------
    molecule_list : list[irmsd.Molecule]
        Structures to analyze.
    heavy: bool
        Consider only heavy atoms

    Returns
    -------
    list[np.ndarray]
        One integer array with the canonical ranks per structure, same order as ``molecule_list``.
    """

    results: List[np.ndarray] = []
    for i, mol in enumerate(molecule_list, start=1):
        rank = mol.get_canonical(heavy=heavy)
        results.append(rank)
    if not run_multiple:
        print_molecule_summary(molecule_list, **{"Canonical ID": results})
    return results


# ------------------------------------------------------
# CMDs for "compare" runtypes
# ------------------------------------------------------


def get_ref_and_align_molecules(
    molecule_list: Sequence["Molecule"],
    idx_ref: int,
    idx_align: int,
) -> Tuple["Molecule", "Molecule"]:
    n_molecules = len(molecule_list)
    if n_molecules < 2:
        raise ValueError("At least two structures are required to compute iRMSD.")
    if n_molecules > 2:
        print(
            f"{n_molecules} structures were provided, comparing only structures {idx_ref} and {idx_align}."
        )
    if (
        idx_ref > n_molecules - 1
        or idx_ref < 0
        or idx_align > n_molecules - 1
        or idx_align < 0
    ):
        raise IndexError(
            f"Reference or align index is out of range. Max index is {n_molecules - 1}, got ref: {idx_ref}, align: {idx_align}."
        )
    if idx_ref == idx_align:
        raise ValueError(
            "Reference and align indices must be different. Both are {idx_ref}."
        )
    mol_ref = molecule_list[idx_ref]
    mol_align = molecule_list[idx_align]
    print("Input structures:")
    print_conformer_structures(
        mol_ref,
        mol_align,
        labels=["Reference", "Probe"],
    )
    print()
    return mol_ref, mol_align


def compute_quaternion_rmsd_and_print(
    molecule_list: Sequence["Molecule"],
    heavy=False,
    outfile=None,
    idx_ref=0,
    idx_align=1,
) -> None:
    """Computes the canonical atom identifiers for a SINGLE PAIR of molecules
    and print the RMSD in Angström between them.

    Parameters
    ----------
    molecule_list : list[irmsd.Molecule]
        Structures to analyze. Must contain exactly two strucutres
    heavy : bool, optional
        If True, only heavy atoms are considered in the RMSD calculation.
    outfile : str or None, optional
        If not None, write the aligned structure to this file.
    idx_ref : int, optional
        Index of the reference structure in molecule_list (default: 0).
    idx_align : int, optional
        Index of the structure to align in molecule_list (default: 1).

    Returns
    -------
    list[np.ndarray]
        One integer array with the canonical ranks per structure, same order as ``molecule_list``.
    """

    mol_ref, mol_align = get_ref_and_align_molecules(molecule_list, idx_ref, idx_align)
    if heavy:
        mask0 = mol_align.get_atomic_numbers() > 1
    else:
        mask0 = None
    rmsd, new_atoms, umat = get_rmsd_molecule(mol_ref, mol_align, mask=mask0)

    if outfile is not None:
        print(f"\nAligned structure written to {outfile}")
        write_structures(outfile, new_atoms)
    else:
        print("Aligned structures:")
        print_conformer_structures(mol_ref, new_atoms, labels=["Reference", "Aligned"])

    print_pretty_array("\nU matrix (Fortran order)", umat)
    print(f"Cartesian RMSD: {rmsd:.10f} Å")


def compute_irmsd_and_print(
    molecule_list: Sequence["Molecule"],
    inversion=None,
    outfile=None,
    idx_ref=0,
    idx_align=1,
) -> None:
    """Computes the iRMSD between a SINGLE PAIR of molecules and print the
    iRMSD value.

    Parameters
    ----------
    molecule_list : list[irmsd.Molecule]
        Structures to analyze. Must contain exactly two strucutres
    inversion :
        parameter to instruct inversion in iRMSD routine
    outfile : str or None, optional
        If not None, write the aligned structures to this file.
    idx_ref : int, optional
        Index of the reference structure in molecule_list (default: 0).
    idx_align : int, optional
        Index of the structure to align in molecule_list (default: 1).

    Returns
    -------
    None
    """
    mol_ref, mol_align = get_ref_and_align_molecules(molecule_list, idx_ref, idx_align)

    if inversion is not None:
        print(f"Inversion check: {inversion}\n")

    if inversion is not None:
        iinversion = {"auto": 0, "on": 1, "off": 2}[inversion]

    irmsd_value, new_atoms_ref, new_atoms_aligned = get_irmsd_molecule(
        mol_ref, mol_align, iinversion=iinversion
    )

    if outfile is not None:
        print(f"\nAligned reference structure written to {outfile}")
        outfile_ref = Path(outfile)
        outfile_ref = outfile_ref.with_stem(outfile_ref.stem + "_ref")
        write_structures(outfile_ref, new_atoms_ref)
        print(f"\nAligned probe structure written to {outfile}")
        outfile_aligned = Path(outfile)
        outfile_aligned = outfile_aligned.with_stem(outfile_aligned.stem + "_aligned")
        write_structures(outfile_aligned, new_atoms_aligned)
    else:
        print("Aligned structures:")
        print_conformer_structures(
            new_atoms_ref, new_atoms_aligned, labels=["Reference", "Aligned"]
        )

    print(f"\niRMSD: {irmsd_value:.10f} Å")


# ------------------------------------------------------
# CMDs for "sort"/"prune" runtypes
# ------------------------------------------------------


def sort_structures_and_print(
    molecule_list: Sequence["Molecule"],
    rthr: float,
    inversion: str = None,
    allcanon: bool = True,
    printlvl: int = 0,
    maxprint: int = 25,
    ethr: float | None = None,
    ewin: float | None = None,
    outfile: str | None = None,
) -> None:
    """
    Convenience wrapper around presorted_sort_structures_and_print:

    - Analyzes the molecule_list to separate them by composition
    - Sorts by energy if applicable.
    - Calls presorted_sort_structures_and_print for each group

    Parameters
    ----------
    molecule_list : sequence of irmsd.Molecule
        Input structures.
    rthr : float | None
        Distance threshold for sorter_irmsd_molecule.
    inversion : str, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.
    maxprint : int, optional
        Max number of lines to print for each structure result table
    ethr : float | None
        Optional inter-conformer energy threshold for more efficient presorting
    ewin : float | None
        Optional energy window to limit ensemble size around lowest energy structure
    outfile : str or None, optional
        If not None, write all resulting structures to this file
        (e.g. 'sorted.xyz') using a write function.
        Gets automatic name appendage if there are more than one
        type of molecule in the molecule_list
    """

    if inversion is not None:
        iinversion = {"auto": 0, "on": 1, "off": 2}[inversion]

    # sort the molecule_list by chemical sum formula
    mol_dict = group_by(
        molecule_list, key=lambda a: a.get_chemical_formula(mode="hill")
    )

    if len(mol_dict) == 1:
        # Exactly one molecule type
        key, molecule_list = next(iter(mol_dict.items()))
        # Sort by energy (if possible)
        energies = get_energies_from_molecule_list(molecule_list)
        molecule_list, energies = sort_by_value(molecule_list, energies)
        print()
        mol_dict[key] = Presorted_sort_structures_and_print(
            molecule_list,
            rthr,
            iinversion,
            allcanon,
            printlvl,
            ethr=ethr,
            ewin=ewin,
            outfile=outfile,
        )
        irmsdvals, _ = delta_irmsd_list_molecule(
            mol_dict[key], iinversion, allcanon=True, printlvl=0
        )
        energies = get_energies_from_molecule_list(mol_dict[key])
        print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)

        # Optionally write all resulting structures to file (e.g. multi-structure XYZ)
        if outfile is not None:
            write_structures(outfile, mol_dict[key])
            repr = len(mol_dict[key])
            if printlvl > 0:
                print(
                    f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile}"
                )

    else:
        # Multiple molecule types
        for key, molecule_list in mol_dict.items():
            if outfile is not None:
                root, ext = os.path.splitext(outfile)
                outfile_key = f"{root}_{key}{ext}"
            else:
                outfile_key = None
            # Sort by energy (if possible)
            energies = get_energies_from_molecule_list(molecule_list)
            molecule_list, energies = sort_by_value(molecule_list, energies)
            print()
            mol_dict[key] = Presorted_sort_structures_and_print(
                molecule_list,
                rthr,
                iinversion,
                allcanon,
                printlvl,
                ethr=ethr,
                ewin=ewin,
                outfile=outfile_key,
            )
            irmsdvals, _ = delta_irmsd_list_molecule(
                mol_dict[key], iinversion, allcanon=True, printlvl=0
            )
            energies = get_energies_from_molecule_list(mol_dict[key])
            print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)

            # Optionally write all resulting structures to file (e.g. multi-structure XYZ)
            if outfile_key is not None:
                write_structures(outfile_key, mol_dict[key])
                repr = len(mol_dict[key])
                if printlvl > 0:
                    print(
                        f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile_key}"
                    )


def Presorted_sort_structures_and_print(
    molecule_list: Sequence["Molecule"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
    outfile: str | None = None,
) -> None:
    """
    Convenience wrapper around sorter_irmsd_molecule:

    - Calls sorter_irmsd_molecule on the given list of ASE Molecule.
    - Prints the resulting groups array.
    - Optionally writes all resulting structures to `outfile` via ASE.

    Parameters
    ----------
    molecule_list : sequence of irmsd.Molecule
        Input structures.
    rthresh : float
        Distance threshold for sorter_irmsd_molecule.
    iinversion : int, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.
    ethr : float | None
        Optional inter-conformer energy threshold for more efficient presorting
    ewin : float | None
        Optional energy window to limit ensemble size around lowest energy structure
    outfile : str or None, optional
        If not None, write all resulting structures to this file
        (e.g. 'sorted.xyz') using a write function.
    """

    # Call the ASE-level sorter
    groups, new_molecule_list = sorter_irmsd_molecule(
        molecule_list=molecule_list,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
        ethr=ethr,
        ewin=ewin,
    )

    new_molecule_list = first_by_assignment(new_molecule_list, groups)

    return new_molecule_list


def sort_get_delta_irmsd_and_print(
    molecule_list: Sequence["Molecule"],
    inversion: str = None,
    allcanon: bool = True,
    printlvl: int = 0,
    maxprint: int = 25,
    outfile: str | None = None,
) -> None:
    """
    Convenience wrapper around presorted_sort_structures_and_print:

    - Analyzes the molecule_list to separate them by composition
    - Sorts by energy if applicable.
    - Calculates iRMSD between structures x_i and x_i-1

    Parameters
    ----------
    molecule_list : sequence of irmsd.Molecule
        Input structures.
    inversion : str, optional
        Inversion symmetry flag, passed through.
    allcanon : bool, optional
        Canonicalization flag, passed through.
    printlvl : int, optional
        Verbosity level, passed through.
    maxprint : int, optional
        Max number of lines to print for each structure result table
    outfile : str or None, optional
        If not None, write all resulting structures to this file
        (e.g. 'sorted.xyz') using a write function.
        Gets automatic name appendage if there are more than one
        type of molecule in the molecule_list
    """

    if inversion is not None:
        iinversion = {"auto": 0, "on": 1, "off": 2}[inversion]

    # sort the molecule_list by chemical sum formula
    mol_dict = group_by(
        molecule_list, key=lambda a: a.get_chemical_formula(mode="hill")
    )

    if len(mol_dict) == 1:
        # Exactly one molecule type
        key, molecule_list = next(iter(mol_dict.items()))
        # Sort by energy (if possible)
        energies = get_energies_from_molecule_list(molecule_list)
        molecule_list, energies = sort_by_value(molecule_list, energies)
        print()
        irmsdvals, mol_dict[key] = delta_irmsd_list_molecule(
            molecule_list, iinversion, allcanon, printlvl
        )
        energies = get_energies_from_molecule_list(mol_dict[key])
        print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)

    else:
        # Multiple molecule types
        for key, molecule_list in mol_dict.items():
            if outfile is not None:
                root, ext = os.path.splitext(outfile)
                outfile_key = f"{root}_{key}{ext}"
            else:
                outfile_key = None
            # Sort by energy (if possible)
            energies = get_energies_from_molecule_list(molecule_list)
            molecule_list, energies = sort_by_value(molecule_list, energies)
            print()
            irmsdvals, mol_dict[key] = delta_irmsd_list_molecule(
                molecule_list, iinversion, allcanon, printlvl
            )
            energies = get_energies_from_molecule_list(mol_dict[key])
            print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)


def run_cregen_and_print(
    molecule_list: Sequence["Molecule"],
    rthr: float,
    ethr: float,
    bthr: float,
    ewin: float | None = None,
    printlvl: int = 0,
    maxprint: int = 25,
    outfile: str | None = None,
) -> None:
    """Convenience wrapper around cregen() from mol_interface. Splits according
    to sum formula, if necessary.

    Parameters
    ----------
    molecule_list : sequence of irmsd.Molecule
        Input structures.
    rthr: float
        RMSD thershold for conformer identification
    ethr: float
        Energy threshold for conformer identification
    bthr: float
        Rotational constant threshold for conformer identification
    printlvl : int, optional
        Verbosity level, passed through.
    maxprint : int, optional
        Max number of lines to print for each structure result table
    outfile : str or None, optional
        If not None, write all resulting structures to this file
        (e.g. 'sorted.xyz') using a write function.
        Gets automatic name appendage if there are more than one
        type of molecule in the molecule_list
    """

    # sort the molecule_list by chemical sum formula
    mol_dict = group_by(
        molecule_list, key=lambda a: a.get_chemical_formula(mode="hill")
    )

    if len(mol_dict) == 1:
        # Exactly one molecule type
        key, molecule_list = next(iter(mol_dict.items()))
        # Sort by energy (if possible)
        print()
        mol_dict[key] = cregen(
            molecule_list, rthr, ethr, bthr, ewin=ewin, printlvl=printlvl
        )

        # allcanon can be False here because CREGEN requires same atom order.
        irmsdvals, _ = delta_irmsd_list_molecule(
            mol_dict[key], iinversion=0, allcanon=False, printlvl=0
        )
        energies = get_energies_from_molecule_list(mol_dict[key])

        print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)

        if outfile is not None:
            write_structures(outfile, mol_dict[key])
            if printlvl > 0:
                repr = len(mol_dict[key])
                print(
                    f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile}"
                )

    else:
        # Multiple molecule types
        for key, molecule_list in mol_dict.items():
            if outfile is not None:
                root, ext = os.path.splitext(outfile)
                outfile_key = f"{root}_{key}{ext}"
            else:
                outfile_key = None
            # Sort by energy (if possible)
            print()
            mol_dict[key] = cregen(
                molecule_list, rthr, ethr, bthr, ewin=ewin, printlvl=printlvl
            )

            # allcanon can be False here because CREGEN requires same atom order.
            irmsdvals, _ = delta_irmsd_list_molecule(
                mol_dict[key], iinversion=0, allcanon=False, printlvl=0
            )
            energies = get_energies_from_molecule_list(mol_dict[key])
            print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)

            if outfile_key is not None:
                write_structures(outfile_key, mol_dict[key])

            if printlvl > 0:
                repr = len(mol_dict[key])
                print(
                    f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile_key}"
                )
