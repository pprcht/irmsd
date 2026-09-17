"""CLI command implementations. With ``run_multiple=True`` nothing is printed."""

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

_INVERSION_CODES = {"auto": 0, "on": 1, "off": 2}

# CMDs for "prop" runtypes


def compute_cn_and_print(
    molecule_list: Sequence["Molecule"],
    run_multiple: bool = False,
) -> List[np.ndarray]:
    """Compute and print coordination numbers; one integer array per structure."""

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
) -> List[dict]:
    """Compute and print rotational constants and principal axes.

    Returns
    -------
    list[dict]
        Per structure, "Rotational constants (MHz)" (3,) and
        "Rotation matrix" (3, 3).
    """

    results: List[dict] = []
    for i, mol in enumerate(molecule_list, start=1):
        axd = dict()
        rot, avmom, evec = mol.get_axis()
        axd["Rotational constants (MHz)"] = rot
        axd["Rotation matrix"] = evec
        results.append(axd)
    if not run_multiple:
        print_molecule_summary(molecule_list, axis=results)
    return results


def _summarize_operations(symbol: str, ops) -> str:
    """Class-style summary, e.g. "E, 8 C3, 3 C2, 6 S4, 6 sigma" for Td.
    Counted per element type and power, not conjugacy class: D4h gives "5 C2"."""
    counts: dict[tuple, int] = {}
    for op in ops:
        rank = {"E": 0, "C": 1, "i": 2, "S": 3, "sigma": 4}[op.kind]
        name, k = op.kind, 0
        if op.kind in ("C", "S"):
            # X_n^k and its inverse X_n^(p-k) are listed as one class
            period = 2 * op.order if op.kind == "S" and op.order % 2 else op.order
            k = min(op.power, period - op.power)
            name = f"{op.kind}{op.order}" + (f"^{k}" if k > 1 else "")
        key = (rank, -op.order, k, name)
        counts[key] = counts.get(key, 0) + 1
    parts = [
        name if n == 1 else f"{n} {name}"
        for (_, _, _, name), n in sorted(counts.items())
    ]
    if symbol in ("Cinfv", "Dinfh"):
        parts.insert(1, "Cinf (all rotations about the molecular axis)")
    return ", ".join(parts)


def compute_symmetry_and_print(
    molecule_list: Sequence["Molecule"],
    run_multiple: bool = False,
    **settings,
) -> List[dict]:
    """Determine and print the Schoenflies point group and operations.

    Parameters
    ----------
    **settings
        threshold, primary_threshold, max_axis_order, max_opt_cycles,
        max_atoms; see :func:`irmsd.get_point_group`.

    Returns
    -------
    list[dict]
        Per structure, the printed "Point group" and "Symmetry operations"
        strings and the ``irmsd.SymmetryOperation`` list under "operations".
        Structures skipped for size carry only a "skipped" point group.
    """
    results: List[dict] = []
    for mol in molecule_list:
        symbol, ops = mol.get_symmetry_operations(**settings)
        if symbol is None:
            results.append({"Point group": "skipped (too many atoms)"})
            continue
        results.append(
            {
                "Point group": symbol,
                "Symmetry operations": _summarize_operations(symbol, ops),
                "operations": ops,
            }
        )
    if not run_multiple:
        print_molecule_summary(molecule_list, symmetry=results)
    return results


def compute_canonical_and_print(
    molecule_list: Sequence["Molecule"],
    heavy: bool = False,
    run_multiple: bool = False,
) -> List[np.ndarray]:
    """Compute and print canonical atom ranks (heavy atoms only if ``heavy``);
    one integer array per structure."""

    results: List[np.ndarray] = []
    for i, mol in enumerate(molecule_list, start=1):
        rank = mol.get_canonical(heavy=heavy)
        results.append(rank)
    if not run_multiple:
        print_molecule_summary(molecule_list, **{"Canonical ID": results})
    return results


# CMDs for "compare" runtypes


def get_ref_and_align_molecules(
    molecule_list: Sequence["Molecule"],
    idx_ref: int,
    idx_align: int,
) -> Tuple["Molecule", "Molecule"]:
    """Select and print the reference/probe pair; raise on bad indices."""
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
            f"Reference and align indices must be different. Both are {idx_ref}."
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
    """Align one pair and print the Cartesian RMSD (Angstrom) and U matrix.

    Parameters
    ----------
    heavy : bool
        Restrict the RMSD to heavy atoms.
    outfile : str or None
        Write the aligned probe here instead of printing it.
    idx_ref, idx_align : int
        Reference and probe indices in ``molecule_list``.
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
    """Align one pair and print the iRMSD (Angstrom).

    Parameters
    ----------
    inversion : {"auto", "on", "off"}
        Inversion handling in the iRMSD routine.
    outfile : str or None
        Write the aligned pair to ``<stem>_ref`` and ``<stem>_aligned``
        instead of printing it.
    idx_ref, idx_align : int
        Reference and probe indices in ``molecule_list``.
    """
    mol_ref, mol_align = get_ref_and_align_molecules(molecule_list, idx_ref, idx_align)

    if inversion is not None:
        print(f"Inversion check: {inversion}\n")

    iinversion = _INVERSION_CODES[inversion or "auto"]

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


# CMDs for "sort"/"prune" runtypes


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
    """Split by sum formula, energy-sort, prune each group with the iRMSD
    sorter, and print a summary per group.

    Parameters
    ----------
    rthr : float
        Distance threshold for the sorter.
    inversion : {"auto", "on", "off"}
    maxprint : int
        Max rows per printed result table.
    ethr : float or None
        Inter-conformer energy threshold for presorting.
    ewin : float or None
        Energy window around the lowest-energy structure.
    outfile : str or None
        Write the representatives here; with several formulas, one
        ``<root>_<formula><ext>`` file each.
    """

    iinversion = _INVERSION_CODES[inversion or "auto"]

    mol_dict = group_by(
        molecule_list, key=lambda a: a.get_chemical_formula(mode="hill")
    )

    if len(mol_dict) == 1:
        key, molecule_list = next(iter(mol_dict.items()))
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

        if outfile is not None:
            write_structures(outfile, mol_dict[key])
            repr = len(mol_dict[key])
            if printlvl > 0:
                print(
                    f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile}"
                )

    else:
        for key, molecule_list in mol_dict.items():
            if outfile is not None:
                root, ext = os.path.splitext(outfile)
                outfile_key = f"{root}_{key}{ext}"
            else:
                outfile_key = None
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
) -> List["Molecule"]:
    """Run the iRMSD sorter on one presorted group; return the first structure
    of each resulting group.

    Parameters as in :func:`sort_structures_and_print`, except ``iinversion``
    is already mapped to 0/1/2 (auto/on/off). ``outfile`` is unused.
    """

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
    """Split by sum formula, energy-sort, and print the iRMSD between
    consecutive structures of each group.

    Parameters
    ----------
    inversion : {"auto", "on", "off"}
    maxprint : int
        Max rows per printed result table.
    outfile : str or None
        Unused; no structures are written.
    """

    iinversion = _INVERSION_CODES[inversion or "auto"]

    mol_dict = group_by(
        molecule_list, key=lambda a: a.get_chemical_formula(mode="hill")
    )

    if len(mol_dict) == 1:
        key, molecule_list = next(iter(mol_dict.items()))
        energies = get_energies_from_molecule_list(molecule_list)
        molecule_list, energies = sort_by_value(molecule_list, energies)
        print()
        irmsdvals, mol_dict[key] = delta_irmsd_list_molecule(
            molecule_list, iinversion, allcanon, printlvl
        )
        energies = get_energies_from_molecule_list(mol_dict[key])
        print_structure_summary(key, energies, irmsdvals, max_rows=maxprint)

    else:
        for key, molecule_list in mol_dict.items():
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
    """Run CREGEN per sum-formula group and print a summary per group.

    Parameters
    ----------
    rthr, ethr, bthr : float
        RMSD, energy, and rotational-constant thresholds for conformer
        identification.
    maxprint : int
        Max rows per printed result table.
    outfile : str or None
        Write the representatives here; with several formulas, one
        ``<root>_<formula><ext>`` file each.
    """

    mol_dict = group_by(
        molecule_list, key=lambda a: a.get_chemical_formula(mode="hill")
    )

    if len(mol_dict) == 1:
        key, molecule_list = next(iter(mol_dict.items()))
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
        for key, molecule_list in mol_dict.items():
            if outfile is not None:
                root, ext = os.path.splitext(outfile)
                outfile_key = f"{root}_{key}{ext}"
            else:
                outfile_key = None
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
                repr = len(mol_dict[key])
                if printlvl > 0:
                    print(
                        f"--> wrote {repr} REPRESENTATIVE structure{'s' if repr != 1 else ''} to: {outfile_key}"
                    )
