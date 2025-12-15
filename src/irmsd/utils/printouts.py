from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from ..core import Molecule

HARTREE_TO_KCAL_MOL = 627.509474

BANNER = r"""
    ██╗██████╗ ███╗   ███╗███████╗██████╗ 
    ╠═╣██╔══██╗████╗ ████║██╔════╝██╔══██╗
    ██║██████╔╝██╔████╔██║███████╗██║  ██║
    ██║██╔══██╗██║╚██╔╝██║╚════██║██║  ██║
    ██║██║  ██║██║ ╚═╝ ██║███████║██████╔╝
    ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝╚═════╝
       A tool for structure comparison 
            and ensemble pruning
   ────────────────────────────────────────
     © 2025 Philipp Pracht, Tobias Kaczun
   https://doi.org/10.1021/acs.jcim.4c02143 
       https://github.com/pprcht/irmsd
"""


def print_pretty_array(title: str, arr: np.ndarray, fmt="{:8.4f}", sep="    ") -> None:
    """Pretty-print a 1D or 2D numpy array with a header.

    Parameters
    ----------
    title : str
        Header title to print before the data.
    arr : numpy.ndarray
        1D or 2D array to print.
    fmt : str, optional
        Format string for each element (default: "{:8.4f}").
    sep : str, optional
        Separator between elements (default: four spaces).

    Raises
    ------
    ValueError
        If the array is not 1D or 2D.
    """
    print(title)
    if arr.ndim == 1:
        print(sep + sep.join(fmt.format(x) for x in arr))

    elif arr.ndim == 2:
        for row in arr:
            print(sep + sep.join(fmt.format(x) for x in row))

    else:
        raise ValueError("Only 1D or 2D arrays are supported.")


def _print_atomwise_table(
    mol,
    properties: Mapping[str, np.ndarray],
) -> None:
    """Pretty-print multiple atom-wise properties for a single Molecule."""
    if not properties:
        return

    def _infer_fmt(arr):
        if np.issubdtype(arr.dtype, np.integer):
            return "{:14d}"
        elif np.issubdtype(arr.dtype, np.floating):
            return "{:14.6f}"
        elif np.issubdtype(arr.dtype, np.bool_):
            return "{:>14}"  # prints True/False
            # or "{:14d}" to print 1/0
        else:
            return "{:>14}"  # fallback for strings or objects

    nat = len(mol)
    # basic checks
    for name, arr in properties.items():
        arr = np.asarray(arr)
        if arr.ndim != 1:
            raise ValueError(f"Property '{name}' is not 1D (shape={arr.shape}).")
        if len(arr) != nat:
            raise ValueError(
                f"Property '{name}' length {len(arr)} != number of atoms {nat}"
            )

    prop_names = list(properties.keys())

    # header
    header = f"{'Atom':>4} {'Symbol':>6}"
    for name in prop_names:
        header += f" {name:>14}"
    print(header)

    # separator
    sep = "---- ------"
    for _ in prop_names:
        sep += " " + "-" * 14
    print(sep)

    symbols = mol.get_chemical_symbols()

    for i in range(nat):
        row = f"{i+1:4d} {symbols[i]:>6}"
        for name in prop_names:
            arr = np.asarray(properties[name])
            fmt_this = _infer_fmt(arr)
            row += " " + fmt_this.format(arr[i])
        print(row)

    print()


def print_molecule_summary(
    molecule_list: Sequence[Any],
    **results_by_name: Sequence[Any],
) -> None:
    """Print a summary for each molecule, plus a combined atom-wise table for
    any results that are 1D per-atom arrays.

    Parameters
    ----------
    molecule_list : sequence of Molecule
        List of molecules.
    **results_by_name :
        Each keyword argument is a sequence aligned with molecule_list,
        e.g. energies=[...], charges=[...], spin=[...].
    """
    n_mol = len(molecule_list)

    # sanity: all result lists must match length of molecule_list
    for name, seq in results_by_name.items():
        if len(seq) != n_mol:
            raise ValueError(
                f"Result '{name}' has length {len(seq)}, expected {n_mol}."
            )

    for idx, mol in enumerate(molecule_list):
        print("\n" + "=" * 60)
        print(f"###  MOLECULE {idx+1:>3}  ###")
        print("=" * 60)
        print()

        # Split per-molecule vs atom-wise for THIS molecule
        per_mol_values: dict[str, Any] = {}
        atomwise_values: dict[str, np.ndarray] = {}

        for name, seq in results_by_name.items():
            value = seq[idx]

            # Detect atom-wise: 1D array, length == number of atoms
            if (
                isinstance(value, np.ndarray)
                and value.ndim == 1
                and len(value) == len(mol)
            ):
                atomwise_values[name] = value
            else:
                per_mol_values[name] = value

        # 1) Print per-molecule values
        for name, value in per_mol_values.items():
            # Case A: the value itself is a dict → iterate through it
            if isinstance(value, dict):
                for subname, subval in value.items():
                    print_pretty_array(f"{subname}:", subval)
                print()

            # Case B: normal scalar or non-dict value
            else:
                print(f"{name}: {value}")
                print()

        # 2) Combined atom-wise table (if any)
        if atomwise_values:
            _print_atomwise_table(mol, atomwise_values)

        print()  # spacing between molecules


def print_conformer_structures(*mols, labels=None) -> None:
    """Print multiple Molecule objects representing different conformers of the
    same molecule in a combined XYZ-like format side-by-side.

    Parameters
    ----------
    *mol : Molecule
        One or more Molecule instances to print.

    Raises
    ------
    TypeError
        If any input is not a Molecule.
    ValueError
        If the Molecule objects do not have the same number of atoms
        or atom ordering.
    """
    assert len(mols) > 0, "At least one Molecule must be provided"
    for i, m in enumerate(mols):
        if not isinstance(m, Molecule):
            raise TypeError(f"Argument {i} is not a Molecule object")

    nat = len(mols[0])
    for m in mols:
        if len(m) != nat:
            raise ValueError("All Molecule objects must have the same number of atoms")

    sep = " │"
    if labels is not None:
        if len(labels) != len(mols):
            raise ValueError("Number of labels must match number of Molecule objects")
        label_line = sep.join(f"{label:^41}" for label in labels)
        print(label_line)
    for i in range(nat):
        fields = []
        for m in mols:
            symbols = m.get_chemical_symbols()
            positions = m.get_positions()
            x, y, z = positions[i]
            fields.append(f"{symbols[i]:>2} {x:>12.6f} {y:>12.6f} {z:>12.6f}")
        print(sep.join(fields))


def print_structure_summary(
    key: str,
    energies_hartree: Sequence[float] | None = None,
    delta_irmsd: Sequence[float] | None = None,
    max_rows: int | None = None,
) -> None:
    """Pretty-print a table summarising structures and associated quantities.

    Parameters
    ----------
    key : str
        A label/title for this block (e.g. method name, run ID, etc.).
    energies_hartree : 1D sequence of float, optional
        Energies in Hartree. If given, an additional 'ΔE / kcal mol⁻¹'
        column is printed relative to the first structure.
    delta_irmsd : 1D sequence of float, optional
        Delta iRMSD values.
    max_rows : int, optional
        Maximum number of data rows to print. If the total number of
        structures is larger, the table is truncated, an extra row
        of "..." is printed, and a message indicates how many entries
        were skipped. If None, all rows are printed.

    Notes
    -----
    - If *all* arrays are None, nothing is printed.
    - All provided arrays must have the same length.
    - First column is always 'structure {i}', i starting at 1.
    """

    if max_rows is not None and max_rows < 1:
        raise ValueError("max_rows must be >= 1 or None.")

    # --- collect numeric columns ---
    columns: list[tuple[str, list[str]]] = []  # (header, cells-as-strings)
    n: int | None = None

    def add_column(
        header: str,
        values: Sequence[float] | None,
        fmt: str,
    ) -> None:
        """Internal helper to add a numeric column."""
        nonlocal n
        if values is None:
            return

        vals = [float(v) for v in values]

        if n is None:
            n = len(vals)
        elif len(vals) != n:
            raise ValueError(
                f"All arrays must have the same length; "
                f"expected {n}, got {len(vals)} for column '{header}'."
            )

        cells = [fmt.format(v) for v in vals]
        columns.append((header, cells))

    # Add the explicit columns requested
    add_column("E / Eh", energies_hartree, "{: .10f}")
    # If we have energies, also add ΔE in kcal/mol relative to first entry
    if energies_hartree is not None:
        e0 = float(energies_hartree[0])
        delta_e_kcal = [(float(e) - e0) * HARTREE_TO_KCAL_MOL for e in energies_hartree]
        add_column("ΔE / kcal mol⁻¹", delta_e_kcal, "{: .3f}")

    add_column("ΔRMSD / Å", delta_irmsd, "{: .4f}")
    # If no arrays were provided at all: do not print anything
    if n is None or n == 0:
        return

    # --- structure labels column ---
    struct_labels = [f" {i+1}" for i in range(n)]
    all_columns = [("Structure", struct_labels)] + columns

    # --- compute column widths ---
    widths: list[int] = []
    for header, cells in all_columns:
        max_cell_len = max(len(c) for c in cells) if cells else 0
        widths.append(max(len(header), max_cell_len))

    # --- determine how many rows to print ---
    if max_rows is None or max_rows >= n:
        rows_to_print = n
        truncated = False
    else:
        rows_to_print = max_rows
        truncated = True

    # --- print the table ---
    print(f"\n=== {key} ===")

    header_line = "  ".join(
        header.ljust(w) for (header, _), w in zip(all_columns, widths)
    )
    sep_line = "  ".join("-" * w for w in widths)
    print(header_line)
    print(sep_line)

    # data rows
    for i in range(rows_to_print):
        row_cells = [col[i] for _, col in all_columns]
        line = "  ".join(cell.ljust(w) for cell, w in zip(row_cells, widths))
        print(line)

    # ellipsis row + summary, if truncated
    if truncated:
        ellipsis_cells = [" (...)" for _ in all_columns]
        ellipsis_line = "  ".join(
            cell.ljust(w) for cell, w in zip(ellipsis_cells, widths)
        )
        print(ellipsis_line)
        remaining = n - rows_to_print
        print(
            f"({remaining} additional entries not shown, use `--maxprint` to increase)"
        )
