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
    """Print ``title``, then a 1D or 2D array one row per line (ValueError otherwise)."""
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
        else:
            return "{:>14}"  # fallback for strings or objects

    nat = len(mol)
    for name, arr in properties.items():
        arr = np.asarray(arr)
        if arr.ndim != 1:
            raise ValueError(f"Property '{name}' is not 1D (shape={arr.shape}).")
        if len(arr) != nat:
            raise ValueError(
                f"Property '{name}' length {len(arr)} != number of atoms {nat}"
            )

    prop_names = list(properties.keys())

    header = f"{'Atom':>4} {'Symbol':>6}"
    for name in prop_names:
        header += f" {name:>14}"
    print(header)

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
    """Print per-molecule results; 1D per-atom arrays go into one atom-wise table.

    Each ``results_by_name`` entry is a sequence aligned with
    ``molecule_list``, e.g. ``energies=[...]``.
    """
    n_mol = len(molecule_list)

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

        per_mol_values: dict[str, Any] = {}
        atomwise_values: dict[str, np.ndarray] = {}

        for name, seq in results_by_name.items():
            value = seq[idx]

            if (
                isinstance(value, np.ndarray)
                and value.ndim == 1
                and len(value) == len(mol)
            ):
                atomwise_values[name] = value
            else:
                per_mol_values[name] = value

        for name, value in per_mol_values.items():
            if isinstance(value, dict):
                for subname, subval in value.items():
                    if isinstance(subval, np.ndarray):
                        print_pretty_array(f"{subname}:", subval)
                    elif isinstance(subval, str):
                        print(f"{subname}: {subval}")
                print()

            else:
                print(f"{name}: {value}")
                print()

        if atomwise_values:
            _print_atomwise_table(mol, atomwise_values)

        print()  # spacing between molecules


def print_conformer_structures(*mols, labels=None) -> None:
    """Print conformers of one molecule side by side in XYZ-like columns.

    ``labels``, if given, head the columns, one per molecule. Raises TypeError
    for non-Molecule input, ValueError on atom-count or label-count mismatch.
    """
    if not mols:
        raise ValueError("At least one Molecule must be provided")
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
    """Print a per-structure table of energies and delta-iRMSD values.

    Parameters
    ----------
    key : str
        Block title.
    energies_hartree : sequence of float, optional
        Energies in Hartree; adds a kcal/mol ΔE column relative to the first.
    delta_irmsd : sequence of float, optional
        Delta iRMSD values in Å.
    max_rows : int, optional
        Truncate after this many rows and report the skipped count.

    Prints nothing if both arrays are None; given arrays must match in length.
    """
    if max_rows is not None and max_rows < 1:
        raise ValueError("max_rows must be >= 1 or None.")

    columns: list[tuple[str, list[str]]] = []  # (header, cells-as-strings)
    n: int | None = None

    def add_column(
        header: str,
        values: Sequence[float] | None,
        fmt: str,
    ) -> None:
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

    add_column("E / Eh", energies_hartree, "{: .10f}")
    if energies_hartree is not None:
        e0 = float(energies_hartree[0])
        delta_e_kcal = [(float(e) - e0) * HARTREE_TO_KCAL_MOL for e in energies_hartree]
        add_column("ΔE / kcal mol⁻¹", delta_e_kcal, "{: .3f}")

    add_column("ΔRMSD / Å", delta_irmsd, "{: .4f}")
    if n is None or n == 0:
        return

    struct_labels = [f" {i+1}" for i in range(n)]
    all_columns = [("Structure", struct_labels)] + columns

    widths: list[int] = []
    for header, cells in all_columns:
        max_cell_len = max(len(c) for c in cells) if cells else 0
        widths.append(max(len(header), max_cell_len))

    if max_rows is None or max_rows >= n:
        rows_to_print = n
        truncated = False
    else:
        rows_to_print = max_rows
        truncated = True

    print(f"\n=== {key} ===")

    header_line = "  ".join(
        header.ljust(w) for (header, _), w in zip(all_columns, widths)
    )
    sep_line = "  ".join("-" * w for w in widths)
    print(header_line)
    print(sep_line)

    for i in range(rows_to_print):
        row_cells = [col[i] for _, col in all_columns]
        line = "  ".join(cell.ljust(w) for cell, w in zip(row_cells, widths))
        print(line)

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
