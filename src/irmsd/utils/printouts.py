from collections.abc import Sequence
import numpy as np
from ..core import Molecule

HARTREE_TO_KCAL_MOL = 627.509474


def print_array(title: str, arr: np.ndarray) -> None:
    """Pretty-print a numpy array with a header and spacing."""
    print(title)
    with np.printoptions(precision=6, suppress=True):
        print(arr)
    print()


def print_structure(mol) -> None:
    """
    Print basic information about a Molecule object in a simple XYZ-like format.

    Parameters
    ----------
    mol : Molecule
        Molecule instance to print.

    Raises
    ------
    TypeError
        If the input is not a Molecule.
    """
    if not isinstance(mol, Molecule):
        raise TypeError("print_structure expects a Molecule object")

    nat = len(mol)
    symbols = mol.get_chemical_symbols()
    positions = mol.get_positions()

    print(f"{nat}\n")
    for sym, (x, y, z) in zip(symbols, positions):
        print(f"{sym:2} {x:12.6f} {y:12.6f} {z:12.6f}")


def print_structure_summary(
    key: str,
    energies_hartree: Sequence[float] | None = None,
    delta_irmsd: Sequence[float] | None = None,
    max_rows: int | None = None,
) -> None:
    """
    Pretty-print a table summarising structures and associated quantities.

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
