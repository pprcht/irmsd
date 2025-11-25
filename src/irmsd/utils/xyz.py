from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence, TextIO
import io
import shlex

import numpy as np

from ..core.molecule import Molecule


# ---------------------------------------------------------------------------
# Helpers for handling file-like vs. path
# ---------------------------------------------------------------------------


def _open_maybe(path_or_file: str | Path | TextIO, mode: str) -> tuple[TextIO, bool]:
    """
    Accept either a path-like object or an open file object.
    Returns (fileobj, should_close).
    """
    if isinstance(path_or_file, (str, Path)):
        f = open(path_or_file, mode, encoding="utf8")
        return f, True
    else:
        # Assume file-like
        return path_or_file, False


# ---------------------------------------------------------------------------
# Helpers for parsing/formatting comment-line key=value pairs
# ---------------------------------------------------------------------------


def _parse_value(s: str) -> Any:
    """Try to interpret a string as int, float, bool, or fall back to str."""
    sl = s.lower()
    if sl in {"true", "t", "yes"}:
        return True
    if sl in {"false", "f", "no"}:
        return False

    # int?
    try:
        return int(s)
    except ValueError:
        pass

    # float?
    try:
        return float(s)
    except ValueError:
        pass

    return s


def _parse_cell_value(val: str) -> np.ndarray | None:
    """
    Parse a cell string into a (3, 3) array.

    Expected formats (very tolerant):
    - 'a11 a12 a13 a21 a22 a23 a31 a32 a33'
    - '[a11, a12, ..., a33]'
    - 3 values → interpreted as diagonal cell lengths
    """
    # Remove brackets/commas to be permissive
    cleaned = val.replace("[", " ").replace("]", " ").replace(",", " ")
    parts = cleaned.split()
    if not parts:
        return None

    try:
        values = [float(x) for x in parts]
    except ValueError:
        return None

    if len(values) == 9:
        arr = np.array(values, dtype=float).reshape(3, 3)
        return arr
    if len(values) == 3:
        # Diagonal cell
        arr = np.diag(values).astype(float)
        return arr

    # Unrecognized cell format
    return None


def _parse_pbc_value(val: str) -> tuple[bool, bool, bool] | None:
    """
    Parse a PBC string into a 3-tuple of bools.

    Typical formats:
    - 'T T T'
    - '1 1 0'
    - '[T, F, T]'
    """
    cleaned = val.replace("[", " ").replace("]", " ").replace(",", " ")
    parts = cleaned.split()
    if len(parts) != 3:
        return None

    def to_bool(x: str) -> bool:
        xl = x.lower()
        if xl in {"t", "true", "1", "yes"}:
            return True
        if xl in {"f", "false", "0", "no"}:
            return False
        return False

    return tuple(to_bool(x) for x in parts)  # type: ignore[return-value]


def _parse_comment_line(
    line: str,
) -> tuple[
    dict[str, Any], float | None, np.ndarray | None, tuple[bool, bool, bool] | None
]:
    """
    Parse an extended-XYZ comment line into:

    - info dict of generic key→value entries
    - energy (if 'energy=' present)
    - cell (if 'cell=' present)
    - pbc  (if 'pbc='  present)

    Remaining key=value pairs go into the info dict.
    """
    info: dict[str, Any] = {}
    energy: float | None = None
    cell: np.ndarray | None = None
    pbc: tuple[bool, bool, bool] | None = None

    # Use shlex to respect quotes in values: key="value with spaces"
    tokens = shlex.split(line, comments=False, posix=True)

    i = 0
    while i < len(tokens):
        token = tokens[i]

        # Must contain '=' or we skip
        if "=" not in token:
            i += 1
            continue

        key, val = token.split("=", 1)
        key = key.strip()
        val = val.strip()

        # CASE: "key=" followed by separate value token
        if val == "" and i + 1 < len(tokens):
            # Accept next token as the value
            nxt = tokens[i + 1].strip()
            val = nxt
            i += 1  # Skip over value token (we consumed it)

        kl = key.lower()

        if kl == "energy":
            try:
                energy = float(val)
            except ValueError:
                info[key] = _parse_value(val)

        elif kl == "cell":
            parsed = _parse_cell_value(val)
            if parsed is not None:
                cell = parsed
            else:
                info[key] = val

        elif kl == "pbc":
            parsed = _parse_pbc_value(val)
            if parsed is not None:
                pbc = parsed
            else:
                info[key] = val

        else:
            info[key] = _parse_value(val)

        i += 1

    return info, energy, cell, pbc



def _format_cell_value(cell: np.ndarray) -> str:
    """Format (3,3) cell into a compact string suitable for key=value."""
    flat = cell.reshape(-1)
    # Keep it simple: space-separated, quoted
    vals = " ".join(f"{x:.10f}" for x in flat)
    return f'"{vals}"'


def _format_pbc_value(pbc: tuple[bool, bool, bool]) -> str:
    """Format (3,) bool PBC into string."""
    chars = ["T" if b else "F" for b in pbc]
    return '"' + " ".join(chars) + '"'


# ---------------------------------------------------------------------------
# Reader: extended XYZ → Molecule / list[Molecule]
# ---------------------------------------------------------------------------


def read_extxyz(path_or_file: str | Path | TextIO) -> Molecule | list[Molecule]:
    """
    Read an extended-XYZ file and return either:

    - a single Molecule (if only one structure is present)
    - a list[Molecule] if there are multiple structures

    The extended-XYZ comment line may contain key=value pairs.
    Special handling:
    - 'energy=' → stored in Molecule.energy
    - 'cell='   → stored in Molecule.cell
    - 'pbc='    → stored in Molecule.pbc
    All other key=value pairs go into Molecule.info.
    """
    f, should_close = _open_maybe(path_or_file, "r")
    molecules: list[Molecule] = []

    try:
        while True:
            # Read number of atoms
            line = f.readline()
            if not line:
                break  # EOF
            line = line.strip()
            if not line:
                # Skip empty lines between frames
                continue

            try:
                natoms = int(line)
            except ValueError as e:
                raise ValueError(
                    f"Failed to parse number of atoms from line: {line!r}"
                ) from e

            # Read comment line (can be empty)
            comment_line = f.readline()
            if comment_line is None:
                raise ValueError(
                    "Unexpected EOF while reading extended-XYZ comment line"
                )

            info, energy, cell, pbc = _parse_comment_line(comment_line.strip())

            # Read natoms atomic lines
            symbols: list[str] = []
            positions = np.zeros((natoms, 3), dtype=float)

            for i in range(natoms):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError("Unexpected EOF while reading atom coordinates")

                parts = atom_line.split()
                if len(parts) < 4:
                    raise ValueError(
                        f"Atom line {i+1} has fewer than 4 fields: {atom_line!r}"
                    )

                sym = parts[0]
                try:
                    x, y, z = map(float, parts[1:4])
                except ValueError as e:
                    raise ValueError(
                        f"Failed to parse coordinates on line: {atom_line!r}"
                    ) from e

                symbols.append(sym)
                positions[i, 0] = x
                positions[i, 1] = y
                positions[i, 2] = z

                # We ignore any per-atom extra fields for now.

            mol = Molecule(
                symbols=symbols,
                positions=positions,
                energy=energy,
                info=info,
                cell=cell,
                pbc=pbc,
            )
            molecules.append(mol)

    finally:
        if should_close:
            f.close()

    if not molecules:
        raise ValueError("No structures found in extended-XYZ file")

    if len(molecules) == 1:
        return molecules[0]
    return molecules


# ---------------------------------------------------------------------------
# Writer: Molecule / sequence[Molecule] → extended XYZ
# ---------------------------------------------------------------------------


def _iter_molecules(
    obj: Molecule | Sequence[Molecule],
) -> Iterable[Molecule]:
    """Normalize to an iterable of Molecule objects."""
    if isinstance(obj, Molecule):
        yield obj
    else:
        for m in obj:
            if not isinstance(m, Molecule):
                raise TypeError(
                    "write_extxyz expects Molecule or a sequence of Molecule"
                )
            yield m


def write_extxyz(
    path_or_file: str | Path | TextIO,
    molecules: Molecule | Sequence[Molecule],
    mode: str = "w",
) -> None:
    """
    Write one or many Molecule objects to an extended-XYZ file.

    Special behavior:
    - If Molecule.energy is not None, writes 'energy=<value>' in the comment line.
    - If Molecule.cell is not None, writes 'cell="<a11 ... a33>"'.
    - If Molecule.pbc is not None, writes 'pbc="T T T"' etc.
    - All entries in Molecule.info are written as additional key=value pairs.

    Parameters
    ----------
    path_or_file : str | Path | TextIO
        Output path or already opened file object.
    molecules : Molecule | Sequence[Molecule]
        One Molecule or a sequence of Molecule objects.
    mode : str
        File open mode, default "w". Use "a" to append.
    """
    f, should_close = _open_maybe(path_or_file, mode)

    try:
        for mol in _iter_molecules(molecules):
            natoms = mol.natoms

            # 1) number of atoms
            f.write(f"{natoms:d}\n")

            # 2) construct comment line with key=value pairs
            parts: list[str] = []

            # energy
            if mol.energy is not None:
                parts.append(f"energy={mol.energy:.12g}")

            # cell
            if mol.cell is not None:
                parts.append(f"cell={_format_cell_value(mol.cell)}")

            # pbc
            if mol.pbc is not None:
                parts.append(f"pbc={_format_pbc_value(mol.pbc)}")

            # info dict (do not overwrite energy/cell/pbc even if present)
            for key, value in mol.info.items():
                kl = key.lower()
                if kl in {"energy", "cell", "pbc"}:
                    continue

                if isinstance(value, bool):
                    sval = "T" if value else "F"
                elif isinstance(value, (int, float)):
                    sval = repr(value)
                else:
                    sval = str(value)
                    # Quote if there is whitespace
                    if any(c.isspace() for c in sval):
                        sval = f'"{sval}"'

                parts.append(f"{key}={sval}")

            comment_line = " ".join(parts) if parts else "generated_by=irmsd"
            f.write(comment_line + "\n")

            # 3) atom lines
            positions = mol.get_positions(copy=False)
            symbols = mol.get_chemical_symbols()
            for sym, (x, y, z) in zip(symbols, positions):
                f.write(f"{sym:2s} {x: .15f} {y: .15f} {z: .15f}\n")

    finally:
        if should_close:
            f.close()
