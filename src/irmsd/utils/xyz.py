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
# Energy units
# ---------------------------------------------------------------------------
#
# Molecule.energy is stored in Hartree throughout the package. The extended-XYZ
# comment line may carry an optional ``energy_units=`` token; when absent the
# energy is assumed to already be in Hartree (the long-standing CREST/QC
# convention). We deliberately avoid depending on ASE here, so the conversion
# factor is hard-coded (CODATA 2018 value, identical to ``ase.units.Hartree``).

#: eV per Hartree (CODATA 2018; matches ``ase.units.Hartree``).
EV_PER_HARTREE = 27.211386245988

#: Values (case-insensitive) flagging an energy already given in Hartree.
_HARTREE_UNIT_ALIASES = frozenset(
    {"hartree", "hartrees", "ha", "au", "a.u.", "eh", "e_h", "atomic"}
)
#: Values (case-insensitive) flagging an energy given in electronvolt.
_EV_UNIT_ALIASES = frozenset({"ev", "electronvolt", "electronvolts"})


def _energy_to_hartree(value: float | None, units: str | None) -> float | None:
    """Convert a comment-line energy to Hartree given an optional unit token.

    - ``units is None`` (no marker) → assumed Hartree, returned unchanged.
    - ``units`` an eV alias → divided by :data:`EV_PER_HARTREE`.
    - ``units`` a Hartree alias (or anything unrecognized) → returned unchanged.
    """
    if value is None:
        return None
    if units is None:
        return value
    u = str(units).strip().lower()
    if u in _EV_UNIT_ALIASES:
        return value / EV_PER_HARTREE
    return value  # Hartree alias or unknown -> treat as already Hartree


# ---------------------------------------------------------------------------
# Per-atom column schema (extended-XYZ ``Properties=`` token)
# ---------------------------------------------------------------------------
#
# The atom block layout is described by a ``Properties=`` token in the comment
# line, e.g. ``Properties=species:S:1:pos:R:3``. We support the standard
# species/pos columns plus an optional integer ``canonical_id`` column carrying
# the per-atom canonical identifiers (``Molecule.ids``). When no ``Properties``
# token is present the default ``species:S:1:pos:R:3`` layout is assumed, which
# reproduces the historical behavior (extra columns ignored).

#: Name of the per-atom column holding canonical atom identifiers.
CANONICAL_ID_COLUMN = "canonical_id"

#: Default per-atom column schema when no ``Properties=`` token is given.
_DEFAULT_PROPERTIES = "species:S:1:pos:R:3"

#: Schema string emitted when a Molecule carries per-atom IDs.
_PROPERTIES_WITH_IDS = f"species:S:1:pos:R:3:{CANONICAL_ID_COLUMN}:I:1"


def _parse_properties_schema(spec: str | None) -> list[tuple[str, str, int]] | None:
    """Parse a ``Properties`` spec into a list of ``(name, type, count)`` triples.

    Returns None if the spec is empty or malformed (caller falls back to the
    default layout).
    """
    if not spec:
        return None
    toks = [t for t in str(spec).split(":") if t != ""]
    if not toks or len(toks) % 3 != 0:
        return None
    fields: list[tuple[str, str, int]] = []
    for i in range(0, len(toks), 3):
        name = toks[i]
        typ = toks[i + 1].upper()
        try:
            count = int(toks[i + 2])
        except ValueError:
            return None
        if count <= 0:
            return None
        fields.append((name, typ, count))
    return fields


def _column_layout(
    fields: list[tuple[str, str, int]] | None,
) -> tuple[int, int | None, int | None, int]:
    """Resolve a schema to column offsets.

    Returns ``(sym_idx, pos_idx, id_idx, ncols)`` where the indices are 0-based
    offsets into the whitespace-split atom line. ``pos_idx``/``id_idx`` are None
    when the corresponding column is absent. The species column is mandatory and
    defaults to 0.
    """
    if fields is None:
        fields = _parse_properties_schema(_DEFAULT_PROPERTIES)

    sym_idx: int | None = None
    pos_idx: int | None = None
    id_idx: int | None = None
    offset = 0
    for name, typ, count in fields:  # type: ignore[union-attr]
        lname = name.lower()
        if sym_idx is None and (lname == "species" or typ == "S"):
            sym_idx = offset
        if pos_idx is None and (lname == "pos" or (typ == "R" and count == 3)):
            pos_idx = offset
        if id_idx is None and lname == CANONICAL_ID_COLUMN:
            id_idx = offset
        offset += count

    if sym_idx is None:
        sym_idx = 0
    return sym_idx, pos_idx, id_idx, offset


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
    dict[str, Any],
    float | None,
    np.ndarray | None,
    tuple[bool, bool, bool] | None,
    str | None,
]:
    """
    Parse an extended-XYZ comment line into:

    - info dict of generic key→value entries
    - energy (if 'energy=' present), returned in Hartree
    - cell (if 'cell=' present)
    - pbc  (if 'pbc='  present)
    - properties (the raw ``Properties=`` schema string, if present)

    Remaining key=value pairs go into the info dict.

    The energy is interpreted according to an optional ``energy_units=`` token:
    absent (or a Hartree alias) means the value is already Hartree, ``eV`` is
    converted to Hartree. The ``energy_units`` token is consumed and not placed
    in the info dict. The ``Properties`` token is likewise consumed (it describes
    the atom-block column layout, not per-frame metadata).
    """
    info: dict[str, Any] = {}
    energy: float | None = None
    energy_units: str | None = None
    cell: np.ndarray | None = None
    pbc: tuple[bool, bool, bool] | None = None
    properties: str | None = None

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

        elif kl == "energy_units":
            # Consume the unit marker; it is applied below and intentionally
            # kept out of Molecule.info so it cannot trigger a double
            # conversion on a later round-trip.
            energy_units = val

        elif kl == "properties":
            # Consume the atom-block column schema; it drives the per-atom
            # parsing below and is not per-frame metadata.
            properties = val

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

    # Normalize the energy to the internal Hartree convention.
    energy = _energy_to_hartree(energy, energy_units)

    return info, energy, cell, pbc, properties



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
    - 'Properties=' → describes the atom-block column layout; a
      ``canonical_id:I:1`` column (if present) is read into Molecule.ids
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

            info, energy, cell, pbc, properties = _parse_comment_line(
                comment_line.strip()
            )

            # Resolve the per-atom column layout from the Properties schema
            # (falls back to the default species:S:1:pos:R:3 layout).
            fields = _parse_properties_schema(properties)
            sym_idx, pos_idx, id_idx, ncols = _column_layout(fields)
            if pos_idx is None:
                pos_idx = 1  # default: coordinates follow the species column

            # Read natoms atomic lines
            symbols: list[str] = []
            positions = np.zeros((natoms, 3), dtype=float)
            ids = (
                np.zeros(natoms, dtype=np.int32) if id_idx is not None else None
            )

            min_cols = max(sym_idx + 1, pos_idx + 3, (id_idx + 1) if id_idx is not None else 0)
            for i in range(natoms):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError("Unexpected EOF while reading atom coordinates")

                parts = atom_line.split()
                if len(parts) < min_cols:
                    raise ValueError(
                        f"Atom line {i+1} has fewer than {min_cols} fields: {atom_line!r}"
                    )

                sym = parts[sym_idx]
                try:
                    x, y, z = map(float, parts[pos_idx : pos_idx + 3])
                except ValueError as e:
                    raise ValueError(
                        f"Failed to parse coordinates on line: {atom_line!r}"
                    ) from e

                symbols.append(sym)
                positions[i, 0] = x
                positions[i, 1] = y
                positions[i, 2] = z

                if ids is not None:
                    try:
                        # Tolerate integers written as floats (e.g. "2.0").
                        ids[i] = int(round(float(parts[id_idx])))
                    except ValueError as e:
                        raise ValueError(
                            f"Failed to parse canonical_id on line: {atom_line!r}"
                        ) from e

                # Any further per-atom columns are ignored.

            mol = Molecule(
                symbols=symbols,
                positions=positions,
                energy=energy,
                info=info,
                cell=cell,
                pbc=pbc,
                ids=ids,
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
    - If Molecule.energy is not None, writes 'energy=<value> energy_units=Hartree'
      in the comment line (Molecule energies are stored in Hartree).
    - If Molecule.cell is not None, writes 'cell="<a11 ... a33>"'.
    - If Molecule.pbc is not None, writes 'pbc="T T T"' etc.
    - If Molecule.ids is not None, writes a per-atom ``canonical_id`` column and
      advertises it via 'Properties=species:S:1:pos:R:3:canonical_id:I:1'.
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

            # energy (stored in Hartree); stamp an explicit unit marker so the
            # file is self-documenting and read back unambiguously by both the
            # native reader and a manual ASE read.
            if mol.energy is not None:
                parts.append(f"energy={mol.energy:.12g}")
                parts.append("energy_units=Hartree")

            # cell
            if mol.cell is not None:
                parts.append(f"cell={_format_cell_value(mol.cell)}")

            # pbc
            if mol.pbc is not None:
                parts.append(f"pbc={_format_pbc_value(mol.pbc)}")

            # per-atom column schema: only advertise the canonical_id column when
            # the Molecule actually carries IDs (otherwise stay byte-compatible
            # with the historical species/pos-only output).
            write_ids = mol.ids is not None
            if write_ids:
                parts.append(f"Properties={_PROPERTIES_WITH_IDS}")

            # info dict (do not overwrite energy/cell/pbc/properties even if present)
            for key, value in mol.info.items():
                kl = key.lower()
                if kl in {"energy", "energy_units", "cell", "pbc", "properties"}:
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
            if write_ids:
                ids = mol.get_ids(copy=False)
                for sym, (x, y, z), aid in zip(symbols, positions, ids):
                    f.write(
                        f"{sym:2s} {x: .15f} {y: .15f} {z: .15f} {int(aid):d}\n"
                    )
            else:
                for sym, (x, y, z) in zip(symbols, positions):
                    f.write(f"{sym:2s} {x: .15f} {y: .15f} {z: .15f}\n")

    finally:
        if should_close:
            f.close()
