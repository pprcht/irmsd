from __future__ import annotations

from pathlib import Path
from typing import Union, List, Sequence, overload

import sys
import numpy as np

from ..core import Molecule
from .xyz import read_extxyz, write_extxyz
from ..interfaces.ase_io import ase_to_molecule, molecule_to_ase
from .utils import require_ase


def check_frames(obj, src: str):
    """
    Check how many frames are in a provided object.

    Parameters
    ----------
    obj : single object or sequence
        Typically a single structure (Molecule or ASE Atoms) or a sequence
        of them as returned by a reader.
    src : str
        Source label, usually the file path. Used only for informational
        messages.

    Returns
    -------
    obj
        The input object unchanged. If a non-empty sequence with more than
        one frame is provided, an informational message is printed.
    """
    if isinstance(obj, (list, tuple)):
        if len(obj) == 0:
            raise ValueError(f"No frames found in '{src}'.")
        if len(obj) > 1:
            print(f"ℹ️  '{src}' has multiple frames; {len(obj)}.")
        return obj
    return obj


def read_structures(paths: Union[str,Sequence[str]]) -> List[Molecule]:
    """
    Read an arbitrary number of structures and return them as Molecule objects.

    For each path, this routine behaves as follows:

    - If the file extension is ``.xyz`` or ``.extxyz``, it uses the internal
      ``read_extxyz`` helper to obtain one or more Molecule objects.

    - For all other file types, it attempts to import ASE via ``require_ase()``,
      uses ``ase.io.read`` to read one or more ASE Atoms objects, and converts
      them into Molecule objects using ``ase_to_molecule``.

    Multi-frame files:
      If a file contains multiple frames, all frames are read and appended to
      the output list. A short informational message is printed indicating the
      number of frames that were found.

    Parameters
    ----------
    paths : Sequence[str]
        File paths to read.

    Returns
    -------
    structures : list[Molecule]
        One Molecule per frame found across all input paths.
    """
    molecules: List[Molecule] = []

    if isinstance(paths, str):
        paths = [paths]

    for p in paths:
        path = str(p)
        ext = Path(path).suffix.lower()

        try:
            # --- Our own extended XYZ reader branch ---
            if ext in {".xyz", ".extxyz"}:
                obj = read_extxyz(path)  # Molecule or list[Molecule]
                obj = check_frames(obj, path)

                if isinstance(obj, list):
                    # obj is list[Molecule]
                    molecules.extend(obj)
                else:
                    # single Molecule
                    molecules.append(obj)

            # --- ASE fallback for other formats ---
            else:
                ase = require_ase()
                from ase.io import read as ase_read  # type: ignore

                frames = ase_read(path, index=":")  # all frames
                frames = check_frames(frames, path)

                if isinstance(frames, list):
                    # list[ase.Atoms] -> list[Molecule]
                    mols = ase_to_molecule(frames)
                    molecules.extend(mols)
                else:
                    # single ase.Atoms -> Molecule
                    mol = ase_to_molecule(frames)
                    molecules.append(mol)

        except Exception as e:  # pragma: no cover
            print(f"❌ Failed to read '{path}': {e}", file=sys.stderr)
            raise

    return molecules


@overload
def write_structures(
    filename: str | Path,
    structures: Molecule,
    mode: str = "w",
) -> None: ...
@overload
def write_structures(
    filename: str | Path,
    structures: Sequence[Molecule],
    mode: str = "w",
) -> None: ...


def write_structures(
    filename: str | Path,
    structures: Molecule | Sequence[Molecule],
    mode: str = "w",
) -> None:
    """
    High-level structure writer for the irmsd Molecule type.

    This routine mirrors the behaviour of `read_structures` on the output side:
    it chooses the appropriate backend based on the file extension and accepts
    either a single Molecule or a sequence of Molecule objects.

    Dispatch rules
    --------------
    - If the filename has *no* extension, '.xyz' is appended and the
      internal extended-XYZ writer is used.

    - If the filename ends with '.xyz', '.extxyz' or '.trj' (case-insensitive),
      the structures are written using the internal extended-XYZ writer
      (`write_extxyz`). The `mode` argument is passed through and controls
      whether the file is overwritten ('w', default) or appended to ('a').

    - For all other filename extensions, ASE is used as a backend. The
      structures are first converted to ASE `Atoms` objects using
      `molecule_to_ase`, and then written via `ase.io.write`. In this case
      the `mode` argument is currently ignored and ASE's default behaviour
      for the chosen format is used.

    Parameters
    ----------
    filename : str or pathlib.Path
        Output filename. Its extension determines the backend.
    structures : Molecule or Sequence[Molecule]
        A single Molecule or a sequence of Molecules to be written.
        For extended XYZ, multiple Molecules are written as consecutive
        frames in one file.
    mode : {"w", "a"}, optional
        File open mode for extended XYZ output. Ignored for non-XYZ formats
        handled via ASE.

    Raises
    ------
    RuntimeError
        If ASE is required (non-XYZ formats) but not installed.
    TypeError
        If `structures` is not a Molecule or a sequence of Molecules.
    """
    path = Path(filename)

    # If no extension is given, default to '.xyz' and use our internal writer
    if path.suffix == "":
        path = path.with_suffix(".xyz")

    ext = path.suffix.lower()

    # Extended XYZ branch (internal writer)
    if ext in {".xyz", ".extxyz", ".trj"}:
        write_extxyz(path, structures, mode=mode)
        return

    if ext in {".pkl"}:
        dump_results_to_pickle(structures,path)
        return

    # ASE branch for all other extensions
    require_ase()
    from ase.io import write as ase_write  # type: ignore[import]

    # molecule_to_ase must itself handle single vs sequence
    ase_obj = molecule_to_ase(structures)
    ase_write(str(path), ase_obj)


# -----------------------------------------------------------
# More general data dump: pickle files!
# -----------------------------------------------------------
import os
import pickle
import gzip
import bz2
import lzma

from collections.abc import Sequence, Mapping
from typing import Any


def dump_results_to_pickle(
    molecules: Sequence[Any],
    outfile: str,
    results: Mapping[str, Any] | None = None,
    compress: str | None = None,
) -> str:
    """
    Dump a sequence of Molecule objects (and optionally a results mapping)
    into a pickle file, optionally compressed.

    Parameters
    ----------
    molecules : Sequence[Any]
        Sequence of Molecule objects.
    outfile : str
        Output filename. If extension is missing or wrong, '.pkl' is used.
    results : Mapping[str, Any] or None
        Optional results dictionary. If provided, its entries become top-level
        keys in the payload and 'molecules' is added. If None, the payload
        is simply `molecules`.
    compress : {None, "gz", "bz2", "xz"}, optional
        Compression type.

    Returns
    -------
    str
        Final filename used.
    """

    # --- Validate compression option ---
    valid_compress = {None, "gz", "bz2", "xz"}
    if compress not in valid_compress:
        raise ValueError(
            f"Invalid compress={compress!r}; " f"expected {valid_compress}."
        )

    # --- Normalize filename ---
    base, ext = os.path.splitext(outfile)

    # Handle cases like "foo.pkl.gz"
    if ext in (".gz", ".bz2", ".xz"):
        base, pkl_ext = os.path.splitext(base)
        ext = pkl_ext

    if ext == "" or ext.lower() != ".pkl":
        if ext not in ("", ".pkl"):
            print(f"Warning: changing extension '{ext}' to '.pkl'")
        ext = ".pkl"

    if compress is None:
        outfile = base + ext
    else:
        outfile = f"{base}{ext}.{compress}"

    # --- Build payload ---
    if results is None:
        # simplest case: payload is just the sequence of molecules
        payload = molecules
    else:
        if "molecules" in results:
            raise ValueError(
                "The results mapping contains a key 'molecules', "
                "which would conflict with the payload key."
            )

        payload = {
            **results,
            "molecules": molecules,
        }

    # --- Select opener ---
    if compress is None:
        open_fn = open
    elif compress == "gz":
        open_fn = gzip.open
    elif compress == "bz2":
        open_fn = bz2.open
    elif compress == "xz":
        open_fn = lzma.open

    # --- Write pickle ---
    with open_fn(outfile, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)

    return outfile
