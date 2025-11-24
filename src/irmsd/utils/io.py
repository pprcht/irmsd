from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

import sys
import numpy as np  # only if you need it elsewhere; not actually used below

from ..core import Molecule
from .xyz import read_extxyz
from ..interfaces.ase_io import ase_to_molecule
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


def read_structures(paths: Sequence[str]) -> List[Molecule]:
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

