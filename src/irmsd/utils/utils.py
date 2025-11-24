from __future__ import annotations

import sys
from typing import TYPE_CHECKING, List, Optional

import numpy as np

if TYPE_CHECKING:  # avoid hard dependency at import time
    from ase import Atoms  # type: ignore

__all__ = [
    "require_ase",
    "print_array",
]


# -----------------------------------------------------------------------------
# Dependency guard
# -----------------------------------------------------------------------------

def require_ase() -> ModuleType:
    """
    Import and return the ASE module, or raise a clear error if it is missing.
    """
    try:
        import ase  # type: ignore[import]
    except ImportError as exc:
        raise RuntimeError(
            "This function requires ASE, but it is not installed. "
            "Install it with `pip install ase`."
        ) from exc
    return ase


def require_rdkit() -> None:
    """Ensure rdkit is importable; raise a helpful ImportError otherwise.

    Use this at the *start* of any function that depends on ASE. We keep this
    separate so that importing `irmsd` does not immediately require ASE.
    """
    try:
        import rdkit  # noqa: F401
    except Exception as e:  # pragma: no cover
        raise ImportError(
            "rdkit is required for this function. Install optional extra: pip install 'irmsd[rdkit]'"
        ) from e


# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------


#def check_frames(obj: "Atoms | List[Atoms]", src: str) -> "Atoms":
#    """Check how many frames are in a provided atoms list."""
#    # Late import to avoid hard dependency unless called
#    from ase import Atoms  # type: ignore
#
#    if isinstance(obj, list) and len(obj) > 1:
#        if len(obj) == 0:
#            raise ValueError(f"No frames found in '{src}'.")
#        print(f"ℹ️  '{src}' has multiple frames; {len(obj)}.")
#        return obj
#    return obj
#
#
#def read_structures(paths: List[str]) -> List["Atoms"]:
#    """Read an arbitrary number of structures using ASE.
#
#    Parameters
#    ----------
#    paths : list of str
#        File paths to read.
#    Returns
#    -------
#    list[ase.Atoms]
#        One Atoms per input path (first frame if multi-frame file).
#    """
#    require_ase()
#    from ase.io import read as ase_read  # type: ignore
#
#    atoms_list: List["Atoms"] = []
#    for p in paths:
#        try:
#            ext = p.lower().endswith(".xyz")
#            if ext :
#               frames = ase_read(p, index=":", format="extxyz")  # force XYZ reader
#            else:
#               frames = ase_read(p, index=":")  # let ASE guess
#            atoms = check_frames(frames, p)
#            if isinstance(atoms, list):
#                atoms_list.extend(atoms)  # add elements individually
#            else:
#                atoms_list.append(atoms)  # add a single Atoms object
#        except Exception as e:  # pragma: no cover
#            print(f"❌ Failed to read '{p}': {e}", file=sys.stderr)
#            raise
#    return atoms_list


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------


def print_array(title: str, arr: np.ndarray) -> None:
    """Pretty-print a numpy array with a header and spacing."""
    print(title)
    with np.printoptions(precision=6, suppress=True):
        print(arr)
    print()


def print_structur(atom) -> None:
    """Print basic information about an ASE Atoms object."""
    require_ase()
    from ase import Atoms  # type: ignore

    if not isinstance(atom, Atoms):
        raise TypeError("print_structur expects an ase.Atoms object")

    print(f"{len(atom)}\n")
    for sym, pos in zip(atom.get_chemical_symbols(), atom.get_positions()):
        print(f"{sym:2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}")
