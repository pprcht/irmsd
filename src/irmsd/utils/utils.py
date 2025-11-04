from __future__ import annotations

"""ASE-dependent utility functions for irmsd.

These helpers deliberately *avoid* importing ASE at module import time so that
`import irmsd` remains robust when ASE is not installed. Functions that need ASE
call `require_ase()` first and then perform local (late) imports.
"""

import sys
from typing import List, Optional, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # avoid hard dependency at import time
    from ase import Atoms  # type: ignore

__all__ = [
    "require_ase",
    "read_structures",
    "print_array",
]


# -----------------------------------------------------------------------------
# Dependency guard
# -----------------------------------------------------------------------------

def require_ase() -> None:
    """Ensure ASE is importable; raise a helpful ImportError otherwise.

    Use this at the *start* of any function that depends on ASE. We keep this
    separate so that importing `irmsd` does not immediately require ASE.
    """
    try:
        import ase  # noqa: F401
    except Exception as e:  # pragma: no cover
        raise ImportError(
            "ASE is required for this function. Install optional extra: pip install 'irmsd[ase]'"
        ) from e


# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------

def check_frames(obj: "Atoms | List[Atoms]", src: str) -> "Atoms":
    """
    Check how many frames are in a provided atoms list
    """
    # Late import to avoid hard dependency unless called
    from ase import Atoms  # type: ignore

    if isinstance(obj, list) and len(obj) > 1:
        if len(obj) == 0:
            raise ValueError(f"No frames found in '{src}'.")
        print(f"ℹ️  '{src}' has multiple frames; {len(obj)}.")
        return obj
    return obj


def read_structures(paths: List[str]) -> List["Atoms"]:
    """Read an arbitrary number of structures using ASE.

    Parameters
    ----------
    paths : list of str
        File paths to read.
    Returns
    -------
    list[ase.Atoms]
        One Atoms per input path (first frame if multi-frame file).
    """
    require_ase()
    from ase.io import read as ase_read  # type: ignore

    atoms_list: List["Atoms"] = []
    for p in paths:
        try:
            frames = ase_read(p, index=":")
            atoms = check_frames(frames, p)
            if isinstance(atoms, list):
                atoms_list.extend(atoms)   # add elements individually
            else:
                atoms_list.append(atoms)   # add a single Atoms object
        except Exception as e:  # pragma: no cover
            print(f"❌ Failed to read '{p}': {e}", file=sys.stderr)
            raise
    return atoms_list


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

def print_array(title: str, arr: np.ndarray) -> None:
    """Pretty-print a numpy array with a header and spacing."""
    print(title)
    with np.printoptions(precision=6, suppress=True):
        print(arr)
    print()


