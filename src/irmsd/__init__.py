# src/irmsd/__init__.py
from __future__ import annotations

from .api.cn_exposed import get_cn_fortran

# ---- Core API (what you already had) ----------------------------------------
from .api.xyz_bridge import xyz_to_fortran, xyz_to_fortran_pair

# Try to expose ase_to_fortran if ASE is present; otherwise export a stub that errors nicely.
try:
    from .utils.ase_io import (
        ase_to_fortran,
        ase_to_fortran_pair,
        get_axis_ase,
        get_cn_ase,
    )
except Exception:

    def ase_to_fortran(*args, **kwargs):  # type: ignore
        raise ImportError(
            "ASE-based helper not available. Install optional extra: pip install 'irmsd[ase]'"
        )

    def ase_to_fortran_pair(*args, **kwargs):  # type: ignore
        raise ImportError(
            "ASE-based helper not available. Install optional extra: pip install 'irmsd[ase]'"
        )


# ---- New: re-export Python utilities ----------------------------------------
from .utils.utils import print_array, read_structures

# ---- Optional: command helpers ----------
try:
    from .utils.cmds import compute_axis_and_print, compute_cn_and_print  # noqa: F401
except Exception:
    # Safe to ignore so `import irmsd` never breaks due to optional pieces.
    pass

__all__ = [
    # core API
    "xyz_to_fortran",
    "xyz_to_fortran_pair",
    "ase_to_fortran",
    "ase_to_fortran_pair",
    "get_cn_fortran",
    "get_cn_ase",
    # utils
    "read_structures",
    "print_array",
    # optional cmds
    "compute_cn_and_print",
    "compute_axis_and_print",
]
