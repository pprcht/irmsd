# src/irmsd/__init__.py
from __future__ import annotations

from .api.axis_exposed import get_axis_fortran
from .api.canonical_exposed import get_canonical_fortran
from .api.cn_exposed import get_cn_fortran

# ---- Core API (what you already had) ----------------------------------------
from .api.xyz_bridge import xyz_to_fortran, xyz_to_fortran_pair

# Try to expose ase_to_fortran if ASE is present; otherwise export a stub that errors nicely.
try:
    from .utils.ase_io import (
        ase_to_fortran,
        ase_to_fortran_pair,
        get_axis_ase,
        get_canonical_ase,
        get_cn_ase,
        get_irmsd_ase,
        get_rmsd_ase,
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
    from .utils.cmds import (
        compute_axis_and_print,
        compute_canonical_and_print,
        compute_cn_and_print,
        compute_irmsd_and_print,
        compute_quaternion_rmsd_and_print,
    )
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
    "get_axis_fortran",
    "get_canonical_fortran",
    "get_quaternion_rmsd_fortran",
    "get_cn_ase",
    "get_axis_ase",
    "get_canonical_ase",
    "get_rmsd_ase" "get_irmsd_ase",
    # utils
    "read_structures",
    "print_array",
    # optional cmds
    "compute_cn_and_print",
    "compute_axis_and_print",
    "compute_canonical_and_print",
    "compute_quaternion_rmsd_and_print",
    "compute_irmsd_and_print",
]
