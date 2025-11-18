# src/irmsd/__init__.py
from __future__ import annotations

from .api.axis_exposed import get_axis_fortran
from .api.canonical_exposed import get_canonical_fortran
from .api.cn_exposed import get_cn_fortran
from .api.irmsd_exposed import get_irmsd
from .api.sorter_exposed import sorter_irmsd

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
        sorter_irmsd_ase,
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


try:
    from .utils.rdkit_io import (
        get_axis_rdkit,
        get_canonical_rdkit,
        get_cn_rdkit,
        get_rmsd_rdkit,
        rdkit_to_fortran,
        rdkit_to_fortran_pair,
    )
except Exception:

    def rdkit_to_fortran(*args, **kwargs):
        raise ImportError(
            "RDKit-based helper not available. Install optional extra: pip install 'irmsd[rdkit]'"
        )

    def rdkit_to_fortran_pair(*args, **kwargs):
        raise ImportError(
            "RDKit-based helper not available. Install optional extra: pip install 'irmsd[rdkit]'"
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
    "get_irmsd",
    "sorter_irmsd",
    # ase utils
    "get_cn_ase",
    "get_axis_ase",
    "get_canonical_ase",
    "get_irmsd_ase",
    "get_rmsd_ase",
    "sorter_irmsd_ase",
    # rdkit utils
    "get_cn_rdkit",
    "get_axis_rdkit",
    "get_canonical_rdkit",
    "get_rmsd_rdkit",
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
