# src/irmsd/__init__.py
from __future__ import annotations

from .core.molecule import Molecule
from .api.axis_exposed import get_axis
from .api.canonical_exposed import get_canonical_fortran
from .api.cn_exposed import get_cn_fortran
from .api.irmsd_exposed import get_irmsd
from .api.sorter_exposed import sorter_irmsd, delta_irmsd_list
from . import sorting

# ---- Core API ----------------------------------------
# Try to expose ase_to_fortran if ASE is present; otherwise export a stub that errors nicely.
try:
    from .interfaces.ase_io import (
        ase_to_molecule,
        get_axis_ase,
        get_canonical_ase,
        get_cn_ase,
        get_irmsd_ase,
        get_rmsd_ase,
        sorter_irmsd_ase,
    )
except Exception:
    pass
# Same for RDkit
try:
    from .interfaces.rdkit_io import (
        get_axis_rdkit,
        get_canonical_rdkit,
        get_cn_rdkit,
        get_irmsd_rdkit,
        get_rmsd_rdkit,
    )
except Exception:
    pass

# ---- New: re-export Python utilities ----------------------------------------
from .utils.utils import print_array, read_structures

# ---- Optional: command helpers ----------
try:
    from .interfaces.cmds import (
        compute_axis_and_print,
        compute_canonical_and_print,
        compute_cn_and_print,
        compute_irmsd_and_print,
        compute_quaternion_rmsd_and_print,
        sort_structures_and_print,
        sort_get_delta_irmsd_and_print,
    )
except Exception:
    # Safe to ignore so `import irmsd` never breaks due to optional pieces.
    pass

__all__ = [
    "Molecule",
    # core API
    "get_cn_fortran",
    "get_axis",
    "get_canonical_fortran",
    "get_quaternion_rmsd_fortran",
    "get_irmsd",
    "sorter_irmsd",
    # ase utils
    "ase_to_molecule",
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
    "get_irmsd_rdkit",
    # sorting
    "sorting",
    # optional cmds
    "compute_cn_and_print",
    "compute_axis_and_print",
    "compute_canonical_and_print",
    "compute_quaternion_rmsd_and_print",
    "compute_irmsd_and_print",
    "sort_structures_and_print",
    "sort_get_delta_irmsd_and_print",
]
