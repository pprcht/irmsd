# src/irmsd/__init__.py
from __future__ import annotations

from ._version import __version__
from . import sorting
from .api.axis_exposed import get_axis
from .api.canonical_exposed import get_canonical_fortran
from .api.cn_exposed import get_cn_fortran
from .api.irmsd_exposed import get_irmsd
from .api.rmsd_exposed import get_quaternion_rmsd_fortran
from .api.sorter_exposed import delta_irmsd_list, sorter_irmsd
from .core.molecule import Molecule

# ---- Core API ----------------------------------------
from .interfaces.mol_interface import (
    delta_irmsd_list_molecule,
    get_energies_from_molecule_list,
    get_irmsd_molecule,
    get_rmsd_molecule,
    sorter_irmsd_molecule,
)

# Try to expose ase_to_fortran if ASE is present; otherwise export a stub that errors nicely.
try:
    from .interfaces.ase_io import (
        ase_to_molecule,
        get_axis_ase,
        get_canonical_ase,
        get_cn_ase,
        get_irmsd_ase,
        get_rmsd_ase,
        molecule_to_ase,
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
        molecule_to_rdkit,
        rdkit_to_molecule,
        sorter_irmsd_rdkit,
    )
except Exception:
    pass

from .interfaces.cmds import (
    compute_axis_and_print,
    compute_canonical_and_print,
    compute_cn_and_print,
    compute_irmsd_and_print,
    compute_quaternion_rmsd_and_print,
    sort_get_delta_irmsd_and_print,
    sort_structures_and_print,
    run_cregen_and_print,
)
from .interfaces.mol_interface import (
    delta_irmsd_list_molecule,
    get_energies_from_molecule_list,
    prune_by_energy_window,
    get_irmsd_molecule,
    get_rmsd_molecule,
    sorter_irmsd_molecule,
    cregen,
    prune,
)

# ---- New: re-export Python utilities ----------------------------------------
from .utils.io import read_structures, write_structures

__all__ = [
    "__version__",
    "Molecule",
    "delta_irmsd_list_molecule",
    "get_energies_from_molecule_list",
    "get_irmsd_molecule",
    "get_rmsd_molecule",
    "sorter_irmsd_molecule",
    # core API
    "get_cn_fortran",
    "get_axis",
    "get_canonical_fortran",
    "get_quaternion_rmsd_fortran",
    "get_irmsd",
    "sorter_irmsd",
    # molecule utils
    "get_rmsd_molecule",
    "get_irmsd_molecule",
    "delta_irmsd_list_molecule",
    "sorter_irmsd_molecule",
    "get_energies_from_molecule_list",
    "cregen",
    "prune",
    # ase utils
    "ase_to_molecule",
    "molecule_to_ase",
    "get_cn_ase",
    "get_axis_ase",
    "get_canonical_ase",
    "get_irmsd_ase",
    "get_rmsd_ase",
    "sorter_irmsd_ase",
    "prune_ase",
    # rdkit utils
    "rdkit_to_molecule",
    "molecule_to_rdkit",
    "get_cn_rdkit",
    "get_axis_rdkit",
    "get_canonical_rdkit",
    "get_rmsd_rdkit",
    "get_irmsd_rdkit",
    "sorter_irmsd_rdkit",
    "prune_rdkit",
    # sorting
    "sorting",
    # optional cmds
    "read_structures",
    "write_structures",
    "compute_cn_and_print",
    "compute_axis_and_print",
    "compute_canonical_and_print",
    "compute_quaternion_rmsd_and_print",
    "compute_irmsd_and_print",
    "sort_structures_and_print",
    "sort_get_delta_irmsd_and_print",
    "run_cregen_and_print",
]
