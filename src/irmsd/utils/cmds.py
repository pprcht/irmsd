from __future__ import annotations

from typing import List, TYPE_CHECKING

import numpy as np

try:
    from .ase_io import get_cn_ase
except Exception:  # pragma: no cover
    get_cn_ase = None  # type: ignore

from .utils import print_array, require_ase

if TYPE_CHECKING:
    from ase import Atoms  # type: ignore

__all__ = ["compute_cn_and_print"]


def compute_cn_and_print(atoms_list: List["Atoms"]) -> List[np.ndarray]:
    """Compute coordination numbers for each structure and print them.

    Parameters
    ----------
    atoms_list : list[ase.Atoms]
        Structures to analyze.

    Returns
    -------
    list[np.ndarray]
        One integer array per structure, same order as ``atoms_list``.
    """
    # Ensure ASE is present only when this command is actually invoked
    require_ase()

    results: List[np.ndarray] = []
    for i, atoms in enumerate(atoms_list, start=1):
        if get_cn_ase is not None:
            cn_vec = get_cn_ase(atoms)
        else:
            cn_vec = None
        results.append(cn_vec)
        print_array(f"CN[ structure {i} ] (n={len(atoms)})", cn_vec)
    return results

