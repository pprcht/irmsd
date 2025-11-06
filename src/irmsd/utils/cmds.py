from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

import numpy as np

try:
    from .ase_io import get_axis_ase, get_canonical_ase, get_cn_ase
except Exception:  # pragma: no cover
    get_cn_ase = None  # type: ignore
    get_axis_ase = None  # type: ignore
    get_canonical_ase = None  # type: ignore

from .utils import print_array, require_ase

if TYPE_CHECKING:
    from ase import Atoms  # type: ignore

__all__ = ["compute_cn_and_print", "compute_axis_and_print"]


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


def compute_axis_and_print(
    atoms_list: List["Atoms"],
) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    # Ensure ASE is present only when this command is actually invoked
    require_ase()

    results: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for i, atoms in enumerate(atoms_list, start=1):
        if get_axis_ase is not None:
            rot, avmom, evec = get_axis_ase(atoms)
        else:
            rot, avmom, evec = None, None, None
        results.append((rot, avmom, evec))
        print_array(f"Rotational constants (MHz) for structure {i}", rot)
        print(f"Average momentum a.u. (10⁻⁴⁷kg m²) for structure {i}: {avmom}")
        print_array(f"Rotation matrix for structure {i}", evec)
    return results


def compute_canonical_and_print(atoms_list: List["Atoms"]):
    # Ensure ASE is present only when this command is actually invoked
    require_ase()

    results: List[Tuple[np.ndarray, np.ndarray]] = []
    for i, atoms in enumerate(atoms_list, start=1):
        if get_canonical_ase is not None:
            rank, invariants = get_canonical_ase(atoms)
        else:
            rank, invariants = None, None
        results.append((rank, invariants))
        print_array(f"Canonical rank for structure {i}", rank)
    return results
