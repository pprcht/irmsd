from __future__ import annotations
from typing import Tuple
import numpy as np

def ase_to_fortran(atoms) -> Tuple["Atoms", np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core xyz_to_fortran, and returns:
      - a NEW Atoms with updated positions
      - the 3×3 matrix.

    Raises a helpful ImportError if ASE is not available.
    """
    try:
        from ase import Atoms  # type: ignore
    except Exception as exc:
        raise ImportError(
            "ASE is required for `ase_to_fortran(atoms)`. "
            "Install the optional extra: pip install 'irmsd[ase]'"
        ) from exc

    from ..api.xyz_bridge import xyz_to_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("ase_to_fortran expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()            # (N,)
    pos = atoms.get_positions()               # (N, 3) float64

    new_pos, mat = xyz_to_fortran(Z, pos)

    new_atoms = atoms.copy()
    new_atoms.set_positions(new_pos, apply_constraint=False)
    return new_atoms, mat

