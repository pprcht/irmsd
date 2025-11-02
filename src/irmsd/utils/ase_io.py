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

###############################################################################

def ase_to_fortran_pair(atoms1, atoms2) -> Tuple["Atoms", np.ndarray, np.ndarray]:
    """
    Optional ASE utility: operate on TWO ASE Atoms. Returns ONLY the modified second Atoms
    plus both 3×3 matrices produced by the Fortran routine.

    Parameters
    ----------
    atoms1 : ase.Atoms
    atoms2 : ase.Atoms

    Returns
    -------
    new_atoms2 : ase.Atoms
        Copy of `atoms2` with updated positions (atoms1 is NOT returned).
    mat1 : (3, 3) float64 ndarray (Fortran-ordered)
        Matrix corresponding to the first structure.
    mat2 : (3, 3) float64 ndarray (Fortran-ordered)
        Matrix corresponding to the second structure.

    Notes
    -----
    - ASE is an optional dependency. If it's not installed, this function raises ImportError
      suggesting to install the extra: `pip install 'irmsd[ase]'`.
    - The corresponding core API is `irmsd.api.xyz_bridge.xyz_to_fortran_pair`, which takes
      plain NumPy arrays and has no ASE dependency.
    """
    # Lazy import to keep ASE optional
    try:
        from ase import Atoms  # type: ignore
    except Exception as exc:
        raise ImportError(
            "ASE is required for `ase_to_fortran_pair(atoms1, atoms2)`. "
            "Install the optional extra: pip install 'irmsd[ase]'"
        ) from exc

    from ..api.xyz_bridge import xyz_to_fortran_pair

    if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
        raise TypeError("ase_to_fortran_pair expects two ase.Atoms objects")

    Z1 = atoms1.get_atomic_numbers()          # (N1,)
    P1 = atoms1.get_positions()               # (N1, 3)
    Z2 = atoms2.get_atomic_numbers()          # (N2,)
    P2 = atoms2.get_positions()               # (N2, 3)

    new_P1, new_P2, M1, M2 = xyz_to_fortran_pair(Z1, P1, Z2, P2)

    # Return ONLY the modified second structure, as requested
    new_atoms2 = atoms2.copy()
    new_atoms2.set_positions(new_P2, apply_constraint=False)

    return new_atoms2, M1, M2

