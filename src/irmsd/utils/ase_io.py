from __future__ import annotations

from typing import Tuple

import numpy as np

from .utils import require_ase


def ase_to_fortran(atoms) -> Tuple["Atoms", np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core xyz_to_fortran, and returns:
      - a NEW Atoms with updated positions
      - the 3×3 matrix.
    """

    require_ase()

    from ..api.xyz_bridge import xyz_to_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("ase_to_fortran expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

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
    - The corresponding core API is `irmsd.api.xyz_bridge.xyz_to_fortran_pair`, which takes
      plain NumPy arrays and has no ASE dependency.
    """

    require_ase()

    from ..api.xyz_bridge import xyz_to_fortran_pair

    if not isinstance(atoms1, Atoms) or not isinstance(atoms2, Atoms):
        raise TypeError("ase_to_fortran_pair expects two ase.Atoms objects")

    Z1 = atoms1.get_atomic_numbers()  # (N1,)
    P1 = atoms1.get_positions()  # (N1, 3)
    Z2 = atoms2.get_atomic_numbers()  # (N2,)
    P2 = atoms2.get_positions()  # (N2, 3)

    new_P1, new_P2, M1, M2 = xyz_to_fortran_pair(Z1, P1, Z2, P2)

    # Return ONLY the modified second structure, as requested
    new_atoms2 = atoms2.copy()
    new_atoms2.set_positions(new_P2, apply_constraint=False)

    return new_atoms2, M1, M2


#####################################################################################


def get_cn_ase(atoms) -> Tuple["Atoms", np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core get_cn_fortran, and
    returns a numpy array with the coordination numbers per atom
    """

    require_ase()

    from ase import Atoms

    from ..api.cn_exposed import get_cn_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("get_cn_ase expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    new_cn = get_cn_fortran(Z, pos)

    return new_cn


def get_axis_ase(atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core get_axis_fortran, and
    returns a numpy array with the rotation constants in Mhz, the average  momentum in a.u.
    and the rot. matrix.
    """
    require_ase()

    from ase import Atoms

    from ..api.axis_exposed import get_axis_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("get_cn_ase expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    rot, avmom, evec = get_axis_fortran(Z, pos)

    return rot, avmom, evec


def get_canonical_ase(
    atoms, wbo=None, invtype="apsp+", heavy: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Optional utility: accepts ASE Atoms, adapts to core get_canonical_fortran, and
    returns rank and invariants arrays.
    """

    require_ase()

    from ase import Atoms

    from ..api.canonical_exposed import get_canonical_fortran

    if not isinstance(atoms, Atoms):
        raise TypeError("get_canonical_ase expects an ase.Atoms object")

    Z = atoms.get_atomic_numbers()  # (N,)
    pos = atoms.get_positions()  # (N, 3) float64

    rank, invariants = get_canonical_fortran(
        Z, pos, wbo=wbo, invtype=invtype, heavy=heavy
    )

    return rank, invariants
