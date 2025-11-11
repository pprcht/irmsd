from __future__ import annotations

from typing import Tuple

import numpy as np

from .utils import require_rdkit


def rdkit_to_fortran():
    raise NotImplementedError("rdkit_to_fortran is not yet implemented")


def rdkit_to_fortran_pair():
    raise NotImplementedError("rdkit_to_fortran_pair is not yet implemented")


def get_cn_rdkit():
    raise NotImplementedError("get_cn_rdkit is not yet implemented")


def get_axis_rdkit():
    raise NotImplementedError("get_axis_rdkit is not yet implemented")


# TODO: Dicsuss. There are two ways how to deal with this, rdkits Conformer objects only contain
#       atom positons, but not atomic numbers. We can either require the user to pass in the Mol object
#       and then specify via indeces the conformers of the Mol to be compared, or we can require the user
#       to pass in two Conformers and then call GetOwningMol() to get the parent Mol object from which we can
#       get the atomic numbers.
def get_rmsd_rdkit(
    molecule, conf_id_ref=0, conf_id_align=1, mask=None
) -> Tuple[np.ndarray, "Mol", np.ndarray]:

    require_rdkit()

    from rdkit import Chem

    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("get_rmsd_rdkit expects rdkit.Chem.Mol objects")

    conformer_ref = molecule.GetConformer(conf_id_ref)
    conformer_align = molecule.GetConformer(conf_id_align)

    if not conformer_ref.Is3D() or not conformer_align.Is3D():
        raise ValueError("get_rmsd_rdkit expects 3D conformers")

    Z1 = Z2 = molecule.GetAtomicNumbers()
    P1 = conformer_ref.GetPositions()

    P2 = conformer_align.GetPositions()

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    # Also what do we want to return? A new conformer object? Or a full new molecule?
    # better a copied molecule
    molecule_ret = molecule.Copy()
    molecule_ret.GetConformer(conf_id_align).SetPositions(new_P2)

    return rmsdval, molecule_ret, umat


def get_rmsd_rdkit(
    conformer1, conformer2, mask=None
) -> Tuple[np.ndarray, "Conformer", np.ndarray]:
    require_rdkit()

    from rdkit import Chem

    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(conformer1, Chem.Conformer) or not isinstance(
        conformer2, Chem.Conformer
    ):
        raise TypeError("get_rmsd_rdkit expects rdkit.Chem.Mol objects")

    if not conformer1.Is3D() or not conformer2.Is3D():
        raise ValueError("get_rmsd_rdkit expects 3D conformers")

    Z1 = conformer1.GetOwningMol().GetAtomicNumbers()
    P1 = conformer1.GetPositions()

    Z2 = conformer2.GetOwningMol().GetAtomicNumbers()
    P2 = conformer2.GetPositions()

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    conformer2_ret = conformer2.Clone()
    conformer2_ret.SetPositions(new_P2)

    return rmsdval, conformer2_ret, umat
