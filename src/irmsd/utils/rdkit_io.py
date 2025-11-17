from __future__ import annotations

import copy
from typing import Sequence, Tuple

import numpy as np

from .utils import require_rdkit


def conformer_iterator(molecule: "Mol", conf_ids: list[int]) -> "Conformer":
    for conf_id in conf_ids:
        yield molecule.GetConformer(conf_id)


def conf_id_to_iterator(molecule: "Mol", conf_id: None | int | Sequence) -> "Conformer":
    if conf_id is None:
        conf_iterator = molecule.GetConformers()
    elif isinstance(conf_id, int):
        conf_iterator = [molecule.GetConformer(conf_id)]
    elif isinstance(conf_id, list):
        conf_iterator = conformer_iterator(molecule, conf_id)
    else:
        raise TypeError("conf_id must be None, int, or list of int")
    return conf_iterator


def get_atom_numbers_rdkit(molecule) -> np.ndarray:
    Z = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    return np.array(Z, dtype=np.int32)


def rdkit_to_fortran(molecule, conf_id=-1) -> Tuple["Mol", np.ndarray]:
    require_rdkit()

    from rdkit import Chem

    from ..api.xyz_bridge import xyz_to_fortran

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran expects rdkit.Chem.Mol objects")

    Z = get_atom_numbers_rdkit(molecule)  # (N,)
    conformer = molecule.GetConformer(id=conf_id)

    if not conformer.Is3D():
        raise ValueError("rdkit_to_fortran expects a 3D conformer")

    pos = conformer.GetPositions()  # (N, 3) float64

    new_pos, mat = xyz_to_fortran(Z, pos)

    new_molecule = molecule.Copy()
    new_molecule.GetConformer(id=conf_id).SetPositions(new_pos)

    return new_molecule, mat


def rdkit_to_fortran_pair(
    molecule1, molecule2, conf_id1=-1, conf_id2=-1
) -> Tuple["Mol", np.ndarray, np.ndarray]:
    require_rdkit()

    from rdkit import Chem

    from ..api.xyz_bridge import xyz_to_fortran_pair

    if not isinstance(molecule1, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    if not isinstance(molecule2, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    Z1 = molecule1.GetAtomicNumbers()  # (N,)
    conformer1 = molecule1.GetConformer(id=conf_id1)

    if not conformer1.Is3D():
        raise ValueError("rdkit_to_fortran_pair expects a 3D conformer")

    pos1 = conformer1.GetPositions()  # (N, 3) float64

    Z2 = molecule2.GetAtomicNumbers()  # (N,)
    conformer2 = molecule2.GetConformer(id=conf_id2)

    if not conformer2.Is3D():
        raise ValueError("rdkit_to_fortran_pair expects a 3D conformer")

    pos2 = conformer2.GetPositions()  # (N, 3) float64

    new_pos1, new_pos2, mat1, mat2 = xyz_to_fortran_pair(Z1, pos1, Z2, pos2)

    new_molecule2 = molecule2.Copy()
    new_molecule2.GetConformer(id=conf_id2).SetPositions(new_pos2)

    return new_molecule2, mat1, mat2


def get_cn_rdkit(molecule, conf_id: None | int | Sequence = None) -> np.ndarray:
    """Optional RDKit utility: compute coordination numbers for one or more conformers of a molecule
    and return as a numpy array with shape (n_conf, n_atoms)."""

    require_rdkit()

    from rdkit import Chem

    from ..api.cn_exposed import get_cn_fortran

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    Z = get_atom_numbers_rdkit(molecule)  # (N,)

    conf_iterator = conf_id_to_iterator(molecule, conf_id)

    all_cn = []
    for conformer in conf_iterator:
        if not conformer.Is3D():
            raise ValueError("get_cn_rdkit expects 3D conformers")

        pos = conformer.GetPositions()
        cn = get_cn_fortran(Z, pos)
        all_cn.append(cn)
    return np.array(all_cn).squeeze()


def get_axis_rdkit(
    molecule, conf_id: None | int | Sequence = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Optional RDKit utility: compute principal axes for one or more conformers of a molecule"""
    require_rdkit()

    from rdkit import Chem

    from ..api.axis_exposed import get_axis_fortran

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    Z = get_atom_numbers_rdkit(molecule)  # (N,)

    conf_iterator = conf_id_to_iterator(molecule, conf_id)

    all_rot = []
    all_avmom = []
    all_evec = []

    for conformer in conf_iterator:
        if not conformer.Is3D():
            raise ValueError("get_axis_rdkit expects 3D conformers")

        pos = conformer.GetPositions()
        rot, avmom, evec = get_axis_fortran(Z, pos)
        all_rot.append(rot)
        all_avmom.append(avmom)
        all_evec.append(evec)

    return (
        np.array(all_rot).squeeze(),
        np.array(all_avmom).squeeze(),
        np.array(all_evec).squeeze(),
    )


def get_canonical_rdkit(
    molecule,
    conf_id: None | int | Sequence = None,
    wbo=None,
    invtype="apsp+",
    heavy: bool = False,
) -> np.ndarray:
    """Optional RDKit utility: compute canonical ranks for one or more conformers of a molecule"""
    require_rdkit()

    from rdkit import Chem

    from ..api.canonical_exposed import get_canonical_fortran

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    Z = get_atom_numbers_rdkit(molecule)  # (N,)

    conf_iterator = conf_id_to_iterator(molecule, conf_id)

    all_rank = []

    for conformer in conf_iterator:
        if not conformer.Is3D():
            raise ValueError("get_canonical_rdkit expects 3D conformers")

        pos = conformer.GetPositions()
        rank = get_canonical_fortran(Z, pos, wbo=wbo, invtype=invtype, heavy=heavy)
        all_rank.append(rank)

    return np.array(all_rank).squeeze()


def get_rmsd_rdkit(
    molecule_ref, molecule_align, conf_id_ref=-1, conf_id_align=-1, mask=None
) -> Tuple[float, "Mol", np.ndarray]:
    """Optional Rdkit utility: operate on two Rdkit Molecules. Returns the RMSD in Angström,
    the molecule object with both Conformers aligned.
    """

    require_rdkit()

    from rdkit import Chem

    from ..api.rmsd_exposed import get_quaternion_rmsd_fortran

    if not isinstance(molecule_ref, Chem.Mol) or not isinstance(
        molecule_align, Chem.Mol
    ):
        raise TypeError("get_rmsd_rdkit expects rdkit.Chem.Mol objects")

    conformer_ref = molecule_ref.GetConformer(conf_id_ref)
    conformer_align = molecule_align.GetConformer(conf_id_align)

    if not conformer_ref.Is3D() or not conformer_align.Is3D():
        raise ValueError("get_rmsd_rdkit expects 3D conformers")

    Z1 = get_atom_numbers_rdkit(molecule_ref)
    P1 = conformer_ref.GetPositions()

    Z2 = get_atom_numbers_rdkit(molecule_align)
    P2 = conformer_align.GetPositions()

    rmsdval, new_P2, umat = get_quaternion_rmsd_fortran(Z1, P1, Z2, P2, mask=mask)

    molecule_ret = Chem.Mol(molecule_ref, confId=conf_id_ref)
    align_id = molecule_ret.AddConformer(conformer_align, assignId=True)
    molecule_ret.GetConformer(align_id).SetPositions(new_P2)

    return rmsdval, molecule_ret, umat


def get_irmsd_rdkit(
    molecule_ref, molecule_align, conf_id_ref=-1, conf_id_align=-1, iinversion: int = 0
) -> Tuple[float, "Mol"]:
    """
    Optional Rdkit utility: operate on TWO Rdkit Molecules. Returns the iRMSD in Angström,
    the molecule object with both Conformers permuted and aligned.
    """
    require_rdkit()

    from rdkit import Chem

    from ..api.irmsd_exposed import get_irmsd_fortran

    if not isinstance(molecule_ref, Chem.Mol) or not isinstance(
        molecule_align, Chem.Mol
    ):
        raise TypeError("get_rmsd_rdkit expects rdkit.Chem.Mol objects")

    conformer_ref = molecule_ref.GetConformer(conf_id_ref)
    conformer_align = molecule_align.GetConformer(conf_id_align)

    if not conformer_ref.Is3D() or not conformer_align.Is3D():
        raise ValueError("get_rmsd_rdkit expects 3D conformers")

    Z1 = get_atom_numbers_rdkit(molecule_ref)
    P1 = conformer_ref.GetPositions()

    Z2 = get_atom_numbers_rdkit(molecule_align)
    P2 = conformer_align.GetPositions()

    irmsdval, new_Z1, new_P1, new_Z2, new_P2 = get_irmsd_fortran(
        Z1, P1, Z2, P2, iinversion=iinversion
    )

    molecule_ret = Chem.Mol(molecule_ref, confId=conf_id_ref)
    molecule_ret.GetConformer(conf_id_ref).SetPositions(new_P1)

    align_id = molecule_ret.AddConformer(conformer_align, assignId=True)
    molecule_ret.GetConformer(align_id).SetPositions(new_P2)

    return irmsdval, molecule_ret
