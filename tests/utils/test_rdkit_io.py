import numpy as np
import pytest

pytest.importorskip("rdkit")

from rdkit import Chem

from irmsd.utils.rdkit_io import (
    get_axis_rdkit,
    get_canonical_rdkit,
    get_cn_rdkit,
    get_irmsd_rdkit,
    get_rmsd_rdkit,
)


def test_get_cn_rdkit(caffeine_cn_test_data):
    caffeine_xyz, expected_cn = caffeine_cn_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    cn = get_cn_rdkit(caffeine_rdkit)
    assert pytest.approx(expected_cn, abs=1e-6) == cn


def test_get_axis_rdkit(caffeine_axis_test_data):
    caffeine_xyz, expected_rot, expected_avmom, expected_evec = caffeine_axis_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    rot, avmom, evec = get_axis_rdkit(caffeine_rdkit)
    assert pytest.approx(expected_rot, abs=1e-6) == rot
    assert pytest.approx(expected_avmom, abs=1e-6) == avmom
    assert pytest.approx(expected_evec, abs=1e-6) == evec


def test_get_canonical_rdkit(caffeine_canonical_test_data):
    caffeine_xyz, heavy, expected_canonical = caffeine_canonical_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    canonical = get_canonical_rdkit(caffeine_rdkit, heavy=heavy)
    assert pytest.approx(expected_canonical, abs=1e-6) == canonical


def test_get_rmsd_rdkit(caffeine_rmsd_test_data):
    (
        caffeine_xyz1,
        caffeine_xyz2,
        heavy,
        expected_rmsd,
        caffeine_aligned,
        expected_Umat,
    ) = caffeine_rmsd_test_data

    caffeine_rdkit_1 = Chem.MolFromXYZBlock(caffeine_xyz1)
    caffeine_rdkit_2 = Chem.MolFromXYZBlock(caffeine_xyz2)
    expected_aligned_mol = Chem.MolFromXYZBlock(caffeine_aligned)

    atomic_numbers = np.asarray(
        [atom.GetAtomicNum() for atom in caffeine_rdkit_1.GetAtoms()]
    )
    mask = atomic_numbers != 1 if heavy else None  # exclude hydrogens
    rmsd, aligned_mol, Umat = get_rmsd_rdkit(
        caffeine_rdkit_1, caffeine_rdkit_2, mask=mask
    )
    assert pytest.approx(expected_rmsd, abs=1e-6) == rmsd
    assert pytest.approx(expected_Umat, abs=1e-6) == Umat
    assert (
        pytest.approx(expected_aligned_mol.GetConformer().GetPositions(), abs=1e-6)
        == aligned_mol.GetConformer(1).GetPositions()
    )


def test_get_irmsd_rdkit(caffeine_irmsd_test_data):
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )

    mol1 = Chem.MolFromXYZBlock(conformer1)
    mol2 = Chem.MolFromXYZBlock(conformer2)
    expected_aligned_mol = Chem.MolFromXYZBlock(expected_aligned_conformer)

    rmsd, aligned_mol = get_irmsd_rdkit(mol1, mol2, iinversion=1)

    assert pytest.approx(expected_irmsd, abs=1e-6) == rmsd
