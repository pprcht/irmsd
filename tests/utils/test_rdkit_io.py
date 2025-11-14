import pytest

pytest.importorskip("rdkit")

from rdkit import Chem

from irmsd.utils.rdkit_io import (
    get_axis_rdkit,
    get_canonical_rdkit,
    get_cn_rdkit,
    get_rmsd_rdkit,
)


def test_get_cn_rdkit(caffeine_cn_test_data):
    caffeine_xyz, expected_cn = caffeine_cn_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    cn = get_cn_rdkit(caffeine_rdkit)
    assert pytest.approx(cn, abs=1e-6) == expected_cn


def test_get_axis_rdkit(caffeine_axis_test_data):
    caffeine_xyz, expected_rot, expected_avmom, expected_evec = caffeine_axis_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    rot, avmom, evec = get_axis_rdkit(caffeine_rdkit)
    assert pytest.approx(rot, abs=1e-6) == expected_rot
    assert pytest.approx(avmom, abs=1e-6) == expected_avmom
    assert pytest.approx(evec, abs=1e-6) == expected_evec


def test_get_canonical_rdkit(caffeine_canonical_test_data):
    caffeine_xyz, heavy, expected_canonical = caffeine_canonical_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    canonical = get_canonical_rdkit(caffeine_rdkit, heavy=heavy)
    assert pytest.approx(canonical, abs=1e-6) == expected_canonical


def test_get_rmsd_rdkit(caffeine_rmsd_test_data):
    caffeine_xyz1, caffeine_xyz2, expected_rmsd, caffeine_aligned, expected_Umat = (
        caffeine_rmsd_test_data
    )

    caffeine_rdkit_1 = Chem.MolFromXYZBlock(caffeine_xyz1)
    caffeine_rdkit_2 = Chem.MolFromXYZBlock(caffeine_xyz2)
    expected_aligned_mol = Chem.MolFromXYZBlock(caffeine_aligned)

    rmsd, aligned_mol, Umat = get_rmsd_rdkit(caffeine_rdkit_1, caffeine_rdkit_2)
    assert pytest.approx(rmsd, abs=1e-6) == expected_rmsd
    assert pytest.approx(Umat, abs=1e-6) == expected_Umat
    assert (
        pytest.approx(aligned_mol.GetConformer().GetPositions(), abs=1e-6)
        == expected_aligned_mol.GetConformer().GetPositions()
    )
