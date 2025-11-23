# tests/test_utils_xyz_io.py

import numpy as np
import pytest

from irmsd.core.molecule import Molecule
from irmsd.utils.xyz import read_extxyz, write_extxyz


def test_read_single_extxyz(tmp_path):
    """Read a simple single-frame extended XYZ and check Molecule fields."""
    content = """3
energy=-40.5 cell="1.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 3.0" pbc="T F T" comment=foo
C 0.0 0.0 0.0
H 0.0 0.0 1.0
H 0.0 1.0 0.0
"""
    path = tmp_path / "single.extxyz"
    path.write_text(content)

    mol = read_extxyz(path)
    assert isinstance(mol, Molecule)
    assert mol.natoms == 3

    # symbols and positions
    assert mol.get_chemical_symbols() == ["C", "H", "H"]
    pos = mol.get_positions()
    assert pos.shape == (3, 3)
    np.testing.assert_allclose(pos[0], [0.0, 0.0, 0.0])

    # energy
    assert mol.energy == pytest.approx(-40.5)

    # cell
    assert mol.cell is not None
    assert mol.cell.shape == (3, 3)
    np.testing.assert_allclose(mol.cell, np.diag([1.0, 2.0, 3.0]))

    # pbc
    assert mol.pbc == (True, False, True)

    # info
    assert "comment" in mol.info
    assert mol.info["comment"] == "foo"


def test_read_multi_extxyz(tmp_path):
    """Read a multi-frame extended XYZ and ensure a list of Molecules is returned."""
    content = """2
energy=-1.0
H 0.0 0.0 0.0
H 0.0 0.0 0.8
2
energy=-2.0
H 0.0 0.0 0.0
H 0.0 0.0 1.0
"""
    path = tmp_path / "multi.extxyz"
    path.write_text(content)

    mols = read_extxyz(path)
    assert isinstance(mols, list)
    assert len(mols) == 2

    assert mols[0].natoms == 2
    assert mols[1].natoms == 2

    assert mols[0].energy == pytest.approx(-1.0)
    assert mols[1].energy == pytest.approx(-2.0)


def test_write_and_read_roundtrip_single(tmp_path):
    """Write a single Molecule to extxyz and read it back."""
    symbols = ["C", "H", "H", "H", "H"]
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.1],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        dtype=float,
    )
    cell = np.diag([5.0, 5.0, 5.0])

    mol = Molecule(
        symbols=symbols,
        positions=positions,
        energy=-40.123,
        info={"tag": "test_molecule", "spin": 1},
        cell=cell,
        pbc=(True, True, True),
    )

    path = tmp_path / "roundtrip_single.extxyz"
    write_extxyz(path, mol)

    mol_rt = read_extxyz(path)
    assert isinstance(mol_rt, Molecule)
    assert mol_rt.natoms == mol.natoms
    assert mol_rt.get_chemical_symbols() == mol.get_chemical_symbols()
    np.testing.assert_allclose(mol_rt.get_positions(), mol.get_positions())
    assert mol_rt.energy == pytest.approx(mol.energy)
    np.testing.assert_allclose(mol_rt.cell, mol.cell)
    assert mol_rt.pbc == mol.pbc

    # info should survive for non-reserved keys
    assert mol_rt.info.get("tag") == "test_molecule"
    assert mol_rt.info.get("spin") == 1


def test_write_and_read_roundtrip_multi(tmp_path):
    """Write a list of Molecule objects and read them back as a list."""
    positions1 = np.zeros((2, 3), float)
    positions2 = np.ones((2, 3), float)

    mol1 = Molecule(symbols=["H", "H"], positions=positions1, energy=-1.0)
    mol2 = Molecule(symbols=["H", "H"], positions=positions2, energy=-2.0)

    path = tmp_path / "roundtrip_multi.extxyz"
    write_extxyz(path, [mol1, mol2])

    mols_rt = read_extxyz(path)
    assert isinstance(mols_rt, list)
    assert len(mols_rt) == 2

    assert mols_rt[0].energy == pytest.approx(-1.0)
    assert mols_rt[1].energy == pytest.approx(-2.0)

    np.testing.assert_allclose(mols_rt[0].get_positions(), positions1)
    np.testing.assert_allclose(mols_rt[1].get_positions(), positions2)
