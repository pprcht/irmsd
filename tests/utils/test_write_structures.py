# tests/test_write_structures.py

import numpy as np
import pytest

from irmsd.core import Molecule
from irmsd.utils.xyz import read_extxyz
from irmsd.utils.io import write_structures


def test_write_structures_extxyz_single(tmp_path):
    """write_structures should use internal write_extxyz for .xyz (single Molecule)."""
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
        info={"tag": "extxyz_single"},
        cell=cell,
        pbc=(True, True, True),
    )

    path = tmp_path / "out.xyz"
    write_structures(path, mol)

    # Read back via our own extended XYZ reader
    mol_rt = read_extxyz(path)
    assert isinstance(mol_rt, Molecule)
    assert mol_rt.natoms == mol.natoms
    assert mol_rt.get_chemical_symbols() == mol.get_chemical_symbols()
    np.testing.assert_allclose(mol_rt.get_positions(), mol.get_positions())
    assert mol_rt.energy == pytest.approx(mol.energy)
    np.testing.assert_allclose(mol_rt.cell, mol.cell)
    assert mol_rt.pbc == mol.pbc
    assert mol_rt.info.get("tag") == "extxyz_single"


def test_write_structures_extxyz_multi(tmp_path):
    """write_structures should write multiple Molecules into one extxyz trajectory."""
    pos1 = np.zeros((2, 3), float)
    pos2 = np.ones((2, 3), float)

    mol1 = Molecule(symbols=["H", "H"], positions=pos1, energy=-1.0, info={"id": 1})
    mol2 = Molecule(symbols=["H", "H"], positions=pos2, energy=-2.0, info={"id": 2})

    path = tmp_path / "traj.extxyz"
    write_structures(path, [mol1, mol2])

    mols_rt = read_extxyz(path)
    assert isinstance(mols_rt, list)
    assert len(mols_rt) == 2

    m1, m2 = mols_rt
    assert m1.energy == pytest.approx(-1.0)
    assert m2.energy == pytest.approx(-2.0)
    np.testing.assert_allclose(m1.get_positions(), pos1)
    np.testing.assert_allclose(m2.get_positions(), pos2)
    assert m1.info.get("id") == 1
    assert m2.info.get("id") == 2


@pytest.mark.skipif(
    pytest.importorskip("ase", reason="ASE not installed") is None,
    reason="ASE not installed",
)
def test_write_structures_ase_backend_single(tmp_path):
    """
    For non-xyz extensions, write_structures should go through ASE backend.

    We test this by writing a Molecule to a '.traj' file (handled by ASE),
    then reading it back using ASE and checking basic consistency.
    """
    import ase.io

    symbols = ["O", "H", "H"]
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.96],
            [0.0, 0.75, -0.24],
        ],
        dtype=float,
    )

    mol = Molecule(
        symbols=symbols,
        positions=positions,
        energy=-76.0,
        info={"tag": "ase_backend"},
    )

    path = tmp_path / "out.traj"  # non-xyz extension → ASE path
    write_structures(path, mol)

    # ASE's read returns a single Atoms object for a single-frame .traj
    ats = ase.io.read(str(path))
    assert ats.get_chemical_symbols() == symbols
    np.testing.assert_allclose(ats.get_positions(), positions)


@pytest.mark.skipif(
    pytest.importorskip("ase", reason="ASE not installed") is None,
    reason="ASE not installed",
)
def test_write_structures_ase_backend_multi(tmp_path):
    """Multiple Molecules with non-xyz extension should be written via ASE as trajectory."""
    import ase.io

    pos1 = np.zeros((2, 3), float)
    pos2 = np.ones((2, 3), float)

    mol1 = Molecule(symbols=["H", "H"], positions=pos1, energy=-1.0)
    mol2 = Molecule(symbols=["H", "H"], positions=pos2, energy=-2.0)

    path = tmp_path / "traj.traj"
    write_structures(path, [mol1, mol2])

    # Read all frames with ASE: ":" slice returns a list
    ats_list = ase.io.read(str(path), index=":")
    assert len(ats_list) == 2

    np.testing.assert_allclose(ats_list[0].get_positions(), pos1)
    np.testing.assert_allclose(ats_list[1].get_positions(), pos2)

