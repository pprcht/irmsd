# tests/test_utils_xyz_io.py

import numpy as np
import pytest

from irmsd.core.molecule import Molecule
from irmsd.utils.xyz import (
    CANONICAL_ID_COLUMN,
    EV_PER_HARTREE,
    read_extxyz,
    write_extxyz,
)


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


# ---------------------------------------------------------------------------
# energy_units handling (Molecule.energy is always Hartree)
# ---------------------------------------------------------------------------


def test_read_energy_no_marker_is_hartree(tmp_path):
    """Without an energy_units marker the value is taken as Hartree (default)."""
    content = "1\nenergy=-76.16\nH 0.0 0.0 0.0\n"
    path = tmp_path / "bare.xyz"
    path.write_text(content)

    mol = read_extxyz(path)
    assert mol.energy == pytest.approx(-76.16)


def test_read_energy_units_hartree_passthrough(tmp_path):
    """energy_units=Hartree (any case) passes through unchanged and is consumed."""
    content = "1\nenergy=-76.16 energy_units=Hartree\nH 0.0 0.0 0.0\n"
    path = tmp_path / "ha.xyz"
    path.write_text(content)

    mol = read_extxyz(path)
    assert mol.energy == pytest.approx(-76.16)
    # marker must not leak into info
    assert "energy_units" not in mol.info


def test_read_energy_units_ev_converted(tmp_path):
    """energy_units=eV is converted to Hartree on read."""
    ev_value = -76.16 * EV_PER_HARTREE
    content = f"1\nenergy={ev_value:.12g} energy_units=eV\nH 0.0 0.0 0.0\n"
    path = tmp_path / "ev.xyz"
    path.write_text(content)

    mol = read_extxyz(path)
    assert mol.energy == pytest.approx(-76.16)
    assert "energy_units" not in mol.info


def test_write_stamps_hartree_marker(tmp_path):
    """The writer emits energy_units=Hartree alongside the energy."""
    mol = Molecule(symbols=["H"], positions=[[0.0, 0.0, 0.0]], energy=-1.234)
    path = tmp_path / "out.xyz"
    write_extxyz(path, mol)

    text = path.read_text()
    assert "energy=-1.234" in text
    assert "energy_units=Hartree" in text


def test_write_no_energy_no_marker(tmp_path):
    """No energy means no energy_units marker is written."""
    mol = Molecule(symbols=["H"], positions=[[0.0, 0.0, 0.0]])
    path = tmp_path / "noe.xyz"
    write_extxyz(path, mol)

    assert "energy_units" not in path.read_text()


def test_write_read_marker_roundtrip_idempotent(tmp_path):
    """write -> read -> write keeps the energy in Hartree and stays clean."""
    mol = Molecule(symbols=["H", "H"], positions=[[0, 0, 0], [0, 0, 0.74]], energy=-1.10584)
    p1 = tmp_path / "rt1.xyz"
    write_extxyz(p1, mol)

    back = read_extxyz(p1)
    assert back.energy == pytest.approx(-1.10584)
    assert "energy_units" not in back.info  # consumed, not carried

    p2 = tmp_path / "rt2.xyz"
    write_extxyz(p2, back)
    # exactly one energy_units marker, no duplication
    assert p2.read_text().count("energy_units=") == 1


# ---------------------------------------------------------------------------
# per-atom canonical_id column (Properties schema)
# ---------------------------------------------------------------------------


def test_read_canonical_id_column(tmp_path):
    """A Properties schema with canonical_id:I:1 populates Molecule.ids."""
    content = (
        "3\n"
        "Properties=species:S:1:pos:R:3:canonical_id:I:1\n"
        "C 0.0 0.0 0.0 1\n"
        "H 0.0 0.0 1.0 2\n"
        "H 0.0 1.0 0.0 2\n"
    )
    path = tmp_path / "ids.extxyz"
    path.write_text(content)

    mol = read_extxyz(path)
    assert mol.ids is not None
    assert mol.ids.dtype == np.int32
    np.testing.assert_array_equal(mol.ids, [1, 2, 2])
    # positions still parsed from the correct columns
    np.testing.assert_allclose(mol.get_positions()[1], [0.0, 0.0, 1.0])


def test_read_no_properties_no_ids(tmp_path):
    """Without a Properties token, ids stays None and extra columns are ignored."""
    content = "2\nenergy=-1.0\nH 0.0 0.0 0.0 7\nH 0.0 0.0 0.8 9\n"
    path = tmp_path / "noprop.xyz"
    path.write_text(content)

    mol = read_extxyz(path)
    assert mol.ids is None
    assert mol.get_chemical_symbols() == ["H", "H"]


def test_write_canonical_id_column(tmp_path):
    """The writer emits the Properties token and the per-atom id column."""
    mol = Molecule(
        symbols=["C", "H", "H"],
        positions=np.zeros((3, 3)),
        ids=[1, 2, 2],
    )
    path = tmp_path / "out.extxyz"
    write_extxyz(path, mol)

    text = path.read_text()
    assert f"Properties=species:S:1:pos:R:3:{CANONICAL_ID_COLUMN}:I:1" in text
    # last token of each atom line is the integer id
    atom_lines = text.splitlines()[2:]
    assert [ln.split()[-1] for ln in atom_lines] == ["1", "2", "2"]


def test_write_without_ids_stays_legacy(tmp_path):
    """No ids => no Properties token (byte-compatible legacy output)."""
    mol = Molecule(symbols=["H"], positions=[[0.0, 0.0, 0.0]])
    path = tmp_path / "legacy.xyz"
    write_extxyz(path, mol)

    text = path.read_text()
    assert "Properties=" not in text
    assert len(text.splitlines()[2].split()) == 4  # sym x y z only


def test_canonical_id_roundtrip(tmp_path):
    """write -> read preserves the per-atom ids."""
    mol = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
        ids=[5, 3, 3],
    )
    path = tmp_path / "rt.extxyz"
    write_extxyz(path, mol)

    back = read_extxyz(path)
    np.testing.assert_array_equal(back.ids, [5, 3, 3])
    assert "Properties" not in back.info  # consumed, not leaked


def test_read_canonical_id_as_float(tmp_path):
    """Integer ids written as floats (e.g. '2.0') are tolerated on read."""
    content = (
        "2\n"
        "Properties=species:S:1:pos:R:3:canonical_id:I:1\n"
        "H 0.0 0.0 0.0 1.0\n"
        "H 0.0 0.0 0.8 2.0\n"
    )
    path = tmp_path / "floatid.extxyz"
    path.write_text(content)

    mol = read_extxyz(path)
    np.testing.assert_array_equal(mol.ids, [1, 2])
