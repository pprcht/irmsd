from io import StringIO

import numpy as np
import pytest

pytest.importorskip("ase")
import ase
from ase import Atoms
from ase.io import read as ase_read

from irmsd import Molecule
from irmsd.interfaces.ase_io import (
    ase_to_molecule,
    delta_irmsd_list_ase,
    get_axis_ase,
    get_canonical_ase,
    get_cn_ase,
    get_irmsd_ase,
    get_rmsd_ase,
    sorter_irmsd_ase,
)


def test_get_axis_ase(caffeine_axis_test_data):
    caffeine_xzy, expected_rot, expected_avmom, expected_evec = caffeine_axis_test_data
    caffeine_file = StringIO(caffeine_xzy)

    atoms = ase_read(caffeine_file, format="xyz")

    rot, avmom, evec = get_axis_ase(atoms)

    assert pytest.approx(rot, abs=1e-6) == expected_rot
    assert pytest.approx(avmom, rel=1e-4) == expected_avmom
    assert pytest.approx(evec, abs=1e-6) == expected_evec


def test_get_cn_ase(caffeine_cn_test_data):
    caffeine_xzy, expected_cn = caffeine_cn_test_data
    caffeine_file = StringIO(caffeine_xzy)

    atoms = ase_read(caffeine_file, format="xyz")

    cn = get_cn_ase(atoms)

    assert pytest.approx(cn, abs=1e-6) == expected_cn


def test_get_canonical_ase(caffeine_canonical_test_data):
    caffeine_xzy, heavy, expected_rank = caffeine_canonical_test_data
    caffeine_file = StringIO(caffeine_xzy)

    atoms = ase_read(caffeine_file, format="xyz")

    rank = get_canonical_ase(atoms, heavy=heavy)

    assert pytest.approx(rank, abs=1e-6) == expected_rank


def test_get_rmsd_ase(caffeine_rmsd_test_data):
    (
        caffeine_xzy_1,
        caffeine_xzy_2,
        heavy,
        expected_rmsd,
        expected_aligned,
        expected_Umat,
    ) = caffeine_rmsd_test_data
    caffeine_file_1 = StringIO(caffeine_xzy_1)
    caffeine_file_2 = StringIO(caffeine_xzy_2)
    expected_aligned_file = StringIO(expected_aligned)

    atoms1 = ase_read(caffeine_file_1, format="xyz")
    atoms2 = ase_read(caffeine_file_2, format="xyz")
    expected_aligned = ase_read(expected_aligned_file, format="xyz")

    mask = atoms1.get_atomic_numbers() != 1 if heavy else None  # exclude hydrogens
    rmsd, atoms_aligned, Umat = get_rmsd_ase(atoms1, atoms2, mask=mask)

    assert pytest.approx(rmsd, abs=1e-6) == expected_rmsd
    assert (
        pytest.approx(atoms_aligned.get_positions(), abs=1e-6)
        == expected_aligned.get_positions()
    )
    assert pytest.approx(Umat, abs=1e-6) == expected_Umat


def test_get_irmsd_ase(caffeine_irmsd_test_data):
    caffeine_xzy_1, caffeine_xzy_2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )
    caffeine_file_1 = StringIO(caffeine_xzy_1)
    caffeine_file_2 = StringIO(caffeine_xzy_2)
    # expected_aligned_file = StringIO(expected_aligned)

    atoms1 = ase_read(caffeine_file_1, format="xyz")
    atoms2 = ase_read(caffeine_file_2, format="xyz")
    # expected_aligned = ase_read(expected_aligned_file, format="xyz")

    irmsd, atoms_aligned1, atoms_aligned2 = get_irmsd_ase(atoms1, atoms2)

    assert pytest.approx(irmsd, abs=1e-6) == expected_irmsd


def test_ase_to_molecule_single_basic():
    """Single ASE Atoms → Molecule: structure, cell, pbc, energy, info."""
    # Build a simple water molecule
    symbols = ["O", "H", "H"]
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
    ]
    cell = np.diag([10.0, 10.0, 10.0])

    a = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=[True, False, True])

    # Attach metadata
    a.info["energy"] = -76.1234
    a.info["tag"] = "test_water"

    mol = ase_to_molecule(a)
    assert isinstance(mol, Molecule)

    # Structure
    assert mol.natoms == 3
    assert mol.get_chemical_symbols() == symbols
    np.testing.assert_allclose(mol.get_positions(), np.array(positions, float))

    # Cell and PBC
    assert mol.cell is not None
    np.testing.assert_allclose(mol.cell, cell)
    assert mol.pbc == (True, False, True)

    # Energy
    assert mol.energy == pytest.approx(-76.1234)

    # Info (at least the extra key survives)
    assert mol.info.get("tag") == "test_water"


def test_ase_to_molecule_single_no_energy():
    """If ASE Atoms has no energy info or calculator, Molecule.energy should be
    None."""
    a = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])

    # No info['energy'], no calculator
    assert "energy" not in a.info
    assert a.calc is None

    mol = ase_to_molecule(a)
    assert isinstance(mol, Molecule)
    assert mol.natoms == 2
    assert mol.energy is None  # nothing to pick up


def test_ase_to_molecule_sequence():
    """Sequence of ASE Atoms → list[Molecule]."""
    a1 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]])
    a2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.8]])

    a1.info["energy"] = -1.0
    a2.info["energy"] = -2.0

    mols = ase_to_molecule([a1, a2])

    assert isinstance(mols, list)
    assert len(mols) == 2
    assert all(isinstance(m, Molecule) for m in mols)

    # Check ordering and basic fields
    assert mols[0].natoms == 2
    assert mols[1].natoms == 2

    assert mols[0].energy == pytest.approx(-1.0)
    assert mols[1].energy == pytest.approx(-2.0)

    np.testing.assert_allclose(
        mols[0].get_positions(),
        np.array([[0, 0, 0], [0, 0, 0.7]], float),
    )
    np.testing.assert_allclose(
        mols[1].get_positions(),
        np.array([[0, 0, 0], [0, 0, 0.8]], float),
    )


def test_sorter_irmsd_ase(caffeine_sorter_irmsd_test_data):
    conformer_list, rthr, expected_groups = caffeine_sorter_irmsd_test_data
    atoms_list = []
    for conf_xyz in conformer_list:
        conf_file = StringIO(conf_xyz)
        atoms = ase_read(conf_file, format="xyz")
        atoms_list.append(atoms)

    groups, new_atoms_list = sorter_irmsd_ase(atoms_list, rthr)
    assert all(groups == expected_groups)


def test_delta_irmsd_list_ase(caffeine_delta_irmsd_list_test_data):
    conformer_list, expected_delta = caffeine_delta_irmsd_list_test_data
    atoms_list = []
    for conf_xyz in conformer_list:
        conf_file = StringIO(conf_xyz)
        atoms = ase_read(conf_file, format="xyz")
        atoms_list.append(atoms)

    delta, new_atoms = delta_irmsd_list_ase(atoms_list)
    assert pytest.approx(delta, abs=1e-6) == expected_delta
