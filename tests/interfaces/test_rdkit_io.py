import numpy as np
import pytest

from irmsd.core.molecule import Molecule

pytest.importorskip("rdkit")

from rdkit import Chem
from rdkit.Chem import AllChem

from irmsd.interfaces.rdkit_io import (
    conf_id_to_iterator,
    conformer_iterator,
    delta_irmsd_list_rdkit,
    get_axis_rdkit,
    get_canonical_rdkit,
    get_cn_rdkit,
    get_irmsd_rdkit,
    get_rmsd_rdkit,
    rdkit_to_molecule,
    sorter_irmsd_rdkit,
    cregen_rdkit,
    prune_rdkit,
)

CAFFEINE_WBO = np.asarray(
    [
        # fmt: off
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0],
        [1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
        [0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        # fmt: on
    ],
    dtype=np.float64,
)

# ----------------------------------------------------------------------
# Fixture: build a simple RDKit Mol with multiple conformers
# ----------------------------------------------------------------------


@pytest.fixture
def mol():
    """Return an RDKit molecule with 3 conformers."""

    m = Chem.AddHs(Chem.MolFromSmiles("CC"))  # simple ethane
    AllChem.EmbedMultipleConfs(m, numConfs=3)
    assert m.GetNumConformers() == 3
    return m


# ======================================================================
# Tests for conformer_iterator
# ======================================================================


def test_conformer_iterator_valid(mol):
    """Iterate over multiple conformer IDs in order."""
    ids = [2, 0, 1]  # deliberately non-sorted
    it = conformer_iterator(mol, ids)

    returned = [conf.GetId() for conf in it]
    assert returned == ids


def test_conformer_iterator_invalid_id(mol):
    """RDKit should raise if conformer ID is invalid."""
    ids = [999]  # does not exist

    it = conformer_iterator(mol, ids)

    # Accessing the conformer should raise:
    with pytest.raises(ValueError):
        list(it)  # triggers GetConformer(999)


# ======================================================================
# Tests for conf_id_to_iterator
# ======================================================================


def test_conf_id_none_returns_all_conformers(mol):
    it = conf_id_to_iterator(mol, None)
    returned = [c.GetId() for c in it]
    assert returned == [0, 1, 2]


def test_conf_id_int_single_conformer(mol):
    it = conf_id_to_iterator(mol, 1)
    returned = [c.GetId() for c in it]
    assert returned == [1]


def test_conf_id_list_multiple(mol):
    it = conf_id_to_iterator(mol, [2, 0])
    returned = [c.GetId() for c in it]
    assert returned == [2, 0]


def test_conf_id_wrong_type_string(mol):
    with pytest.raises(TypeError, match="conf_id must be None, int, or list of int"):
        conf_id_to_iterator(mol, "abc")


def test_conf_id_wrong_type_float(mol):
    with pytest.raises(TypeError, match="conf_id must be None, int, or list of int"):
        conf_id_to_iterator(mol, 1.23)


def test_conf_id_wrong_type_tuple(mol):
    with pytest.raises(TypeError, match="conf_id must be None, int, or list of int"):
        conf_id_to_iterator(mol, (0, 1))


def test_conf_id_list_empty(mol):
    # Empty list yields an empty generator (allowed behavior)
    it = conf_id_to_iterator(mol, [])
    assert list(it) == []


# =============================================
# Tests for properties calculation functions
# =============================================


def test_get_cn_rdkit(caffeine_cn_test_data):
    caffeine_xyz, expected_cn = caffeine_cn_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    cn = get_cn_rdkit(caffeine_rdkit)
    assert pytest.approx(expected_cn, abs=1e-6) == cn


def test_get_cn_rdkit_multiple(caffeine_cn_test_data):
    caffeine_xyz, expected_cn = caffeine_cn_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)
    # Add multiple conformers (identical for this test)
    caffeine_rdkit.AddConformer(caffeine_rdkit.GetConformer(), assignId=True)
    caffeine_rdkit.AddConformer(caffeine_rdkit.GetConformer(), assignId=True)

    cn = get_cn_rdkit(caffeine_rdkit)

    assert cn.shape == (3, expected_cn.shape[0])
    for i, conf_cn in enumerate(cn):
        assert pytest.approx(expected_cn, abs=1e-6) == conf_cn


def test_get_axis_rdkit(caffeine_axis_test_data):
    caffeine_xyz, expected_rot, expected_avmom, expected_evec = caffeine_axis_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    rot, avmom, evec = get_axis_rdkit(caffeine_rdkit)
    assert pytest.approx(expected_rot, abs=1e-6) == rot
    assert pytest.approx(expected_avmom, abs=1e-6) == avmom
    assert pytest.approx(expected_evec, abs=1e-6) == evec


def test_get_axis_rdkit_multiple(caffeine_axis_test_data):
    caffeine_xyz, expected_rot, expected_avmom, expected_evec = caffeine_axis_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)
    # Add multiple conformers (identical for this test)
    caffeine_rdkit.AddConformer(caffeine_rdkit.GetConformer(), assignId=True)
    caffeine_rdkit.AddConformer(caffeine_rdkit.GetConformer(), assignId=True)

    results = get_axis_rdkit(caffeine_rdkit)

    for i, (rot, avmom, evec) in enumerate(results):
        assert pytest.approx(expected_rot, abs=1e-6) == rot
        assert pytest.approx(expected_avmom, abs=1e-6) == avmom
        assert pytest.approx(expected_evec, abs=1e-6) == evec


def test_get_canonical_rdkit(caffeine_canonical_test_data):
    caffeine_xyz, heavy, expected_canonical = caffeine_canonical_test_data

    caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)

    canonical = get_canonical_rdkit(caffeine_rdkit, heavy=heavy)
    assert pytest.approx(expected_canonical, abs=1e-6) == canonical


@pytest.mark.parametrize(
    "conformer_xyz,wbo,expected_canonical",
    [
        (
            [
                "caffeine_molecule_xyz",
                "caffeine_molecule_xyz_obabel",
                "caffeine_molecule_xyz_gff",
            ],
            None,
            np.array(
                # fmt: off
        [[ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15, 15, 18, 17, 17, 17, 16, 16, 16],
       [ 1,  7,  5, 10, 13, 14, 12,  6, 11,  8,  2,  9,  3,  4, 15, 15, 15, 18, 16, 16, 16, 17, 17, 17],
       [ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15, 15, 18, 17, 17, 17, 16, 16, 16]],
                # fmt: on
            ),
        ),
        (
            [
                "caffeine_molecule_xyz",
                "caffeine_molecule_xyz_obabel",
                "caffeine_molecule_xyz_gff",
            ],
            CAFFEINE_WBO,
            # fmt: off
        np.asarray(
            # fmt: off
                    [[ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15,
        15, 18, 17, 17, 17, 16, 16, 16],
       [ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15,
        15, 18, 17, 17, 17, 16, 16, 16],
       [ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15,
        15, 18, 17, 17, 17, 16, 16, 16]],
            # fmt: on
            dtype=int,
        ),
        ),
        (
            [
                "caffeine_molecule_xyz",
                "caffeine_molecule_xyz_obabel",
                "caffeine_molecule_xyz_gff",
            ],
            np.stack([CAFFEINE_WBO, CAFFEINE_WBO, CAFFEINE_WBO]),
            # fmt: off
        np.asarray(
            # fmt: off
                    [[ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15,
        15, 18, 17, 17, 17, 16, 16, 16],
       [ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15,
        15, 18, 17, 17, 17, 16, 16, 16],
       [ 1,  8,  5,  7, 13, 14, 12,  6, 10,  9,  2, 11,  4,  3, 15, 15,
        15, 18, 17, 17, 17, 16, 16, 16]],
            # fmt: on
            dtype=int,
        ),
        ),
    ],
)
def test_get_canonical_rdkit_specific(conformer_xyz, wbo, expected_canonical, request):
    caffeine_rdkit = None
    for xyz_fixture in conformer_xyz:
        caffeine_xyz = request.getfixturevalue(xyz_fixture)
        if caffeine_rdkit is None:
            caffeine_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)
        else:
            tmp_rdkit = Chem.MolFromXYZBlock(caffeine_xyz)
            caffeine_rdkit.AddConformer(tmp_rdkit.GetConformer(), assignId=True)

    canonical = get_canonical_rdkit(caffeine_rdkit, wbo=wbo)
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
        == aligned_mol.GetConformer().GetPositions()
    )


def test_get_irmsd_rdkit(caffeine_irmsd_test_data):
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )

    mol1 = Chem.MolFromXYZBlock(conformer1)
    mol2 = Chem.MolFromXYZBlock(conformer2)
    expected_aligned_mol = Chem.MolFromXYZBlock(expected_aligned_conformer)

    rmsd, aligned_mol_1, aligned_mol_2 = get_irmsd_rdkit(mol1, mol2, iinversion=1)

    assert pytest.approx(expected_irmsd, abs=1e-6) == rmsd


def test_rdkit_to_molecule_single_basic():
    """Single ASE Atoms → Molecule: structure, cell, pbc, energy, info."""
    # Build a simple water molecule
    symbols = ["O", "H", "H"]
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
    ]

    xyz_block = "3\n"
    xyz_block += "Water molecule\n"
    for sym, pos in zip(symbols, positions):
        xyz_block += f"{sym} {pos[0]} {pos[1]} {pos[2]}\n"

    rdkit_mol = Chem.MolFromXYZBlock(xyz_block)
    # Attach metadata
    rdkit_mol.SetProp("tag", "test_water")
    rdkit_mol.GetConformer().SetProp("energy", "-76.1234")

    mol = rdkit_to_molecule(rdkit_mol)
    assert isinstance(mol, Molecule)

    # Structure
    assert mol.natoms == 3
    assert mol.get_chemical_symbols() == symbols
    np.testing.assert_allclose(mol.get_positions(), np.array(positions, float))

    # Cell and PBC
    assert mol.cell is None

    # Energy
    assert mol.energy == pytest.approx(-76.1234)

    # Info (at least the extra key survives)
    assert mol.info.get("tag") == "test_water"


def test_rdkit_to_molecule_single_no_energy():
    """If ASE Atoms has no energy info or calculator, Molecule.energy should be
    None."""
    symbols = ["O", "H", "H"]
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
    ]
    xyz_block = "3\n"
    xyz_block += "Water molecule\n"
    for sym, pos in zip(symbols, positions):
        xyz_block += f"{sym} {pos[0]} {pos[1]} {pos[2]}\n"

    rdkit_mol = Chem.MolFromXYZBlock(xyz_block)

    mol = rdkit_to_molecule(rdkit_mol)
    assert isinstance(mol, Molecule)
    assert mol.natoms == 3
    assert mol.energy is None  # nothing to pick up
    assert mol.info == {}


def test_rdkit_to_molecule_sequence():
    """Sequence of ASE Atoms → list[Molecule]."""
    symbols = ["H", "H"]
    positions1 = np.asarray([[0, 0, 0], [0, 0, 0.7]])
    positions2 = np.asarray([[0, 0, 0], [0, 0, 0.8]])
    energies = [-1.0, -2.0]

    rdkti_mol = Chem.RWMol()
    for i, sym in enumerate(symbols):
        atom = Chem.Atom(sym)
        rdkti_mol.AddAtom(atom)

    for energy, pos in zip(energies, [positions1, positions2]):
        conformer = Chem.Conformer(len(symbols))
        conformer.SetPositions(pos)
        conformer.SetDoubleProp("energy", energy)
        rdkti_mol.AddConformer(conformer, assignId=True)

    mols = rdkit_to_molecule(rdkti_mol)

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


def test_rdkit_sequence_to_molecule_sequence():
    symbols = ["H", "H"]
    positions1 = np.asarray([[0, 0, 0], [0, 0, 0.7]])
    positions2 = np.asarray([[0, 0, 0], [0, 0, 0.8]])
    energies = [-1.0, -2.0]

    rdkit_mols = []
    for energy, pos in zip(energies, [positions1, positions2]):
        rdkti_mol = Chem.RWMol()
        for i, sym in enumerate(symbols):
            atom = Chem.Atom(sym)
            rdkti_mol.AddAtom(atom)

        conformer = Chem.Conformer(len(symbols))
        conformer.SetPositions(pos)
        conformer.SetDoubleProp("energy", energy)
        rdkti_mol.AddConformer(conformer, assignId=True)
        rdkit_mols.append(rdkti_mol)

    mols = rdkit_to_molecule(rdkit_mols)

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


def test_sorter_irmsd_rdkit(caffeine_sorter_irmsd_test_data):
    conformer_list, rthr, expected_groups = caffeine_sorter_irmsd_test_data
    rdkit_mol_list = []
    for conf_xyz in conformer_list:
        rdkit_mol = Chem.MolFromXYZBlock(conf_xyz)
        rdkit_mol_list.append(rdkit_mol)

    groups, new_rdkit_mol_list = sorter_irmsd_rdkit(rdkit_mol_list, rthr)
    assert all(groups == expected_groups)


def test_delta_irmsd_list_rdkit(caffeine_delta_irmsd_list_test_data):
    conformer_list, expected_delta = caffeine_delta_irmsd_list_test_data
    rdkit_mol_list = []
    for conf_xyz in conformer_list:
        rdkit_mol = Chem.MolFromXYZBlock(conf_xyz)
        rdkit_mol_list.append(rdkit_mol)

    delta, new_rdkit_mol_list = delta_irmsd_list_rdkit(rdkit_mol_list)
    assert pytest.approx(delta, abs=1e-4) == expected_delta


def test_cregen_rdkit(caffeine_sorter_irmsd_test_data):
    conformer_list, rthr, expected_groups = caffeine_sorter_irmsd_test_data
    rdkit_mol_list = []
    for conf_xyz in conformer_list:
        rdkit_mol = Chem.MolFromXYZBlock(conf_xyz)
        rdkit_mol_list.append(rdkit_mol)

    new_rdkit_mol_list = cregen_rdkit(rdkit_mol_list, rthr)
    assert len(new_rdkit_mol_list) == np.max(expected_groups)


def test_prune_rdkit(caffeine_sorter_irmsd_test_data):
    conformer_list, rthr, expected_groups = caffeine_sorter_irmsd_test_data
    rdkit_mol_list = []
    for conf_xyz in conformer_list:
        rdkit_mol = Chem.MolFromXYZBlock(conf_xyz)
        rdkit_mol_list.append(rdkit_mol)

    new_rdkit_mol_list = prune_rdkit(rdkit_mol_list, rthr)
    assert len(new_rdkit_mol_list) == np.max(expected_groups)
