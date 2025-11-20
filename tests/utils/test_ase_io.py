from io import StringIO

import pytest

pytest.importorskip("ase")

import ase
from ase.io import read as ase_read

from irmsd.utils.ase_io import (
    get_axis_ase,
    get_canonical_ase,
    get_cn_ase,
    get_irmsd_ase,
    get_rmsd_ase,
)
from irmsd.utils.utils import read_structures


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
