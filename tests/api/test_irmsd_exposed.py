import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.irmsd_exposed import get_irmsd


def test_get_irmsd(caffeine_irmsd_test_data):
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = (
        caffeine_irmsd_test_data
    )
    atom_numbers1, positions1 = get_atom_num_and_pos_from_xyz(conformer1)
    atom_numbers2, positions2 = get_atom_num_and_pos_from_xyz(conformer2)
    expected_atom_numbers, expected_positions = get_atom_num_and_pos_from_xyz(
        expected_aligned_conformer
    )

    (
        rmsd,
        atom_numbers_aligned1,
        positions_aligned1,
        atom_numbers_aligned2,
        positions_aligned2,
    ) = get_irmsd(
        atom_numbers1, positions1, atom_numbers2, positions2, iinversion=1
    )

    assert pytest.approx(expected_irmsd, abs=1e-6) == rmsd


import numpy as np

from irmsd.api.canonical_exposed import get_canonical_fortran


def test_get_irmsd_with_external_ranks_matches_internal(caffeine_irmsd_test_data):
    """Passing valid canonical ranks reproduces the internally-computed iRMSD."""
    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    # baseline: let the backend compute the ranks itself
    rmsd_internal, *_ = get_irmsd(Z1, P1, Z2, P2, iinversion=1)

    # compute the same canonical ranks up front and feed them back in
    rank1 = get_canonical_fortran(Z1, P1)
    rank2 = get_canonical_fortran(Z2, P2)
    rmsd_external, *_ = get_irmsd(
        Z1, P1, Z2, P2, iinversion=1, ranks1=rank1, ranks2=rank2
    )

    assert pytest.approx(rmsd_internal, abs=1e-9) == rmsd_external


def test_get_irmsd_zero_ranks_is_sentinel(caffeine_irmsd_test_data):
    """All-zero ranks (the sentinel) behave exactly like passing nothing."""
    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    rmsd_default, *_ = get_irmsd(Z1, P1, Z2, P2, iinversion=1)
    rmsd_zeros, *_ = get_irmsd(
        Z1,
        P1,
        Z2,
        P2,
        iinversion=1,
        ranks1=np.zeros(len(Z1), dtype=np.int32),
        ranks2=np.zeros(len(Z2), dtype=np.int32),
    )

    assert pytest.approx(rmsd_default, abs=1e-12) == rmsd_zeros


def test_get_irmsd_inconsistent_ranks_fall_back(caffeine_irmsd_test_data):
    """Ranks rejected by the consistency check fall back to recomputation."""
    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    rmsd_default, *_ = get_irmsd(Z1, P1, Z2, P2, iinversion=1)

    # all-distinct vs all-equal -> different max rank -> checkranks() rejects
    bad1 = np.arange(1, len(Z1) + 1, dtype=np.int32)
    bad2 = np.ones(len(Z2), dtype=np.int32)
    rmsd_fallback, *_ = get_irmsd(
        Z1, P1, Z2, P2, iinversion=1, ranks1=bad1, ranks2=bad2
    )

    assert pytest.approx(rmsd_default, abs=1e-9) == rmsd_fallback


def test_get_irmsd_molecule_use_ids(caffeine_irmsd_test_data):
    """Molecule.ids are wired through and reproduce the no-ids result."""
    from irmsd.core.molecule import Molecule
    from irmsd.interfaces.mol_interface import get_irmsd_molecule

    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    syms1 = [int(z) for z in Z1]
    syms2 = [int(z) for z in Z2]
    from irmsd.core.molecule import _INV_PERIODIC

    mol1 = Molecule(symbols=[_INV_PERIODIC[z] for z in syms1], positions=P1)
    mol2 = Molecule(symbols=[_INV_PERIODIC[z] for z in syms2], positions=P2)

    rmsd_no_ids, *_ = get_irmsd_molecule(mol1, mol2, iinversion=1, use_ids=False)

    mol1.set_ids(get_canonical_fortran(Z1, P1))
    mol2.set_ids(get_canonical_fortran(Z2, P2))
    rmsd_ids, *_ = get_irmsd_molecule(mol1, mol2, iinversion=1, use_ids=True)

    assert pytest.approx(rmsd_no_ids, abs=1e-9) == rmsd_ids


def test_external_ranks_match_internal_auto_inversion(caffeine_irmsd_test_data):
    """iinversion=0: self-filling stereo detection on provided ranks must
    reproduce the fully-recomputed result (exercises has_stereo self-fill)."""
    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    rmsd_internal, *_ = get_irmsd(Z1, P1, Z2, P2, iinversion=0)

    rank1 = get_canonical_fortran(Z1, P1)
    rank2 = get_canonical_fortran(Z2, P2)
    rmsd_external, *_ = get_irmsd(
        Z1, P1, Z2, P2, iinversion=0, ranks1=rank1, ranks2=rank2
    )

    assert pytest.approx(rmsd_internal, abs=1e-9) == rmsd_external


def test_external_ranks_one_side_only(caffeine_irmsd_test_data):
    """Providing ranks for only one structure recomputes the other side and
    still reconciles to the internal result (per-structure granularity)."""
    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    rmsd_internal, *_ = get_irmsd(Z1, P1, Z2, P2, iinversion=0)

    rank1 = get_canonical_fortran(Z1, P1)  # only structure 1 provided
    rmsd_mixed, *_ = get_irmsd(
        Z1, P1, Z2, P2, iinversion=0, ranks1=rank1, ranks2=None
    )

    assert pytest.approx(rmsd_internal, abs=1e-9) == rmsd_mixed
