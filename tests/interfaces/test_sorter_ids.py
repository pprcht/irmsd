"""End-to-end tests for wiring per-atom canonical ids through the ensemble
(sorter / delta) path. The provided ids should reproduce the result obtained
when the backend recomputes the canonical ranks itself."""

import numpy as np
import pytest

from irmsd import read_structures
from irmsd.interfaces.mol_interface import (
    delta_irmsd_list_molecule,
    sorter_irmsd_molecule,
)


@pytest.fixture(scope="module")
def pentane_subset():
    """A small slice of the shuffled pentane ensemble (all one conformer)."""
    root = __import__("pathlib").Path(__file__).resolve().parents[1]
    path = root / "data" / "ensembles_shuffled" / "pentane_100_shuffled.xyz"
    mols = read_structures([str(path)])
    return mols[:12]


def _with_ids(mols):
    """Return copies carrying their computed canonical ids."""
    out = []
    for m in mols:
        c = m.copy()
        c.set_ids(m.get_canonical())
        out.append(c)
    return out


def test_sorter_ids_match_recomputed(pentane_subset):
    """Attaching canonical ids to every structure reproduces the id-less sort."""
    base_groups, _ = sorter_irmsd_molecule(
        pentane_subset, rthr=0.125, iinversion=1, use_ids=False
    )
    id_groups, _ = sorter_irmsd_molecule(
        _with_ids(pentane_subset), rthr=0.125, iinversion=1, use_ids=True
    )
    np.testing.assert_array_equal(base_groups, id_groups)
    # the shuffled pentane ensemble is a single conformer
    assert len(np.unique(id_groups)) == 1


def test_sorter_ids_mixed_some_missing(pentane_subset):
    """A mix of structures with and without ids still matches the baseline."""
    base_groups, _ = sorter_irmsd_molecule(
        pentane_subset, rthr=0.125, iinversion=1, use_ids=False
    )

    mixed = []
    for i, m in enumerate(pentane_subset):
        c = m.copy()
        if i % 2 == 0:  # only every other structure carries ids
            c.set_ids(m.get_canonical())
        mixed.append(c)

    mixed_groups, _ = sorter_irmsd_molecule(
        mixed, rthr=0.125, iinversion=1, use_ids=True
    )
    np.testing.assert_array_equal(base_groups, mixed_groups)


def test_sorter_inconsistent_ids_fall_back(pentane_subset):
    """Garbage ids that fail the consistency check fall back gracefully."""
    base_groups, _ = sorter_irmsd_molecule(
        pentane_subset, rthr=0.125, iinversion=1, use_ids=False
    )

    nat = pentane_subset[0].natoms
    bad = []
    for i, m in enumerate(pentane_subset):
        c = m.copy()
        # structure 0 keeps a valid (all-distinct) ranking, others all-equal:
        # checkranks rejects the pair -> fallback to atom types
        if i == 0:
            c.set_ids(np.arange(1, nat + 1, dtype=np.int32))
        else:
            c.set_ids(np.ones(nat, dtype=np.int32))
        bad.append(c)

    fb_groups, _ = sorter_irmsd_molecule(
        bad, rthr=0.125, iinversion=1, use_ids=True
    )
    # falling back to atom-type ranks still recognizes the single conformer
    np.testing.assert_array_equal(base_groups, fb_groups)


def test_delta_ids_match_recomputed(pentane_subset):
    """Provided ids reproduce the id-less delta-iRMSD list."""
    base_delta, _ = delta_irmsd_list_molecule(
        pentane_subset, iinversion=1, use_ids=False
    )
    id_delta, _ = delta_irmsd_list_molecule(
        _with_ids(pentane_subset), iinversion=1, use_ids=True
    )
    np.testing.assert_allclose(base_delta, id_delta, atol=1e-9)
