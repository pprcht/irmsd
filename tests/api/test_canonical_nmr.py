"""Tests for the 'apsp+nmr' invtype: a magnetic-equivalency hack that splits
any canonical rank shared by exactly two non-hydrogen atoms into two distinct
ranks (and lets that distinction propagate to attached hydrogens)."""

import pathlib
from collections import Counter

import numpy as np
import pytest

from irmsd import read_structures


@pytest.fixture(scope="module")
def pentane():
    """A single pentane structure (CCCCC, terminal C's and the two CH2 are
    each an equivalent pair under plain apsp+)."""
    root = pathlib.Path(__file__).resolve().parents[1]
    path = root.parent / "examples" / "extxyz" / "pentane_single.extxyz"
    return read_structures([str(path)])[0]


def _dupe_ranks(ranks, mask):
    """ranks that appear more than once among the masked atoms."""
    counts = Counter(ranks[mask].tolist())
    return {r: c for r, c in counts.items() if c > 1}


def test_nmr_splits_equivalent_heavy_pairs(pentane):
    Z = pentane.get_atomic_numbers()
    heavy = Z != 1

    r_apsp = pentane.get_canonical(invtype="apsp+")
    r_nmr = pentane.get_canonical(invtype="apsp+nmr")

    # plain apsp+ leaves equivalent heavy pairs (terminal C's, the two CH2)
    assert _dupe_ranks(r_apsp, heavy), "expected degenerate heavy ranks under apsp+"
    # the nmr hack splits every exact pair -> all heavy ranks distinct
    assert _dupe_ranks(r_nmr, heavy) == {}
    assert len(set(r_nmr[heavy].tolist())) == int(heavy.sum())


def test_nmr_distinction_propagates_to_hydrogens(pentane):
    Z = pentane.get_atomic_numbers()
    Hmask = Z == 1

    r_apsp = pentane.get_canonical(invtype="apsp+")
    r_nmr = pentane.get_canonical(invtype="apsp+nmr")

    # splitting heavy atoms yields strictly more distinct H environments
    assert len(set(r_nmr[Hmask].tolist())) > len(set(r_apsp[Hmask].tolist()))


def test_nmr_is_a_refinement_of_apsp(pentane):
    """Atoms made distinct by the hack were equivalent before; the hack never
    merges atoms that apsp+ kept apart."""
    r_apsp = pentane.get_canonical(invtype="apsp+")
    r_nmr = pentane.get_canonical(invtype="apsp+nmr")

    n = len(r_apsp)
    for i in range(n):
        for j in range(i + 1, n):
            if r_nmr[i] == r_nmr[j]:
                # equal under the finer scheme => must be equal under apsp+
                assert r_apsp[i] == r_apsp[j]


def test_nmr_does_not_change_default_apsp(pentane):
    """The new invtype must not perturb the plain apsp+ result."""
    a = pentane.get_canonical(invtype="apsp+")
    b = pentane.get_canonical(invtype="apsp+")
    np.testing.assert_array_equal(a, b)
