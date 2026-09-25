"""Regression test for the equal-mass axis fallback (see this folder's
README.md and the main README's "Equal-Mass Axis Fallback" section).

Fe4Cage.xyz is a T-symmetric Fe4 cage consisting of 4 heavy metal centers and organic linkers connected in a tetrahedron framework.
The file pointgroup_isomers.json holds 200 atom permutations, each a point-group operation (simple rotation) of this structure.
Each of these isomers must score 0.00 A against the original, unpermuted structure.
Before the equal-mass fallback, get_irmsd scored only 170/200 of these correctly (worst case ~2 A).
With this new change, all 200/200 are correctly identified as identical (worst case ~2e-7 A).
"""

import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ase")
from ase.io import read as ase_read

from irmsd.api.irmsd_exposed import get_irmsd

HERE = Path(__file__).parent


@pytest.fixture(scope="module")
def fe4cage():
    atoms = ase_read(HERE / "Fe4Cage.xyz")
    return atoms.get_atomic_numbers(), atoms.get_positions()

@pytest.fixture(scope="module")
def pointgroup_isomers():
    with open(HERE / "pointgroup_isomers.json") as f:
        return json.load(f)


def test_all_pointgroup_isomers_score_near_zero(fe4cage, pointgroup_isomers):
    Z, pos = fe4cage
    worst = 0.0
    n_near_zero = 0
    for rec in pointgroup_isomers:
        pos_iso = pos[rec["perm"]]
        rmsdval, *_ = get_irmsd(Z, pos, Z, pos_iso)
        worst = max(worst, float(rmsdval))
        n_near_zero += rmsdval < 0.01

    assert n_near_zero == len(pointgroup_isomers), (
        f"only {n_near_zero}/{len(pointgroup_isomers)} known-exact point-group "
        f"operations scored near-zero (worst case {worst:.3f} A) -- the "
        f"equal-mass axis fallback may be missing or broken"
    )
    assert worst < 1e-4
