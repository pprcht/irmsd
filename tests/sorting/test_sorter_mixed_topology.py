"""Regression tests for the iRMSD sorter on ensembles with more than one
rank pattern, and for its OpenMP parallelization.

The sorter reuses one rmsd_cache for all pairs. min_rmsd used to count the
atoms per rank only while that count was still empty, so every later pair
reused the counts of the first pair. With an ensemble mixing isopropanol
(two equivalent methyl groups) and 1-propanol (same atoms, different classes),
the 1-propanol near-duplicates were then never merged.
"""

import json
import os
import subprocess
import sys

import numpy as np
import pytest

from irmsd.api.sorter_exposed import sorter_irmsd

ISOPROPANOL_Z = [6, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1]
ISOPROPANOL_XYZ = [
    [1.3863, 0.0106, 0.1074],
    [-0.0415, 0.3720, -0.2693],
    [-0.9934, -0.7911, -0.0353],
    [-0.4466, 1.4803, 0.5278],
    [1.7543, -0.8306, -0.4877],
    [1.4588, -0.2511, 1.1689],
    [2.0509, 0.8670, -0.0504],
    [-0.0762, 0.6795, -1.3200],
    [-1.0136, -1.0751, 1.0227],
    [-2.0169, -0.5082, -0.3038],
    [-0.7074, -1.6648, -0.6287],
    [-1.3548, 1.7117, 0.2684],
]
PROPANOL_Z = [6, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1]
PROPANOL_XYZ = [
    [-1.3115, -0.2035, 0.4439],
    [-0.1645, 0.6789, -0.0273],
    [0.9366, -0.1202, -0.7131],
    [1.5247, -1.0487, 0.1872],
    [-1.7576, -0.7466, -0.3953],
    [-0.9733, -0.9337, 1.1857],
    [-2.0929, 0.4083, 0.9059],
    [-0.5506, 1.4340, -0.7207],
    [0.2480, 1.2047, 0.8416],
    [1.7228, 0.5476, -1.0792],
    [0.5406, -0.6772, -1.5680],
    [1.8777, -0.5436, 0.9394],
]


def copies(Z, xyz, n, rng, sigma=0.01):
    """n near-duplicates: noise of sigma (A), random rotation, shuffled atoms"""
    Z = np.asarray(Z, dtype=np.int32)
    P = np.asarray(xyz, dtype=float)
    out = []
    for _ in range(n):
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        perm = rng.permutation(len(Z))
        Q = (P + rng.normal(scale=sigma, size=P.shape)) @ q.T
        out.append((Z[perm], Q[perm]))
    return out


def sort_groups(ensemble):
    groups, _, _ = sorter_irmsd(
        [z for z, _ in ensemble],
        [p for _, p in ensemble],
        len(ensemble[0][0]),
        rthr=0.125,
        iinversion=0,
        allcanon=True,
        printlvl=0,
    )
    return [int(g) for g in groups]


@pytest.mark.parametrize("order", ["isopropanol_first", "propanol_first", "interleaved"])
def test_mixed_rank_patterns_merge_duplicates(order):
    rng = np.random.default_rng(0)
    iso = [("iso", s) for s in copies(ISOPROPANOL_Z, ISOPROPANOL_XYZ, 2, rng)]
    prop = [("prop", s) for s in copies(PROPANOL_Z, PROPANOL_XYZ, 4, rng)]
    if order == "isopropanol_first":
        ens = iso + prop
    elif order == "propanol_first":
        ens = prop + iso
    else:
        ens = [prop[0], iso[0], prop[1], prop[2], iso[1], prop[3]]

    groups = sort_groups([s for _, s in ens])
    by_name = {}
    for (name, _), g in zip(ens, groups):
        by_name.setdefault(name, set()).add(g)

    assert len(by_name["iso"]) == 1, groups
    assert len(by_name["prop"]) == 1, groups
    assert by_name["iso"] != by_name["prop"], groups


THREAD_SCRIPT = """
import json, sys
import numpy as np
sys.path.insert(0, {testdir!r})
from test_sorter_mixed_topology import (ISOPROPANOL_XYZ, ISOPROPANOL_Z, PROPANOL_XYZ,
                                        PROPANOL_Z, copies, sort_groups)
rng = np.random.default_rng(7)
ens = []
for k in range(12):
    sigma = 0.01 if k % 3 else 0.3
    ens += copies(PROPANOL_Z, PROPANOL_XYZ, 2, rng, sigma)
    ens += copies(ISOPROPANOL_Z, ISOPROPANOL_XYZ, 1, rng, sigma)
print(json.dumps(sort_groups(ens)))
"""


def test_thread_count_does_not_change_grouping():
    """Each pair is compared independently, so the grouping must not depend on
    the number of OpenMP threads (a build without OpenMP passes trivially)."""
    script = THREAD_SCRIPT.format(testdir=os.path.dirname(os.path.abspath(__file__)))
    results = []
    for threads in ("1", "3"):
        env = dict(os.environ, OMP_NUM_THREADS=threads)
        out = subprocess.run(
            [sys.executable, "-c", script], env=env, capture_output=True, text=True, check=True
        )
        results.append(json.loads(out.stdout.strip().splitlines()[-1]))
    assert results[0] == results[1]
    assert len(set(results[0])) > 1  # the ensemble really has several groups
