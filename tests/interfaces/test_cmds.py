from io import StringIO

import numpy as np
import pytest

from irmsd.interfaces.cmds import (
    compute_axis_and_print,
    compute_canonical_and_print,
    compute_cn_and_print,
    compute_irmsd_and_print,
    compute_quaternion_rmsd_and_print,
)
from irmsd.utils.xyz import read_extxyz


@pytest.fixture(scope="module")
def axis_test_data(caffeine_axis_test_data_all):
    """Caffeine test data for axis computation tests of cmds."""
    CAFFEINE_AXIS_TEST_DATA = caffeine_axis_test_data_all
    atoms_list = []
    expected_rot = []
    expected_avmom = []
    expected_evecs = []
    for xyz_string, rot, avmom, evecs in CAFFEINE_AXIS_TEST_DATA:
        xyz_file = StringIO(xyz_string)
        atoms_list.append(read_extxyz(xyz_file))
        expected_rot.append(rot)
        expected_avmom.append(avmom)
        expected_evecs.append(evecs)

    return atoms_list, expected_rot, expected_avmom, expected_evecs


def test_compute_axis_and_print(axis_test_data):
    atoms_list, expected_rot, expected_avmom, expected_evecs = axis_test_data
    results = compute_axis_and_print(atoms_list)
    for i, (rot, avmom, evecs) in enumerate(results):
        assert pytest.approx(rot, abs=1e-6) == expected_rot[i]
        assert pytest.approx(avmom, rel=1e-4) == expected_avmom[i]
        assert pytest.approx(evecs, abs=1e-6) == expected_evecs[i]


@pytest.fixture(scope="module")
def cn_test_data(caffeine_cn_test_data_all):
    """Caffeine test data for CN computation tests of cmds."""
    CAFFEINE_CN_TEST_DATA = caffeine_cn_test_data_all
    atoms_list = []
    expected_cn = []
    for xyz_string, cn in CAFFEINE_CN_TEST_DATA:
        xyz_file = StringIO(xyz_string)
        atoms_list.append(read_extxyz(xyz_file))
        expected_cn.append(cn)

    return atoms_list, expected_cn


def test_compute_cn_and_print(cn_test_data):
    atoms_list, expected_cn = cn_test_data
    results = compute_cn_and_print(atoms_list)
    for i, cn in enumerate(results):
        assert pytest.approx(cn, abs=1e-6) == expected_cn[i]


@pytest.fixture(scope="module")
def canonical_test_data(caffeine_canonical_test_data_all):
    """Caffeine test data for canonical rank computation tests of cmds."""
    CAFFEINE_CANONICAL_TEST_DATA = caffeine_canonical_test_data_all
    atoms_list = [[], []]
    expected_ranks = [[], []]
    heavies = [False, True]
    for xyz_string, heavy, rank in CAFFEINE_CANONICAL_TEST_DATA:
        xyz_file = StringIO(xyz_string)
        k = 1 if heavy else 0
        atoms_list[k].append(read_extxyz(xyz_file))
        expected_ranks[k].append(rank)

    return atoms_list, heavies, expected_ranks


def test_compute_canonical_and_print(canonical_test_data):
    atoms_list, heavies, expected_ranks = canonical_test_data
    for k, heavy in enumerate(heavies):
        results = compute_canonical_and_print(atoms_list[k], heavy=heavy)
        for i, rank in enumerate(results):
            if heavies[i] == heavy:
                assert pytest.approx(rank, abs=1e-6) == expected_ranks[k][i]
