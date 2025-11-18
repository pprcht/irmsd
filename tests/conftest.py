from random import shuffle

import data
import numpy as np
import pytest


@pytest.fixture(scope="session")
def caffeine_molecule_xyz():
    """Fixture that returns a caffeine molecule."""
    return data.CAFFEINE_XTB


@pytest.fixture(scope="session")
def caffeine_molecule_xyz_rotated():
    """Fixture that returns a caffeine molecule."""
    return data.CAFFEINE_XTB_ROTATED


@pytest.fixture(scope="session")
def caffeine_molecule_xyz_obabel():
    """Fixture that returns a caffeine molecule."""
    return data.CAFFEINE_OBABEL


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, Conformer2, heavy, rmsd, aligned_conformer, Umat
        (
            data.CAFFEINE_XTB,
            data.CAFFEINE_XTB_ROTATED,
            False,
            0.00000042,
            data.CAFFEINE_XTB,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    [0.000586, 0.025178, 0.999683],
                    [0.025178, 0.999366, -0.025185],
                    [-0.999683, 0.025185, -0.000049],
                ],
                dtype=np.float64,
            ),
        ),
        (
            data.CAFFEINE_XTB,
            data.CAFFEINE_OBABEL,
            False,
            0.13939435,
            data.CAFFEINE_OBABEL_ALIGNED,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    [0.165823, -0.986155, -0.000907],
                    [0.986155, 0.165823, 0.000503],
                    [-0.000346, -0.000978, 0.999999],
                ],
                dtype=np.float64,
            ),
        ),
        (
            data.CAFFEINE_XTB,
            data.CAFFEINE_XTB_ROTATED,
            True,
            0.00000042,
            data.CAFFEINE_XTB,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    [0.000586, 0.025178, 0.999683],
                    [0.025178, 0.999366, -0.025185],
                    [-0.999683, 0.025185, -0.000049],
                ],
                dtype=np.float64,
            ),
        ),
        (
            data.CAFFEINE_XTB,
            data.CAFFEINE_OBABEL,
            True,
            0.1725815540,
            data.CAFFEINE_OBABEL_HEAVY_ALIGNED,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    [0.171404, -0.9852, -0.001332],
                    [0.985201, 0.171405, -0.00012],
                    [0.000347, -0.001292, 0.999999],
                ],
                dtype=np.float64,
            ),
        ),
    ],
)
def caffeine_rmsd_test_data(request):
    """Fixture that returns atom numbers and positions for caffeine molecule in
    different conformations."""
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, rot constants in MHz, avmom in a.u., evec
        (
            data.CAFFEINE_XTB,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                [
                    1068.073127,
                    710.711798,
                    430.26121,
                ],
                dtype=np.float64,
            ),
            1.30564546e-44,
            np.asarray(
                [
                    [0.671158, -0.741313, 0.001206],
                    [0.741314, 0.671157, -0.001149],
                    [0.000042, 0.001665, 0.999999],
                ],
                dtype=np.float64,
            ),
        ),
    ],
)
def caffeine_axis_test_data(request):
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, cn
        (
            data.CAFFEINE_OBABEL,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                # fmt: off
            [3.689705,3.866219,2.934731,3.626334,3.144874,3.791101,2.768159,0.858364,2.743387,
             2.716063,0.860503,2.747386,3.698299,3.698564,0.924137,0.924083,0.924084,0.925654,
             0.924202,0.92413,0.924141,0.924245,0.924105,0.924098],
                # fmt: on
                dtype=np.float64,
            ),
        ),
    ],
)
def caffeine_cn_test_data(request):
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, heavy, canonical_atom_id
        (
            data.CAFFEINE_XTB,
            False,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                # fmt: off
                [1,8,5,7,13,14,12,6,10,9,2,11,4,3,15,15,15,18,17,17,17,16,16,16],
                # fmt: on
                dtype=int,
            ),
        ),
        (
            data.CAFFEINE_XTB,
            True,
            # FIXME: check with crest, this is actually taken from this code itself but the rmsd is the same as crest
            np.asarray(
                # fmt: off
                [1,8,5,7,13,14,12,6,10,9,2,11,4,3,0,0,0,0,0,0,0,0,0,0],
                # fmt: on
                dtype=int,
            ),
        ),
    ],
)
def caffeine_canonical_test_data(request):
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        # Conformer1, Conformer2, mask, rmsd, aligned_conformer, Umat
        (
            data.CAFFEINE_XTB,
            data.CAFFEINE_XTB,
            0.00000000,
            data.CAFFEINE_XTB,
        ),
        (
            data.CAFFEINE_XTB,
            data.CAFFEINE_XTB_ROTATED,
            0.00000042,
            data.CAFFEINE_XTB,
        ),
        (data.CAFFEINE_XTB, data.CAFFEINE_GFF, 0.06160473, data.CAFFEINE_GFF_ALIGNED),
    ],
)
def caffeine_irmsd_test_data(request):
    conformer1, conformer2, expected_irmsd, expected_aligned_conformer = request.param
    conformer2_list = conformer2.splitlines()
    conformer2_header = conformer2_list[:2]
    conformer2_molecule = conformer2_list[2:]
    shuffle(conformer2_molecule)

    # conformer2_shuffled = conformer2_header + conformer2_molecule[::-1]
    conformer2_shuffled = conformer2_header + conformer2_molecule

    conformer2 = "\n".join(conformer2_shuffled)
    return conformer1, conformer2, expected_irmsd, expected_aligned_conformer
