import numpy as np
import pytest
from irmsd import read_structures
from irmsd.interfaces.mol_interface import (
    delta_irmsd_list_molecule,
    get_irmsd_molecule,
    get_rmsd_molecule,
    sorter_irmsd_molecule,
)


def pytest_generate_tests(metafunc):
    if "ensembles_shuffled" in metafunc.fixturenames:
        root_dir = metafunc.config.rootpath
        folder = root_dir / "tests" / "data" / "ensembles_shuffled"
        files = list(folder.glob("*.xyz"))
        metafunc.parametrize("ensembles_shuffled", files, indirect=True)


@pytest.fixture
def ensembles_shuffled(request):
    file = request.param
    return read_structures([file])


@pytest.mark.slow
def test_ensembles(ensembles_shuffled):
    ensembles = ensembles_shuffled
    n_molecules = len(ensembles)
    rmsd_mat = np.zeros((n_molecules, n_molecules))
    for i, mol_ref in enumerate(ensembles):
        for j, mol_cmp in enumerate(ensembles[i + 1 :]):
            rmsd, mol_aligned, rotmat = get_rmsd_molecule(mol_ref, mol_cmp)

            rmsd_mat[i, j + i + 1] = rmsd
            rmsd_mat[j + i + 1, i] = rmsd

    # turn inversion on
    iinversion = 1
    irmsd_threshold = 0.125
    # group the molecules into unique conformers
    group_ids, new_molecules = sorter_irmsd_molecule(
        ensembles, irmsd_threshold, iinversion=iinversion
    )
    unique_group_ids = np.unique(group_ids)
    assert len(unique_group_ids) == 1

    # compute irmsd matrix
    irmsd_mat = np.zeros((len(new_molecules), len(new_molecules)))
    for i, mol_ref in enumerate(new_molecules):
        for j, mol_cmp in enumerate(new_molecules[i + 1 :]):
            irmsd, _, _ = get_irmsd_molecule(mol_ref, mol_cmp, iinversion=iinversion)
            irmsd_mat[i, j + i + 1] = irmsd
            irmsd_mat[j + i + 1, i] = irmsd

    # compute averages
    # upper tril indeces
    triu_indeces = np.triu_indices(n_molecules, k=1)
    rmsd_avg = np.mean(rmsd_mat[triu_indeces])
    irmsd_avg = np.mean(irmsd_mat[triu_indeces])

    assert not np.isclose(rmsd_avg, 0.0, atol=1e-6)
    assert np.isclose(irmsd_avg, 0.0, atol=1e-6)
