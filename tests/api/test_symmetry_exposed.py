import numpy as np
import pytest
from helpers.utils import symmetry_references as _reference_structures

from irmsd.api.symmetry_exposed import (
    get_point_group,
    get_symmetry_elements,
    get_symmetry_operations,
)


@pytest.mark.parametrize(
    "label,mol", _reference_structures(), ids=lambda x: x if isinstance(x, str) else ""
)
def test_get_point_group_references(label, mol):
    pg = get_point_group(mol.get_atomic_numbers(), mol.get_positions())
    assert pg.lower() == label


def test_get_point_group_invariant_to_rigid_motion():
    mol = dict(_reference_structures())["c3v"]
    rng = np.random.default_rng(42)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    pos = mol.get_positions() @ q.T + np.array([3.0, -1.5, 7.0])
    assert get_point_group(mol.get_atomic_numbers(), pos) == "C3v"


def test_get_point_group_threshold():
    mol = dict(_reference_structures())["c2v"]
    pos = mol.get_positions().copy()
    pos[1, 0] += 0.02  # breaks the water C2 axis by ~0.04 Bohr
    Z = mol.get_atomic_numbers()
    assert get_point_group(Z, pos, threshold=0.001) == "Cs"
    assert get_point_group(Z, pos, threshold=0.1) == "C2v"


def test_get_point_group_max_atoms():
    mol = dict(_reference_structures())["c2v"]
    Z, pos = mol.get_atomic_numbers(), mol.get_positions()
    assert get_point_group(Z, pos, max_atoms=2) is None
    assert get_point_group(Z, pos, max_atoms=None) == "C2v"


def test_get_point_group_bad_shape():
    with pytest.raises(ValueError):
        get_point_group(np.array([1, 1]), np.zeros((2, 2)))


def test_molecule_get_point_group():
    mol = dict(_reference_structures())["td"]
    assert mol.get_point_group() == "Td"


GROUP_ORDERS = {
    "c1": 1, "ci": 2, "cs": 2, "c2": 2, "c2h": 4, "c2v": 4, "c3": 3,
    "c3v": 6, "c4v": 8, "c5": 5, "c5v": 10, "d2": 4, "d2d": 8, "d2h": 8,
    "d3": 6, "d3d": 12, "d3h": 12, "d4": 8, "d4h": 16, "d5d": 20,
    "d5h": 20, "d6h": 24, "d7h": 28, "d8h": 32, "s4": 4, "td": 24,
    "oh": 48, "ih": 120,
}  # fmt: skip

# Elements are least-squares fits, so the worst atom of a slightly
# asymmetric reference may exceed the 0.1 Bohr detection tolerance a bit.
TOL_ANG = 2 * 0.1 * 0.52917726


def _find(ops, R):
    return [k for k, o in enumerate(ops) if np.abs(o.matrix - R).max() < 0.1]


@pytest.mark.parametrize(
    "label,mol",
    [r for r in _reference_structures() if r[0] in GROUP_ORDERS],
    ids=lambda x: x if isinstance(x, str) else "",
)
def test_symmetry_operations_references(label, mol):
    Z, pos = mol.get_atomic_numbers(), mol.get_positions()
    symbol, ops = get_symmetry_operations(Z, pos)
    assert symbol.lower() == label
    assert len(ops) == GROUP_ORDERS[label]
    assert ops[0].label == "E"

    for op in ops:
        np.testing.assert_allclose(op.matrix @ op.matrix.T, np.eye(3), atol=1e-8)
        assert np.all(Z[op.permutation] == Z)
        assert sorted(op.permutation) == list(range(len(Z)))
        assert np.abs(op.apply(pos) - pos[op.permutation]).max() < TOL_ANG

    # closed under composition, with permutations composing alongside
    for a in ops:
        for b in ops:
            hit = _find(ops, a.matrix @ b.matrix)
            assert len(hit) == 1
            assert np.all(ops[hit[0]].permutation == a.permutation[b.permutation])


def test_symmetry_operation_labels():
    mol = dict(_reference_structures())["td"]
    _, ops = get_symmetry_operations(mol.get_atomic_numbers(), mol.get_positions())
    counts = {}
    for op in ops:
        key = op.label.split("^")[0]
        counts[key] = counts.get(key, 0) + 1
    assert counts == {"E": 1, "C3": 8, "C2": 3, "S4": 6, "sigma": 6}


def test_symmetry_elements_linear():
    mol = dict(_reference_structures())["dinfh"]
    Z, pos = mol.get_atomic_numbers(), mol.get_positions()
    symbol, els = get_symmetry_elements(Z, pos)
    assert symbol == "Dinfh"
    cinf = [e for e in els if e.label == "Cinf"]
    assert len(cinf) == 1 and cinf[0].order == 0
    # the Cinf axis runs along the molecule
    bond = pos[1] - pos[0]
    assert abs(abs(np.dot(cinf[0].axis, bond)) - np.linalg.norm(bond)) < 1e-3

    _, ops = get_symmetry_operations(Z, pos)
    assert all(op.label != "Cinf" for op in ops)


def test_symmetry_operations_off_origin():
    mol = dict(_reference_structures())["c3v"]
    pos = mol.get_positions() + np.array([10.0, -4.0, 2.5])
    _, ops = get_symmetry_operations(mol.get_atomic_numbers(), pos)
    assert len(ops) == 6
    for op in ops:
        assert np.abs(op.apply(pos) - pos[op.permutation]).max() < TOL_ANG


def test_symmetry_operations_skipped():
    mol = dict(_reference_structures())["c2v"]
    Z, pos = mol.get_atomic_numbers(), mol.get_positions()
    assert get_symmetry_operations(Z, pos, max_atoms=2) == (None, [])
    assert get_symmetry_elements(Z, pos, max_atoms=2) == (None, [])


def test_molecule_get_symmetry_operations():
    mol = dict(_reference_structures())["c2v"]
    symbol, ops = mol.get_symmetry_operations()
    assert symbol == "C2v"
    assert sorted(op.label for op in ops) == ["C2", "E", "sigma", "sigma"]


def test_symmetry_settings_max_axis_order():
    mol = dict(_reference_structures())["ih"]
    Z, pos = mol.get_atomic_numbers(), mol.get_positions()
    assert get_point_group(Z, pos, max_axis_order=10) == "Ih"
    # S10 axes are no longer searched, so Ih is not recognized
    assert get_point_group(Z, pos, max_axis_order=5) != "Ih"
    with pytest.raises(ValueError):
        get_point_group(Z, pos, max_axis_order=1)


def test_symmetry_settings_primary_threshold():
    mol = dict(_reference_structures())["td"]
    Z = mol.get_atomic_numbers()
    rng = np.random.default_rng(1)
    pos = mol.get_positions() + rng.normal(scale=0.03, size=(len(Z), 3))
    assert get_point_group(Z, pos, threshold=0.3) == "Td"
    # candidates are no longer paired up at all
    assert get_point_group(Z, pos, threshold=0.3, primary_threshold=1e-4) == "C1"


def test_symmetry_settings_forwarded():
    mol = dict(_reference_structures())["c2v"]
    assert mol.get_point_group(max_atoms=2) is None
    assert mol.get_symmetry_operations(max_atoms=2) == (None, [])
    assert mol.get_point_group(threshold=0.2, primary_threshold=0.5) == "C2v"
