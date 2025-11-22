import pytest

from irmsd.sorting import group_by


def test_group_by_integers_parity():
    nums = [1, 2, 3, 4, 5, 6]
    groups = group_by(nums, key=lambda x: x % 2)

    # We expect two groups: odd (1) and even (0)
    assert set(groups.keys()) == {0, 1}
    assert groups[0] == [2, 4, 6]
    assert groups[1] == [1, 3, 5]


def test_group_by_strings_first_letter():
    words = ["apple", "ape", "bat", "boat", "banana"]
    groups = group_by(words, key=lambda w: w[0])

    assert set(groups.keys()) == {"a", "b"}
    assert groups["a"] == ["apple", "ape"]
    # Order must be preserved within each group
    assert groups["b"] == ["bat", "boat", "banana"]


def test_group_by_strings_length():
    words = ["a", "bb", "ccc", "dd", "eee", "f"]
    groups = group_by(words, key=len)

    assert set(groups.keys()) == {1, 2, 3}
    assert groups[1] == ["a", "f"]
    assert groups[2] == ["bb", "dd"]
    assert groups[3] == ["ccc", "eee"]


def test_group_by_atoms_formula():
    ase = pytest.importorskip("ase")
    molecule = pytest.importorskip("ase.build").molecule

    atoms_list = [
        molecule("H2O"),
        molecule("CH4"),
        molecule("H2O"),
        molecule("NH3"),
        molecule("CH4"),
    ]

    groups = group_by(
        atoms_list,
        key=lambda a: a.get_chemical_formula(mode="hill"),
    )

    # We expect three distinct formulas
    assert set(groups.keys()) == {"H2O", "CH4", "H3N"} or set(groups.keys()) == {"H2O", "CH4", "NH3"}

    # Check counts & identity via formula
    assert len(groups["H2O"]) == 2
    assert len(groups["CH4"]) == 2

    # The third group is the ammonia variant (name may depend on ASE version)
    other_keys = [k for k in groups.keys() if k not in {"H2O", "CH4"}]
    assert len(other_keys) == 1
    nh3_key = other_keys[0]
    assert len(groups[nh3_key]) == 1

