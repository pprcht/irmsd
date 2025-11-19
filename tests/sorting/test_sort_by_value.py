import pytest
import numpy as np
from irmsd.sorting import sort_by_value


def test_sort_simple():
    items = ["a", "b", "c"]
    values = [3.0, 1.0, 2.0]

    items_s, values_s = sort_by_value(items, values)

    assert items_s == ["b", "c", "a"]
    assert values_s == [1.0, 2.0, 3.0]


def test_sort_descending():
    items = ["a", "b", "c"]
    values = [3.0, 1.0, 2.0]

    items_s, values_s = sort_by_value(items, values, descending=True)

    assert items_s == ["a", "c", "b"]
    assert values_s == [3.0, 2.0, 1.0]


def test_numpy_values():
    items = list("abcd")
    values = np.array([0.2, -5.0, 10.0, 3.0], dtype=np.float64)

    items_s, values_s = sort_by_value(items, values)

    assert items_s == ["b", "a", "d", "c"]
    assert values_s == [-5.0, 0.2, 3.0, 10.0]


def test_mismatched_length():
    with pytest.raises(ValueError):
        sort_by_value([1, 2, 3], [1.0, 2.0])  # mismatch

