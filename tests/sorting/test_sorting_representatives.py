import pytest
import numpy as np
from irmsd.sorting import first_by_assignment


def test_simple_python_objects():
    items = list("abcdefgh")
    assignments = [1, 1, 1, 2, 3, 3, 3, 3]

    result = first_by_assignment(items, assignments)

    assert result == ["a", "d", "e"]


def test_numpy_assignments():
    items = list(range(8))
    assignments = np.array([1, 1, 1, 2, 3, 3, 3, 3], dtype=np.int32)

    result = first_by_assignment(items, assignments)

    assert result == [0, 3, 4]


def test_error_on_mismatched_lengths():
    items = [1, 2, 3]
    assignments = [1, 1]

    with pytest.raises(ValueError):
        first_by_assignment(items, assignments)


def test_generic_object_sequence():
    items = [np.array([i, i, i]) for i in range(8)]
    assignments = [0, 0, 1, 1, 2, 2, 2, 2]

    result = first_by_assignment(items, assignments)

    assert len(result) == 3
    assert np.array_equal(result[0], np.array([0, 0, 0]))
    assert np.array_equal(result[1], np.array([2, 2, 2]))
    assert np.array_equal(result[2], np.array([4, 4, 4]))

