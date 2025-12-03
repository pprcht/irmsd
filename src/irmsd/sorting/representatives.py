from __future__ import annotations

from typing import List, Sequence, Set, SupportsInt, TypeVar

T = TypeVar("T")


def first_by_assignment(
    items: Sequence[T],
    assignments: Sequence[SupportsInt],
) -> list[T]:
    """Select the first representative of each assignment label.

    Parameters
    ----------
    items
        Sequence of arbitrary objects (e.g. ASE Atoms, numpy arrays, ints, ...).
    assignments
        Sequence of integer-like labels (same length as `items`).
        The i-th entry is the assignment for items[i].

    Returns
    -------
    list[T]
        A new list containing the first occurrence of each assignment,
        in the order the *labels* first appear.

    Example
    -------
    >>> items = list("abcdefgh")
    >>> assignments = [1, 1, 1, 2, 3, 3, 3, 3]
    >>> first_by_assignment(items, assignments)
    ['a', 'd', 'e']
    """
    if len(items) != len(assignments):
        raise ValueError(
            f"items and assignments must have the same length "
            f"(got {len(items)} and {len(assignments)})"
        )

    seen: Set[int] = set()
    selected: List[T] = []

    for item, label in zip(items, assignments):
        ilabel = int(label)  # works for int, numpy.int32, etc.
        if ilabel not in seen:
            seen.add(ilabel)
            selected.append(item)

    return selected
