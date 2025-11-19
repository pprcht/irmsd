from __future__ import annotations

from typing import TypeVar, SupportsFloat, SupportsInt, List, Tuple, Sequence

T = TypeVar("T")
V = TypeVar("V", bound=SupportsFloat | SupportsInt)


def sort_by_value(
    items: Sequence[T],
    values: Sequence[V],
    *,
    descending: bool = False,
) -> tuple[list[T], list[V]]:
    """
    Sort items according to the corresponding numeric values.

    Parameters
    ----------
    items
        Sequence of arbitrary objects (ASE Atoms, numpy arrays, ints, ...).
    values
        Sequence of numeric values (same length as `items`).
        Must contain int-, float-, or numpy-number-like entries.
    descending
        If True, sort from largest to smallest.

    Returns
    -------
    (sorted_items, sorted_values) : (list[T], list[V])
        The items and their values sorted together.

    Example
    -------
    >>> items = ["a", "b", "c"]
    >>> values = [3.0, 1.0, 2.0]
    >>> sort_by_value(items, values)
    (['b', 'c', 'a'], [1.0, 2.0, 3.0])

    >>> sort_by_value(items, values, descending=True)
    (['a', 'c', 'b'], [3.0, 2.0, 1.0])
    """
    if len(items) != len(values):
        raise ValueError(
            f"items and values must have same length (got {len(items)} and {len(values)})"
        )

    # Convert to list of (value, index) pairs
    # Using float(...) ensures numpy scalars behave consistently.
    paired = [(float(v), i) for i, v in enumerate(values)]

    # Sort by the float-converted value
    paired.sort(reverse=descending)

    # Rebuild sorted items and values
    sorted_items: List[T] = [items[i] for _, i in paired]
    sorted_values: List[V] = [values[i] for _, i in paired]

    return sorted_items, sorted_values

