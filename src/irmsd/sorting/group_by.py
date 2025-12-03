from collections import defaultdict
from typing import Callable, Dict, Hashable, Iterable, List, TypeVar

T = TypeVar("T")  # item type
K = TypeVar("K", bound=Hashable)  # key type


def group_by(items: Iterable[T], key: Callable[[T], K]) -> Dict[K, List[T]]:
    """Group arbitrary items according to a key function.

    This is a fully general utility: it works for integers, strings, objects,
    or domain-specific types such as ASE Atoms.

    Parameters
    ----------
    items : iterable of T
        The objects to be grouped.
    key : callable
        Function mapping an item -> hashable key (str, int, tuple, ...).

    Returns
    -------
    groups : dict
        Mapping key -> list of items having that key.

    Examples
    --------
    **1. Grouping ASE Atoms by chemical formula**

    >>> from ase.build import molecule
    >>> atoms_list = [molecule("H2O"), molecule("CH4"), molecule("H2O")]
    >>> groups = group_by(atoms_list, key=lambda a: a.get_chemical_formula())
    >>> groups.keys()
    dict_keys(['H2O', 'CH4'])

    **2. Grouping integers by parity**

    >>> nums = [1, 2, 3, 4, 5]
    >>> group_by(nums, key=lambda x: x % 2)
    {1: [1, 3, 5], 0: [2, 4]}

    **3. Grouping strings by first character**

    >>> words = ["apple", "ape", "bat", "boat"]
    >>> group_by(words, key=lambda w: w[0])
    {'a': ['apple', 'ape'], 'b': ['bat', 'boat']}
    """
    groups: Dict[K, List[T]] = defaultdict(list)
    for item in items:
        groups[key(item)].append(item)
    return dict(groups)
