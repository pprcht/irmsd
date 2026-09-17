from __future__ import annotations

from dataclasses import dataclass
from math import gcd

import numpy as np

from ..bindings import symmetry_exposed as _F

# transform_type codes of the Fortran analyzer
_MIRROR, _INVERT, _ROTATE, _ROTREFLECT = 1, 2, 3, 4

# Max-abs matrix difference; distinct operations differ by far more.
_SAME_OP_TOL = 0.1


@dataclass
class SymmetryOperation:
    """Symmetry operation ``x' = matrix @ x + translation``.

    Atom ``i`` maps onto ``permutation[i]``, so ``op.apply(pos)[i]``
    matches ``pos[op.permutation[i]]`` within the symmetry tolerance.

    Attributes
    ----------
    label : str
        e.g. "E", "C3", "C3^2", "S6^5", "i", "sigma", "Cinf".
    kind : str
        One of "E", "C", "S", "i", "sigma".
    order : int
        n of C_n^k / S_n^k; 0 for Cinf, 1 for E, 2 for i and sigma.
    power : int
        k of C_n^k / S_n^k; 1 for E, i, sigma.
    matrix : (3, 3) ndarray of float64
        Orthogonal part.
    translation : (3,) ndarray of float64
        In Å; nonzero when the element misses the origin.
    axis : (3,) ndarray of float64 or None
        Rotation axis, or plane normal for sigma; None for E and i.
    point : (3,) ndarray of float64 or None
        A point on the element in Å; None for E.
    permutation : (N,) ndarray of int32
        0-based image atom of each atom.
    maxdev : float or None
        Largest atom displacement in Å; None for derived operations.
    """

    label: str
    kind: str
    order: int
    power: int
    matrix: np.ndarray
    translation: np.ndarray
    axis: np.ndarray | None
    point: np.ndarray | None
    permutation: np.ndarray
    maxdev: float | None = None

    def apply(self, positions: np.ndarray) -> np.ndarray:
        """Apply the operation to (N, 3) positions in Å."""
        return np.asarray(positions) @ self.matrix.T + self.translation


def _classify(proper: bool, n: int, k: int) -> tuple[str, int, int, str]:
    """Reduced (kind, order, power, label) of C_n^k (proper) or S_n^k."""
    if proper or k % 2 == 0:
        # S_n^k with even k is C_n^k
        k %= n
        if k == 0:
            return "E", 1, 1, "E"
        g = gcd(k, n)
        n, k = n // g, k // g
        return "C", n, k, f"C{n}" if k == 1 else f"C{n}^{k}"
    g = gcd(k, n)
    n, k = n // g, k // g
    if n == 1:
        return "sigma", 2, 1, "sigma"
    if n == 2:
        return "i", 2, 1, "i"
    return "S", n, k, f"S{n}" if k == 1 else f"S{n}^{k}"


def _element_to_operation(el: dict, i: int) -> SymmetryOperation:
    etype, n = int(el["type"][i]), int(el["order"][i])
    if etype == _MIRROR:
        kind, order, label = "sigma", 2, "sigma"
    elif etype == _INVERT:
        kind, order, label = "i", 2, "i"
    elif n == 0:
        kind, order, label = "C", 0, "Cinf"
    else:
        kind, order, _, label = _classify(etype == _ROTATE, n, 1)
    return SymmetryOperation(
        label=label,
        kind=kind,
        order=order,
        power=1,
        matrix=el["matrix"][i].copy(),
        translation=el["translation"][i].copy(),
        axis=None if etype == _INVERT else el["axis"][i].copy(),
        point=el["point"][i].copy(),
        permutation=el["permutation"][i].copy(),
        maxdev=float(el["maxdev"][i]),
    )


_SETTINGS_DOC = """threshold : float, optional
        Final symmetry tolerance in Bohr: the largest atom displacement an
        element may cause to be accepted (default: 0.1, as in CREST).
        Increase for noisy or loosely optimized geometries; large systems
        need this more, since the worst atom deviation grows with size.
        Too tight a value on such structures is also slow, as every
        near-miss candidate is fully optimized before being rejected.
    primary_threshold : float, optional
        Tolerance in Bohr for the initial pairing of symmetry-related atoms
        and the distance prescreening of candidates (default: 0.5). Too small
        values miss elements of distorted structures; larger values test
        many more candidates, which dominates the cost for large systems.
    max_axis_order : int, optional
        Highest rotation axis order searched, for proper and improper axes
        alike (default: 10). An S2n axis needs at least 2n, so values below
        10 lose e.g. the S10 axes of Ih and D5d and with them the group.
    max_opt_cycles : int, optional
        Maximum optimization cycles per candidate element (default: 100).
    max_atoms : int or None, optional
        Skip the analysis for structures with more atoms than this, since the
        search scales steeply with system size (default: 200, as in CREST).
        None disables the limit."""


def _with_settings_doc(func):
    func.__doc__ = func.__doc__.replace("{settings}", _SETTINGS_DOC)
    return func


def _run(
    raw,
    atom_numbers,
    positions,
    threshold,
    primary_threshold,
    max_axis_order,
    max_opt_cycles,
    max_atoms,
):
    """Validate input and call a raw binding; None if skipped for size."""
    atom_numbers = np.ascontiguousarray(atom_numbers, dtype=np.int32)
    pos = np.ascontiguousarray(positions, dtype=np.float64)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = int(pos.shape[0])
    if max_atoms is not None and n > max_atoms:
        return None
    if max_axis_order < 2:
        raise ValueError("max_axis_order must be at least 2")
    coords_flat = pos.reshape(-1).copy(order="C")
    return raw(
        n,
        atom_numbers,
        coords_flat,
        threshold,
        primary_threshold,
        max_axis_order,
        max_opt_cycles,
    )


@_with_settings_doc
def get_point_group(
    atom_numbers: np.ndarray,
    positions: np.ndarray,
    *,
    threshold: float = 0.1,
    primary_threshold: float = 0.5,
    max_axis_order: int = 10,
    max_opt_cycles: int = 100,
    max_atoms: int | None = 200,
) -> str | None:
    """Schoenflies point group from the brute-force analyzer ported from CREST.

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.
    {settings}

    Returns
    -------
    str or None
        Schoenflies symbol, e.g. "C2v", "Dinfh"; the highest proper axis
        (e.g. "C3") if no tabulated group matches; None if skipped by
        ``max_atoms``.

    Raises
    ------
    ValueError
        If positions is not (N, 3) or max_axis_order < 2.
    """
    return _run(
        _F.get_symmetry_fortran_raw,
        atom_numbers,
        positions,
        threshold,
        primary_threshold,
        max_axis_order,
        max_opt_cycles,
        max_atoms,
    )


@_with_settings_doc
def get_symmetry_elements(
    atom_numbers: np.ndarray,
    positions: np.ndarray,
    *,
    threshold: float = 0.1,
    primary_threshold: float = 0.5,
    max_axis_order: int = 10,
    max_opt_cycles: int = 100,
    max_atoms: int | None = 200,
) -> tuple[str | None, list[SymmetryOperation]]:
    """Point group and located symmetry elements, each as its generating operation.

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.
    {settings}

    Returns
    -------
    symbol : str or None
        Schoenflies symbol; None if skipped by ``max_atoms``.
    elements : list[SymmetryOperation]
        Empty if skipped. A Cinf axis has ``order == 0``.

    Raises
    ------
    ValueError
        If positions is not (N, 3) or max_axis_order < 2.
    """
    res = _run(
        _F.get_symmetry_elements_fortran_raw,
        atom_numbers,
        positions,
        threshold,
        primary_threshold,
        max_axis_order,
        max_opt_cycles,
        max_atoms,
    )
    if res is None:
        return None, []
    symbol, el = res
    return symbol, [_element_to_operation(el, i) for i in range(len(el["type"]))]


@_with_settings_doc
def get_symmetry_operations(
    atom_numbers: np.ndarray,
    positions: np.ndarray,
    *,
    threshold: float = 0.1,
    primary_threshold: float = 0.5,
    max_axis_order: int = 10,
    max_opt_cycles: int = 100,
    max_atoms: int | None = 200,
) -> tuple[str | None, list[SymmetryOperation]]:
    """Point group and all distinct powers of the located elements.

    Order: E, i, mirror planes, proper, then improper rotations. For finite
    groups the length equals the group order; for Cinfv/Dinfh only the
    operations generated by i, sigma and C2 are returned (the Cinf axis
    comes from :func:`get_symmetry_elements`).

    Parameters
    ----------
    atom_numbers : (N,) int32-like
        Atomic numbers (or types).
    positions : (N, 3) float64-like
        Cartesian coordinates in Å.
    {settings}

    Returns
    -------
    symbol : str or None
        Schoenflies symbol; None if skipped by ``max_atoms``.
    operations : list[SymmetryOperation]
        Empty if skipped.

    Raises
    ------
    ValueError
        If positions is not (N, 3) or max_axis_order < 2.
    """
    symbol, elements = get_symmetry_elements(
        atom_numbers,
        positions,
        threshold=threshold,
        primary_threshold=primary_threshold,
        max_axis_order=max_axis_order,
        max_opt_cycles=max_opt_cycles,
        max_atoms=max_atoms,
    )
    if symbol is None:
        return None, []

    nat = np.asarray(atom_numbers).size
    ops = [
        SymmetryOperation(
            label="E",
            kind="E",
            order=1,
            power=1,
            matrix=np.eye(3),
            translation=np.zeros(3),
            axis=None,
            point=None,
            permutation=np.arange(nat, dtype=np.int32),
        )
    ]

    def add(op: SymmetryOperation) -> None:
        for other in ops:
            if np.abs(other.matrix - op.matrix).max() < _SAME_OP_TOL:
                return
        ops.append(op)

    # lower orders first, so C4^2 is labelled by its own C2 axis
    rank = {"i": 0, "sigma": 1, "C": 2, "S": 3}
    for el in sorted(elements, key=lambda e: (rank[e.kind], e.order)):
        if el.kind in ("i", "sigma"):
            add(el)
            continue
        n = el.order
        if n == 0:
            continue
        proper = el.kind == "C"
        kmax = n if proper or n % 2 == 0 else 2 * n
        R, t, perm = el.matrix, el.translation, el.permutation
        for k in range(1, kmax):
            if k > 1:
                R = el.matrix @ R
                t = el.matrix @ t + el.translation
                perm = el.permutation[perm]
            kind, order, power, label = _classify(proper, n, k)
            if kind == "E":
                continue
            add(
                SymmetryOperation(
                    label=label,
                    kind=kind,
                    order=order,
                    power=power,
                    matrix=R.copy(),
                    translation=t.copy(),
                    axis=el.axis,
                    point=el.point,
                    permutation=perm.copy(),
                    maxdev=el.maxdev if k == 1 else None,
                )
            )
    return symbol, ops
