from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

# Symbol -> Z map, kept local to avoid an ASE dependency.
# fmt: off
_PERIODIC = {
    "H": 1,   "He": 2,
    "Li": 3,  "Be": 4,  "B": 5,   "C": 6,   "N": 7,   "O": 8,   "F": 9,   "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,  "S": 16,  "Cl": 17, "Ar": 18,
    "K": 19,  "Ca": 20, "Sc": 21, "Ti": 22, "V": 23,  "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39,  "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53,  "Xe": 54,
    "Cs": 55, "Ba": 56,
    "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74,  "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79,
    "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88,
    "Ac": 89, "Th": 90, "Pa": 91, "U": 92,  "Np": 93, "Pu": 94, "Am": 95, "Cm": 96,
    "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
    "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
    "Rg": 111, "Cn": 112,
    "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
}

_INV_PERIODIC = {Z: sym for sym, Z in _PERIODIC.items()}
# fmt: on


@dataclass
class Molecule:
    """Dependency-free replacement for ase.Atoms.

    ``energy`` is in Hartree. ``ids`` holds optional per-atom canonical atom
    identifiers, one int32 per atom. Construction raises ValueError on bad
    shapes or unknown symbols.
    """

    symbols: list[str]
    positions: np.ndarray
    energy: float | None = None
    info: dict[str, Any] = field(default_factory=dict)
    cell: np.ndarray | None = None
    pbc: tuple[bool, bool, bool] | None = None
    ids: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.symbols = [str(s) for s in self.symbols]
        n = len(self.symbols)

        try:
            self.numbers = np.ascontiguousarray(
                [_PERIODIC[s] for s in self.symbols],
                dtype=np.int32,
            )
        except KeyError as e:
            raise ValueError(f"Unknown chemical symbol: {e.args[0]!r}")

        self.positions = np.ascontiguousarray(self.positions, dtype=np.float64)
        if self.positions.shape != (n, 3):
            raise ValueError(
                f"positions must have shape ({n}, 3), got {self.positions.shape}"
            )

        if self.cell is not None:
            self.cell = np.asarray(self.cell, dtype=np.float64)
            if self.cell.shape != (3, 3):
                raise ValueError("cell must be (3,3)")

        if self.pbc is not None:
            if len(self.pbc) != 3:
                raise ValueError("pbc must be length-3")
            self.pbc = tuple(bool(x) for x in self.pbc)

        if self.ids is not None:
            self.ids = np.ascontiguousarray(self.ids, dtype=np.int32)
            if self.ids.shape != (n,):
                raise ValueError(
                    f"ids must have shape ({n},), got {self.ids.shape}"
                )

        self.info = dict(self.info)

    @property
    def natoms(self) -> int:
        return len(self.symbols)

    def __len__(self) -> int:
        return self.natoms

    def get_chemical_symbols(self) -> list[str]:
        return list(self.symbols)

    def get_atomic_numbers(self) -> np.ndarray:
        """Return atomic numbers as int32 array."""
        return self.numbers.copy()

    def get_positions(self, copy: bool = True) -> np.ndarray:
        """Return (N, 3) positions; ``copy=False`` returns the internal array."""
        return self.positions.copy() if copy else self.positions

    def get_ids(self, copy: bool = True) -> np.ndarray | None:
        """Return (N,) canonical atom IDs or None; ``copy=False`` returns the internal array."""
        if self.ids is None:
            return None
        return self.ids.copy() if copy else self.ids

    def set_ids(self, ids: Sequence[int] | np.ndarray | None) -> None:
        """Set canonical atom IDs, one per atom (ValueError otherwise); None clears them."""
        if ids is None:
            self.ids = None
            return
        arr = np.ascontiguousarray(ids, dtype=np.int32)
        if arr.shape != (self.natoms,):
            raise ValueError(
                f"ids must have shape ({self.natoms},), got {arr.shape}"
            )
        self.ids = arr

    def get_potential_energy(self) -> float:
        if self.energy is None:
            raise AttributeError("Potential energy not set.")
        return float(self.energy)

    def get_chemical_formula(self, mode: str = "hill") -> str:
        """Return the chemical formula, omitting counts of 1 as ASE does.

        ``mode="hill"`` puts C and H first, then the rest alphabetically; any
        other mode sorts all elements alphabetically, matching ASE.
        """
        from collections import Counter

        counts = Counter(self.symbols)

        if mode.lower() == "hill":
            order = []
            if "C" in counts:
                order.append("C")
            if "H" in counts:
                order.append("H")

            others = sorted(sym for sym in counts if sym not in ("C", "H"))
            order.extend(others)
        else:
            order = sorted(counts.keys())

        fragments = []
        for sym in order:
            n = counts[sym]
            if n == 1:
                fragments.append(sym)
            else:
                fragments.append(f"{sym}{n}")

        return "".join(fragments)

    def copy(self) -> "Molecule":
        """Return a deep copy; no array, ``info`` or ``symbols`` is shared."""
        return Molecule(
            symbols=list(self.symbols),
            positions=self.positions.copy(),
            energy=self.energy,
            info=dict(self.info),
            cell=None if self.cell is None else self.cell.copy(),
            pbc=None if self.pbc is None else tuple(self.pbc),
            ids=None if self.ids is None else self.ids.copy(),
        )

    def set_positions(self, positions: Sequence[Sequence[float]]) -> None:
        new_pos = np.ascontiguousarray(positions, dtype=np.float64)
        if new_pos.shape != self.positions.shape:
            raise ValueError("New positions have wrong shape.")
        self.positions[...] = new_pos

    def set_potential_energy(self, energy: float | None):
        self.energy = None if energy is None else float(energy)

    def set_atomic_numbers(self, numbers: Sequence[int]) -> None:
        if len(numbers) != self.natoms:
            raise ValueError(
                f"Expected {self.natoms} atomic numbers, got {len(numbers)}"
            )
        Z = np.ascontiguousarray(numbers, dtype=np.int32)
        self.numbers = Z
        try:
            self.symbols = [_INV_PERIODIC[int(z)] for z in Z]
        except KeyError as e:
            raise KeyError(f"Unknown atomic number: {e.args[0]}") from e

    def get_cn(self) -> np.ndarray:
        """Return (N,) coordination numbers from the Fortran core."""
        from ..api.cn_exposed import get_cn_fortran

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        new_cn = get_cn_fortran(Z, pos)
        return new_cn

    def get_axis(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return rotation constants, average momentum and rotation matrix.

        Returns
        -------
        rot : (3,) ndarray
            Rotation constants in MHz.
        avmom : (1,) ndarray
            Average momentum in a.u.
        evec : (3, 3) ndarray
        """
        from ..api.axis_exposed import get_axis

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        rot, avmom, evec = get_axis(Z, pos)
        return rot, avmom, evec

    def get_point_group(self, **settings) -> str | None:
        """Return the Schoenflies symbol (e.g. "C2v"), or None if skipped.

        ``settings`` are forwarded; see :func:`irmsd.get_point_group`.
        """
        from ..api.symmetry_exposed import get_point_group

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        return get_point_group(Z, pos, **settings)

    def get_symmetry_operations(self, **settings) -> tuple[str | None, list]:
        """Return the Schoenflies symbol and the symmetry operations.

        ``settings`` are forwarded; see :func:`irmsd.get_point_group`.

        Returns
        -------
        symbol : str or None
            None if skipped.
        operations : list[irmsd.SymmetryOperation]
            Identity first; empty if skipped.
        """
        from ..api.symmetry_exposed import get_symmetry_operations

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        return get_symmetry_operations(Z, pos, **settings)

    def get_canonical(
        self,
        wbo: np.ndarray | None = None,
        invtype: str = "apsp+",
        heavy: bool = False,
    ) -> np.ndarray:
        """Return (N,) int32 canonical ranks from the Fortran core.

        Parameters
        ----------
        wbo : (N, N) ndarray, optional
            Wiberg bond orders; required for ``invtype="cangen"``.
        invtype : {"apsp+", "cangen", "apsp+nmr"}
            ``"apsp+nmr"`` runs apsp+, then splits any rank shared by exactly
            two non-hydrogen atoms (propagated to attached hydrogens), a hack
            for NMR magnetic (in)equivalencies. Rank requests only, not iRMSD.
        heavy : bool
            Rank heavy atoms only.
        """
        from ..api.canonical_exposed import get_canonical_fortran

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        rank = get_canonical_fortran(Z, pos, wbo=wbo, invtype=invtype, heavy=heavy)
        return rank
