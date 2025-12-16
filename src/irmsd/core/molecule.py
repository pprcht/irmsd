from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

# Very small symbol → Z map so we don't have to rely on ASE at all.
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

# Precomputed inverse table: Z → symbol
_INV_PERIODIC = {Z: sym for sym, Z in _PERIODIC.items()}
# fmt: on


@dataclass
class Molecule:
    """Lightweight replacement for ase.Atoms.

    Raises
    ------
    ValueError
        If input data has incorrect shape or contains unknown symbols.

    Notes
    -----
    This class is designed to be a lightweight alternative to ASE's Atoms
    class to minimize dependencies.
    """

    symbols: list[str]
    positions: np.ndarray
    energy: float | None = None
    info: dict[str, Any] = field(default_factory=dict)
    cell: np.ndarray | None = None
    pbc: tuple[bool, bool, bool] | None = None

    def __post_init__(self) -> None:
        # Normalize symbols
        self.symbols = [str(s) for s in self.symbols]
        n = len(self.symbols)

        # Convert symbols → atomic numbers
        try:
            self.numbers = np.ascontiguousarray(
                [_PERIODIC[s] for s in self.symbols],
                dtype=np.int32,
            )
        except KeyError as e:
            raise ValueError(f"Unknown chemical symbol: {e.args[0]!r}")

        # Normalize positions
        self.positions = np.ascontiguousarray(self.positions, dtype=np.float64)
        if self.positions.shape != (n, 3):
            raise ValueError(
                f"positions must have shape ({n}, 3), got {self.positions.shape}"
            )

        # Normalize cell
        if self.cell is not None:
            self.cell = np.asarray(self.cell, dtype=np.float64)
            if self.cell.shape != (3, 3):
                raise ValueError("cell must be (3,3)")

        # Normalize pbc
        if self.pbc is not None:
            if len(self.pbc) != 3:
                raise ValueError("pbc must be length-3")
            self.pbc = tuple(bool(x) for x in self.pbc)

        self.info = dict(self.info)

    # --- Basic info ------------------------------------------------------------
    @property
    def natoms(self) -> int:
        return len(self.symbols)

    def __len__(self) -> int:
        return self.natoms

    # --- Minimal ASE-like API -------------------------------------------------------
    def get_chemical_symbols(self) -> list[str]:
        return list(self.symbols)

    def get_atomic_numbers(self) -> np.ndarray:
        """Return atomic numbers as int32 array."""
        return self.numbers.copy()

    def get_positions(self, copy: bool = True) -> np.ndarray:
        """Return atomic positions.

        Parameters
        ----------
        copy : bool
            If True, return a copy of the positions array.
            If False, return the internal array (may be modified).
        Returns
        -------
        positions : (N, 3) ndarray of float64
        """
        return self.positions.copy() if copy else self.positions

    def get_potential_energy(self) -> float:
        if self.energy is None:
            raise AttributeError("Potential energy not set.")
        return float(self.energy)

    def get_chemical_formula(self, mode: str = "hill") -> str:
        """Return a chemical formula string.

        Parameters
        ----------
        mode : str
            "hill"  → C, H first, then alphabetical (standard Hill formula)
            others → alphabetical order of all elements

        Returns
        -------
        formula : str
        """
        # Count symbols
        from collections import Counter

        counts = Counter(self.symbols)

        # Sorting rules
        if mode.lower() == "hill":
            order = []
            # Hill system: C first, H second
            if "C" in counts:
                order.append("C")
            if "H" in counts:
                order.append("H")

            # Then all others alphabetically
            others = sorted(sym for sym in counts if sym not in ("C", "H"))
            order.extend(others)
        else:
            # Alphabetical mode (same as ASE if not "hill")
            order = sorted(counts.keys())

        # Build formula: omit "1" as ASE does
        fragments = []
        for sym in order:
            n = counts[sym]
            if n == 1:
                fragments.append(sym)
            else:
                fragments.append(f"{sym}{n}")

        return "".join(fragments)

    def copy(self) -> "Molecule":
        """Return a deep copy of the Molecule.

        All arrays (positions, numbers, cell) are copied, and both `info` and
        `symbols` are duplicated so that modifying the copy has no effect on
        the original.
        """
        return Molecule(
            symbols=list(self.symbols),
            positions=self.positions.copy(),
            energy=self.energy,
            info=dict(self.info),
            cell=None if self.cell is None else self.cell.copy(),
            pbc=None if self.pbc is None else tuple(self.pbc),
        )

    # --- Optional setters ------------------------------------------------------
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

    # --- Larger Functions interfaced to Fortran ---------------------------------
    def get_cn(self) -> np.ndarray:
        """Optional utility: calls core get_cn_fortran and returns
        a numpy array with the coordination numbers per atom.

        Returns
        -------
        cn : (N,) ndarray of float64
        """
        from ..api.cn_exposed import get_cn_fortran

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        new_cn = get_cn_fortran(Z, pos)
        return new_cn

    def get_axis(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Optional utility: calls core get_axis and returns
        rotation constants in MHz, average momentum in a.u.,
        and the rotation matrix.

        Returns
        -------
        rot : (3,) ndarray of float64
            Rotation constants in MHz.
        avmom : (1,) ndarray of float64
            Average momentum in a.u.
        evec : (3, 3) ndarray of float64
        """
        from ..api.axis_exposed import get_axis

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        rot, avmom, evec = get_axis(Z, pos)
        return rot, avmom, evec

    def get_canonical(
        self,
        wbo: np.ndarray | None = None,
        invtype: str = "apsp+",
        heavy: bool = False,
    ) -> np.ndarray:
        """Optional utility: calls core get_canonical_fortran and
        returns the rank (and/or invariants, depending on backend).

        Parameters
        ----------
        wbo : (N, N) ndarray of float64, optional
            Wiberg bond order matrix, required if invtype is 'cangen'.
        invtype : str, optional
            Algorithm type for invariants calculation (default: 'apsp+'),
            alternatively 'cangen'.
        heavy : bool, optional
            Whether to consider only heavy atoms (default: False).

        Returns
        -------
        rank : (N,) ndarray of int32
            Rank array.
        """
        from ..api.canonical_exposed import get_canonical_fortran

        Z = self.get_atomic_numbers()  # (N,)
        pos = self.get_positions()  # (N, 3) float64

        rank = get_canonical_fortran(Z, pos, wbo=wbo, invtype=invtype, heavy=heavy)
        return rank
