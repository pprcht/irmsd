from pathlib import Path

import numpy as np

from irmsd import Molecule

ATOM_SYMBOLS2NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
}


def get_atom_num_and_pos_from_xyz(xyzstring):
    """Helper function to extract atom numbers and positions from an XYZ
    string."""
    lines = xyzstring.strip().splitlines()
    num_atoms = int(lines[0])
    atom_numbers = np.zeros(num_atoms, dtype=np.int32)
    positions = np.zeros((num_atoms, 3), dtype=np.float64)
    for i in range(num_atoms):
        parts = lines[i + 2].split()
        atom_symbol = parts[0].upper()
        x, y, z = map(float, parts[1:4])
        atom_numbers[i] = ATOM_SYMBOLS2NUMBERS[atom_symbol]
        positions[i] = [x, y, z]

    return atom_numbers, positions


# Reference set with the expected point group (lowercase) as comment line.
SYMMETRIES_XYZ = Path(__file__).parents[1] / "data" / "symmetries.xyz"


def symmetry_references():
    with open(SYMMETRIES_XYZ) as f:
        lines = f.read().splitlines()
    refs = []
    i = 0
    while i < len(lines):
        nat = int(lines[i])
        label = lines[i + 1].strip()
        block = [ln.split() for ln in lines[i + 2 : i + 2 + nat]]
        symbols = [b[0] for b in block]
        positions = np.array([b[1:4] for b in block], dtype=float)
        refs.append((label, Molecule(symbols=symbols, positions=positions)))
        i += nat + 2
    return refs
