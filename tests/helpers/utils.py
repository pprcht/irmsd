import numpy as np

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
