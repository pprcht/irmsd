import numpy as np
import pytest

from irmsd.core.molecule import Molecule

# =====================================================================
# SYMBOL VALIDATION
# =====================================================================


def test_unknown_symbol():
    with pytest.raises(ValueError, match="Unknown chemical symbol"):
        Molecule(symbols=["C", "Xx", "H"], positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])


def test_numbers_correspond_to_symbols():
    m = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    assert m.numbers[0] == 6
    assert m.numbers[1] == 1


# =====================================================================
# POSITIONS VALIDATION
# =====================================================================


def test_positions_wrong_shape():
    with pytest.raises(ValueError, match="positions must have shape"):
        Molecule(symbols=["H", "H"], positions=[[0, 0, 0]])


def test_positions_wrong_dtype_coerced():
    m = Molecule(symbols=["H"], positions=[[1, 2, 3]])
    assert m.positions.dtype == np.float64


def test_set_positions_wrong_shape():
    m = Molecule(symbols=["H", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    with pytest.raises(ValueError, match="wrong shape"):
        m.set_positions([[0, 0, 0]])  # wrong length


def test_set_positions_valid():
    m = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    new = np.array([[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]])
    m.set_positions(new)
    assert np.allclose(m.positions, new)


# =====================================================================
# CELL VALIDATION
# =====================================================================


def test_cell_wrong_shape():
    with pytest.raises(ValueError, match="cell must be \\(3,3\\)"):
        Molecule(symbols=["H"], positions=[[0, 0, 0]], cell=np.zeros((2, 3)))


def test_cell_correct_shape():
    cell = np.eye(3)
    m = Molecule(symbols=["H"], positions=[[0, 0, 0]], cell=cell)
    assert np.allclose(m.cell, cell)


# =====================================================================
# PBC VALIDATION
# =====================================================================


def test_pbc_wrong_length():
    with pytest.raises(ValueError, match="pbc must be length-3"):
        Molecule(symbols=["H"], positions=[[0, 0, 0]], pbc=(True, False))


def test_pbc_converted_to_tuple_of_bool():
    m = Molecule(symbols=["H"], positions=[[0, 0, 0]], pbc=[1, 0, 1])
    assert m.pbc == (True, False, True)


# =====================================================================
# set_atomic_numbers VALIDATION
# =====================================================================


def test_set_atomic_numbers_wrong_length():
    m = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    with pytest.raises(ValueError, match="Expected 2 atomic numbers"):
        m.set_atomic_numbers([6])


def test_set_atomic_numbers_unknown_Z():
    m = Molecule(symbols=["C"], positions=[[0, 0, 0]])
    with pytest.raises(KeyError, match="Unknown atomic number"):
        m.set_atomic_numbers([999])


def test_set_atomic_numbers_valid():
    m = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]])
    m.set_atomic_numbers([1, 6])  # swap symbols
    assert list(m.get_chemical_symbols()) == ["H", "C"]
    assert list(m.get_atomic_numbers()) == [1, 6]


# =====================================================================
# ENERGY VALIDATION
# =====================================================================


def test_get_potential_energy_missing():
    m = Molecule(symbols=["H"], positions=[[0, 0, 0]])
    with pytest.raises(AttributeError):
        m.get_potential_energy()


def test_set_potential_energy_accepts_none_and_float():
    m = Molecule(symbols=["H"], positions=[[0, 0, 0]])
    m.set_potential_energy(12.3)
    assert m.get_potential_energy() == 12.3
    m.set_potential_energy(None)
    with pytest.raises(AttributeError):
        m.get_potential_energy()


# =====================================================================
# COPY VALIDATION
# =====================================================================


def test_copy_is_deep():
    m = Molecule(symbols=["C", "H"], positions=[[0, 0, 0], [1, 0, 0]], info={"a": 1})
    c = m.copy()

    # Must be independent
    c.positions[0, 0] = 99
    assert not np.allclose(m.positions, c.positions)

    c.info["a"] = 2
    assert m.info["a"] == 1

    c.symbols[0] = "H"
    assert m.symbols[0] == "C"
