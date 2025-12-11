import numpy as np
import pytest

from irmsd.api.sorter_exposed import delta_irmsd_list, sorter_irmsd, cregen_raw

# -------------------------------------------------------------------
# Fake Fortran backend to intercept calls
# -------------------------------------------------------------------


class FakeFBackend:
    def __init__(self):
        self.sorter_called = False
        self.delta_called = False
        self.sorter_args = None
        self.delta_args = None

    def sorter_exposed_xyz_fortran_raw(
        self,
        nat,
        nall,
        xyzall,
        atall,
        groups,
        rthr,
        iinv,
        allcanon,
        printlvl,
        ethr,
        energies_list,
    ):
        self.sorter_called = True
        self.sorter_args = (
            nat,
            nall,
            xyzall.copy(),
            atall.copy(),
            groups.copy(),
            rthr,
            iinv,
            allcanon,
            printlvl,
            ethr,
            energies_list.copy(),
        )
        # Mutate output arrays to verify propagation
        groups[:] = np.arange(groups.size, dtype=np.int32)

    def delta_irmsd_list_fortran_raw(
        self, nat, nall, xyzall, atall, iinv, delta, allcanon, printlvl
    ):
        self.delta_called = True
        self.delta_args = (
            nat,
            nall,
            xyzall.copy(),
            atall.copy(),
            iinv,
            delta.copy(),
            allcanon,
            printlvl,
        )
        delta[:] = np.linspace(0, 1, len(delta))

    def cregen_fortran_raw(
        self, nat, nall, xyzall, atall, groups, rthr, ethr, bthr, printlvl, energies
    ):
        self.cregen_called = True
        self.cregen_args = (
            nat,
            nall,
            xyzall.copy(),
            atall.copy(),
            groups.copy(),
            rthr,
            ethr,
            bthr,
            printlvl,
            energies.copy(),
        )
        groups[:] = np.arange(groups.size, dtype=np.int32)


@pytest.fixture
def fake_backend(monkeypatch):

    import irmsd.api.sorter_exposed as mod

    fake = FakeFBackend()
    monkeypatch.setattr(mod, "_F", fake)
    return fake


# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------


def make_valid_structs(nall=3, N=5):
    atom_numbers = [np.arange(N, dtype=np.int32) for _ in range(nall)]
    positions = [np.zeros((N, 3), dtype=np.float64) for _ in range(nall)]
    return atom_numbers, positions


# ===================================================================
# sorter_irmsd — VALID CASE
# ===================================================================


def test_sorter_valid(fake_backend):
    atoms, pos = make_valid_structs(nall=3, N=4)
    groups, xyz_out, Z_out = sorter_irmsd(
        atoms, pos, nat=3, rthr=1.5, iinversion=1, allcanon=False, printlvl=0
    )

    assert fake_backend.sorter_called is True
    assert groups.shape == (3,)
    # groups mutated by fake backend
    assert np.all(groups == np.arange(3))


# ===================================================================
# sorter_irmsd — ERROR CASES
# ===================================================================


def test_sorter_empty_atom_list():
    with pytest.raises(ValueError, match="must contain at least one structure"):
        sorter_irmsd([], [], 1, 1.0)


def test_sorter_length_mismatch():
    atoms, pos = make_valid_structs()
    with pytest.raises(ValueError, match="same length"):
        sorter_irmsd(atoms, pos[:-1], 1, 1.0)


def test_sorter_positions_wrong_shape():
    atoms, pos = make_valid_structs()
    pos[0] = np.zeros((5,), dtype=np.float64)  # wrong shape
    with pytest.raises(ValueError, match="positions array must have shape"):
        sorter_irmsd(atoms, pos, 1, 1.0)


def test_sorter_Z_wrong_shape():
    atoms, pos = make_valid_structs()
    atoms[0] = atoms[0].reshape(5, 1)  # wrong shape
    with pytest.raises(ValueError, match="atom_numbers array must be shape"):
        sorter_irmsd(atoms, pos, 1, 1.0)


def test_sorter_nat_out_of_range():
    atoms, pos = make_valid_structs(nall=1, N=4)
    with pytest.raises(ValueError, match="nat must satisfy"):
        sorter_irmsd(atoms, pos, nat=0, rthr=1.0)
    with pytest.raises(ValueError, match="nat must satisfy"):
        sorter_irmsd(atoms, pos, nat=10, rthr=1.0)


def test_sorter_inconsistent_positions():
    atoms, pos = make_valid_structs(nall=3, N=5)
    pos[1] = np.zeros((4, 3), dtype=np.float64)
    with pytest.raises(ValueError, match="positions_list\\[1\\]"):
        sorter_irmsd(atoms, pos, 2, 1.0)


def test_sorter_inconsistent_Z():
    atoms, pos = make_valid_structs(nall=3, N=5)
    atoms[2] = np.arange(4, dtype=np.int32)
    with pytest.raises(ValueError, match="atom_numbers_list\\[2\\]"):
        sorter_irmsd(atoms, pos, 2, 1.0)


# ===================================================================
# delta_irmsd_list — VALID CASE
# ===================================================================


def test_delta_valid(fake_backend):
    atoms, pos = make_valid_structs(nall=3, N=4)
    delta, xyz_out, Z_out = delta_irmsd_list(
        atoms, pos, nat=2, iinversion=1, allcanon=False, printlvl=1
    )

    assert fake_backend.delta_called is True
    assert delta.shape == (3,)
    # fake backend writes linspace output
    assert np.allclose(delta, np.linspace(0, 1, 3))


# ===================================================================
# delta_irmsd_list — ERROR CASES
# ===================================================================


def test_delta_empty_list():
    with pytest.raises(ValueError):
        delta_irmsd_list([], [], 1)


def test_delta_length_mismatch():
    atoms, pos = make_valid_structs()
    with pytest.raises(ValueError):
        delta_irmsd_list(atoms, pos[:-1], 1)


def test_delta_positions_wrong_shape():
    atoms, pos = make_valid_structs()
    pos[0] = np.zeros((5,), dtype=np.float64)
    with pytest.raises(ValueError, match="positions array must have shape"):
        delta_irmsd_list(atoms, pos, 1)


def test_delta_Z_wrong_shape():
    atoms, pos = make_valid_structs()
    atoms[0] = atoms[0].reshape(5, 1)
    with pytest.raises(ValueError, match="atom_numbers array must be shape"):
        delta_irmsd_list(atoms, pos, 1)


def test_delta_nat_out_of_range():
    atoms, pos = make_valid_structs(nall=1, N=4)
    with pytest.raises(ValueError):
        delta_irmsd_list(atoms, pos, nat=0)
    with pytest.raises(ValueError):
        delta_irmsd_list(atoms, pos, nat=6)


def test_delta_inconsistent_positions():
    atoms, pos = make_valid_structs(nall=3, N=5)
    pos[2] = np.zeros((4, 3), dtype=np.float64)
    with pytest.raises(ValueError):
        delta_irmsd_list(atoms, pos, nat=2)


def test_delta_inconsistent_Z():
    atoms, pos = make_valid_structs(nall=3, N=5)
    atoms[1] = np.arange(4, dtype=np.int32)
    with pytest.raises(ValueError):
        delta_irmsd_list(atoms, pos, nat=2)


# ===================================================================
# cregen — VALID CASE
# ===================================================================


def test_cregen_valid(fake_backend):
    atoms, pos = make_valid_structs(nall=3, N=4)
    energies = np.zeros((3,), dtype=np.float64)
    groups, xyz_out, energies_out = cregen_raw(
        atoms, pos, energies, nat=3, rthr=0.125, ethr=8.0e-5, bthr=0.01, printlvl=0
    )
    assert fake_backend.cregen_called is True
    assert groups.shape == (3,)
    # groups mutated by fake backend
    assert np.all(groups == np.arange(3))


# ===================================================================
# cregen — ERROR CASES
# ===================================================================


def test_cregen_empty_atom_list():
    with pytest.raises(ValueError, match="must contain at least one structure"):
        cregen_raw([], [], [], 1, 1.0)


def test_cregen_length_mismatch():
    atoms, pos = make_valid_structs()
    energies = np.zeros((1,), dtype=np.float64)
    with pytest.raises(ValueError, match="same length"):
        cregen_raw(atoms, pos[:-1], energies, 1, 1.0)
    with pytest.raises(ValueError, match="same length"):
        cregen_raw(atoms, pos, energies[:-1], 1, 1.0)


def test_cregen_positions_wrong_shape():
    atoms, pos = make_valid_structs()
    pos[0] = np.zeros((5,), dtype=np.float64)  # wrong shape
    energies = np.zeros((3,), dtype=np.float64)
    with pytest.raises(ValueError, match="positions array must have shape"):
        cregen_raw(atoms, pos, energies, 1, 1.0)


def test_cregen_Z_wrong_shape():
    atoms, pos = make_valid_structs()
    atoms = [a.reshape(5, 1) for a in atoms]
    energies = np.zeros((3,), dtype=np.float64)
    with pytest.raises(ValueError, match="atom_numbers array must be shape"):
        cregen_raw(atoms, pos, energies, 1, 1.0)


def test_cregen_atoms_wrong_order():
    atoms, pos = make_valid_structs()
    atoms[0] = np.full((3,), 1, dtype=np.int32)
    energies = np.zeros((3,), dtype=np.float64)
    with pytest.raises(ValueError, match="differs from the first array"):
        cregen_raw(atoms, pos, energies, 1, 1.0)


def test_cregen_nat_out_of_range():
    atoms, pos = make_valid_structs(nall=1, N=4)
    energies = np.zeros((1,), dtype=np.float64)
    with pytest.raises(ValueError, match="nat must satisfy"):
        cregen_raw(atoms, pos, energies, nat=0, rthr=1.0)
    with pytest.raises(ValueError, match="nat must satisfy"):
        cregen_raw(atoms, pos, energies, nat=10, rthr=1.0)


def test_cregen_inconsistent_positions():
    atoms, pos = make_valid_structs(nall=3, N=5)
    pos[1] = np.zeros((4, 3), dtype=np.float64)
    energies = np.zeros((3,), dtype=np.float64)
    with pytest.raises(ValueError, match="positions_list\\[1\\]"):
        cregen_raw(atoms, pos, energies, 2, 1.0)


def test_cregen_inconsistent_Z():
    atoms, pos = make_valid_structs(nall=3, N=5)
    atoms[2] = np.arange(4, dtype=np.int32)
    energies = np.zeros((3,), dtype=np.float64)
    with pytest.raises(ValueError, match="atom_numbers_list\\[2\\]"):
        cregen_raw(atoms, pos, energies, 2, 1.0)
