from __future__ import annotations

from typing import Generator, List, Sequence, Tuple, overload

import numpy as np

from irmsd.interfaces.mol_interface import (
    delta_irmsd_list_molecule,
    get_irmsd_molecule,
    get_rmsd_molecule,
    sorter_irmsd_molecule,
    cregen,
    prune,
)

from ..core.molecule import Molecule
from ..utils.utils import require_rdkit


def conformer_iterator(molecule: "Mol", conf_ids: list[int]) -> "Conformer":
    """Yield the conformers of `molecule` listed in `conf_ids`."""
    for conf_id in conf_ids:
        yield molecule.GetConformer(conf_id)


def conf_id_to_iterator(
    molecule: "Mol", conf_id: None | int | Sequence
) -> Generator | List["Mol"]:
    """Iterate over the conformers selected by `conf_id` (None selects all).

    Raises TypeError unless `conf_id` is None, an int, or a list.
    """
    if conf_id is None:
        conf_iterator = molecule.GetConformers()
    elif isinstance(conf_id, int):
        conf_iterator = [molecule.GetConformer(conf_id)]
    elif isinstance(conf_id, list):
        conf_iterator = conformer_iterator(molecule, conf_id)
    else:
        raise TypeError("conf_id must be None, int, or list of int")
    return conf_iterator


def get_atom_symbols_rdkit(molecule) -> list[str]:
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    return symbols


def get_energy_rdkit(conformer) -> float | None:
    """Float of the conformer's "energy" property, or None if unset."""
    if conformer.HasProp("energy"):
        energy = float(conformer.GetProp("energy"))
    else:
        energy = None
    return energy


@overload
def rdkit_to_molecule(
    molecules: "Mol", conf_id: int | Sequence[int] | None = None
) -> Molecule | list[Molecule]: ...
@overload
def rdkit_to_molecule(
    molecules: Sequence["Mol"], conf_id: int | Sequence[int] | None = None
) -> list[Molecule]: ...


def rdkit_to_molecule(
    molecules, conf_id: int | Sequence[int] | None = None
) -> Molecule | list[Molecule]:
    """Convert RDKit Mol(s) to irmsd Molecule(s), one per selected conformer.

    Mol-level and conformer-level properties are merged into ``info``;
    conformer properties win on key clashes.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
    conf_id : int, list of int, or None, optional
        Conformer ID(s); None selects all conformers.

    Returns
    -------
    Molecule or list of Molecule
        A single Molecule if exactly one conformer was converted, else a list.

    Raises
    ------
    TypeError
        If the input is not a Mol or list of Mols, or a conformer is not 3D.
    """

    require_rdkit()

    from rdkit import Chem

    if isinstance(molecules, Chem.Mol):
        molecules = [molecules]

    if not isinstance(molecules, Sequence):
        raise TypeError(
            "rdkit_to_molecule expects rdkit.Chem.Mol objects or a list of them"
        )

    all_mols = []
    for mol in molecules:
        if not isinstance(mol, Chem.Mol):
            raise TypeError("rdkit_to_molecule expects rdkit.Chem.Mol objects")

        mol_info = {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}
        conf_iterator = conf_id_to_iterator(mol, conf_id)

        for conformer in conf_iterator:
            if not conformer.Is3D():
                raise TypeError("rdkit_to_molecule expects 3D conformers")

            symbols = get_atom_symbols_rdkit(mol)
            pos = conformer.GetPositions()  # (N, 3)
            conf_info = {
                prop: conformer.GetProp(prop) for prop in conformer.GetPropNames()
            }
            energy = get_energy_rdkit(conformer)

            info = mol_info | conf_info

            new_mol = Molecule(symbols=symbols, positions=pos, info=info, energy=energy)
            all_mols.append(new_mol)

    if len(all_mols) == 1:
        return all_mols[0]
    else:
        return all_mols


@overload
def molecule_to_rdkit(molecule: Molecule) -> "Mol": ...
@overload
def molecule_to_rdkit(molecules: Sequence[Molecule]) -> list["Mol"]: ...


def molecule_to_rdkit(molecule: Molecule | Sequence[Molecule]) -> "Mol" | list["Mol"]:
    """Convert irmsd Molecule(s) to RDKit Mol(s) carrying atoms and one conformer.

    No bonds or properties are transferred. Returns a single Mol if exactly
    one Molecule results, else a list.

    Raises
    ------
    TypeError
        If the input is not a Molecule or a list of them.
    """

    require_rdkit()

    from rdkit import Chem
    from rdkit.Chem import AllChem

    if isinstance(molecule, Molecule):
        molecule = [molecule]

    if not isinstance(molecule, Sequence):
        raise TypeError(
            "molecule_to_rdkit expects irmsd.core.Molecule objects or a list of them"
        )

    all_mols = []
    for mol in molecule:
        if not isinstance(mol, Molecule):
            raise TypeError("molecule_to_rdkit expects irmsd.core.Molecule objects")

        rdkit_mol = Chem.RWMol()
        atom_indices = []
        for symbol in mol.symbols:
            atom = Chem.Atom(symbol)
            atom_idx = rdkit_mol.AddAtom(atom)
            atom_indices.append(atom_idx)

        conformer = Chem.Conformer(len(mol.symbols))
        conformer.SetPositions(mol.positions)

        rdkit_mol.AddConformer(conformer, assignId=True)
        all_mols.append(rdkit_mol.GetMol())

    if len(all_mols) == 1:
        return all_mols[0]
    else:
        return all_mols


def get_cn_rdkit(molecule, conf_id: None | int | Sequence = None) -> np.ndarray:
    """Coordination numbers for the selected conformers.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
    conf_id : int, list of int, or None
        Conformer ID(s); None selects all conformers.

    Returns
    -------
    np.ndarray
        Shape (n_atoms,) for one conformer, (n_conf, n_atoms) for several.

    Raises
    ------
    TypeError
        If `molecule` is not an RDKit Mol.
    """

    require_rdkit()

    from rdkit import Chem

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    core_mol = rdkit_to_molecule(molecule, conf_id=conf_id)
    if isinstance(core_mol, Molecule):
        return core_mol.get_cn()
    else:
        all_cn = [mol.get_cn() for mol in core_mol]
        return np.array(all_cn)


def get_axis_rdkit(
    molecule, conf_id: None | int | Sequence = None
) -> (
    Tuple[np.ndarray, np.ndarray, np.ndarray]
    | List[Tuple[np.ndarray, np.ndarray, np.ndarray]]
):
    """Principal axes for the selected conformers.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
    conf_id : int, list of int, or None
        Conformer ID(s); None selects all conformers.

    Returns
    -------
    tuple of np.ndarray, or list of such tuples
        (rotational constants, average moments, eigenvectors) per conformer.

    Raises
    ------
    TypeError
        If `molecule` is not an RDKit Mol.
    """
    require_rdkit()

    from rdkit import Chem

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    core_mol = rdkit_to_molecule(molecule, conf_id=conf_id)
    if isinstance(core_mol, Molecule):
        return core_mol.get_axis()
    else:
        all_results = [mol.get_axis() for mol in core_mol]
        return all_results


def get_point_group_rdkit(
    molecule,
    conf_id: None | int | Sequence = None,
    **settings,
) -> str | None | List[str | None]:
    """Schoenflies point group for the selected conformers.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
    conf_id : int, list of int, or None
        Conformer ID(s); None selects all conformers.
    **settings
        Analyzer settings, see :func:`irmsd.get_point_group`.

    Returns
    -------
    str or None, or list of those
        Symbol per conformer; None if the analyzer skipped it.

    Raises
    ------
    TypeError
        If `molecule` is not an RDKit Mol.
    """
    require_rdkit()

    from rdkit import Chem

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("get_point_group_rdkit expects rdkit.Chem.Mol objects")

    core_mol = rdkit_to_molecule(molecule, conf_id=conf_id)
    if isinstance(core_mol, Molecule):
        return core_mol.get_point_group(**settings)
    return [mol.get_point_group(**settings) for mol in core_mol]


def get_canonical_rdkit(
    molecule,
    conf_id: None | int | Sequence = None,
    wbo: None | np.ndarray = None,
    invtype="apsp+",
    heavy: bool = False,
) -> np.ndarray:
    """Canonical atom ranks for the selected conformers.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
    conf_id : int, list of int, or None
        Conformer ID(s); None selects all conformers.
    wbo : np.ndarray, optional
        Wiberg bond orders, required for ``invtype='cangen'``. Shape
        (n_atoms, n_atoms) to share across conformers, or
        (n_conf, n_atoms, n_atoms) for one per conformer.
    invtype : str
        Invariant type.
    heavy : bool
        Rank heavy atoms only.

    Returns
    -------
    np.ndarray
        Shape (n_atoms,) for one conformer, (n_conf, n_atoms) for several.

    Raises
    ------
    TypeError
        If `molecule` is not an RDKit Mol.
    """

    require_rdkit()

    from rdkit import Chem

    if not isinstance(molecule, Chem.Mol):
        raise TypeError("rdkit_to_fortran_pair expects rdkit.Chem.Mol objects")

    core_mol = rdkit_to_molecule(molecule, conf_id=conf_id)
    if isinstance(core_mol, Molecule):
        return core_mol.get_canonical(invtype=invtype, wbo=wbo, heavy=heavy)
    else:
        if wbo is None:
            all_canonical = [
                mol.get_canonical(invtype=invtype, heavy=heavy) for mol in core_mol
            ]
        elif wbo.ndim == 2:
            all_canonical = [
                mol.get_canonical(invtype=invtype, wbo=wbo, heavy=heavy)
                for mol in core_mol
            ]
        else:
            all_canonical = [
                mol.get_canonical(invtype=invtype, wbo=wbo_mol, heavy=heavy)
                for mol, wbo_mol in zip(core_mol, wbo)
            ]
        return np.array(all_canonical)


def get_rmsd_rdkit(
    molecule_ref, molecule_align, conf_id_ref=-1, conf_id_align=-1, mask=None
) -> Tuple[float, "Mol", np.ndarray]:
    """RMSD between two conformers after alignment, without permutation.

    Parameters
    ----------
    molecule_ref, molecule_align : rdkit.Chem.Mol
    conf_id_ref, conf_id_align : int
        Conformer IDs; -1 is the RDKit default conformer.
    mask : array-like of bool, optional
        Atoms to include in the RMSD.

    Returns
    -------
    rmsd : float
        In Angstrom.
    aligned : rdkit.Chem.Mol
    rotmat : np.ndarray
        Rotation matrix.

    Raises
    ------
    TypeError
        If either input is not an RDKit Mol.
    """

    require_rdkit()

    from rdkit import Chem

    if not isinstance(molecule_ref, Chem.Mol) or not isinstance(
        molecule_align, Chem.Mol
    ):
        raise TypeError("get_rmsd_rdkit expects rdkit.Chem.Mol objects")

    molecule_ref_core = rdkit_to_molecule(molecule_ref, conf_id=conf_id_ref)
    molecule_align_core = rdkit_to_molecule(molecule_align, conf_id=conf_id_align)

    rmsd, molecule_new_core, rotmat = get_rmsd_molecule(
        molecule_ref_core, molecule_align_core, mask=mask
    )

    molecule_ret = molecule_to_rdkit(molecule_new_core)
    return rmsd, molecule_ret, rotmat


def get_irmsd_rdkit(
    molecule_ref, molecule_align, conf_id_ref=-1, conf_id_align=-1, iinversion: int = 0
) -> Tuple[float, "Mol", "Mol"]:
    """iRMSD between two conformers after permutation and alignment.

    Parameters
    ----------
    molecule_ref, molecule_align : rdkit.Chem.Mol
    conf_id_ref, conf_id_align : int
        Conformer IDs; -1 is the RDKit default conformer.
    iinversion : int
        Inversion handling: 0 auto, 1 on, 2 off.

    Returns
    -------
    irmsd : float
        In Angstrom.
    ref, aligned : rdkit.Chem.Mol
        Processed reference and aligned molecules.

    Raises
    ------
    TypeError
        If either input is not an RDKit Mol.
    """
    require_rdkit()

    from rdkit import Chem

    if not isinstance(molecule_ref, Chem.Mol) or not isinstance(
        molecule_align, Chem.Mol
    ):
        raise TypeError("get_rmsd_rdkit expects rdkit.Chem.Mol objects")

    molecule_ref_core = rdkit_to_molecule(molecule_ref, conf_id=conf_id_ref)
    molecule_align_core = rdkit_to_molecule(molecule_align, conf_id=conf_id_align)

    irmsd, molecule_1_core, molecule_2_core = get_irmsd_molecule(
        molecule_ref_core, molecule_align_core, iinversion
    )
    molecule_1_ret = molecule_to_rdkit(molecule_1_core)
    molecule_2_ret = molecule_to_rdkit(molecule_2_core)
    return irmsd, molecule_1_ret, molecule_2_ret


def sorter_irmsd_rdkit(
    molecules: "Mol" | Sequence["Mol"],
    rthr: float = 0.125,  # aligned with '--rthr' in src/irmsd/cli.py
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
) -> Tuple[np.ndarray, List["Mol"]]:
    """Group an ensemble by iRMSD.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        A single Mol must carry more than one conformer.
    rthr : float
        iRMSD threshold in Angstrom.
    iinversion : int
        Inversion handling: 0 auto, 1 on, 2 off.
    allcanon, printlvl
        Canonicalization flag and verbosity, passed to the backend.
    ethr : float or None
        Energy pre-sorting threshold.
    ewin : float or None
        Energy window in Hartree above the lowest-energy structure.

    Returns
    -------
    groups : np.ndarray
        Integer group index per structure.
    new_molecules_list : list of rdkit.Chem.Mol
        Sorted structures.

    Raises
    ------
    TypeError
        If the input is not an RDKit Mol or a list of them.
    """
    require_rdkit()

    from rdkit import Chem

    if isinstance(molecules, Chem.Mol):
        assert (
            molecules.GetNumConformers() > 1
        ), "Molecule must have multiple conformers"
    else:
        for mol in molecules:
            if not isinstance(mol, Chem.Mol):
                raise TypeError("sorter_irmsd_rdkit expects rdkit.Chem.Mol objects")

    mols = rdkit_to_molecule(molecules)

    groups, new_mols = sorter_irmsd_molecule(
        molecule_list=mols,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
        ewin=ewin,
    )

    new_molecules_list = molecule_to_rdkit(new_mols)

    return groups, new_molecules_list


def delta_irmsd_list_rdkit(
    molecules: "Mol" | Sequence["Mol"],
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
) -> Tuple[np.ndarray, List["Mol"]]:
    """iRMSD deltas over an ensemble.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        A single Mol must carry more than one conformer.
    iinversion : int
        Inversion handling: 0 auto, 1 on, 2 off.
    allcanon, printlvl
        Canonicalization flag and verbosity, passed to the backend.

    Returns
    -------
    delta : np.ndarray
        See ``delta_irmsd_list`` for semantics.
    new_molecules_list : list of rdkit.Chem.Mol
        Processed structures.

    Raises
    ------
    TypeError
        If the input is not an RDKit Mol or a list of them.
    """
    require_rdkit()

    from rdkit import Chem

    if isinstance(molecules, Chem.Mol):
        assert (
            molecules.GetNumConformers() > 1
        ), "Molecule must have multiple conformers"
    else:
        for mol in molecules:
            if not isinstance(mol, Chem.Mol):
                raise TypeError("sorter_irmsd_rdkit expects rdkit.Chem.Mol objects")
    mols = rdkit_to_molecule(molecules)

    delta, new_mols = delta_irmsd_list_molecule(
        molecule_list=mols,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
    )

    new_molecules_list = molecule_to_rdkit(new_mols)

    return delta, new_molecules_list


def cregen_rdkit(
    molecules: "Mol" | Sequence["Mol"],
    rthr: float = 0.125,
    ethr: float = 8.0e-5,
    bthr: float = 0.01,
    printlvl: int = 0,
    ewin: float | None = None,
) -> List["Mol"]:
    """CREGEN-style pruning: unique structures by iRMSD, energy and rotational constants.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        A single Mol must carry more than one conformer.
    rthr : float
        iRMSD threshold in Angstrom.
    ethr : float
        Energy threshold in Hartree.
    bthr : float
        Relative rotational-constant threshold.
    printlvl : int
        Verbosity, passed to the backend.
    ewin : float or None
        Energy window in Hartree above the lowest-energy structure.

    Returns
    -------
    list of rdkit.Chem.Mol
        Unique structures.

    Raises
    ------
    TypeError
        If the input is not an RDKit Mol or a list of them.
    """
    require_rdkit()

    from rdkit import Chem

    if isinstance(molecules, Chem.Mol):
        assert (
            molecules.GetNumConformers() > 1
        ), "Molecule must have multiple conformers"
    else:
        for mol in molecules:
            if not isinstance(mol, Chem.Mol):
                raise TypeError("sorter_irmsd_rdkit expects rdkit.Chem.Mol objects")

    mols = rdkit_to_molecule(molecules)

    new_mols = cregen(
        molecule_list=mols,
        rthr=rthr,
        ethr=ethr,
        bthr=bthr,
        printlvl=printlvl,
        ewin=ewin,
    )

    new_molecules_list = molecule_to_rdkit(new_mols)

    return new_molecules_list


def prune_rdkit(
    molecules: "Mol" | Sequence["Mol"],
    rthr: float,
    iinversion: int = 0,
    allcanon: bool = True,
    printlvl: int = 0,
    ethr: float | None = None,
    ewin: float | None = None,
) -> List["Mol"]:
    """Prune an ensemble by iRMSD, keeping the first structure of each group.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        A single Mol must carry more than one conformer.
    rthr : float
        iRMSD threshold in Angstrom.
    iinversion : int
        Inversion handling: 0 auto, 1 on, 2 off.
    allcanon, printlvl
        Canonicalization flag and verbosity, passed to the backend.
    ethr : float or None
        Energy pre-sorting threshold.
    ewin : float or None
        Energy window in Hartree above the lowest-energy structure.

    Returns
    -------
    list of rdkit.Chem.Mol
        Retained structures.

    Raises
    ------
    TypeError
        If the input is not an RDKit Mol or a list of them.
    """
    require_rdkit()

    from rdkit import Chem

    if isinstance(molecules, Chem.Mol):
        assert (
            molecules.GetNumConformers() > 1
        ), "Molecule must have multiple conformers"
    else:
        for mol in molecules:
            if not isinstance(mol, Chem.Mol):
                raise TypeError("sorter_irmsd_rdkit expects rdkit.Chem.Mol objects")

    mols = rdkit_to_molecule(molecules)

    new_mols = prune(
        molecule_list=mols,
        rthr=rthr,
        iinversion=iinversion,
        allcanon=allcanon,
        printlvl=printlvl,
        ewin=ewin,
    )

    new_molecules_list = molecule_to_rdkit(new_mols)

    return new_molecules_list
