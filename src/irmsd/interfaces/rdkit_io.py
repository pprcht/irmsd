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
    """Generator that yields conformers from a molecule given a list of
    conformer IDs.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        RDKit Molecule object containing conformers.
    conf_ids : list of int
        List of conformer IDs to yield.

    Yields
    ------
    rdkit.Chem.Conformer
        The conformer corresponding to each ID in conf_ids.
    """
    for conf_id in conf_ids:
        yield molecule.GetConformer(conf_id)


def conf_id_to_iterator(
    molecule: "Mol", conf_id: None | int | Sequence
) -> Generator | List["Mol"]:
    """Convert conf_id input to an iterator over conformers.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        RDKit Molecule object containing conformers.
    conf_id : None, int, or list of int
        Conformer ID(s) to select. If None, all conformers are selected.

    Returns
    -------
    generator or list of rdkit.Chem.Conformer
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
    """Get atomic symbols from an RDKit Molecule object."""
    symbols = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    return symbols


def get_energy_rdkit(conformer) -> float | None:
    """Get energy from an RDKit Conformer object, if available."""
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
    """Convert one or more RDKit Molecule objects to one or more irmsd Molecule
    objects.

    If conf_id is None, all conformers are converted. If conf_id is an
    int, only that conformer is converted. If conf_id is a list of int,
    only those conformers are converted.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        RDKit Molecule object(s) to convert.
    conf_id : int, list of int, or None, optional
        Conformer ID(s) to convert. If None, all conformers are converted.

    Returns
    -------
    irmsd.core.Molecule or list of irmsd.core.Molecule
        Converted irmsd Molecule object(s).

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule or a list of them. Also if any
        conformer is not 3D.
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

            symbols = get_atom_symbols_rdkit(mol)  # list of str
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
    """Convert one or more irmsd Molecule objects to one or more RDKit Molecule
    objects.

    Parameters
    ----------
    molecule : irmsd.core.Molecule or list of irmsd.core.Molecule
        irmsd Molecule object(s) to convert.

    Returns
    -------
    rdkit.Chem.Mol or list of rdkit.Chem.Mol
        Converted RDKit Molecule object(s).

    Raises
    ------
    TypeError
        If the input is not an irmsd Molecule or a list of them.
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
    """Optional RDKit utility: compute coordination numbers for one or more conformers of a molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        RDKit Molecule object containing conformers.
    conf_id : int, list of int, or None, optional
        Conformer ID(s) to compute coordination numbers for. If None, all conformers are used.

    Returns
    -------
    np.ndarray
        Coordination numbers for each atom in the specified conformers. If multiple
        conformers are specified, returns an array of shape (n_conf, n_atoms).

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule.
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
    """Optional RDKit utility: compute principal axes for one or more conformers of a molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        RDKit Molecule object containing conformers.
    conf_id : int, list of int, or None, optional
        Conformer ID(s) to compute principal axes for. If None, all conformers are used.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray] or list of such tuples.
        (Rotational constants, average moments, eigenvectors) for each specified conformer.

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule.
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


def get_canonical_rdkit(
    molecule,
    conf_id: None | int | Sequence = None,
    wbo: None | np.ndarray = None,
    invtype="apsp+",
    heavy: bool = False,
) -> np.ndarray:
    """Optional RDKit utility: compute coordination numbers for one or more conformers of a molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        RDKit Molecule object containing conformers.
    conf_id : int, list of int, or None, optional
        Conformer ID(s) to compute canonical representations for. If None, all conformers are used.
    wbo : np.ndarray, optional
        Optional weight bond order matrix/matrices for canonicalization, required for the 'cangen' invtype. If given either one per conformer with shape (n_conf, n_atoms, n_atoms) or use the same for all (n_atoms, n_atoms).
    invtype : str, optional
        Type of invariant representation to compute. Default is 'apsp+'.
    heavy : bool, optional
        Whether to consider only heavy atoms in the canonicalization. Default is False.

    Returns
    -------
    np.ndarray
        Canonical ranks for each atom in the specified conformers. If multiple
        conformers are specified, returns an array of shape (n_conf, n_atoms).

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule.
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
    """Optional Rdkit utility: operate on two Rdkit Molecules. Returns the RMSD in Angström,
    the molecule object with both Conformers aligned.

    Parameters
    ----------
    molecule_ref : rdkit.Chem.Mol
        Reference RDKit Molecule object.
    molecule_align : rdkit.Chem.Mol
        RDKit Molecule object to be aligned.
    conf_id_ref : int, optional
        Conformer ID for the reference molecule. Default is -1 (rdkit default).
    conf_id_align : int, optional
        Conformer ID for the molecule to be aligned. Default is -1 (rdkit default).
    mask : array-like of bool, optional

    Returns
    -------
    Tuple[float, rdkit.Chem.Mol, np.ndarray]
        RMSD in Angström, aligned RDKit Molecule object, and rotation matrix.

    Raises
    ------
    TypeError
        If the inputs are not RDKit Molecule objects.
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
    """
    Optional Rdkit utility: operate on TWO Rdkit Molecules. Returns the iRMSD in Angström,
    the molecule object with both Conformers permuted and aligned.

    Parameters
    ----------
    molecule_ref : rdkit.Chem.Mol
        Reference RDKit Molecule object.
    molecule_align : rdkit.Chem.Mol
        RDKit Molecule object to be aligned.
    conf_id_ref : int, optional
        Conformer ID for the reference molecule. Default is -1 (rdkit default).
    conf_id_align : int, optional
        Conformer ID for the molecule to be aligned. Default is -1 (rdkit default).
    iinversion : int, optional
        Inversion type for iRMSD calculation. Default is 0. ( 0: 'auto', 1: 'on', 2: 'off' )

    Returns
    -------
    Tuple[float, rdkit.Chem.Mol, rdkit.Chem.Mol]
        iRMSD in Angström, aligned RDKit Molecule object for reference, aligned RDKit Molecule object for alignment.

    Raises
    ------
    TypeError
        If the inputs are not RDKit Molecule objects.
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
    """
    Optional Rdkit utility: operate on a list of Rdkit Molecules.
    Returns a list of indices corresponding to the sorted molecules based on iRMSD.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        RDKit Molecule object(s) containing multiple conformers.
    rthr : float
        iRMSD threshold for grouping.
    iinversion : int, optional
        Inversion type for iRMSD calculation. Default is 0. ( 0: 'auto', 1: 'on', 2: 'off' )
    allcanon : bool, optional
        Canonicalization flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.
    ethr : float | None
        Optional energy threshold to accelerate by pre-sorting
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    groups : np.ndarray
        Integer array of shape (nat,) with group indices for the first
        ``nat`` atoms (as defined by the backend).
    new_molecules_list : list of rdkit.Chem.Mol
        List of RDKit Molecule objects corresponding to the sorted molecules.

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule or a list of them.
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
    """
    Optional Rdkit utility: operate on a list of Rdkit Molecules.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        RDKit Molecule object(s) containing multiple conformers.
    iinversion : int, optional
        Inversion type for iRMSD calculation. Default is 0. ( 0: 'auto', 1: 'on', 2: 'off' )
    allcanon : bool, optional
        Canonicalization flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.

    Returns
    -------
    delta : np.ndarray
        Float array returned by the backend (see ``delta_irmsd_list`` for
        detailed semantics).
    new_molecules_list : list of rdkit.Chem.Mol
        List of RDKit Molecule objects corresponding to the processed molecules.

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule or a list of them.
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
    """
    Optional Rdkit utility: operate on a list of Rdkit Molecules.
    Returns a pruned list of molecules based on iRMSD.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        RDKit Molecule object(s) containing multiple conformers.
    rthr : float
        iRMSD threshold for grouping.
    ethr : float                                      
        Energy threshold to accelerate by pre-sorting. In Hartree.
    bthr: float
        Rotational constant comparison threshold. Relative value (default: 0.01)
    printlvl : int, optional
        Verbosity level, passed through to the backend.
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    groups : np.ndarray
        Integer array of shape (nat,) with group indices for the first
        ``nat`` atoms (as defined by the backend).
    new_molecules_list : list of rdkit.Chem.Mol
        List of RDKit Molecule objects corresponding to the sorted molecules.

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule or a list of them.
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
    """
    Optional Rdkit utility: operate on a list of Rdkit Molecules.
    Returns a pruned list of molecules based on iRMSD.

    Parameters
    ----------
    molecules : rdkit.Chem.Mol or list of rdkit.Chem.Mol
        RDKit Molecule object(s) containing multiple conformers.
    rthr : float
        iRMSD threshold for grouping.
    iinversion : int, optional
        Inversion type for iRMSD calculation. Default is 0. ( 0: 'auto', 1: 'on', 2: 'off' )
    allcanon : bool, optional
        Canonicalization flag, passed through to the backend.
    printlvl : int, optional
        Verbosity level, passed through to the backend.
    ethr : float | None
        Optional energy threshold to accelerate by pre-sorting
    ewin : float | None
        Optional energy window to limit ensembe size around lowest energy structure.
        In Hartree.

    Returns
    -------
    new_molecules_list : list of rdkit.Chem.Mol
        List of RDKit Molecule objects corresponding to the sorted molecules.

    Raises
    ------
    TypeError
        If the input is not an RDKit Molecule or a list of them.
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
