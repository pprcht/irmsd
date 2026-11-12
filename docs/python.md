# Python interface

The package is organized in three layers:

- a {py:class}`irmsd.Molecule` class with readers and writers for (extended) XYZ files,
- high-level functions that operate on lists of `Molecule` objects, and on ASE `Atoms` or RDKit `Mol` objects,
- thin wrappers around the Fortran routines that take NumPy arrays (`get_irmsd`, `sorter_irmsd`, `get_cn_fortran`, ...).

The full list of functions is given in the [package reference](package_reference.rst).
Coordinates are in Ångström and energies in Hartree throughout.

## Reading and writing structures

{py:func}`irmsd.read_structures` takes one or several file names and returns a flat list of `Molecule` objects, one per frame:

```python
import irmsd

molecules = irmsd.read_structures("ensemble.xyz")
molecules = irmsd.read_structures(["run1.xyz", "run2.xyz"])

mol = molecules[0]
mol.get_chemical_formula()       # 'C5H12'
mol.get_positions()              # (natoms, 3) array, Å
mol.get_potential_energy()       # Hartree, from the comment line
mol.get_point_group()            # e.g. 'C2v'

irmsd.write_structures("out.xyz", molecules)
```

The input formats and the parsing of energies and canonical atom identifiers are the same as for the [command-line tool](cli.md#input-structures).

## Comparing two structures

```python
rmsd, mol1, mol2 = irmsd.get_irmsd_molecule(molecules[0], molecules[1])
```

{py:func}`irmsd.get_irmsd_molecule` returns the iRMSD together with copies of both structures, aligned and with the atom order of the second one matched to the first.
The quaternion RMSD without permutation is available as {py:func}`irmsd.get_rmsd_molecule`.

The `iinversion` argument controls the false-enantiomer check: `0` (auto, default), `1` (on), `2` (off).

## Sorting and pruning ensembles

```python
# group index for each structure, and the processed structures
groups, processed = irmsd.sorter_irmsd_molecule(molecules, rthr=0.125)

# only the unique structures
unique = irmsd.prune(molecules, rthr=0.125)

# with an energy threshold and an energy window, both in Hartree
unique = irmsd.prune(molecules, rthr=0.125, ethr=8.0e-5, ewin=0.01)
```

All structures passed to one call must have the same number of atoms.
With `ethr` set, structures whose energies differ by more than `ethr` are treated as different without computing their iRMSD, which saves time on large ensembles.
`ewin` removes all structures more than `ewin` above the lowest one before sorting.

{py:func}`irmsd.cregen` is the cheaper CREGEN-style sorting via energies, quaternion RMSD, and rotational constants.
It compares the structures in their given atom order and therefore requires the same atom order in all of them.

The iRMSD comparisons run in parallel with OpenMP (see [](installation.md#parallel-sorting-openmp)), and {py:func}`irmsd.openmp_info` reports the thread count.

## ASE and RDKit

With ASE or RDKit installed, the same functionality is available directly on their structure types:

```python
from ase.io import read
import irmsd

frames = read("ensemble.xyz", index=":")
unique = irmsd.prune_ase(frames, rthr=0.125)          # list of ase.Atoms

rmsd, a1, a2 = irmsd.get_irmsd_ase(frames[0], frames[1])
```

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import irmsd

mol = Chem.AddHs(Chem.MolFromSmiles("CCCCC"))
AllChem.EmbedMultipleConfs(mol, numConfs=20, randomSeed=42)
AllChem.MMFFOptimizeMoleculeConfs(mol)

unique = irmsd.prune_rdkit(mol, rthr=0.125)           # 20 conformers -> 4 rdkit Mol
```

{py:func}`irmsd.ase_to_molecule`, {py:func}`irmsd.molecule_to_ase`, {py:func}`irmsd.rdkit_to_molecule` and {py:func}`irmsd.molecule_to_rdkit` convert between the types.
ASE energies are converted from eV to Hartree, unless `atoms.info["energy_units"]` is set to `"Hartree"`.

## Canonical atom identifiers

The canonical identifiers that underlie the iRMSD can be computed and stored with the structure:

```python
mol.set_ids(irmsd.get_canonical_fortran(mol.get_atomic_numbers(), mol.get_positions()))
irmsd.write_structures("with_ids.extxyz", mol)   # adds a canonical_id column
```

When identifiers are present on the input (`Molecule.ids`), `get_irmsd_molecule`, `sorter_irmsd_molecule` and `prune` pass them to the backend, which uses them if they are consistent between the structures and recomputes them otherwise.
`use_ids=False` ignores them.
