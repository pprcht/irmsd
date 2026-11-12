# Command-line tool

Installing the package places an `irmsd` executable on the `PATH`.
It reads an arbitrary number of structures and runs one of three subcommands on them:

| Subcommand | Purpose |
|---|---|
| `irmsd prop` | Structural properties: coordination numbers, rotational constants, point group, canonical atom identifiers |
| `irmsd compare` | iRMSD (default) or quaternion RMSD between two structures, with the aligned structures |
| `irmsd sort` (alias `prune`) | Pruning of an ensemble to its unique structures |

## Input structures

`.xyz` and `.extxyz` files are read by the package itself, and every frame of a multi-frame file is used.
All other formats (e.g. `.pdb`, `.cif`, `POSCAR`) are read via ASE, which requires `pip install "irmsd[ase]"`.

Energies are taken from the comment line of each frame and are used by `irmsd sort`:

- an extended-XYZ `energy=` key, in Hartree unless the line also carries `energy_units=eV`,
- otherwise the first bare number on the comment line, in Hartree (the format written by CREST and xtb).

An extended-XYZ file may carry a per-atom `canonical_id` column (written e.g. by `irmsd prop --canonical -o out.extxyz`).
`compare` and `sort` use these identifiers when they are consistent between the structures and recompute them otherwise; `--recompute-ids` discards them altogether.

## Examples

Properties of a single structure:

```bash
irmsd prop --all molecule.xyz
irmsd prop --canonical --heavy -o ids.extxyz molecule.xyz
```

Compare two structures; the iRMSD restores the atom order, the quaternion RMSD does not:

```bash
irmsd compare struc1.xyz struc2.xyz -o aligned.xyz
irmsd compare --quaternion struc1.xyz struc2.xyz
```

Prune a conformer ensemble and write one representative per unique structure:

```bash
irmsd sort ensemble.xyz -o unique.xyz
irmsd sort ensemble.xyz --rthr 0.2 --ethr          # looser RMSD threshold, energy presorting
irmsd sort ensemble.xyz --classic                  # CREGEN-style sorting (cheaper, more empirical)
```

If the input contains different molecules, `sort` first splits it into separate ensembles by sum formula.
The default RMSD threshold of 0.125 Å is discussed in [Example 4](examples.md#example-4-conformerrotamer-ensemble-pruning).
The iRMSD comparisons inside `sort` run in parallel with OpenMP; see [](installation.md#parallel-sorting-openmp).

## Reference

```{eval-rst}
.. argparse::
   :module: irmsd.cli
   :func: build_parser
   :prog: irmsd
```
