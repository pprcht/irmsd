# Command Line Interface

The iRSMD package comes with an CLI tool `irmsd`. This tool allows you to read multiple structures (e.g. from an extended `xyz`-format file), an perform operations on them. There are three subcommands that can be chosen: `prop`, `compare`, and `sort`:
```text
irmsd {prop,compare,sort} ...

positional arguments:   
{prop,compare,sort}  Subcommand to run.     
prop               Compute structural properties (CN, rotational constants, canonical IDs).     
compare            Compare structures via iRMSD (default) or quaternion RMSD.     
sort               Sort or cluster structures based on inter-structure RMSD.
```

The `sort` functionality exists as an utility function and can be used to determine some atomic properties for each structure provided:
```text
irmsd prop [-h] [--cn] [--rot] [--canonical] [--heavy] structures [structures ...]


positional arguments:
  structures   Paths to structure files (e.g. .xyz, .pdb, .cif).

options:
  -h, --help   show this help message and exit
  --cn         Calculate coordination numbers for each structure and print them as numpy arrays.
  --rot        Calculate the rotational constants.
  --canonical  Calculate the canonical identifiers.
  --heavy      When calculating canonical atom identifiers, consider only heavy atoms.
```

The `compare` subcommand performs a quaternion RMSD or an iRMSD comparison of the *first two* structures provided. It will return the alisgned structures.
```text
usage: irmsd compare [-h] [--quaternion] [--inversion {on,off,auto}] [--heavy] [-o OUTPUT] structures [structures ...]

positional arguments:
  structures            Paths to structure files (e.g. .xyz, .pdb, .cif).

options:
  -h, --help            show this help message and exit
  --quaternion          Use the quaternion-based Cartesian RMSD instead of the invariant RMSD.
  --inversion {on,off,auto}
                        Control coordinate inversion in iRMSD runtypes: 'on', 'off', or 'auto' (default: auto). Used only for iRMSD.
  --heavy               When comparing structures, consider only heavy atoms.
  -o OUTPUT, --output OUTPUT
                        Output file name (optional). If not provided, results are only printed.
```

Finally, the `sort` runtype performs ensemble pruning to remove redundant structures from a given structure list. It also splits the structure list into chemically distinct ensembles, should different molecules be included (currently decided via the sum formula).
```text
usage: irmsd sort [-h] [--rthr RTHR] [--inversion {on,off,auto}] [--align] [--heavy] [--maxprint MAXPRINT] [-o OUTPUT]
                  structures [structures ...]

positional arguments:
  structures            Paths to structure files (e.g. .xyz, .pdb, .cif).

options:
  -h, --help            show this help message and exit
  --rthr RTHR           Inter-structure RMSD threshold for sorting in Angström. Structures closer than this threshold are treated as
                        similar.
  --inversion {on,off,auto}
                        Control coordinate inversion when evaluating RMSDs during sorting: 'on', 'off', or 'auto' (default: auto).
  --align               Just sort by energy and align.
  --maxprint MAXPRINT   Printout option; determine how man rows are printed for each sorted ensemble.
  -o OUTPUT, --output OUTPUT
                        Optional output file for sorted / clustered results.
```
