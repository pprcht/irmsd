import argparse
import sys
import numpy as np
from ase.io import read as ase_read, write as ase_write
import irmsd  # exposes ase_to_fortran and ase_to_fortran_pair

def _print_array(title: str, arr: np.ndarray) -> None:
    print(title)
    with np.printoptions(precision=6, suppress=True):
        print(arr)
    print()  # blank line

def main(argv=None):
    p = argparse.ArgumentParser(
        prog="irmsd",
        description="Read one or two structures with ASE, call the ASE→Fortran test routine, and optionally save the result.",
    )
    p.add_argument(
        "structure",
        help="Path to first structure file (e.g. .xyz, .pdb, .cif)",
    )
    p.add_argument(
        "structure2",
        nargs="?",
        help="Optional second structure file. If provided with --apply, the pair routine is used.",
    )
    p.add_argument(
        "--format",
        default=None,
        help="ASE format hint (optional, e.g. 'xyz'). If omitted, ASE guesses from extension.",
    )
    p.add_argument(
        "--apply",
        action="store_true",
        help="Execute the Fortran routine and print the results.",
    )
    p.add_argument(
        "--output",
        metavar="FILE",
        help="Write the modified structure to this file. Format inferred from extension. "
             "Single-structure: writes modified first structure. Two-structure: writes modified second structure.",
    )

    args = p.parse_args(argv)

    # --- read structure(s) ---
    atoms1 = ase_read(args.structure, format=args.format)
    pos1 = atoms1.get_positions()
    _print_array("INPUT #1 COORDINATES (Å):", pos1)

    atoms2 = None
    if args.structure2:
        atoms2 = ase_read(args.structure2, format=args.format)
        pos2 = atoms2.get_positions()
        _print_array("INPUT #2 COORDINATES (Å):", pos2)

    # --- optional Fortran call ---
    if args.apply:
        if atoms2 is None:
            # single-structure mode
            new_atoms1, mat1 = irmsd.ase_to_fortran(atoms1)
            _print_array("NEW #1 COORDINATES (Å):", new_atoms1.get_positions())
            _print_array("OUTPUT 3×3 MATRIX #1:", mat1)

            if args.output:
                ase_write(args.output, new_atoms1)
                print(f"✅ Modified structure written to: {args.output}")
        else:
            # two-structure mode → return only modified second + both matrices
            new_atoms2, mat1, mat2 = irmsd.ase_to_fortran_pair(atoms1, atoms2)
            _print_array("OUTPUT 3×3 MATRIX #1:", mat1)
            _print_array("OUTPUT 3×3 MATRIX #2:", mat2)
            _print_array("NEW #2 COORDINATES (Å):", new_atoms2.get_positions())

            if args.output:
                ase_write(args.output, new_atoms2)
                print(f"✅ Modified structure (second) written to: {args.output}")

    elif args.output:
        print("⚠️  The --output option has no effect without --apply.", file=sys.stderr)

    return 0

if __name__ == "__main__":
    sys.exit(main())
    
