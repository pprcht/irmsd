import argparse
import sys
import numpy as np
from ase.io import read as ase_read, write as ase_write
import irmsd  # your package exposing api.ase_to_fortran


def _print_array(title: str, arr: np.ndarray) -> None:
    print(title)
    with np.printoptions(precision=6, suppress=True):
        print(arr)
    print()  # blank line


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="ase-bridge",
        description="Read a structure with ASE, call the ASE→Fortran test routine, and optionally save the result.",
    )
    p.add_argument(
        "structure",
        help="Path to structure file (e.g. .xyz, .pdb, .cif)",
    )
    p.add_argument(
        "--format",
        default=None,
        help="ASE format hint (optional, e.g. 'xyz'). If omitted, ASE guesses from extension.",
    )
    p.add_argument(
        "--apply",
        action="store_true",
        help="Execute the Fortran routine on the structure and print the results.",
    )
    p.add_argument(
        "--output",
        metavar="FILE",
        help="Write the modified structure to this file. Format is inferred from extension.",
    )

    args = p.parse_args(argv)

    # --- read structure ---
    atoms = ase_read(args.structure, format=args.format)
    pos = atoms.get_positions()

    _print_array("INPUT COORDINATES (Å):", pos)

    # --- optional Fortran call ---
    if args.apply:
        new_atoms, mat = irmsd.ase_to_fortran(atoms)

        _print_array("NEW COORDINATES (Å):", new_atoms.get_positions())
        _print_array("OUTPUT 3×3 MATRIX:", mat)

        if args.output:
            ase_write(args.output, new_atoms)
            print(f"✅ Modified structure written to: {args.output}")

    elif args.output:
        print("⚠️  The --output option has no effect without --apply.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())

