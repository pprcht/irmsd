import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

import irmsd

# --- CLI ---------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="irmsd",
        description=(
            "CLI to read an arbitrary number of structures with ASE and run "
            "selected analysis commands on them."
        ),
    )
    p.add_argument(
        "structures",
        nargs="+",
        help="Paths to structure files (e.g. .xyz, .pdb, .cif). You can pass many.",
    )
    # Commands (flags). Multiple can be combined; they run in the order defined here.
    p.add_argument(
        "--cn",
        action="store_true",
        help=(
            "Calculate coordination numbers for each structure and print them as numpy arrays. "
        ),
    )
    p.add_argument(
        "--rot",
        action="store_true",
        help=("Calculate the rotational constants. "),
    )
    p.add_argument(
        "--canonical",
        action="store_true",
        help=("Calculate the canonical identifiers. "),
    )
    p.add_argument(
        "--rmsd",
        action="store_true",
        help=(
            "Calculate the Cartesian RMSD between two given structures via a quaternion algorithm. "
        ),
    )
    p.add_argument(
        "--irmsd",
        action="store_true",
        help=("Calculate the invariant Cartesian RMSD between two given structures. "),
    )
    p.add_argument(
        "--heavy",
        action="store_true",
        help=(
            "When calculating RMSD or canonical atom identifier, consider only heavy atoms. "
        ),
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file name (optional). If not provided, nothing is written.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    atoms_list = irmsd.read_structures(args.structures)

    ran_any = False

    heavy = args.heavy

    if args.cn:
        irmsd.compute_cn_and_print(atoms_list)
        ran_any = True

    if args.rot:
        irmsd.compute_axis_and_print(atoms_list)
        ran_any = True

    if args.canonical:
        irmsd.compute_canonical_and_print(atoms_list, heavy=heavy)
        ran_any = True

    if args.rmsd:
        irmsd.compute_quaternion_rmsd_and_print(
            atoms_list, heavy=heavy, outfile=args.output
        )
        ran_any = True

    if args.irmsd:
        irmsd.compute_irmsd_and_print(atoms_list, outfile=args.output)
        ran_any = True

    if not ran_any:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
