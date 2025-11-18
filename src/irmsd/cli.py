import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

import irmsd

# --- CLI ---------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    structures_parser = argparse.ArgumentParser(add_help=False)
    structures_parser.add_argument(
        "structures",
        nargs="+",
        help="Paths to structure files (e.g. .xyz, .pdb, .cif). You can pass many.",
    )

    heavy_parser = argparse.ArgumentParser(add_help=False)
    heavy_parser.add_argument(
        "--heavy",
        action="store_true",
        help="When calculating RMSD or canonical atom identifier, consider only heavy atoms. ",
    )
    outfile_parser = argparse.ArgumentParser(add_help=False)
    outfile_parser.add_argument(
        "-o",
        "--output",
        dest="outfile",
        type=Path,
        default=None,
        help="Output file name (optional). If not provided, nothing is written.",
    )

    main_parser = argparse.ArgumentParser(
        prog="irmsd",
        description=(
            "CLI to read an arbitrary number of structures with ASE and run "
            "selected analysis commands on them."
        ),
    )

    subparsers = main_parser.add_subparsers()
    rmsd_parser = subparsers.add_parser(
        "rmsd",
        help="Calculate RMSD between structures.",
        parents=[structures_parser, outfile_parser, heavy_parser],
    )
    rmsd_parser.set_defaults(func=irmsd.compute_quaternion_rmsd_and_print)
    irmsd_parser = subparsers.add_parser(
        "irmsd",
        help="Calculate invariant RMSD between structures.",
        parents=[structures_parser, outfile_parser],
    )
    irmsd_parser.set_defaults(func=irmsd.compute_irmsd_and_print)
    cn_parser = subparsers.add_parser(
        "cn",
        help="Calculate coordination numbers for structures.",
        parents=[structures_parser],
    )
    cn_parser.set_defaults(func=irmsd.compute_cn_and_print)
    canonical_parser = subparsers.add_parser(
        "canonical",
        help="Calculate canonical identifiers for structures.",
        parents=[structures_parser, heavy_parser],
    )
    canonical_parser.set_defaults(func=irmsd.compute_canonical_and_print)
    rot_parser = subparsers.add_parser(
        "rot",
        help="Calculate rotational constants for structures.",
        parents=[structures_parser],
    )
    rot_parser.set_defaults(func=irmsd.compute_axis_and_print)

    irmsd_parser.add_argument(
        "--inversion",
        choices=["on", "off", "auto"],
        default="auto",
        help="Control coordinate inversion in irmsd runtypes: 'on', 'off', or 'auto' (default: auto).",
    )
    return main_parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    def wrapper_func(args):
        atoms_list = irmsd.read_structures(args.structures)
        tmp_args = vars(args).copy()
        tmp_args = {
            key: val
            for key, val in tmp_args.items()
            if key not in ["structures", "func"]
        }
        return args.func(
            atoms_list,
            **tmp_args,
        )

    wrapper_func(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
