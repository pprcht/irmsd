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

    # Global arguments

    subparsers = p.add_subparsers(
        dest="command",
        required=True,
        help="Subcommand to run.",
    )

    # -------------------------------------------------------------------------
    # prop subparser: structural properties (CN, rotational constants, canonical IDs)
    # -------------------------------------------------------------------------
    p_prop = subparsers.add_parser(
        "prop",
        help="Compute structural properties (CN, rotational constants, canonical IDs).",
    )
    p_prop.add_argument(
        "structures",
        nargs="+",
        help="Paths to structure files (e.g. .xyz, .pdb, .cif).",
    )
    p_prop.add_argument(
        "--cn",
        action="store_true",
        help=(
            "Calculate coordination numbers for each structure and print them as numpy arrays."
        ),
    )
    p_prop.add_argument(
        "--rot",
        action="store_true",
        help="Calculate the rotational constants.",
    )
    p_prop.add_argument(
        "--canonical",
        action="store_true",
        help="Calculate the canonical identifiers.",
    )
    p_prop.add_argument(
        "--heavy",
        action="store_true",
        help=(
            "When calculating canonical atom identifiers, consider only heavy atoms."
        ),
    )

    # -------------------------------------------------------------------------
    # compare subparser: compare (exactly) two structures
    # -------------------------------------------------------------------------
    p_compare = subparsers.add_parser(
        "compare",
        help="Compare structures via iRMSD (default) or quaternion RMSD.",
    )
    p_compare.add_argument(
        "structures",
        nargs="+",
        help="Paths to structure files (e.g. .xyz, .pdb, .cif).",
    )
    p_compare.add_argument(
        "--quaternion",
        action="store_true",
        help=("Use the quaternion-based Cartesian RMSD instead of the invariant RMSD."),
    )
    p_compare.add_argument(
        "--inversion",
        choices=["on", "off", "auto"],
        default="auto",
        help=(
            "Control coordinate inversion in iRMSD runtypes: 'on', 'off', or 'auto' "
            "(default: auto). Used only for iRMSD."
        ),
    )
    p_compare.add_argument(
        "--heavy",
        action="store_true",
        help=("When comparing structures, consider only heavy atoms."),
    )
    p_compare.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file name (optional). If not provided, results are only printed.",
    )

    # -------------------------------------------------------------------------
    # sort subparser: sort / cluster structures based on RMSD threshold
    # -------------------------------------------------------------------------
    p_sort = subparsers.add_parser(
        "sort",
        help="Sort or cluster structures based on inter-structure RMSD.",
    )
    p_sort.add_argument(
        "structures",
        nargs="+",
        help="Paths to structure files (e.g. .xyz, .pdb, .cif).",
    )
    p_sort.add_argument(
        "--rthr",
        type=float,
        required=False,
        default=0.125,  # empirical defualt for typical molecules
        help=(
            "Inter-structure RMSD threshold for sorting in Angström. "
            "Structures closer than this threshold are treated as similar."
        ),
    )
    p_sort.add_argument(
        "--ethr",
        nargs="?",                # 0 or 1 values allowed
        type=float,               # user value is interpreted as Hartree
        default=None,             # if --ethr is not given at all
        const=1.5e-5, 
        help=(
            "Optional inter-structure energy threshold in Hartree. "
            "If set, the default is 1.5e-5 Ha or a user-specified value. "
        ),
    )
    p_sort.add_argument(
        "--inversion",
        choices=["on", "off", "auto"],
        default="auto",
        help=(
            "Control coordinate inversion when evaluating RMSDs during sorting: "
            "'on', 'off', or 'auto' (default: auto)."
        ),
    )
    p_sort.add_argument(
        "--align",
        action="store_true",
        help=("Just sort by energy and align."),
    )
    p_sort.add_argument(
        "--heavy",
        action="store_true",
        # help=("When sorting structures, consider only heavy atoms."),
        help=("TODO for sorting routines"),
    )
    p_sort.add_argument(
        "--maxprint",
        type=int,
        default=25,
        help=(
            "Printout option; determine how man rows are printed for each sorted ensemble."
        ),
    )
    p_sort.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional output file for sorted / clustered results.",
    )

    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    heavy = args.heavy  # exists in all subparsers

    # -------------------------------------------------------------------------
    # prop
    # -------------------------------------------------------------------------
    if args.command == "prop":
        molecule_list = irmsd.read_structures(args.structures)

        ran_any = False

        if args.cn:
            irmsd.compute_cn_and_print(molecule_list)
            ran_any = True

        if args.rot:
            irmsd.compute_axis_and_print(molecule_list)
            ran_any = True

        if args.canonical:
            irmsd.compute_canonical_and_print(molecule_list, heavy=heavy)
            ran_any = True

        if not ran_any:
            # No specific property selected: show help for the whole CLI
            parser.print_help()
            return 1

        return 0

    # -------------------------------------------------------------------------
    # compare
    # -------------------------------------------------------------------------
    if args.command == "compare":
        molecule_list = irmsd.read_structures(args.structures)

        if args.quaternion:
            # Quaternion RMSD (old --rmsd behavior)
            irmsd.compute_quaternion_rmsd_and_print(
                molecule_list,
                heavy=heavy,
                outfile=args.output,
            )
        else:
            # Default: iRMSD (old --irmsd behavior)
            irmsd.compute_irmsd_and_print(
                molecule_list,
                inversion=args.inversion,
                outfile=args.output,
            )

        return 0

    # -------------------------------------------------------------------------
    # sort
    # -------------------------------------------------------------------------
    if args.command == "sort":
        molecule_list = irmsd.read_structures(args.structures)

        if args.heavy:
            print("Heavy-atom mapping in sorting functionality is TODO. Sorry.")
            return 1

        if args.align:
            irmsd.sort_get_delta_irmsd_and_print(
                molecule_list,
                inversion=args.inversion,
                printlvl=1,
                maxprint=args.maxprint,
                outfile=args.output,
            )

        else:
            irmsd.sort_structures_and_print(
                molecule_list,
                rthr=args.rthr,
                inversion=args.inversion,
                printlvl=1,
                maxprint=args.maxprint,
                outfile=args.output,
                ethr=args.ethr,
            )

        return 0

    # Fallback: should not be reached due to required=True on subparsers
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
