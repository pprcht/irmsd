import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

import irmsd
from ._version import __version__

# --- CLI ---------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="irmsd",
        description=(
            "CLI to read an arbitrary number of structures and run "
            "selected analysis commands on them."
        ),
    )

    # Global arguments

    subparsers = p.add_subparsers(
        dest="command",
        required=True,
        help="Subcommand to run.",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
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
            "Calculate coordination numbers according to the covalent D4 scheme."
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
        "--all", action="store_true", help="Calculate all of the above."
    )
    p_prop.add_argument(
        "--heavy",
        action="store_true",
        help=(
            "When calculating canonical atom identifiers, consider only heavy atoms."
        ),
    )
    p_prop.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file name (optional). For properties only pickle (pkl) files are allowed.",
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
    p_compare.add_argument(
        "--ref-idx",
        type=int,
        default=0,
        help=(
            "Index of the reference structure in the provided structure list "
            "(default: 0, i.e., the first structure)."
        ),
    )
    p_compare.add_argument(
        "--align-idx",
        type=int,
        default=1,
        help=(
            "Index of the structure to align to the reference structure "
            "(default: 1, i.e., the second structure). Used only for quaternion RMSD."
        ),
    )

    # -------------------------------------------------------------------------
    # sort subparser: sort / cluster structures based on RMSD threshold
    # -------------------------------------------------------------------------
    p_sort = subparsers.add_parser(
        "sort",
        aliases=["prune"],
        help=(
            "Sort, prune or cluster structures based on inter-structure measures."
            " By default, the more expensive iRMSD version is used. The use of the"
            " molecules' energies is optional (--ethr) but recommended."
            " To fall back to the quicker, but more empirical CREGEN workflow for"
            " ensemble sorting (using energies, quaternion RMSDs and rotational"
            " constants), use --classic"
        ),
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
        nargs="?",  # 0 or 1 values allowed
        type=float,  # user value is interpreted as Hartree
        default=None,  # if --ethr is not given at all
        const=8.0e-5,
        help=(
            "Inter-structure energy threshold in Hartree."
            " If set, the default is 8.0e-5 Ha (≈0.05 kcal/mol)"
            " or a user-specified value."
            " Optional for iRMSD-based runtypes."
        ),
    )
    p_sort.add_argument(
        "--bthr",
        type=float,
        required=False,
        default=0.01,
        help=(
            "Inter-structure rotational threshold used in"
            " the classical CREGEN sorting procedure."
            " The default is 0.01."
        ),
    )
    p_sort.add_argument(
        "--ewin",
        type=float,
        required=False,
        default=None,
        help=(
            "Energy window specification for CREGEN. Structures higher in energy"
            " than this threshold (relative to the lowest energy structure in"
            " the ensemble) will be removed."
            " There is no default (all conformers are considered)."
        ),
    )

    p_sort.add_argument(
        "--inversion",
        choices=["on", "off", "auto"],
        default="auto",
        help=(
            "Control coordinate inversion when evaluating RMSDs during sorting: "
            "'on', 'off', or 'auto' (default: auto)."
            " Only for iRMSD-based runtypes."
        ),
    )
    p_sort.add_argument(
        "--align",
        action="store_true",
        help=("Just sort by energy and align."),
    )
    p_sort.add_argument(
        "--classic",
        "--cregen",
        action="store_true",
        help=(
            "Perform conformer classification with the CREGEN workflow"
            " based on a comparison of quaternion RMSD, energy,"
            " interatomic distances, and rotational constants."
            " This routine is cheaper but more empirical than iRMSD-based sorting."
            " Does NOT restore mismatching atom order."
            " Does not keep individual rotamers."
        ),
    )
    p_sort.add_argument(
        "--heavy",
        action="store_true",
        # help=("When sorting structures, consider only heavy atoms."),
        help=("TODO for sorting routines."),
    )
    p_sort.add_argument(
        "--maxprint",
        type=int,
        default=15,
        help=(
            "Printout option; determine how many rows are printed for each sorted ensemble."
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

    printlevel = 2 if sys.stdout.isatty() else 1

    from .utils.printouts import BANNER
    print(BANNER)

    heavy = args.heavy  # exists in all subparsers

    print(f"Reading structures from: {args.structures}")
    molecule_list = irmsd.read_structures(args.structures)
    print(f"Done! {len(molecule_list)} read in total.")
    sys.stdout.flush()
    sys.stderr.flush()

    # -------------------------------------------------------------------------
    # prop
    # -------------------------------------------------------------------------
    if args.command == "prop":
        print()
        from .utils.printouts import print_molecule_summary

        ran_any = False
        flags = [args.cn, args.rot, args.canonical]
        run_multiple = sum(flags) >= 2 or args.all

        results = dict()
        if args.cn or args.all:
            results["CN"] = irmsd.compute_cn_and_print(molecule_list, run_multiple)
            ran_any = True

        if args.rot or args.all:
            results["axis"] = irmsd.compute_axis_and_print(molecule_list, run_multiple)
            ran_any = True

        if args.canonical or args.all:
            results["Canonical ID"] = irmsd.compute_canonical_and_print(
                molecule_list, heavy=heavy, run_multiple=run_multiple
            )
            ran_any = True

        if not ran_any:
            # No specific property selected: show help for the whole CLI
            parser.print_help()
            return 1

        if run_multiple:
            print_molecule_summary(molecule_list, **results)

        if args.output is not None:
            from .utils.io import dump_results_to_pickle

            outfile = dump_results_to_pickle(molecule_list, args.output, results=results)
            print(f"--> WROTE OUTPUT FILE {outfile}\n") 


        return 0

    # -------------------------------------------------------------------------
    # compare
    # -------------------------------------------------------------------------
    if args.command == "compare":

        if args.quaternion:
            # Quaternion RMSD (old --rmsd behavior)
            irmsd.compute_quaternion_rmsd_and_print(
                molecule_list,
                heavy=heavy,
                outfile=args.output,
                idx_ref=args.ref_idx,
                idx_align=args.align_idx,
            )
        else:
            # Default: iRMSD (old --irmsd behavior)
            irmsd.compute_irmsd_and_print(
                molecule_list,
                inversion=args.inversion,
                outfile=args.output,
                idx_ref=args.ref_idx,
                idx_align=args.align_idx,
            )

        return 0

    # -------------------------------------------------------------------------
    # sort/prune
    # -------------------------------------------------------------------------
    if args.command in ("sort", "prune"):

        if args.heavy:
            print("Heavy-atom mapping in sorting functionality is TODO. Sorry.")
            return 1

        if args.align:
            irmsd.sort_get_delta_irmsd_and_print(
                molecule_list,
                inversion=args.inversion,
                printlvl=printlevel,
                maxprint=args.maxprint,
                outfile=args.output,
            )

        elif args.classic:
            if args.ethr is None:
                ethr = 8.0e-5
            else:
                ethr = args.ethr
            irmsd.run_cregen_and_print(
                molecule_list,
                rthr=args.rthr,
                ethr=ethr,
                bthr=args.bthr,
                ewin=args.ewin,
                maxprint=args.maxprint,
                printlvl=printlevel,
                outfile=args.output,
            )

        else:
            irmsd.sort_structures_and_print(
                molecule_list,
                rthr=args.rthr,
                inversion=args.inversion,
                printlvl=printlevel,
                maxprint=args.maxprint,
                outfile=args.output,
                ethr=args.ethr,
                ewin=args.ewin,
            )

        return 0

    # Fallback: should not be reached due to required=True on subparsers
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
