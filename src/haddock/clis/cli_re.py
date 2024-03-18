"""
HADDOCK 3.0 CLI for interactive tasks.
"""
import argparse
import sys
import os

from haddock.clis.re.clustfcc import add_clustfcc_arguments, reclustfcc
from haddock.clis.re.clustrmsd import add_clustrmsd_arguments, reclustrmsd
from haddock.clis.re.score import add_rescore_arguments, rescore

# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-re",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

subparsers = ap.add_subparsers(title='subcommands',
                               help='valid subcommands',
                              )

# score subcommand
rescore_subcommand = subparsers.add_parser('score')
rescore_subcommand.set_defaults(func=rescore)
rescore_subcommand = add_rescore_arguments(rescore_subcommand)
# clustrmsd subcommand
clustrmsd_subcommand = subparsers.add_parser('clustrmsd')
clustrmsd_subcommand.set_defaults(func=reclustrmsd)
clustrmsd_subcommand = add_clustrmsd_arguments(clustrmsd_subcommand)
# clustfcc subcommand
clustfcc_subcommand = subparsers.add_parser('clustfcc')
clustfcc_subcommand.set_defaults(func=reclustfcc)
clustfcc_subcommand = add_clustfcc_arguments(clustfcc_subcommand)


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def maincli():
    """Execute main client."""
    args = ap.parse_args()
    cmd = vars(load_args(ap))
    if "func" not in cmd.keys():
        sys.stderr.write(f"Error: no subcommand given.{os.linesep}")
        ap.print_help()
        sys.exit(1)
    else:
        cmd.pop("func")
    args.func(**cmd)


if __name__ == "__main__":
    sys.exit(maincli())
