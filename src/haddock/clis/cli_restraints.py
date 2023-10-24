"""
HADDOCK 3.0 CLI for restraints-related tasks.

USAGE::

    haddock3-restraints TASK NAME
"""
import argparse
import sys

from haddock.restraints.restrain_bodies import add_restrain_bodies_arguments, restrain_bodies
from haddock.restraints.passive_from_active import add_pass_from_act_arguments, passive_from_active
from haddock.restraints.validate_tbl import add_validate_tbl_arguments, validate_tbl
from haddock.restraints.active_passive_to_ambig import add_actpass_to_ambig_arguments, actpass_to_ambig


ANA_FOLDER = "interactive"  # name of the analysis folder

# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-restraints",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

subparsers = ap.add_subparsers(title='subcommands',
                               description='valid subcommands',
                               help='additional help')

# score subcommand
restrain_bodies_subcommand = subparsers.add_parser('restrain_bodies')
restrain_bodies_subcommand.set_defaults(func=restrain_bodies)
restrain_bodies_subcommand = add_restrain_bodies_arguments(restrain_bodies_subcommand)

# passive_from_active subcommand
pass_from_act_subcommand = subparsers.add_parser('passive_from_active')
pass_from_act_subcommand.set_defaults(func=passive_from_active)
pass_from_act_subcommand = add_pass_from_act_arguments(pass_from_act_subcommand)

# validate_tbl subcommand
validate_tbl_subcommand = subparsers.add_parser('validate_tbl')
validate_tbl_subcommand.set_defaults(func=validate_tbl)
validate_tbl_subcommand = add_validate_tbl_arguments(validate_tbl_subcommand)

# active_passive_to_ambig subcommand
actpass_to_ambig_subcommand = subparsers.add_parser('active_passive_to_ambig')
actpass_to_ambig_subcommand.set_defaults(func=actpass_to_ambig)
actpass_to_ambig_subcommand = add_actpass_to_ambig_arguments(actpass_to_ambig_subcommand)


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def maincli():
    """Execute main client."""
    args = ap.parse_args()
    cmd = vars(load_args(ap))
    cmd.pop("func")
    args.func(**cmd)


if __name__ == "__main__":
    sys.exit(maincli())