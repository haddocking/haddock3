"""
HADDOCK3 CLI for restraints-related tasks.

DISCLAIMER: these scripts have been ported from old code and are not
optimized for code quality and performance. They are provided as a convenience
for the user.

USAGE::

    haddock3-restraints <TASK_NAME> <TASK_ARGS>

For the list of available tasks, run::

    haddock3-restraints -h

For the list of arguments for a given task, run::

    haddock3-restraints <TASK_NAME> -h
"""

import argparse
import sys

from haddock import log
from haddock.clis.restraints.active_passive_to_ambig import (
    actpass_to_ambig,
    add_actpass_to_ambig_arguments,
    )
from haddock.clis.restraints.calc_accessibility import (
    add_calc_accessibility_arguments,
    calc_accessibility,
    )
from haddock.clis.restraints.passive_from_active import (
    add_pass_from_act_arguments,
    passive_from_active,
    )
from haddock.clis.restraints.restrain_bodies import (
    add_restrain_bodies_arguments,
    restrain_bodies,
    )
from haddock.clis.restraints.validate_tbl import (
    add_validate_tbl_arguments,
    validate_tbl,
    )
from haddock.clis.restraints.z_surface_restraints import (
    add_z_surf_restraints_arguments,
    gen_z_surface_restraints,
    )


# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-restraints",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

subparsers = ap.add_subparsers(
    title="subcommands", description="valid subcommands", help="additional help"
)

# restrain_bodies subcommand
restrain_bodies_subcommand = subparsers.add_parser("restrain_bodies")
restrain_bodies_subcommand.set_defaults(func=restrain_bodies)
restrain_bodies_subcommand = add_restrain_bodies_arguments(restrain_bodies_subcommand)

# passive_from_active subcommand
pass_from_act_subcommand = subparsers.add_parser("passive_from_active")
pass_from_act_subcommand.set_defaults(func=passive_from_active)
pass_from_act_subcommand = add_pass_from_act_arguments(pass_from_act_subcommand)

# validate_tbl subcommand
validate_tbl_subcommand = subparsers.add_parser("validate_tbl")
validate_tbl_subcommand.set_defaults(func=validate_tbl)
validate_tbl_subcommand = add_validate_tbl_arguments(validate_tbl_subcommand)

# active_passive_to_ambig subcommand
actpass_to_ambig_subcommand = subparsers.add_parser("active_passive_to_ambig")
actpass_to_ambig_subcommand.set_defaults(func=actpass_to_ambig)
actpass_to_ambig_subcommand = add_actpass_to_ambig_arguments(
    actpass_to_ambig_subcommand
)

# calc_accessibility subcommand
calc_accessibility_subcommand = subparsers.add_parser("calc_accessibility")
calc_accessibility_subcommand.set_defaults(func=calc_accessibility)
calc_accessibility_subcommand = add_calc_accessibility_arguments(
    calc_accessibility_subcommand
)

# z_surface_restraints subcommand
z_surface_restraints_subcommand = subparsers.add_parser("z_surface_restraints")
z_surface_restraints_subcommand.set_defaults(func=gen_z_surface_restraints)
z_surface_restraints_subcommand = add_z_surf_restraints_arguments(
    z_surface_restraints_subcommand
)


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
