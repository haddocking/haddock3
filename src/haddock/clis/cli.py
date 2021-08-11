#!/usr/bin/env python3

import argparse
import logging
import sys
from haddock.version import CURRENT_VERSION
from haddock.workflow import WorkflowManager
from haddock.gear.greetings import get_adieu, get_initial_greeting
from haddock.gear.prepare_run import setup_run
from haddock.error import HaddockError, ConfigurationError


def main(args=None):

    def positive_int(n):
        n = int(n)
        if n < 0:
            raise argparse.ArgumentTypeError("Minimum value is 0")
        return n

    # Command line interface parser
    parser = argparse.ArgumentParser()
    # Add logging to CLI parser
    levels = ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
    parser.add_argument("--log-level", default="INFO", choices=levels)
    # Restart option
    parser.add_argument("--restart", type=positive_int, default=0,
                        help="Restart the recipe from this course")
    # The recipe to be used
    parser.add_argument("recipe", type=argparse.FileType("r"),
                        help="The input recipe file name")
    # Version
    parser.add_argument("-V", "-v", "--version", help="show version",
                        action="version",
                        version="%s %s" % (parser.prog, CURRENT_VERSION))

    # Special case only using print instead of logging
    options = parser.parse_args()
    if not hasattr(options, "version"):
        print(get_initial_greeting())

    # Configuring logging
    logging.basicConfig(level=options.log_level,
                        format=("[%(asctime)s] %(name)s:L%(lineno)d"
                                " %(levelname)s - %(message)s"))

    try:
        params, other_params = setup_run(options.recipe.name)

    except ConfigurationError as se:
        logging.error(se)
        sys.exit()

    try:
        workflow = WorkflowManager(
            workflow_params=params,
            start=options.restart,
            **other_params,
            )

        # Main loop of execution
        workflow.run()

    except HaddockError as he:
        logging.error(he)

    # Finish
    logging.info(get_adieu())


if __name__ == "__main__":
    sys.exit(main())
