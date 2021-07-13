#!/usr/bin/env python3

import argparse
import logging
import sys
from haddock.version import CURRENT_VERSION
from haddock.cli import greeting, adieu
from haddock.workflow import WorkflowManager
from haddock.error import HaddockError


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
        print(greeting())

    # Configuring logging
    logging.basicConfig(level=options.log_level,
                        format=("[%(asctime)s] %(levelname)s - "
                                "%(name)s: %(message)s"),
                        datefmt="%d/%m/%Y %H:%M:%S")

    try:
        # Let the chef work
        workflow = WorkflowManager(recipe_path=options.recipe.name,
                                   start=options.restart)

        # Main loop of execution
        workflow.run()

    except HaddockError as he:
        logging.error(he)

    # Finish
    logging.info(adieu())


if __name__ == "__main__":
    sys.exit(main())
