"""
Extract the default configuration file for each module.

By default presents the parameters from the three expertise levels:
    `basic`, `intermediate`, and `guru`.

USAGE:
    $ haddock3-cfg -m MODULE
    $ haddock3-cfg -m MODULE -l LEVEL
"""
import argparse
import collections
import importlib
import os
import sys
from pathlib import Path

from haddock import config_expert_levels
from haddock.lib.libio import read_from_yaml
from haddock.modules import modules_category


ap = argparse.ArgumentParser(
    prog="HADDOCK3 config retriever",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    "-m",
    dest="module",
    help="The module for which you want to retrieve the default configuration.",
    required=True,
    choices=sorted(modules_category.keys()),
    )

ap.add_argument(
    "-l",
    dest="explevel",
    required=False,
    help="The expertise level of the parameters. Defaults to \"all\".",
    default="all",
    choices=("basic", "intermediate", "guru", "all"),
    )


# command-line client helper functions
# load_args, cli, maincli
def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(module, explevel):
    """Extract the default configuration file for a given module."""
    module_name = ".".join((
        'haddock',
        'modules',
        modules_category[module],
        module,
        ))
    module_lib = importlib.import_module(module_name)
    cfg = module_lib.DEFAULT_CONFIG

    ycfg = read_from_yaml(cfg)

    new_config = []
    new_config.append(f"[{module}]")

    if explevel == "all":
        for level in config_expert_levels:
            new_config.append(write_params(ycfg[level], module, level))

    else:
        new_config.append(write_params(ycfg[explevel], module, explevel))

    Path(f"haddock3_{module}.cfg").write_text(os.linesep.join(new_config))

    return


def write_params(cfg, module, subtitle=None):
    """Write config parameters with comments."""
    params = []

    if subtitle:
        params.append(f"# {subtitle}")

    for param_name, param in cfg.items():

        if isinstance(param, collections.abc.Mapping) \
                and "default" not in param:

            params.append("")  # give extra space
            curr_module = f"{module}.{param_name}"
            params.append(f"[{curr_module}]")
            params.append(write_params(param, module=curr_module))

        else:

            comment = []
            comment.append(param["hoover"])
            if "min" in param:
                comment.append(f"$min {param['min']} $max {param['max']}")
            if "length" in param:
                comment.append(f"$maxlen {param['length']}")

            params.append("{} = {!r}  # {}".format(
                param_name,
                param["default"],
                " ".join(comment),
                ))

            if param["type"] == "list":
                params.append(os.linesep)

    return os.linesep.join(params)


if __name__ == "__main__":
    sys.exit(maincli())
