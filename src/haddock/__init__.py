"""HADDOCK3 library."""
import logging
import sys
from pathlib import Path


log = logging.getLogger(__name__)
log.handlers.clear()
log.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter(
    "[%(asctime)s %(module)s %(levelname)s] %(message)s"))
log.addHandler(handler)

haddock3_source_path = Path(__file__).resolve().parent
haddock3_repository_path = haddock3_source_path.parents[1]
core_path = Path(haddock3_source_path, "core")
toppar_path = Path(haddock3_source_path, "cns", "toppar")
modules_defaults_path = Path(haddock3_source_path, "modules", "defaults.yaml")

FCC_path = Path(haddock3_source_path.parent, 'fcc')

config_expert_levels = ("easy", "expert", "guru")
# yaml parameters with this `explevel` should be ignored when reading the yaml
_hidden_level = "hidden"


class EmptyPath:
    """Define the type EmptyPath."""

    def __str__(self):
        return ""

    def __repr__(self):
        return ""

    def __bool__(self):
        return False


def get_package_version():
    """Use package egg info to obain build version."""
    # point package path
    ppath = __path__[0]
    # point setup egg data
    pkg_info_fpath = f'{ppath}3.egg-info/PKG-INFO'
    # read this file
    with open(pkg_info_fpath, 'r') as filin:
        # loop over lines
        for _ in filin:
            # find the line where version is written
            if _.startswith('Version'):
                # return the version
                return _.split(':')[-1].strip()


def split_version(version: str) -> tuple:
    """Split version into major, minor and patch."""
    s_version = version.split('.')
    v_major = s_version[0]
    v_minor = s_version[1]
    v_patch = '.'.join(s_version[2:])
    return (v_major, v_minor, v_patch)


# version
version = get_package_version()
v_major, v_minor, v_patch = split_version(version)

contact_us = 'https://github.com/haddocking/haddock3/issues'
