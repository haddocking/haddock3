"""HADDOCK3 library."""

import logging
import sys
from pathlib import Path
from importlib.metadata import version as package_version

log = logging.getLogger(__name__)
log.handlers.clear()
log.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
# The logger is kept at DEBUG so that `--log-level DEBUG` can be honoured
# once `add_log_for_CLI` reconfigures the handlers. This default handler,
# however, must not leak DEBUG records to stdout for code paths that never
# go through `add_log_for_CLI` (e.g. parallel worker subprocesses that
# re-import `haddock`). Default it to INFO.
handler.setLevel(logging.INFO)
handler.setFormatter(
    logging.Formatter("[%(asctime)s %(module)s %(levelname)s] %(message)s")
)
log.addHandler(handler)

haddock3_source_path = Path(__file__).resolve().parent
haddock3_repository_path = haddock3_source_path.parents[1]
core_path = Path(haddock3_source_path, "core")
toppar_path = Path(haddock3_source_path, "cns", "toppar")
modules_defaults_path = Path(haddock3_source_path, "modules", "defaults.yaml")

FCC_path = Path(haddock3_source_path.parent, "fcc")
RMSD_path = Path(haddock3_source_path.parent, "fast-rmsdmatrix")

config_expert_levels = ("easy", "expert", "guru")
# yaml parameters with this `explevel` should be ignored when reading the yaml
_hidden_level = "hidden"


class EmptyPath:
    """Define the type EmptyPath."""

    def __str__(self) -> str:
        return ""

    def __repr__(self) -> str:
        return ""

    def __bool__(self) -> bool:
        return False


# version
version = package_version("haddock3")

v_major, v_minor, v_patch = version.split(".")[:3]

contact_us = "https://github.com/haddocking/haddock3/issues"
