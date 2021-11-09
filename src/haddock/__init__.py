"""HADDOCK3 library."""
import logging
from os import get_terminal_size
from pathlib import Path

from haddock.libs.liblog import add_sysout_handler, add_syserr_handler


log = logging.getLogger(__name__)
log.handlers.clear()
log.setLevel(logging.DEBUG)
add_sysout_handler(log)
add_syserr_handler(log)

haddock3_source_path = Path(__file__).resolve().parent
haddock3_repository_path = haddock3_source_path.parents[1]
toppar_path = Path(haddock3_source_path, "cns", "toppar")

FCC_path = Path(
    Path(__file__).resolve().parents[1],
    'fcc',
    )


# version
version = "3.0.0"
v_major, v_minor, v_patch = version.split('.')

contact_us = 'https://github.com/haddocking/haddock3/issues'
