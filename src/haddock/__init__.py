"""HADDOCK3 library."""
import logging
from pathlib import Path

from haddock.libs.liblog import add_syserr_handler, add_sysout_handler


class PrePath:
    """Predefined Path."""
    def __init__(self, path):
        self.path = path

    def __call__(self, *args, **kwargs):
        return Path(self.path, *args, **kwargs)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Path({str(self.path)!r})"

    @property
    def path(self):
        """Base path."""
        return self._path

    @path.setter
    def path(self, path):
        p = Path(path)
        if not p.exists():
            emsg = f"Path {str(p.resolve())!r} does not exist."
            raise FileNotFoundError(emsg)
        self._path = path


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
