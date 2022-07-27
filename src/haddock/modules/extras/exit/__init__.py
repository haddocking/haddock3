"""
Exit module
===========

Stop the workflow when this module is reached. This allows users to execute
only a certain initial part of the workflow without having to comment/uncomment
the unwanted lines.

Examples
--------

Consider the following config file example::

    [topoaa]
    (...)

    [rigidbody]
    (...)

    [exit]  # <- workflow will stop here

    [flexref]
    (...)

The workflow will stop at ``[exit]`` and ``[flexref]`` will not be performed.

You can also use this option combined with ``--restart`` and ``--extend-run``.
See examples in ``examples/docking-protein-protein/*-exit-test.cfg`` files.
"""
import shutil
from pathlib import Path

from haddock.core.exceptions import HaddockTermination
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")


class HaddockModule(BaseHaddockModule):
    """Stop the workflow when this module is reached."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order,
            path,
            *ignore,
            init_params=DEFAULT_CONFIG,
            **everything,
            ):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        return

    def _run(self):
        # removes the `exit` step folder
        self.log(self.params["message"])
        shutil.rmtree(Path.cwd())
        error = HaddockTermination()
        raise error
