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

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.exceptions import HaddockTermination
from haddock.core.typing import Any, FilePath
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """Stop the workflow when this module is reached."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order: int,
            path: Path,
            *ignore: Any,
            init_params: FilePath = DEFAULT_CONFIG,
            **everything: Any,
            ) -> None:
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if contact executable is compiled."""
        return

    def _run(self) -> None:
        # removes the `exit` step folder
        self.log(self.params["message"])
        shutil.rmtree(Path.cwd())
        error = HaddockTermination()
        raise error
