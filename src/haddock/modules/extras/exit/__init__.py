"""
Exit module.

Stop the workflow when this module is reached.

Examples
--------

```
[topoaa]
(...)

[rigidbody]
(...)

[exit]  # <- workflow will stop here

[flexref]
(...)
"""
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
        error = HaddockTermination(self.params["message"])
        raise error
