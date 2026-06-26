"""Molecular data structures."""
from functools import partial
from pathlib import Path
from typing import Any, Iterable, Optional


class Molecule:
    """
    Input molecule, usually a PDB file.

    Parameters
    ----------
    file_name : :external:py:class:`pathlib.Path`
        The path to the molecule file.

    segid : int, optional
        The ID of the segment. Defaults to ``None``.

    no_parent : boolean
        Whether to add the parent path ``..`` to the
        :py:attr:`haddock.libs.libstructure.Molecule.with_parent`.
        When set to true, the ``with_parent`` attribute returns the same
        as ``file_name``.
    """

    def __init__(self,
                 file_name: Path,
                 segid: Optional[int] = None,
                 no_parent: bool = False) -> None:
        # the rest of the code is too dependent on the Path API
        assert isinstance(file_name, Path), \
            f"`file_name` must be pathlib.Path: {type(file_name)} given"

        self.file_name = file_name
        self.segid = segid
        if no_parent:
            self.with_parent = file_name
        else:
            self.with_parent = Path('..', file_name)


def make_molecules(paths: Iterable[Path], **kwargs: Any) -> list[Molecule]:
    """Get input molecules from the data stream."""
    return list(map(partial(Molecule, **kwargs), paths))
