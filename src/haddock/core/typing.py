"""
Useful type aliases and type variables for type-hints support.

Visit [PEP 484 Type Hints](https://peps.python.org/pep-0484/) for more help.

Variables endswith `T` are TypeVars for generic function signatures.
They are defined to

    - accept constraint types(Any type when no type is given)
    - accept any subclass by using keyword `bound`
    - be parameter of Generic and subclass of Generic,
      espacially for invariant Collections(list,sequence and so on).
    - show more than once in generic functions to indicates that
      these types are exactly same.

Don't import other classes, type aliases or type variables defined
in `haddock` into this module.
This may lead to circular import problem.
"""


from argparse import ArgumentParser, Namespace
from logging import FileHandler, Handler, StreamHandler
from pathlib import Path
from typing import (
    Any,
    Callable,
    Container,
    Generator,
    Generic,
    Iterable,
    Iterator,
    Literal,
    Mapping,
    MutableMapping,
    Optional,
    Protocol,
    Sequence,
    Sized,
    SupportsInt,
    TextIO,
    Type,
    TypedDict,
    TypeVar,
    Union,
    runtime_checkable,
    )

from numpy import float64
from numpy.typing import NDArray
from pandas import DataFrame
from pandas.core.groupby.generic import DataFrameGroupBy
from plotly.graph_objects import Figure


AnyT = TypeVar('AnyT')
"""Arbitrary type for generic function signatures."""

PT = TypeVar('PT')
"""
Arbitrary parameter type for generic `Callable` signatures.

For py3.11+, `ParamSpec` and `Concatenate` may be more explict
while making signatures longer.
"""

# Literals
ImgFormat = Literal["png", "pdf", "svg", "jpeg", "webp"]
"""Supported image formats by plotly."""

ExpertLevel = Literal["all", "easy", "expert", "guru"]
"""The expertise level of the parameters.`hidden` not inclued."""

LogLevel = Literal["INFO", "DEBUG", "ERROR", "WARNING", "CRITICAL"]


# Protocols
@runtime_checkable
class SupportsAdd(Protocol):
    """
    An ABC with one abstract method __add__.
    
    Can be used in `isinstance` and `issubclass`.

    Assuming that the operands and the result's types
    are same for most common situations.

    Examples
    --------

    >>> isinstance(int, SupportsAdd)
    True
    """

    __slots__ = ()

    def __add__(self: AnyT, other: AnyT, /) -> AnyT:
        raise NotImplementedError()


SupportsAddT = TypeVar("SupportsAddT", bound=SupportsAdd)
"""TypeVar of SupportsAdd for Generic type check."""


@runtime_checkable
class SupportsRun(Protocol):
    """
    An ABC with one abstract method run().
    
    Can be used in `isinstance` and `issubclass`.

    Examples
    --------

    >>> from multiprocessing import Process
    >>> isinstance(Process, SupportsRun)
    True
    """

    __slots__ = ()

    def run(self) -> Any:
        """Do something."""
        raise NotImplementedError()


SupportsRunT = TypeVar("SupportsRunT", bound=SupportsRun)
"""TypeVar of SupportsRun for Generic type check."""

# other type alias and type variables
NDFloat = NDArray[float64]

FilePath = Union[str, Path]

FilePathT = TypeVar("FilePathT", bound=FilePath)
"""
Generic type variable for file paths.

If the first annotated variable is str,
the second annotated variable will be str instead of Path,vice versa.
"""

LineIterSource = Union[Iterable[str], FilePath]
"""An Object who can generate lines in some way."""

ParamMap = MutableMapping[str, Any]
"""
A mutable mapping object contains parameters generated from config
files or functions.

The keys are always `str`.
"""

ParamMapT = TypeVar('ParamMapT', bound=ParamMap)
"""
Generic type of `MutableMapping[str, Any]` for generic function sigunature.
"""

ParamDict = dict[str, Any]
"""
A dict contains parameters generated from config files or functions.

The keys are always `str`.

Can be called like `dict` and returns a python `dict`.

Examples
--------

>>> ParamDict(a=1) == {'a':1}
True

>>> ParamDict({'a':1}) == {'a':1}
True

>>> ParamDict().setdefault('a', 1)
1
"""

ParamDictT = TypeVar('ParamDictT', bound=ParamDict)
"""
Generic type of `dict[str, Any]` for generic function sigunature.
"""

ModuleParams = dict[str, ParamDict]
"""A dict contains modules' parameters."""

AtomsDict = dict[str, list[str]]
"""
A dict contains all kinds of residues' atom type infomations
(including non-standard residue and atom type).

Key is residue name, value is a list of the residue's atom types.
"""

StreamHandlerT = TypeVar('StreamHandlerT',
                         bound=Union[StreamHandler, FileHandler])
