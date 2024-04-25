"""HADDOCK3 modules related to model analysis."""

from typing import Iterable


modules_using_resdic = ("caprieval", "rmsdmatrix", "alascan", "sasascore")


def confirm_resdic_chainid_length(params: Iterable[str]) -> None:
    """
    Confirm resdic params have chain IDs of length one.

    Parameters
    ----------
    params : dict, or list of keys
        The parameters.

    Raises
    ------
    ValueError
        If a `resdic_*` parameter has chain IDs with more than one letter.
        For example, `resdic_AB`.
    """
    resdic_params = (p for p in params if p.startswith('resdic_'))
    for param in resdic_params:
        chainid = param.split('_')[-1]
        if len(chainid) > 1:
            raise ValueError(
                f"We found the parameter {param!r} which has "
                "more than one character in the chain "
                "identifier. Chain IDs should have only one character."
                )

def get_analysis_exec_mode(mode: str) -> str:
    """
    Get the execution mode for analysis modules.

    Parameters
    ----------
    exec_mode : str
        The execution mode to use.

    Returns
    -------
    str
        The execution mode to use for the analysis modules. 
        If it's "batch", it will be changed to "local".
    """
    if mode != "batch":
        exec_mode = mode
    else:
        exec_mode = "local"
    return exec_mode