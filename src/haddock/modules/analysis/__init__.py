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
    resdic_params = (p for p in params if p.startswith("resdic_"))
    for param in resdic_params:
        chainid = param.split("_")[-1]
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
    # =====================================================================#
    # EXPLANATION #
    # =====================================================================#
    # This function is called in the analysis-based modules to OVERWRITE
    #  the mode selected by the user.
    # This must be done because the analysis modules do not support being
    #  executed over the HPC system or the GRID
    #
    # The function below returns "local" if the user selected either
    #  `batch` or `grid` as execution mode. Thus effectively forcing the
    #  execution to be local.
    return "local" if mode in ("batch", "grid") else mode
