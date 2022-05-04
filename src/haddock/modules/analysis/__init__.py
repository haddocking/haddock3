"""HADDOCK3 modules related to model analysis."""


modules_using_resdic = ("caprieval", "rmsdmatrix")


def confirm_resdic_chainid_length(params):
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
