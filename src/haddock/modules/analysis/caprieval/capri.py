"""CAPRIeval-specific helpers: weight dumping and CNS step detection.

All shared CAPRI computation that are being used both in caprieval and
in caprifilter (namely CAPRI class and metric functions) are in haddock.libs.libcapri.
"""

import json
from pathlib import Path

from haddock import log
from haddock.core.defaults import CNS_MODULES
from haddock.core.typing import Union
from haddock.gear.config import load as read_config
from haddock.modules import get_module_steps_folders


WEIGHTS = ["w_elec", "w_vdw", "w_desolv", "w_bsa", "w_air"]


def get_previous_cns_step(sel_steps: list, st_order: int) -> Union[str, None]:
    """Get the previous CNS step.

    Parameters
    ----------
    sel_steps : list
        List of folder names that are modules.
    st_order : int
        Index of the current module.

    Returns
    -------
    cns_step : str or None
        Name of the most recent CNS step.
    """
    cns_step = None
    # keep only well-formed "<order>_<module>" step folders; the module name
    # itself may contain underscores, so do not assume a single underscore
    sel_steps = [step for step in sel_steps if "_" in step]
    mod = min(st_order - 1, len(sel_steps) - 1)
    while mod > -1:
        st_name = sel_steps[mod].split("_", 1)[1]
        if st_name in CNS_MODULES:
            cns_step = sel_steps[mod]
            break
        mod -= 1
    return cns_step


def save_scoring_weights(cns_step: str) -> Path:
    """Save the scoring weights in a json file.

    Parameters
    ----------
    cns_step : str
        Name of the CNS step.

    Returns
    -------
    scoring_params_fname : Path
        Path to the json file.
    """
    # `params.cfg` is a mix of toml and cfg format, set strict to false
    cns_params = read_config(Path("..", cns_step, "params.cfg"), strict=False)
    key = list(cns_params["final_cfg"].keys())[0]
    scoring_pars = {kv: cns_params["final_cfg"][key][kv] for kv in WEIGHTS}

    scoring_params_fname = Path("weights_params.json")
    with open(scoring_params_fname, "w", encoding="utf-8") as jsonf:
        json.dump(scoring_pars, jsonf, indent=4)
    return scoring_params_fname


def dump_weights(order: int) -> None:
    sel_steps = get_module_steps_folders(Path(".."))
    cns_step = get_previous_cns_step(sel_steps=sel_steps, st_order=order)
    if cns_step:
        log.info(f"Found previous CNS step: {cns_step}")
        scoring_params_fname = save_scoring_weights(cns_step)
        log.info(f"Saved scoring weights to: {scoring_params_fname}")
    else:
        log.info("No previous CNS step found. Cannot save scoring weights.")
