"""haddock3-restraints calc_accessibility subcommand.

Given a pdb file, it will calculate the relative accessibility of
the side chains and return a list of surface exposed residues.

Nucleic acids bases are considered to be always accessible.

Usage:
    haddock3-restraints calc_accessibility
        <input_pdb_file>            Path to a PDB file to analyse
        [-c <cutoff>]               Relative side-chain accessibility cutoff
        [--log_level <log_level>]   DEBUG, INFO, WARNING, ERROR, or CRITICAL
        [--export_to_actpass]       Flag to export accessible resiudes
        [--probe_radius <probe_radius>]   Probe radius for the accessibility calculation
"""


import logging
import os
from pathlib import Path
from freesasa import Structure

from haddock.core.typing import Union
from haddock.libs.librestraints import DEFAULT_PROBE_RADIUS, REL_ASA


# As this script is a subcommand of the `haddock3-restraints` cli,
# it requires its own options and arguments that are managed here.
def add_calc_accessibility_arguments(calc_accessibility_subcommand):
    """Add arguments to the calc_accessibility subcommand."""
    calc_accessibility_subcommand.add_argument(
        "input_pdb_file",
        type=str,
        help="input PDB structure.",
        )

    calc_accessibility_subcommand.add_argument(
        "-c",
        "--cutoff",
        help="Relative cutoff for sidechain accessibility",
        required=False,
        default=0.4,
        type=float,
        )
    
    calc_accessibility_subcommand.add_argument(
        "--log_level",
        default='INFO',
        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'),
        help="Logging level",
        required=False,
        )
    
    calc_accessibility_subcommand.add_argument(
        "--export_to_actpass",
        default=False,
        action="store_true",
        help="Export the exposed residues as passive to an actpass file",
        required=False,
        )
    
    calc_accessibility_subcommand.add_argument(
        "--probe_radius",
        default=DEFAULT_PROBE_RADIUS,
        type=float,
        help="Probe radius for the accessibility calculation",
        required=False,
        )

    return calc_accessibility_subcommand


##########################################################
# Define set of resiudes / bases handeled by this script #
##########################################################
DEFAULT_BASES = ["DA", "DC", "DG", "DT", "A", "C", "G", "U"]
MODIFIED_BASES = ["DJ"]
STANDARD_AA = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    ]
MODIFIED_AA = [
    'ACE', 'ALY', 'ASH', 'CFE', 'CHX', 'CSP', 'CTN', 'CYC', 'CYF', 'CYM',
    'DDZ', 'DUM', 'GLH', 'HLY', 'HYP', 'M3L', 'MLY', 'MLZ', 'MSE', 'NEP',
    'PNS', 'PTR', 'SEP', 'TOP', 'TYP', 'TYS',
    ]  # , 'THP'='TOP'?
VALID_AA = STANDARD_AA + MODIFIED_AA
VALID_BASES = DEFAULT_BASES + MODIFIED_BASES
VALID_IONS = [
    "LI", "F", "NA", "MG", "AL", "CL", "K", "CA", "V", "CR", "MN", "FE", "NI",
    "CO", "CU", "ZN", "BR", "KR", "SR", "MO", "AG", "CD", "I", "CS", "HO",
    "YB", "OS", "IR", "PT", "AU", "HG", "PB",
    ]


def get_accessibility(
        pdb_f: Union[Path, str],
        probe_radius: float = DEFAULT_PROBE_RADIUS,
        ) -> dict[str, dict[int, dict[str, float]]]:
    """Compute per-residue accessibility values.
    
    Calls `FreeSASA <https://freesasa.github.io/>`_ using its Python API
    and returns per-residue accessibility values.

    Parameters
    ----------
    pdb_f : Union[Path, str]
        Path to the PDB file of interest.
    
    Return
    ------
    resid_access : dict[str, dict[int, dict[str, float]]]
        Dictionary containing a list of accessible residues for each chain(s).
    """
    naccess_unsupported_aa = ['HEC', 'TIP', 'ACE', 'THP', 'HEB', 'CTN']
    logging.info("Calculate accessibility...")
    try:
        from freesasa import Classifier, calc, Parameters
    except ImportError as err:
        logging.error("calc_accessibility requires the 'freesasa' Python API")
        raise ImportError(err)

    # Initiate data holders
    asa_data: dict[tuple[str, str, int, str], float] = {}
    rsa_data: dict[tuple[str, str, int], float] = {}
    rel_main_chain: dict[tuple[str, str, int], float] = {}
    rel_side_chain: dict[tuple[str, str, int], float] = {}
    # Point relative asa data
    _rsa = REL_ASA['total']
    _rsa_bb = REL_ASA['bb']
    _rsa_sc = REL_ASA['sc']

    # Workaround to get relative accessibility values for Nucleic Acids
    #  -> will always be accessible !
    for na in VALID_BASES:
        _rsa[na] = 1
        _rsa_bb[na] = 1
        _rsa_sc[na] = 1

    _script_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    config_path = _script_path + '/naccess.config'
    classifier = Classifier(config_path)

    struct = Structure(pdb_f, classifier, options={})
    struct.setRadiiWithClassifier(classifier)
    # if the probe_radius is different from the default value
    # we need to redefine the parameters
    if probe_radius != DEFAULT_PROBE_RADIUS:
        new_parameters = Parameters(
            {
                'algorithm': 'LeeRichards',
                'probe-radius': probe_radius,
                'n-points': Parameters.defaultParameters['n-points'],
                'n-slices': Parameters.defaultParameters['n-slices'],
                'n-threads': Parameters.defaultParameters['n-threads'],
                }
            )
        result = calc(struct, parameters=new_parameters)
    else:
        result = calc(struct)

    # iterate over all atoms to get SASA and residue name
    for idx in range(struct.nAtoms()):
        atname = struct.atomName(idx).strip()
        resname = struct.residueName(idx).strip()
        resid = int(struct.residueNumber(idx))
        chain = struct.chainLabel(idx)
        at_uid = (chain, resname, resid, atname)
        res_uid = (chain, resname, resid)

        asa = result.atomArea(idx)
        asa_data[at_uid] = asa
        # add asa to residue
        rsa_data[res_uid] = rsa_data.get(res_uid, 0) + asa

        # If residue name is unknown, this is probably an ion,
        #  set a relative accessibility of -1
        if resname not in _rsa and (
                resname not in VALID_AA + VALID_BASES
                or resname[:2] in VALID_IONS
                or resname in naccess_unsupported_aa
                ):
            logging.warning(f"UPDATED RSA for {resname}")
            _rsa[resname] = -1
            _rsa_bb[resname] = -1
            _rsa_sc[resname] = -1

        # 3 cases: Regular amino-acid, regular nucleic acid, ion
        if (resname not in VALID_BASES and atname in ('C', 'N', 'O')) or \
                (resname in VALID_BASES and atname in ('P', 'C1', 'C9')):
            rel_main_chain[res_uid] = rel_main_chain.get(res_uid, 0) + asa
        elif all([
                resname not in VALID_AA + VALID_BASES,
                resname[:2] in VALID_IONS,
                ]):
            rel_main_chain[res_uid] = rel_main_chain.get(res_uid, 0) + asa
            rel_side_chain[res_uid] = rel_side_chain.get(res_uid, 0) + asa
        else:
            rel_side_chain[res_uid] = rel_side_chain.get(res_uid, 0) + asa
    # convert total asa to relative asa
    rsa_data.update(
        (res_uid, asa / _rsa[res_uid[1]])
        for res_uid, asa in rsa_data.items()
        )
    rel_main_chain.update(
        (res_uid, asa / _rsa_bb[res_uid[1]] * 100)
        for res_uid, asa in rel_main_chain.items()
        )
    rel_side_chain.update(
        (res_uid, asa / _rsa_sc[res_uid[1]] * 100)
        for res_uid, asa in rel_side_chain.items()
        )
    # We format to fit the pipeline
    resid_access: dict[str, dict[int, dict[str, float]]] = {}
    for res_uid, access in rel_main_chain.items():
        chain = res_uid[0]
        # resname = res_uid[1]
        resnum = res_uid[2]
        if chain not in resid_access:
            resid_access[chain] = {}
        resid_access[chain][resnum] = {
            'side_chain_rel': rel_side_chain.get(res_uid, 0),
            'main_chain_rel': access,
            }
    # Display accessible residues
    for chain in resid_access:
        logging.info(f"Chain: {chain} - {len(resid_access[chain])} residues")
    return resid_access


def apply_cutoff(
        access_data: dict[str, dict[int, dict[str, float]]],
        cutoff: float,
        ) -> dict[str, list[int]]:
    """Apply a cutoff to the sidechain relative accessibility."""
    logging.info(f'Applying cutoff to side_chain_rel - {cutoff}')
    # saving the results in a dictionary
    result_dict: dict[str, list[int]] = {}
    for chain in access_data:
        result_list: list[int] = []
        for res in access_data[chain]:
            sc_rel_accessibility = access_data[chain][res]['side_chain_rel']
            if sc_rel_accessibility >= cutoff * 100:
                result_list.append(res)
        result_list = list(set(result_list))
        result_list.sort()
        result_str = ','.join(map(str, result_list))
        logging.info(f'Chain {chain} - {result_str}')
        # putting the data in the dictionary
        result_dict[chain] = result_list
    return result_dict


def export_passive(
        result_dict: dict[str, list[int]],
        prefix: str = '',
        ) -> None:
    """Export the exposed residues as passive to an actpass file."""
    for chain in result_dict:
        filename = Path(f"{prefix}passive_{chain}.actpass")
        logging.info(f"Writing exposed residues as passive in: {filename}")
        if filename.exists():
            logging.error(
                f"File {filename} already exists, nothing performed!"
                )
        else:
            with open(filename, "w") as f:
                f.write(os.linesep)  # Add empty line as no active residues
                f.write(f"{' '.join(map(str, result_dict[chain]))}")


def calc_accessibility(
        input_pdb_file: Union[Path, str],
        cutoff: float = 0.4,
        log_level: str = "INFO",
        export_to_actpass: bool = False,
        probe_radius: float = DEFAULT_PROBE_RADIUS,
        ) -> None:
    """Calculate the accessibility of the side chains and apply a cutoff.
    
    Parameters
    ----------
    input_pdb_file : str
        Path to the PDB file.
    cutoff : float
        Relative cutoff for sidechain accessibility, default = 0.4
    log_level : str
        Logging level.
    export_to_actpass : bool
        Export the exposed residues as passive to an actpass file.
    probe_radius : float
        Probe radius for the accessibility calculation.
    """
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s L%(lineno)d %(levelname)s - %(message)s',
        datefmt='%d/%m/%Y %H:%M:%S',
        )
    # Compute per-residues accessibilities
    access_dic = get_accessibility(input_pdb_file, probe_radius)
    # Filter residues based on accessibility cutoff
    result_dict = apply_cutoff(access_dic, cutoff)

    # export residues as passive to an actpass file
    if export_to_actpass:
        export_passive(
            result_dict,
            prefix=f'{Path(input_pdb_file).stem}_',
            )
