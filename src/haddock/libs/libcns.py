"""CNS scripts util functions."""
import itertools
from os import linesep
from pathlib import Path

from haddock import log
from haddock.core import cns_paths
from haddock.libs import libpdb
from haddock.libs.libfunc import false, true
from haddock.libs.libmath import RandomNumberGenerator
from haddock.libs.libontology import PDBFile


RND = RandomNumberGenerator()


def generate_default_header(protonation=None):
    """Generate CNS default header."""
    link = load_link(cns_paths.link_file)
    topology_protonation = load_protonation_state(protonation)
    trans_vec = load_trans_vectors(cns_paths.translation_vectors)
    tensor = load_tensor(cns_paths.tensors)
    scatter = load_scatter(cns_paths.scatter_lib)
    axis = load_axis(cns_paths.axis)
    water_box = load_waterbox(cns_paths.water_box["boxtyp20"])

    return (
        link,
        topology_protonation,
        trans_vec,
        tensor,
        scatter,
        axis,
        water_box,
        )


def filter_empty_vars(v):
    """
    Filter empty variables.

    See: https://github.com/haddocking/haddock3/issues/162

    Returns
    -------
    bool
        Returns `True` if the variable is not empty, and `False` if
        the variable is empty. That is, `False` reflects those variables
        that should not be written in CNS.

    Raises
    ------
    TypeError
        If the type of `value` is not supported by CNS.
    """
    cases = (
        (lambda x: isinstance(x, str), true),
        (lambda x: isinstance(x, bool), true),  # it should return True
        (lambda x: isinstance(x, Path), true),
        (lambda x: type(x) in (int, float), true),
        (lambda x: x is None, false),
        )

    for detect, give in cases:
        if detect(v):
            return give(v)
    else:
        emsg = f"Value {v!r} has a unknown type for CNS: {type(v)}."
        log.error(emsg)
        raise TypeError(emsg)


def load_workflow_params(
        params,
        param_header=f"{linesep}! Parameters{linesep}",
        ):
    """
    Write the values at the header section.

    "Empty variables" are ignored. These are defined accoring to
    :func:`filter_empty_vars`.

    Parameters
    ----------
    params : dict
        Dictionary containing the key:value pars for the parameters to
        be written to CNS. Values cannot be of dictionary type.

    Returns
    -------
    str
        The string with the CNS parameters defined.
    """
    non_empty_parameters = (
        (k, v) for k, v in params.items() if filter_empty_vars(v)
        )

    # types besides the ones in the if-statements should not enter this loop
    for param, v in non_empty_parameters:
        param_header += write_eval_line(param, v)

    return param_header


def write_eval_line(param, value, eval_line="eval (${}={})"):
    """Write the CNS eval line depending on the type of `value`."""
    eval_line += linesep

    if isinstance(value, bool):
        value = str(value).lower()
        return eval_line.format(param, value)

    elif isinstance(value, str):
        value = '"' + value + '"'
        return eval_line.format(param, value)

    elif isinstance(value, Path):
        value = '"' + str(value) + '"'
        return eval_line.format(param, value)

    elif isinstance(value, (int, float)):
        return eval_line.format(param, value)

    else:
        emsg = f"Unexpected type when writing CNS header: {type(value)}"
        log.error(emsg)
        raise TypeError(emsg)


def load_link(mol_link):
    """Add the link header."""
    link_header = f"{linesep}! Link file{linesep}"
    link_header += f'eval ($link_file = "{mol_link}" ){linesep}'

    return link_header


def load_trans_vectors(trans_vectors):
    """Add translation vectors."""
    trans_header = f"{linesep}! Translation vectors{linesep}"
    i = 0
    for vector_id in trans_vectors:
        vector_file = trans_vectors[vector_id]
        trans_header += f'eval ($trans_vector_{i} = "{vector_file}" ){linesep}'
        i += 1

    return trans_header


def load_tensor(tensor):
    """Add tensor information."""
    tensor_header = f"{linesep}! Tensors{linesep}"
    for tensor_id in tensor:
        tensor_file = tensor[tensor_id]
        tensor_header += f'eval (${tensor_id} = "{tensor_file}" ){linesep}'

    return tensor_header


def load_axis(axis):
    """Add axis."""
    axis_header = f"{linesep}! Axis{linesep}"
    for axis_id in axis:
        axis_file = axis[axis_id]
        axis_header += f'eval (${axis_id} = "{axis_file}" ){linesep}'

    return axis_header


def load_scatter(scatter_lib):
    """Add scatter library."""
    scatter_header = f"{linesep}! Scatter lib{linesep}"
    scatter_header += f'eval ($scatter_lib = "{scatter_lib}" ){linesep}'

    return scatter_header


def load_waterbox(waterbox_param):
    """Add waterbox information."""
    water_header = f"{linesep}! Water box{linesep}"
    water_header += f'eval ($boxtyp20 = "{waterbox_param}" ){linesep}'

    return water_header


def load_ambig(ambig_f):
    """Add ambig file."""
    ambig_str = f'eval ($ambig_fname="{str(ambig_f)}"){linesep}'
    return ambig_str


def load_unambig(unambig_f):
    """Add unambig file."""
    unambig_str = f'eval ($unambig_fname="{str(unambig_f)}"){linesep}'
    return unambig_str


def load_hbond(hbond_f):
    """Add hbond file."""
    hbond_str = f'eval ($hbond_fname="{hbond_f}"){linesep}'
    return hbond_str


def load_dihe(dihe_f):
    """Add dihedral file."""
    dihe_str = f'eval ($dihe_fname="{dihe_f}"){linesep}'
    return dihe_str


def load_tensor_tbl(tensor_f):
    """Add tensor tbl file."""
    tensor_str = f'eval ($tensor_tbl="{tensor_f}"){linesep}'
    return tensor_str


def prepare_output(output_psf_filename, output_pdb_filename):
    """Output of the CNS file."""
    output = f"{linesep}! Output structure{linesep}"
    output += (
        "eval ($output_psf_filename=" f' "{output_psf_filename}"){linesep}'
        )
    output += (
        "eval ($output_pdb_filename=" f' "{output_pdb_filename}"){linesep}'
        )
    return output


def load_protonation_state(protononation):
    """Prepare the CNS protononation."""
    protonation_header = ""
    if protononation and isinstance(protononation, dict):
        protonation_header += f"{linesep}! Protonation states{linesep}"

        for i, chain in enumerate(protononation):
            hise_l = [0] * 10
            hisd_l = [0] * 10
            hisd_counter = 0
            hise_counter = 0
            for res in protononation[chain]:
                state = protononation[chain][res].lower()
                if state == "hise":
                    hise_l[hise_counter] = res
                    hise_counter += 1
                if state == "hisd":
                    hisd_l[hisd_counter] = res
                    hisd_counter += 1

            hise_str = ""
            for e in [(i + 1, c + 1, r) for c, r in enumerate(hise_l)]:
                hise_str += (
                    f"eval ($toppar.hise_resid_{e[0]}_{e[1]}"
                    f" = {e[2]}){linesep}"
                    )
            hisd_str = ""
            for e in [(i + 1, c + 1, r) for c, r in enumerate(hisd_l)]:
                hisd_str += (
                    f"eval ($toppar.hisd_resid_{e[0]}_{e[1]}"
                    f" = {e[2]}){linesep}"
                    )

            protonation_header += hise_str
            protonation_header += hisd_str

    return protonation_header


# This is used by docking
def prepare_multiple_input(pdb_input_list, psf_input_list):
    """Prepare multiple input files."""
    input_str = f"{linesep}! Input structure{linesep}"
    for psf in psf_input_list:
        input_str += f"structure{linesep}"
        input_str += f"  @@{psf}{linesep}"
        input_str += f"end{linesep}"

    ncount = 1
    for pdb in pdb_input_list:
        input_str += f"coor @@{pdb}{linesep}"
        input_str += (
            f"eval ($input_pdb_filename_{ncount}=" f' "{pdb}"){linesep}'
            )
        ncount += 1

    # check how many chains there are across all the PDBs
    chain_l = []
    for pdb in pdb_input_list:
        for element in libpdb.identify_chainseg(pdb):
            chain_l.append(element)
    ncomponents = len(set(itertools.chain(*chain_l)))
    input_str += f"eval ($ncomponents={ncomponents}){linesep}"

    seed = RND.randint(100, 999)
    input_str += f"eval ($seed={seed}){linesep}"

    return input_str


# This is used by Topology and Scoring
def prepare_single_input(pdb_input, psf_input=None):
    """Input of the CNS file.

    This section will be written for any recipe even if some CNS variables
    are not used, it should not be an issue.
    """
    input_str = f"{linesep}! Input structure{linesep}"

    if psf_input:
        # if isinstance(psf_input, str):
        input_str += f"structure{linesep}"
        input_str += f"  @@{psf_input}{linesep}"
        input_str += f"end{linesep}"
        input_str += f"coor @@{pdb_input}{linesep}"
        if isinstance(psf_input, list):
            input_str += f"structure{linesep}"
            for psf in psf_input:
                input_str += f"  @@{psf}{linesep}"
            input_str += f"end{linesep}"
    # $file variable is still used by some CNS recipes, need refactoring!
    input_str += f'eval ($file="{pdb_input}"){linesep}'
    segids, chains = libpdb.identify_chainseg(pdb_input)
    chainsegs = sorted(list(set(segids) | set(chains)))

    ncomponents = len(chainsegs)

    input_str += f"eval ($ncomponents={ncomponents}){linesep}"

    for i, segid in enumerate(chainsegs):
        input_str += f'eval ($prot_segid_{i+1}="{segid}"){linesep}'

    seed = RND.randint(100, 99999)
    input_str += f"eval ($seed={seed}){linesep}"

    return input_str


def prepare_cns_input(
        model_number,
        input_element,
        step_path,
        recipe_str,
        defaults,
        identifier,
        ambig_fname=None,
        native_segid=False,
        ):
    """Generate the .inp file needed by the CNS engine."""
    # read the default parameters
    default_params = load_workflow_params(defaults)
    (
        link,  # unused
        topology_protonation,
        trans_vec,  # unused
        tensor,  # unused
        scatter,  # unused
        axis,  # unused
        water_box,  # unused
        ) = generate_default_header()

    # write the PDBs
    pdb_list = []
    if isinstance(input_element, (list, tuple)):
        for pdb in input_element:
            pdb_list.append(pdb.full_name)
    else:
        pdb_list.append(input_element.full_name)

    # write the PSFs
    psf_list = []
    if isinstance(input_element, (list, tuple)):
        for pdb in input_element:
            if isinstance(pdb.topology, (list, tuple)):
                for psf in pdb.topology:
                    psf_fname = psf.full_name
                    psf_list.append(psf_fname)
            else:
                psf_fname = pdb.topology.full_name
                psf_list.append(psf_fname)
    elif isinstance(input_element.topology, (list, tuple)):
        pdb = input_element  # for clarity
        for psf in pdb.topology:
            psf_fname = psf.full_name
            psf_list.append(psf_fname)
    else:
        pdb = input_element  # for clarity
        psf_fname = pdb.topology.full_name
        psf_list.append(psf_fname)

    input_str = prepare_multiple_input(pdb_list, psf_list)

    if ambig_fname:
        ambig_str = load_ambig(ambig_fname)
    else:
        ambig_str = ""

    output_pdb_filename = step_path / f"{identifier}_{model_number}.pdb"
    output = f"{linesep}! Output structure{linesep}"
    output += (
        f"eval ($output_pdb_filename=" f' "{output_pdb_filename}"){linesep}'
        )

    segid_str = ""
    if native_segid:
        pdb_list = []
        if isinstance(input_element, (list, tuple)):
            id_counter = 0
            for pdb in input_element:
                segids, chains = libpdb.identify_chainseg(
                    pdb.full_name, sort=False
                    )
                chainsegs = sorted(list(set(segids) | set(chains)))
                for i, _ in enumerate(chainsegs, start=1):
                    segid_str += (f"eval ($prot_segid_{i}=\"{id_counter}\")"
                                  f"{linesep}")
                    id_counter += 1
        else:
            segids, chains = libpdb.identify_chainseg(
                input_element.full_name, sort=False
                )
            chainsegs = sorted(list(set(segids) | set(chains)))
            for i, id in enumerate(chainsegs, start=1):
                segid_str += f"eval ($prot_segid_{i}=\"{id}\"){linesep}"

    output += f"eval ($count=" f" {model_number}){linesep}"
    inp = (
        default_params
        + input_str
        + output
        + topology_protonation
        + ambig_str
        + segid_str
        + recipe_str
        )

    inp_file = Path(step_path, f"{identifier}_{model_number}.inp")
    with open(inp_file, "w") as fh:
        fh.write(inp)

    return inp_file


def prepare_expected_pdb(model_obj, model_nb, path, identifier):
    """Prepare a PDBobject."""
    expected_pdb_fname = Path(path, f"{identifier}_{model_nb}.pdb")
    pdb = PDBFile(expected_pdb_fname, path=path)
    if type(model_obj) == tuple:
        pdb.topology = [p.topology for p in model_obj]
    else:
        pdb.topology = model_obj.topology
    return pdb
