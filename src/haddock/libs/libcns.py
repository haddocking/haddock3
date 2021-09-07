"""CNS scripts util functions"""
from os import linesep

from haddock.libs.libmath import RandomNumberGenerator
from haddock.libs import libpdb
from haddock.core import cns_paths


RND = RandomNumberGenerator()


def generate_default_header(protonation=None):
    param = load_ff_parameters(cns_paths.parameters_file)
    top = load_ff_topology(cns_paths.topology_file)
    link = load_link(cns_paths.link_file)
    topology_protonation = load_protonation_state(protonation)
    trans_vec = load_trans_vectors(cns_paths.translation_vectors)
    tensor = load_tensor(cns_paths.tensors)
    scatter = load_scatter(cns_paths.scatter_lib)
    axis = load_axis(cns_paths.axis)
    water_box = load_waterbox(cns_paths.water_box['boxtyp20'])

    return (param, top, link, topology_protonation, trans_vec, tensor, scatter,
            axis, water_box)


def load_workflow_params(default_params):
    """Writes the values at the header section"""
    param_header = f'{linesep}! Parameters{linesep}'

    for param, v in default_params.items():

        if isinstance(v, bool):
            v = str(v).lower()
            param_header += f'eval (${param}={v}){linesep}'

        elif isinstance(v, str):
            param_header += f'eval (${param}="{v}"){linesep}'

        elif isinstance(v, int):
            param_header += f'eval (${param}={v}){linesep}'

        elif isinstance(v, float):
            param_header += f'eval (${param}={v}){linesep}'

        elif not v:
            # either 0 or empty string
            if isinstance(v, str):
                v = '\"\"'
                param_header += f'eval (${param}={v}){linesep}'
            if isinstance(v, int):
                v = 0.0
                param_header += f'eval (${param}={v}){linesep}'

    if 'chain' in default_params:
        # load molecule specific things
        for mol in default_params['chain']:
            for param in default_params['chain'][mol]:
                v = default_params['chain'][mol][param]
                # this are LOGICAL, which means no quotes
                param_header += f'eval (${param}_{mol}={v}){linesep}'

    return param_header


def load_ff_parameters(forcefield_parameters):
    """Add force-field specific parameters to its appropriate places"""
    ff_param_header = f'{linesep}! FF parameters{linesep}'
    ff_param_header += f'parameter{linesep}'
    ff_param_header += f'  @@{forcefield_parameters}{linesep}'
    ff_param_header += f'end{linesep}'

    return ff_param_header


def load_ff_topology(forcefield_topology):
    """Add force-field specific topology to its appropriate places"""
    ff_top_header = f'{linesep}! Toplogy{linesep}'
    ff_top_header += f'topology{linesep}'
    ff_top_header += f'  @@{forcefield_topology}{linesep}'
    ff_top_header += f'end{linesep}'

    return ff_top_header


def load_link(mol_link):
    """Add the link header"""
    link_header = f'{linesep}! Link file{linesep}'
    link_header += f'eval ($link_file = "{mol_link}" ){linesep}'

    return link_header


def load_trans_vectors(trans_vectors):
    """Add translation vectors"""
    trans_header = f'{linesep}! Translation vectors{linesep}'
    i = 0
    for vector_id in trans_vectors:
        vector_file = trans_vectors[vector_id]
        trans_header += f'eval ($trans_vector_{i} = "{vector_file}" ){linesep}'
        i += 1

    return trans_header


def load_tensor(tensor):
    """Add tensor information"""
    tensor_header = f'{linesep}! Tensors{linesep}'
    for tensor_id in tensor:
        tensor_file = tensor[tensor_id]
        tensor_header += f'eval (${tensor_id} = "{tensor_file}" ){linesep}'

    return tensor_header


def load_axis(axis):
    """Add axis"""
    axis_header = f'{linesep}! Axis{linesep}'
    for axis_id in axis:
        axis_file = axis[axis_id]
        axis_header += f'eval (${axis_id} = "{axis_file}" ){linesep}'

    return axis_header


def load_scatter(scatter_lib):
    """Add scatter library"""
    scatter_header = f'{linesep}! Scatter lib{linesep}'
    scatter_header += f'eval ($scatter_lib = "{scatter_lib}" ){linesep}'

    return scatter_header


def load_waterbox(waterbox_param):
    """Add waterbox information"""
    water_header = f'{linesep}! Water box{linesep}'
    water_header += f'eval ($boxtyp20 = "{waterbox_param}" ){linesep}'

    return water_header


def load_ambig(ambig_f):
    """Add ambig file"""
    ambig_str = f'eval ($ambig_fname="{str(ambig_f)}"){linesep}'
    return ambig_str


def load_unambig(unambig_f):
    """Add unambig file"""
    unambig_str = f'eval ($unambig_fname="{str(unambig_f)}"){linesep}'
    return unambig_str


def load_hbond(hbond_f):
    """Add hbond file"""
    hbond_str = f'eval ($hbond_fname="{hbond_f}"){linesep}'
    return hbond_str


def load_dihe(dihe_f):
    """Add dihedral file"""
    dihe_str = f'eval ($dihe_fname="{dihe_f}"){linesep}'
    return dihe_str


def load_tensor_tbl(tensor_f):
    """Add tensor tbl file"""
    tensor_str = f'eval ($tensor_tbl="{tensor_f}"){linesep}'
    return tensor_str


def prepare_output(output_psf_filename, output_pdb_filename):
    """Output of the CNS file"""
    output = f'{linesep}! Output structure{linesep}'
    output += ("eval ($output_psf_filename="
               f" \"{output_psf_filename}\"){linesep}")
    output += ("eval ($output_pdb_filename="
               f" \"{output_pdb_filename}\"){linesep}")
    return output


def load_protonation_state(protononation):
    """Prepare the CNS protononation"""
    protonation_header = ''
    if protononation and isinstance(protononation, dict):
        protonation_header += f'{linesep}! Protonation states{linesep}'

        for i, chain in enumerate(protononation):
            hise_l = [0] * 10
            hisd_l = [0] * 10
            hisd_counter = 0
            hise_counter = 0
            for res in protononation[chain]:
                state = protononation[chain][res].lower()
                if state == 'hise':
                    hise_l[hise_counter] = res
                    hise_counter += 1
                if state == 'hisd':
                    hisd_l[hisd_counter] = res
                    hisd_counter += 1

            hise_str = ''
            for e in [(i + 1, c + 1, r) for c, r in enumerate(hise_l)]:
                hise_str += (f'eval ($toppar.hise_resid_{e[0]}_{e[1]}'
                             f' = {e[2]}){linesep}')
            hisd_str = ''
            for e in [(i + 1, c + 1, r) for c, r in enumerate(hisd_l)]:
                hisd_str += (f'eval ($toppar.hisd_resid_{e[0]}_{e[1]}'
                             f' = {e[2]}){linesep}')

            protonation_header += hise_str
            protonation_header += hisd_str

    return protonation_header


# This is used by docking
def prepare_multiple_input(pdb_input_list, psf_input_list):
    input_str = f'{linesep}! Input structure{linesep}'
    for psf in psf_input_list:
        input_str += f'structure{linesep}'
        input_str += f'  @@{psf}{linesep}'
        input_str += f'end{linesep}'

    for pdb in pdb_input_list:
        input_str += f'coor @@{pdb}{linesep}'

    ncomponents = len(psf_input_list)
    input_str += f'eval ($ncomponents={ncomponents}){linesep}'

    seed = RND.randint(100, 999)
    input_str += f'eval ($seed={seed}){linesep}'

    return input_str


# This is used by Topology and Scoring
def prepare_single_input(pdb_input, psf_input=None):
    """Input of the CNS file.

    This section will be written for any recipe even if some CNS variables
    are not used, it should not be an issue.
    """
    input_str = f'{linesep}! Input structure{linesep}'

    if psf_input:
        # if isinstance(psf_input, str):
        input_str += f'structure{linesep}'
        input_str += f'  @@{psf_input}{linesep}'
        input_str += f'end{linesep}'
        input_str += f'coor @@{pdb_input}{linesep}'
        if isinstance(psf_input, list):
            input_str += f'structure{linesep}'
            for psf in psf_input:
                input_str += f'  @@{psf}{linesep}'
            input_str += f'end{linesep}'
    # $file variable is still used by some CNS recipes, need refactoring!
    input_str += f'eval ($file=\"{pdb_input}\"){linesep}'
    segids, chains = libpdb.identify_chainseg(pdb_input)
    chainsegs = sorted(list(set(segids) | set(chains)))

    ncomponents = len(chainsegs)

    input_str += f'eval ($ncomponents={ncomponents}){linesep}'

    for i, segid in enumerate(chainsegs):
        input_str += f'eval ($prot_segid_{i+1}="{segid}"){linesep}'

    seed = RND.randint(100, 999)
    input_str += f'eval ($seed={seed}){linesep}'

    return input_str
