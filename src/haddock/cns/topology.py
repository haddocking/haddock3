"""CNS language translation for topology"""
from os import linesep
from haddock.defaults import Default
from haddock.ontology import Format
from haddock.cns.util import load_recipe_params, prepare_input


def get_topology_header(protonation=None):
    param = load_ff_parameters(Default.PARAMETERS_FILE)
    top = load_ff_topology(Default.TOPOLOGY_FILE)
    link = load_link(Default.LINK_FILE)
    topology_protonation = load_protonation_state(protonation)
    trans_vec = load_trans_vectors(Default.TRANSLATION_VECTORS)
    tensor = load_tensor(Default.TENSORS)
    scatter = load_scatter(Default.SCATTER_LIB)
    axis = load_axis(Default.AXIS)
    water_box = load_waterbox(Default.WATER_BOX['boxtyp20'])

    return (param, top, link, topology_protonation, trans_vec, tensor, scatter,
            axis, water_box)


def generate_topology(input_pdb, course_path, recipe_str, defaults,
                      protonation=None):
    """Generate a HADDOCK topology file from input_pdb"""
    general_param = load_recipe_params(defaults)

    param, top, link, topology_protonation, \
        trans_vec, tensor, scatter, \
        axis, water_box = get_topology_header(protonation)

    abs_path = input_pdb.resolve().parent.absolute()
    output_pdb_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock{input_pdb.suffix}')
    output_psf_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock.{Format.TOPOLOGY}')
    output = prepare_output(output_psf_filename, output_pdb_filename)

    input_str = prepare_input(str(input_pdb.resolve().absolute()), course_path)

    inp = general_param + param + top + input_str + output + link \
        + topology_protonation + trans_vec + tensor + scatter + axis \
        + water_box + recipe_str

    output_inp_filename = abs_path / f'{input_pdb.stem}.{Format.CNS_INPUT}'
    with open(output_inp_filename, 'w') as output_handler:
        output_handler.write(inp)

    return output_inp_filename


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
