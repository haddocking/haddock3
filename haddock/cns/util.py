"""CNS scripts util functions"""
from os import linesep
from haddock.pdbutil import PDBFactory
from haddock.mathutil import RandomNumberGenerator


def load_recipe_params(default_params):
    """Writes the values at the header section"""
    param_header = f'{linesep}! Parameters{linesep}'

    for param in default_params['params']:

        v = default_params['params'][param]

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


def prepare_input(pdb_input, course_path, psf_input=None):
    """Input of the CNS file.

    This section will be written for any recipe even if some CNS variables
    are not used, it should not be an issue.
    """
    input_str = f'{linesep}! Input structure{linesep}'

    if psf_input:
        if isinstance(psf_input, str):
            input_str += f'structure{linesep}'
            input_str += f'  @@{psf_input}{linesep}'
            input_str += f'end{linesep}'
        if isinstance(psf_input, list):
            input_str += f'structure{linesep}'
            for psf in psf_input:
                input_str += f'  @@{psf}{linesep}'
            input_str += f'end{linesep}'

    if isinstance(pdb_input, str):
        if psf_input:
            input_str += f'coor @@{pdb_input}{linesep}'

        # $file variable is still used by some CNS recipes, need refactoring!
        input_str += f'eval ($file=\"{pdb_input}\"){linesep}'

    if isinstance(pdb_input, (list, tuple)):
        for pdb in pdb_input:
            input_str += f'coor @@{pdb}{linesep}'

    segids, chains = PDBFactory.identify_chainseg(pdb_input)
    chainsegs = sorted(list(set(segids) | set(chains)))

    ncomponents = len(chainsegs)

    input_str += f'eval ($ncomponents={ncomponents}){linesep}'

    for i, segid in enumerate(chainsegs):
        input_str += f'eval ($prot_segid_{i+1}="{segid}"){linesep}'

    try:
        ambig_fname = list(course_path.glob('ambig.tbl'))[0]
        input_str += f'eval ($ambig_fname="{ambig_fname}"){linesep}'
    except IndexError:
        input_str += f'eval ($ambig_fname=""){linesep}'

    try:
        unambig_fname = list(course_path.glob('unambig.tbl'))[0]
        input_str += f'eval ($unambig_fname="{unambig_fname}"){linesep}'
    except IndexError:
        input_str += f'eval ($unambig_fname=""){linesep}'

    try:
        hbond_fname = list(course_path.glob('hbond.tbl'))[0]
        input_str += f'eval ($hbond_fname="{hbond_fname}"){linesep}'
    except IndexError:
        input_str += f'eval ($hbond_fname=""){linesep}'

    try:
        dihe_fname = list(course_path.glob('dihe.tbl'))[0]
        input_str += f'eval ($dihe_fname="{dihe_fname}"){linesep}'
    except IndexError:
        input_str += f'eval ($dihe_fname=""){linesep}'

    try:
        tensor_fname = list(course_path.glob('tensor.tbl'))[0]
        input_str += f'eval ($tensor_tbl="{tensor_fname}"){linesep}'
    except IndexError:
        input_str += f'eval ($tensor_fname=""){linesep}'

    rnd = RandomNumberGenerator()
    seed = rnd.randint(100, 999)
    input_str += f'eval ($seed={seed}){linesep}'

    return input_str
