import re
import os
from haddock.workflows.scoring.config import load_parameters

param_dic = load_parameters()
input_file_name = param_dic['input']['ensemble']


def load_structure(pdb_f):
    prot_dic = {}
    for line in open(pdb_f):
        if line.startswith(('ATOM', 'HETATM', 'ANISOU')):
            if 'HN' in line:
                new_line = line[:12] + ' H  ' + line[16:]
                line = new_line
            segid = line[72:76].strip()
            if segid in prot_dic:
                prot_dic[segid].append(line)
            else:
                prot_dic[segid] = [line]
    return prot_dic


def identify_chains(pdb):
    chain_l = []
    for l in open(pdb):
        if 'ATOM' in l[:4]:
            chain = l[21]
            if chain not in chain_l:
                chain_l.append(chain)
    return chain_l


def extract_md5():
    """ Read the whole scoring set and return a dictionary with model number and MD5 """
    regex = r"(\d*)\sMD5\s(.*)"
    md5_dictionary = {}
    with open(input_file_name) as f:
        data = f.readlines()
        for line in data:
            matches = re.finditer(regex, line, re.MULTILINE)
            if matches:
                for match in matches:
                    try:
                        model_number = int(match.group(1))
                        md5 = match.group(2)
                    except ValueError:
                        # this is not the match you are looking for
                        continue
                    md5_dictionary[int(model_number)] = md5
            if 'ATOM' in line[:4]:
                break
    f.close()

    return md5_dictionary


def split_models():
    # shamelessly copied and adapted from joao's pdb-tools (:

    """Splits the contents of the PDB file into new files, each containing a
    MODEL in the original file
    """
    wd = os.getcwd()
    if not os.path.isdir(f'{wd}/structures'):
        os.system(f'mkdir {wd}/structures')

    model_list = []
    # model_lines = []
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    with open(input_file_name) as f:
        for line in f.readlines():
            if line.startswith('MODEL'):
                model_no = int(line[10:14].strip())
                model_str = '0' * (6 - len(str(model_no))) + str(model_no)
                model_name = f'structures/{model_str}.pdb'
                model_list.append(model_name)
                fh = open(model_name, 'w')
                model_lines = []

            elif line.startswith('ENDMDL'):
                fh.write(''.join(model_lines))
                fh.close()

            elif line.startswith(records):
                model_lines.append(line)
    f.close()
    return model_list


def chain2segid(pdbf):
    """ """
    path = os.path.dirname(__file__).replace('pdb', 'src')
    os.system(f'{path}/pdb_chain-to-segid {pdbf} > oo')
    os.system(f'mv oo {pdbf}')


_to_remove = ['REMAR', 'CTERB', 'CTERA', 'NTERA', 'NTERB', 'CONECT']
_to_rename = {'WAT ': 'TIP3', 'HSD': 'HIS', 'HSE': 'HIS', 'HID': 'HIS', 'HIE': 'HIS',
              ' 0.00969': ' 0.00   '}


def sanitize(model_l):
    """Remove problematic portions of a PDB file."""

    for pdb in model_l:
        out_l = []
        # with open(output_file_name, 'w') as output_handler:
        with open(pdb) as input_handler:
            for line in input_handler:
                # Ignoring lines containing any tag from _to_remove
                if not any([tag in line for tag in _to_remove]):
                    for tag, new_tag in _to_rename.items():
                        line = line.replace(tag, new_tag)
                    out_l.append(line)

        out_l = [e for e in out_l if 'TER' not in e]
        out_l.append('END\n')

        with open(pdb, 'w') as output_handler:
            for line in out_l:
                output_handler.write(line)
        output_handler.close()

        # add segid, expect models to have CHAIN
        chain2segid(pdb)
