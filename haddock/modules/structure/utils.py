import os
import string
import configparser
from itertools import chain
from haddock.utils.files import get_full_path

# from haddock.modules.structure.reduce import analyze_protonation_state

etc_folder = get_full_path('haddock', 'etc')
config_file = os.path.join(etc_folder, 'haddock3.ini')

ini = configparser.ConfigParser(os.environ)
ini.read(config_file, encoding='utf-8')

src_path = get_full_path('haddock', 'src')


def histogram(items):
    for c, n in items:
        if n/10 > 1:
            output = int(n/10) * '-'
            print(f' {c}\t{output} ({n})')


def pad_line(line):
    # borrowed from from pdb-tools (:
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character


def identify_known(top_f):
    """ Read the topology file and identify which residues are known """
    known_l = []
    with open(top_f) as f:
        for l in f.readlines():
            if 'resi' in l[:4].casefold():
                res = l.split()[1]
                known_l.append(res)
    return known_l


class PDB:

    def __init__(self):
        self.to_remove = ['REMAR', 'CTERB', 'CTERA', 'NTERA', 'NTERB', 'CONECT']
        self.to_keep = identify_known(ini.get('topology', 'top_file'))
        self.to_rename = {'WAT ': 'TIP3', 'HSD': 'HIS', 'HSE': 'HIS', 'HID': 'HIS', 'HIE': 'HIS',
                          ' 0.00969': ' 0.00   '}
        self.model_list = []

    def treat_ensemble(self, pdb_dic):
        """" Separate a multimodel PDB file and add it to data structure """
        check = False
        new_d = {}
        for mol in pdb_dic:
            new_d[mol] = []
            pdb = pdb_dic[mol]
            with open(pdb) as f:
                # Reverse it to avoid reading the whole file
                reversed_list = reversed(f.readlines())
                for line in reversed_list:
                    if 'ENDMDL' in line:
                        check = True
                        break
                if check:
                    split_models = self.split_models(pdb)
                    for model in split_models:
                        new_d[mol].append(model)
                    continue
                else:
                    new_d[mol].append(pdb)

            f.close()
        return new_d

    @staticmethod
    def load_structure(pdb_f):
        """ Load a PDB structure into a dictionary """
        prot_dic = {}
        with open(pdb_f) as f:
            for line in f.readlines():
                if line.startswith(('ATOM', 'HETATM', 'ANISOU')):
                    if 'HN' in line:
                        new_line = line[:12] + ' H  ' + line[16:]
                        line = new_line
                    ident = line[72:76].strip()  # segid
                    # ident = line[21]  # chain
                    if ident in prot_dic:
                        prot_dic[ident].append(line)
                    else:
                        prot_dic[ident] = [line]
        f.close()
        return prot_dic

    @staticmethod
    def identify_chains(pdb):
        """ Read PDB structure and return chainIDs """
        chain_l = []
        with open(pdb) as f:
            for l in f.readlines():
                if 'ATOM' in l[:4]:
                    chain_id = l[21]
                    if chain_id not in chain_l:
                        chain_l.append(chain_id)
        f.close()
        chain_l = list(set(chain_l))
        chain_l.sort()
        return chain_l

    @staticmethod
    def identify_segids(pdb):
        """ Read PDB structure and return chainIDs """
        segid_l = []
        with open(pdb) as f:
            for l in f.readlines():
                if 'ATOM' in l[:4]:
                    segid = l[72:76].strip()[:1]
                    if segid not in segid_l:
                        segid_l.append(segid)
        f.close()
        segid_l = list(set(segid_l))
        segid_l.sort()
        return segid_l


    @staticmethod
    def split_models(ensamble_f):
        # borrowed from from pdb-tools (:
        model_list = []
        records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
        path = '/'.join(ensamble_f.split('/')[:-1])
        with open(ensamble_f) as f:
            for line in f.readlines():
                if line.startswith('MODEL'):
                    model_no = int(line.split()[1])
                    name = ensamble_f.split('/')[-1].split('_')[0]
                    model_name = f'{path}/{name}_{model_no}.pdb'
                    model_list.append(model_name)
                    fh = open(model_name, 'w')
                    model_lines = []

                elif line.startswith('ENDMDL'):
                    fh.write(''.join(model_lines))
                    fh.close()

                elif line.startswith(records):
                    model_lines.append(line)
        f.close()
        # fh.close()
        return model_list

    @staticmethod
    def fix_id(pdbf, priority='seg', overwrite=True):
        """ Replaces the chainID with segID or vice-versa, based on the priority param """
        # if priority = seg, replace chainid with segid
        # if priority = chain, replace segid with chaid
        new_pdb = pdbf + '_'
        with open(new_pdb, 'w') as out_fh:
            with open(pdbf) as in_fh:
                for line in in_fh.readlines():
                    if 'ATOM' in line[:4]:
                        line = pad_line(line)
                        chainid = line[21].strip()[:1]  # chainID
                        segid = line[72:76].strip()[:1]  # segID
                        if priority == 'seg':
                            if not segid:
                                print('\n+ ERROR: Fix id failed, no segID')
                                exit()
                            line = line[:21] + segid + line[22:]
                        elif priority == 'chain':
                            if not chainid:
                                print('\n+ ERROR: Fix id failed, no chaiID')
                                exit()
                            line = line[:72] + chainid.ljust(4) + line[76:]
                            if not '\n' in line:
                                line += '\n'

                        else:
                            # option not supported
                            exit()
                    out_fh.write(line)

            in_fh.close()
        out_fh.close()

        if overwrite:
            os.rename(new_pdb, pdbf)
            return pdbf
        else:
            return new_pdb

    @staticmethod
    def add_chainseg(pdbf, ident, overwrite=True):
        """" Add ONE chainID and segID to a PDB structure """
        new_pdb = pdbf + '_'
        with open(new_pdb, 'w') as out_fh:
            with open(pdbf) as in_fh:
                for line in in_fh.readlines():
                    if 'ATOM' in line[:4]:
                        c = line[21].strip()[:1]  # chainID
                        s = line[72:76].strip()[:1]  # segID

                        if c != ident:
                            line = line[:21] + ident + line[22:]

                        if s != ident:
                            line = line[:72] + ident.ljust(4) + line[76:]

                    out_fh.write(line)
            in_fh.close()
        out_fh.close()

        if overwrite:
            os.rename(new_pdb, pdbf)
            return pdbf
        else:
            return new_pdb

    @staticmethod
    def identify_chainseg(pdb_inp):
        """" Read PDB structure and return segID OR chainID """
        if type(pdb_inp) == str:
            pdb_inp = [pdb_inp]
        segid_l = []
        for pdb in pdb_inp:
            with open(pdb) as fh:
                for line in fh.readlines():
                    if line.startswith('ATOM'):
                        line = pad_line(line)
                        segid = line[72:76].strip()[:1]
                        chainid = line[21].strip()[:1]
                        if segid:
                            segid_l.append(segid)
                            # return segid
                        elif chainid:
                            segid_l.append(segid)
                            # return chainid

        segid_l = list(set(segid_l))
        # WARNING: This expects chains to be sequential!
        # chain1 = A, chain2 = B or chain1 = D, chain = E, etc
        segid_l.sort()
        return segid_l

    def fix_chainseg(self, pdb_dic):
        """" Separate by molecule id and assign segids accordingly """
        chainseg_check = []
        segid_dic = dict(
                [(int(e.split('mol')[1]), {'mol': None, 'segid': None}) for e in pdb_dic if 'mol' in e])

        for e in pdb_dic:
            if 'mol' in e:
                ident = int(e.split('mol')[1])
                molecule = pdb_dic[e]
                segid_dic[ident]['mol'] = molecule
            if 'segid' in e:
                ident = int(e.split('segid')[1])
                segid = pdb_dic[e]
                segid_dic[ident]['segid'] = segid

        # If segid has not been assigned, check if PDB already has one
        for e in segid_dic:
            structure = segid_dic[e]['mol']
            custom_segid = segid_dic[e]['segid']

            if not custom_segid:
                # Segid not defined in setup, check if it is already present
                molecule_chainseg_list = self.identify_chainseg(structure)

                if molecule_chainseg_list:
                    # is ok
                    pass

                # if molecule_chainseg_list:
                #     # result = self.add_chainseg(structure, molecule_chainseg_list)
                #     chainseg_check.append(molecule_chainseg_list)

                if not molecule_chainseg_list:
                    # define sequentially
                    molecule_chainseg = string.ascii_uppercase[e - 1]
                    result = self.add_chainseg(structure, molecule_chainseg)
                    chainseg_check.append(result)

            if custom_segid:
                result = self.add_chainseg(structure, custom_segid)
                chainseg_check.append(result)

        if all(chainseg_check):
            # remove segid key from dictionary and return
            chainseg_dic = dict([(e, pdb_dic[e]) for e in pdb_dic if 'mol' in e])
            return chainseg_dic
        else:
            print('+ ERROR: Could not edit ChainID/SegID')
            exit()

    def sanitize(self, pdb_dic, overwrite=True):
        """ Remove problematic portions of a PDB file """
        clean_model_list = []
        removed_l = []
        # Get just the pdbs, the data structure does not matter
        pdb_list = list(chain.from_iterable([pdb_dic[e] for e in pdb_dic]))
        for pdb in pdb_list:
            out_l = []
            with open(pdb) as input_handler:
                for line in input_handler:
                    # Ignoring lines containing any tag from _to_remove and what is in to_keep
                    if not any([tag in line for tag in self.to_remove]):
                        for tag, new_tag in self.to_rename.items():
                            line = line.replace(tag, new_tag)
                        # check if this residue is known
                        res = line[17:20].strip()
                        if res and res in self.to_keep:
                            out_l.append(line)
                        elif not res:
                            out_l.append(line)
                        else:
                            removed_l.append(res)

            out_l = [e for e in out_l if 'TER' not in e]
            out_l.append('END\n')

            sanitized_pdb = pdb + '_'
            with open(sanitized_pdb, 'w') as output_handler:
                for line in out_l:
                    output_handler.write(line)
            output_handler.close()

            if overwrite:
                os.rename(sanitized_pdb, pdb)
                clean_model_list.append(pdb)
            else:
                clean_model_list.append(sanitized_pdb)

        removed_l = list(set(removed_l))
        if removed_l:
            if len(removed_l) == 1:
                print('\n+ WARNING: The following residue was removed because it is present in the input but not '
                      'in the topology file')
            else:
                print('\n+ WARNING: The following residues were removed because they are present in the input but '
                      'not in the topology file')
            removed_str = ' ,'.join(removed_l)
            print(f'++ {removed_str}')

        return clean_model_list

    @staticmethod
    def count_atoms(pdbf):
        counter = 0
        with open(pdbf) as fh:
            for l in fh.readlines():
                if l.startswith('ATOM'):
                    counter += 1
        return counter

    def organize_chains(self, pdb_dic):

        # Assume all models have CHAIN
        new_pdb_dic = {}
        for mol in pdb_dic:
            new_pdb_dic[mol] = []
            chain_count_dic = {}
            for pdb in pdb_dic[mol]:
                chains = self.identify_chains(pdb)
                size = self.count_atoms(pdb)
                if len(chains) > 1:
                    new_pdb_dic[mol].append(pdb)
                    chain_count_dic[pdb] = [size, len(chains), chains]
                else:
                    print(f'+ WARNING: {pdb} has {len(chains)} chain, removing')
                    os.remove(pdb)

            # count how many times each one appeared
            sizes = [chain_count_dic[e][0] for e in chain_count_dic]
            sizes_data = list(set([(a, sizes.count(a)) for a in sizes]))

            chains = [chain_count_dic[e][1] for e in chain_count_dic]
            chains_data = list(set([(a, chains.count(a)) for a in chains]))

            chain_names = list(chain.from_iterable([chain_count_dic[e][2] for e in chain_count_dic]))
            chain_names = list(set(chain_names))
            chain_names.sort()
            chain_names_str = ' '.join(chain_names)

            sizes_data.sort()
            chains_data.sort()

            print('\n+ Size distribution')
            histogram(sizes_data)

            print(f'\n+ Chain names found: {chain_names_str}')

        return new_pdb_dic

    @staticmethod
    def replace_chain(pdbf, old_chain, new_chain, overwrite=True):
        new_pdb = pdbf + '_'
        with open(new_pdb, 'w') as out_fh:
            with open(pdbf) as in_fh:
                for line in in_fh.readlines():
                    if line.startswith('ATOM'):
                        current_chain = line[21]
                        if current_chain == old_chain:
                            line = line[:21] + new_chain + line[22:]
                    out_fh.write(line)
            in_fh.close()
        out_fh.close()

        if overwrite:
            os.rename(new_pdb, pdbf)
            return pdbf
        else:
            return new_pdb

    @staticmethod
    def renumber(pdbf, renumber_dic, target_chain, overwrite=True):
        new_pdb = pdbf + '_'
        ignored_res = []
        with open(new_pdb, 'w') as out_fh:
            with open(pdbf) as in_fh:
                for line in in_fh.readlines():
                    if line.startswith('ATOM'):
                        current_res = int(line[22:26])
                        current_chain = line[21]
                        if current_chain == target_chain:
                            try:
                                new_res = renumber_dic[current_res]
                                line = line[:22] + '{:>4}'.format(new_res) + line[26:]
                                out_fh.write(line)
                            except:
                                # Residue not found in reference, IGNORE
                                ignored_res.append(current_res)
                        else:
                            out_fh.write(line)
                    else:
                        out_fh.write(line)
            in_fh.close()
        out_fh.close()

        if ignored_res:
            ignored_res_str = ', '.join(map(str, list(set(ignored_res))))
            # print(f'+ WARNING: {pdbf} Chain {target_chain} Res {ignored_res_str} not found in reference, discarded.')

        if overwrite:
            os.rename(new_pdb, pdbf)
            return pdbf
        else:
            return new_pdb

    @staticmethod
    def load_seq(pdb_f):
        aa_dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                  'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                  'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',  'DC': 'C',  'DA': 'A',  'DG': 'G',  'DT': 'T',
                  'ADE': 'A', 'THY': 'T', 'GUA': 'G', 'CYT': 'C'}
        seq_dic = {}
        for l in open(pdb_f):
            if 'ATOM' in l[:4]:
                segment_id = l[72:76].strip()
                resnum = int(l[22:26])
                resname = l[17:20].split()[0]
                try:
                    _ = seq_dic[segment_id]
                except KeyError:
                    seq_dic[segment_id] = {}
                try:
                    name = aa_dic[resname]
                except KeyError:
                    name = 'X'
                seq_dic[segment_id][resnum] = name
        return seq_dic
