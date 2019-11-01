import os
import string
from itertools import chain
# from haddock.modules.structure.reduce import analyze_protonation_state
from utils.files import get_full_path

src_path = get_full_path('haddock', 'src')


class PDB:

    def __init__(self):
        self.to_remove = ['REMAR', 'CTERB', 'CTERA', 'NTERA', 'NTERB', 'CONECT']
        self.to_rename = {'WAT ': 'TIP3', 'HSD': 'HIS', 'HSE': 'HIS', 'HID': 'HIS', 'HIE': 'HIS',
                          ' 0.00969': ' 0.00   '}
        self.model_list = []

    # @staticmethod
    # def get_protonation(pdbf):
    # 	return analyze_protonation_state(pdbf)

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

    # @staticmethod
    # def load_structure(pdb_f):
    #     """ Load a PDB structure into a dictionary """
    #     prot_dic = {}
    #     with open(pdb_f) as f:
    #         for line in f.readlines():
    #             if line.startswith(('ATOM', 'HETATM', 'ANISOU')):
    #                 if 'HN' in line:
    #                     new_line = line[:12] + ' H  ' + line[16:]
    #                     line = new_line
    #                 # segid = line[72:76].strip()
    #                 chainid = line[21]
    #                 if chain in prot_dic:
    #                     prot_dic[chainid].append(line)
    #                 else:
    #                     prot_dic[chainid] = [line]
    #     f.close()
    #     return prot_dic

    # @staticmethod
    # def identify_chains(pdb):
    #     """ Read PDB structure and return chainIDs """
    #     chain_l = []
    #     with open(pdb) as f:
    #         for l in f.readlines():
    #             if 'ATOM' in l[:4]:
    #                 chain_id = l[21]
    #                 if chain_id not in chain_l:
    #                     chain_l.append(chain_id)
    #     f.close()
    #     return chain_l

    # @staticmethod
    # def extract_md5(pdb_ens):
    #     """ Read the whole scoring set and return a dictionary with model number and MD5 """
    #     regex = r"(\d*)\sMD5\s(.*)"
    #     md5_dictionary = {}
    #     with open(pdb_ens) as f:
    #         data = f.readlines()
    #         for line in data:
    #             matches = re.finditer(regex, line, re.MULTILINE)
    #             if matches:
    #                 for match in matches:
    #                     try:
    #                         model_number = int(match.group(1))
    #                         md5 = match.group(2)
    #                     except ValueError:
    #                         # this is not the match you are looking for
    #                         continue
    #                     md5_dictionary[int(model_number)] = md5
    #             if 'ATOM' in line[:4]:
    #                 break
    #     f.close()
    #
    #     return md5_dictionary

    @staticmethod
    def split_models(ensamble_f):
        # """" adapted """"  from pdb-tools (:
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
    def add_chainseg(pdbf, ident, overwrite=True):
        """" Add chainID and segID to a PDB structure """
        new_pdb = pdbf + '_'
        with open(new_pdb, 'w') as out_fh:
            with open(pdbf) as in_fh:
                for line in in_fh.readlines():
                    if 'ATOM' in line[:4]:
                        c = line[21].strip()[:1]
                        s = line[72:76].strip()[:1]

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
                    if 'ATOM' in line[:4]:
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
        # Get just the pdbs, the data structure does not matter
        pdb_list = list(chain.from_iterable([pdb_dic[e] for e in pdb_dic]))
        for pdb in pdb_list:
            out_l = []
            with open(pdb) as input_handler:
                for line in input_handler:
                    # Ignoring lines containing any tag from _to_remove
                    if not any([tag in line for tag in self.to_remove]):
                        for tag, new_tag in self.to_rename.items():
                            line = line.replace(tag, new_tag)
                        out_l.append(line)

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

        return clean_model_list
