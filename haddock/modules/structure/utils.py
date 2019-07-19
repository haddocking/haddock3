import os

src_path = os.path.join(os.path.dirname(__file__), '../../src')


class PDB:

    def __init__(self):
        # wd = os.getcwd()
        # if not os.path.isdir(f'{wd}/structures'):
        #     os.system(f'mkdir {wd}/structures')

        self.to_remove = ['REMAR', 'CTERB', 'CTERA', 'NTERA', 'NTERB', 'CONECT']
        self.to_rename = {'WAT ': 'TIP3', 'HSD': 'HIS', 'HSE': 'HIS', 'HID': 'HIS', 'HIE': 'HIS',
                      ' 0.00969': ' 0.00   '}
        self.model_list = []

    def prepare(self, ensamble_f):
        split_model_list = self.split_models(ensamble_f)
        clean_model_list = self.sanitize(split_model_list)

        self.model_list = clean_model_list

    @staticmethod
    def load_structure(pdb_f):
        prot_dic = {}
        with open(pdb_f) as f:
            for line in f.readlines():
                if line.startswith(('ATOM', 'HETATM', 'ANISOU')):
                    if 'HN' in line:
                        new_line = line[:12] + ' H  ' + line[16:]
                        line = new_line
                    # segid = line[72:76].strip()
                    chain = line[21]
                    if chain in prot_dic:
                        prot_dic[chain].append(line)
                    else:
                        prot_dic[chain] = [line]
        f.close()
        return prot_dic

    @staticmethod
    def identify_chains(pdb):
        chain_l = []
        with open(pdb) as f:
            for l in f.readlines():
                if 'ATOM' in l[:4]:
                    chain = l[21]
                    if chain not in chain_l:
                        chain_l.append(chain)
        f.close()
        return chain_l

    # @staticmethod
    # def extract_md5():
    #     """ Read the whole scoring set and return a dictionary with model number and MD5 """
    #     regex = r"(\d*)\sMD5\s(.*)"
    #     md5_dictionary = {}
    #     with open(input_file_name) as f:
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
        # copied and adapted from pdb-tools (:

        """Splits the contents of the PDB file into new files, each containing a
        MODEL in the original file
        """
        path = '/'.join(ensamble_f.split('/')[:-1]) + '/structures'

        if not os.path.isdir(path):
            os.system(f'mkdir {path}')

        # model_lines = []
        model_list = []
        records = ('ATOM', 'HETATM', 'ANISOU', 'TER')

        with open(ensamble_f) as f:
            for line in f.readlines():
                if line.startswith('MODEL'):
                    model_no = int(line[10:14].strip())
                    model_str = '0' * (6 - len(str(model_no))) + str(model_no)
                    model_name = f'{path}/{model_str}.pdb'
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
    def chain2segid(pdbf):
        """ Save a backup of the original file and move chainID to segID in a new file """
        new_pdb = pdbf.replace('.pdb', '_+sid.pdb')
        # os.system(f'cp {pdbf} {ori}')
        os.system(f'{src_path}/pdb_chain-to-segid {pdbf} > oo')
        os.system(f'mv oo {new_pdb}')
        return new_pdb

    @staticmethod
    def segid2chain(pdbf):
        """ Save a backup of the original file and move segID to chainID in a new file """
        new_pdb = pdbf.replace('.pdb', '_+cid.pdb')
        # os.system(f'cp {pdbf} {ori}')
        os.system(f'{src_path}/pdb_chain-segid {pdbf} > oo')
        os.system(f'mv oo {new_pdb}')
        return new_pdb

    def sanitize(self, model_list):
        """ Remove problematic portions of a PDB file """
        clean_model_list = []
        for pdb in model_list:
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

            new_pdb = pdb.replace('.pdb', '_c.pdb')
            with open(new_pdb, 'w') as output_handler:
                for line in out_l:
                    output_handler.write(line)
            output_handler.close()

            # add segid, expect models to have CHAIN
            clean_pdb = self.chain2segid(new_pdb)
            clean_model_list.append(clean_pdb)

        return clean_model_list


