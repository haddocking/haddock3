import multiprocessing
import os
import subprocess
import toml
from haddock.modules.functions import load_ini

ini = load_ini('haddock3.ini')


class Contacts:

    def __init__(self):
        # self.nproc = config.param_dic['input']['nproc']
        self.con_exec = ini.get('scripts', 'contacts_exe')
        self.run_param = toml.load('data/run.toml')
        self.nproc = self.run_param['execution_parameters']['nproc']
        self.arg_list = []
        self.contact_file_list = []

        # self.path = os.getcwd() + '/contacts'
        #
        # if not os.path.isdir(self.path):
        #     os.system(f'mkdir {self.path}')

    def calculate_contacts(self, pdb_list, cutoff):

        for input_pdb in pdb_list:
            # model_id = input_pdb.split('/')[-1].split('_')[0]
            # output_name = f'{self.path}/{model_id}.con'
            output_name = input_pdb.replace('.pdb', '.con')
            self.contact_file_list.append(output_name)

            if not os.path.isfile(output_name):
                self.arg_list.append((self.con_exec, input_pdb, output_name, cutoff))

        with multiprocessing.Pool(processes=self.nproc) as pool:
            _ = pool.starmap(self.execute, self.arg_list)

        return self.contact_file_list

    @staticmethod
    def execute(contact_exec, input_f, output_f, d_cutoff):
        with open(output_f, 'w+') as out:
            _ = subprocess.run([contact_exec, input_f, str(d_cutoff)], stdout=out)
        out.close()
