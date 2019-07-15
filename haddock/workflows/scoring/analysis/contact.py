import multiprocessing
import os
import subprocess
import haddock.workflows.scoring.config as config


class Contacts:

    def __init__(self):
        self.nproc = config.param_dic['input']['nproc']
        self.con_exec = config.param_dic['third-party']['contacts_exe']
        self.arg_list = []
        self.contact_file_list = []

        if not os.path.isdir('contacts'):
            os.system('mkdir contacts')

    def calculate_contacts(self, pdb_list, cutoff):

        for input_pdb in pdb_list:
            model_name = input_pdb.split('_')[0].split('/')[1]
            output_name = f'contacts/{model_name}.con'
            self.contact_file_list.append(output_name)

            if not os.path.isfile(output_name):
                self.arg_list.append((self.con_exec, input_pdb, output_name, cutoff))

        with multiprocessing.Pool(processes=self.nproc) as pool:
            _ = pool.starmap(self.execute, self.arg_list)

        # return self.contact_file_list

    @staticmethod
    def execute(contact_exec, input_f, output_f, d_cutoff):
        with open(output_f, 'w+') as out:
            _ = subprocess.run([contact_exec, input_f, str(d_cutoff)], stdout=out)
        out.close()
