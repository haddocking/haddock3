import os
import subprocess
from haddock.workflows.scoring.config import load_parameters

param_dic = load_parameters()
contact_exe = param_dic['third-party']['contacts_exe']


def run_contacts(pdb, d_cutoff=5.0):

    model_name = pdb.split('_')[0].split('/')[1]
    output_name = f'contacts/{model_name}.con'

    if not os.path.isfile(output_name):
        with open(output_name, 'w+') as out:
            _ = subprocess.run([contact_exe, pdb, str(d_cutoff)], stdout=out)
        out.close()

    return output_name
