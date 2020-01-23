import os
import itertools
import subprocess
from haddock.modules.structure.utils import PDB
from haddock.utils.files import get_full_path


# TODO: Implement parallelism?
def fastcontact(pdb_f, fastcontact_exe):

    felecv = .0
    fdesolv = .0
    freee = .0

    prot_dic = PDB.load_structure(pdb_f)

    for segid_x, segid_y in itertools.combinations(prot_dic, 2):

        prot_x = '{}.pdb'.format(segid_x)
        prot_x_str = ''.join(prot_dic[segid_x])
        with open(prot_x, 'w') as fh_prot_x:
            fh_prot_x.write(prot_x_str)
        fh_prot_x.close()

        prot_y = '{}.pdb'.format(segid_y)
        prot_y_str = ''.join(prot_dic[segid_y])
        with open(prot_y, 'w') as fh_prot_y:
            fh_prot_y.write(prot_y_str)
        fh_prot_y.close()

        p = subprocess.Popen([fastcontact_exe, prot_x, prot_y], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.communicate()

        returncode = p.wait()

        if returncode != 0:
            print('+ ERROR: Something wrong with fastcontact')
            exit()

        i, j, k = map(float, str(out[0]).split('\\n')[1].split())
        felecv += i
        fdesolv += j
        freee += k

        os.remove(prot_x)
        os.remove(prot_y)

    return felecv, fdesolv
