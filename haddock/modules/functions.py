import copy
import re
import random
import configparser
import os
import glob
from haddock.modules.structure.utils import PDB
from haddock.utils.files import get_full_path


def treat_ensemble(pdb_dic):
    check = False
    new_d = {}
    for mol in pdb_dic:
        new_d[mol] = []
        pdb = pdb_dic[mol]
        with open(pdb) as f:
            reversed_list = reversed(f.readlines())
            for line in reversed_list:
                if 'ENDMDL' in line:
                    check = True
                    break
            if check:
                split_models = PDB.split_models(pdb)
                for model in split_models:
                    new_d[mol].append(model)
                continue
            else:
                new_d[mol].append(pdb)
        f.close()
    return new_d


def molecularize(dic):
    out_dic = {}
    for e in dic:
        try:
            root = dic[e][0].split('/')[1].split('_')[0]
        except IndexError:
            print(f'+ WARNING: No output generated for job {e}')
            continue
        try:
            _ = out_dic[root]
        except:
            out_dic[root] = []

        out_dic[root].append(tuple(dic[e]))

    return out_dic

def check_failures(jobs):
    fjobs = copy.deepcopy(jobs)
    failed_jobs = []
    for j in jobs.dic:
        inp, out = jobs.dic[j]
        complete_check = None
        for i, line in enumerate(reversed(list(open(out)))):
            if 'CNSsolve>stop' in line:
                complete_check = True
            if not complete_check and 'ABORT' in line:
                failed_jobs.append(j)
            if i == 50:
                break
    # create deletion list
    if failed_jobs:
        deletion_list = [k for k in jobs.dic.keys() if k not in failed_jobs]
        for delk in deletion_list:
            del fjobs.dic[delk]
        return fjobs, True
    else:
        return None, False


def retrieve_output(jobs):
    """" Read the output file and look for < OUTPUT: > to identify if the job was sucessful """
    output_dic = {}
    for j in jobs.dic:
        inp, out = jobs.dic[j]
        output_dic[j] = []
        complete_check = None
        for i, line in enumerate(reversed(list(open(out)))):
            if 'OUTPUT:' in line and 'CNS' not in line:
                v = line.strip().split(':')[1][1:]
                output_dic[j].append(v)
            if 'CNSsolve>stop' in line:
                # passed ok
                complete_check = True
            if not complete_check and 'ABORT' in line:
                # job has failed!
                print(f'+ ERROR: Job {inp} has failed, please check the output {out}')
                exit()
            if i == 50:
                break
    return output_dic


def extract_energies(pdb_file):
    """ Extract energies from the header of the PDB file, according to HADDOCK formatting """
    vdw = .0
    elec = .0
    desolv = .0
    air = .0
    bsa = .0

    vdw_elec_air_regex = r"\s(\-?\d*\.?\d{1,}(E-\d{1,}|)|0\b|\$\w{1,})"
    desolv_regex = r"(\-?\d*\.?\d*)$"
    bsa_regex = r"(\-?\d*\.?\d*)$"

    f = open(pdb_file, 'r')
    for line in f:

        if 'REMARK energies' in line:
            # dirty account for 8.754077E-02 notation and the eventual $DANI or $NOE
            energy_v = re.findall(vdw_elec_air_regex, line)

            temp_v = []
            for v in energy_v:
                v = v[0]
                try:
                    v = float(v)
                except ValueError:
                    v = .0
                temp_v.append(v)
            energy_v = temp_v
            # WARNING: This is highly dependent on the values printed by the recipe, refer to print_coorheader.cns
            total, bonds, angles, improper, dihe, vdw, elec, air, cdih, coup, rdcs, vean, dani, xpcs, rg = energy_v
            vdw = float(vdw)
            elec = float(elec)
            air = float(air)

        if 'REMARK Desolvation' in line:
            desolv = float(re.findall(desolv_regex, line)[0])

        if 'REMARK buried surface area' in line:
            bsa = float(re.findall(bsa_regex, line)[0])
            break

    f.close()

    return vdw, elec, desolv, air, bsa


def calculate_haddock_score(pdb_file, stage):
    """ Calculate the HADDOCK Score of the PDB file using its appropriate weight """

    weight_dic = dict(it0=(0.01, 1.0, 1.0, 0.01, 0.01),
                      it1=(1.0, 1.0, 1.0, 0.1, 0.01),
                      itw=(1.0, 0.2, 1.0, 0.1, 0.0))

    vdw, elec, desolv, air, bsa = extract_energies(pdb_file)
    w = weight_dic[stage]
    haddock_score = w[0] * vdw + w[1] * elec + w[2] * desolv + w[3] * air - w[4] * bsa

    return haddock_score


def bye():
    salutations = ['Tot ziens!', 'Good bye!', 'Até logo!', 'Ciao!', 'Au revoir mille tonnerres!', 'Adéu-siau!', 'Agur!', 'Dovidenia!']
    return '\n'.join(random.sample(salutations, k=3))


def get_begin_molecules(folder):
    # {'mol1': ['run-debug/data/mol1_1.pdb'], 'mol2': ['run-debug/data/mol2_1.pdb']}
    mol_dic = {}
    for pdb in glob.glob(f'{folder}/*pdb'):
        root = pdb.split('/')[1].split('_')[0]
        try:
            mol_dic[root].append(pdb)
        except KeyError:
            mol_dic[root] = [pdb]

    return mol_dic


def load_ini(ini_file):
    etc_folder = get_full_path('haddock', 'etc')
    config_file = os.path.join(etc_folder, ini_file)

    ini = configparser.ConfigParser(os.environ)
    ini.read(config_file, encoding='utf-8')

    return ini
