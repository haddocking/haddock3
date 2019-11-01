import itertools
import math
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import InputGenerator
from haddock.modules.functions import retrieve_output, calculate_haddock_score
from haddock.modules.worker.distribution import JobCreator


class RigidBody:

    def __init__(self):
        pass

    @staticmethod
    def init(recipe, molecule_dic, run_params):

        rigid_body_gen = InputGenerator(recipe_file=recipe,
                                        input_folder='rigid_body')

        jobs = JobCreator(job_id='complex_it0',
                          job_folder='rigid_body')

        combinations = []
        pdb_list = []
        psf_list = []
        for mol in molecule_dic:
            psf = f'topology/{mol}_1.psf'
            psf_list.append(psf)

            for e in molecule_dic[mol]:
                pdb = e[1]
                pdb_list.append(pdb)

        for comb in itertools.combinations(pdb_list, len(psf_list)):
            root_list = []
            for e in comb:
                root = e.split('/')[1].split('_')[0]
                root_list.append(root)
            if len(root_list) == len(set(root_list)):
                combinations.append(comb)

        sampling = run_params['stage']['rigid_body']['sampling']
        if sampling % len(combinations):
            # This will result in unbalanced sampling
            old_sampling = sampling
            sampling = len(combinations) * math.ceil(sampling / len(combinations))
            print(
                f'+ WARNING: it0 sampling was increased to {sampling} to balance ensamble composition, previously {old_sampling}')

        elif sampling != len(combinations):
            sampling = int(sampling / len(combinations))

        models = combinations * int(sampling / len(combinations))
        for i, pdb_list in enumerate(models):
            i_str = '0' * (6 - len(str(i))) + str(i)
            input_f = rigid_body_gen.generate(protonation_dic={},
                                              output_pdb=True,
                                              output_psf=False,
                                              input_pdb=pdb_list,
                                              input_psf=psf_list,
                                              output_fname=f'complex_it0_{i_str}')
            jobs.delegate(job_num=i,
                          input_file_str=input_f)

        return jobs

    @staticmethod
    def run(jobs):
        cns = CNS()
        cns.run(jobs)
        output = retrieve_output(jobs)
        return output

    @staticmethod
    def output(pdb_dic):
    # def output(pdb_dic, folder, stage):
        pdb_list = [pdb_dic[e][0] for e in pdb_dic]
        file_list = []
        for pdb in pdb_list:
            hs = calculate_haddock_score(pdb, 'it0')
            file_list.append((pdb, hs))

        sorted_file_list = sorted(file_list, key=lambda x: x[1])

        # FIXME: dynamically assign the folder (?)
        folder = 'rigid_body'

        with open(f'{folder}/file.list', 'w') as f:
            for e in sorted_file_list:
                pdb_name, haddock_score = e
                pdb_name = pdb_name.split('/')[1]
                tbw = f'"PREVIT:{pdb_name}  {{ {haddock_score:.4f} }}\n'
                f.write(tbw)
        f.close()

        with open(f'{folder}/file.nam', 'w') as f:
            for e in sorted_file_list:
                pdb_name, _ = e
                pdb_name = pdb_name.split('/')[1]
                tbw = f'{pdb_name}\n'
                f.write(tbw)
        f.close()

        with open(f'{folder}/result.dat', 'w') as f:
            for i, e in enumerate(sorted_file_list):
                pdb_name, haddock_score = e
                tbw = f'{i} {pdb_name} {haddock_score}\n'
                f.write(tbw)
        f.close()

        return [e[0] for e in sorted_file_list]
