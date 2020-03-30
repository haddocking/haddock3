import itertools
import math
import glob
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import InputGenerator
from haddock.modules.functions import retrieve_output, calculate_haddock_score, check_failures
from haddock.modules.worker.distribution import JobCreator


class SemiFlexible:

    def __init__(self):
        pass

    @staticmethod
    def init(recipe, models):
        semi_flex_gen = InputGenerator(recipe_file=recipe, input_folder='semi_flexible')

        jobs = JobCreator(job_id='complex_it1', job_folder='semi_flexible')

        psf_list = glob.glob('topology/*psf')
        for i, input_complex in enumerate(models):
            i_str = '0' * (6 - len(str(i))) + str(i)
            input_f = semi_flex_gen.generate(protonation_dic={},
                                             output_pdb=True,
                                             output_psf=False,
                                             input_pdb=input_complex,
                                             input_psf=psf_list,
                                             output_fname=f'complex_it1_{i_str}')
            jobs.delegate(job_num=i,
                          input_file_str=input_f)

        return jobs

    @staticmethod
    def run(jobs, retry_counter=5):
        cns = CNS()
        cns.run(jobs)
        failed_jobs, fail_check = check_failures(jobs)
        if fail_check:
            while retry_counter != 0:
                cns = CNS()
                cns.run(failed_jobs, retry=True)
                failed_jobs, fail_check = check_failures(failed_jobs)
                if not fail_check:
                    retry_counter = 0
                else:
                    retry_counter -= 1
            if fail_check and retry_counter == 0:
                print('+ ERROR: Jobs failed in it1')
                exit()

        output = retrieve_output(jobs)
        return output

    @staticmethod
    def output(pdb_dic):
    # def output(pdb_dic, folder, stage):
        pdb_list = [pdb_dic[e][0] for e in pdb_dic]
        file_list = []
        for pdb in pdb_list:
            hs = calculate_haddock_score(pdb, 'it1')
            file_list.append((pdb, hs))

        sorted_file_list = sorted(file_list, key=lambda x: x[1])

        # FIXME: dynamically assign the folder (?)
        folder = 'semi_flexible'

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
