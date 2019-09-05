import itertools
import math
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.functions import retrieve_output, calculate_haddock_score
from haddock.modules.worker.distribution import JobCreator
import config


class RigidBody:

	def __init__(self):
		pass

	@staticmethod
	def init(molecule_dic):
		jobs = JobCreator()
		recipe_gen = RecipeGenerator()

		combinations = []
		pdb_list = []
		psf_list = []
		for mol in molecule_dic:
			psf = f'begin/{mol}_1.psf'
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

		sampling = config.param_dic['it0']['sampling']
		if sampling % len(combinations):
			# increase sampling to make sure all models are sampled equally
			sampling = len(combinations) * math.ceil(sampling / len(combinations))
			print(f'+ WARNING: it0 sampling was increased to {sampling} to balance ensamble composition')
		else:
			# divide it between combinations
			sampling = int(sampling / len(combinations))

		# TODO: get folder from recipe parameters
		models = combinations * int(sampling/len(combinations))
		for i, pdb_list in enumerate(models):
			i_str = '0' * (6 - len(str(i))) + str(i)
			recipe = recipe_gen.generate(recipe_file=config.param_dic['it0']['recipe'],
			                             molecule_id=None,
			                             protonation_dic={},
			                             prefix_folder='it0',
			                             output_name=f'complex_it0_{i_str}',
			                             pdb=pdb_list,
			                             psf=psf_list)
			jobs.delegate(job_num=i,
			              job_id='complex_it0',
			              recipe_str=recipe,
			              job_folder='it0')

		return jobs

	@staticmethod
	def run(jobs):
		cns = CNS()
		cns.run(jobs)
		output = retrieve_output(jobs)

		return output

	@staticmethod
	def output(pdb_dic):
		pdb_list = [pdb_dic[e][0] for e in pdb_dic]
		file_list = []
		for pdb in pdb_list:
			hs = calculate_haddock_score(pdb, 'it0')
			file_list.append((pdb, hs))

		sorted_file_list = sorted(file_list, key=lambda x: x[1])

		# FIXME: dynamically assign the folder (?)
		folder = 'it0'

		with open(f'{folder}/file.list', 'w') as f:
			for e in sorted_file_list:
				pdb_name, haddock_score = e
				pdb_name = pdb_name.split('/')[1]
				tbw = f'"PREVIT:{pdb_name}  {{ {haddock_score} }}\n'
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

