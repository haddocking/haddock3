import itertools
import math
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.functions import retrieve_output
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
