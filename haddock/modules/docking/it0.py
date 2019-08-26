import itertools
import math

from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.worker.distribution import JobCreator
import config


class RigidBody:

	def __init__(self, mol_dic):
		self.molecule_dic = mol_dic
		self.psf_list = []
		self.combinations = []

	def init(self):
		pdb_list = []
		for mol in self.molecule_dic:
			psf = f'begin/{mol}_1.psf'
			self.psf_list.append(psf)

			for e in self.molecule_dic[mol]:
				pdb = e[1]
				pdb_list.append(pdb)

		for comb in itertools.combinations(pdb_list, len(self.psf_list)):
			root_list = []
			for e in comb:
				root = e.split('/')[1].split('_')[0]
				root_list.append(root)
			if len(root_list) == len(set(root_list)):
				self.combinations.append(comb)

		sampling = config.param_dic['it0']['sampling']
		if sampling % len(self.combinations):
			# increase sampling to make sure all models are sampled equally
			sampling = len(self.combinations) * math.ceil(sampling/len(self.combinations))
			print(f'+ WARNING: it0 sampling was increased to {sampling} to balance ensamble composition')
		else:
			# divide it between combinations
			sampling = int(sampling / len(self.combinations))

		recipe_gen = RecipeGenerator()
		jobs = JobCreator()

		models = self.combinations * sampling
		for i, pdb_list in enumerate(models):
			i_str = '0' * (6 - len(str(i))) + str(i)
			recipe = recipe_gen.generate(recipe_file=config.param_dic['it0']['recipe'],
			                             molecule_id=None,
			                             protonation_dic={},
			                             prefix_folder='it0',
			                             output_name=f'complex_{i_str}')
			jobs.delegate_it0(job_num=i,
			                  job_id='complex_it0',
			                  recipe_str=recipe,
			                  input_psf=self.psf_list,
			                  input_pdb=pdb_list)
			# print(jobs)

	def run(self):
		pass

	def retrieve(self):
		pass
