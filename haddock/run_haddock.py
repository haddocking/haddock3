from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.structure.utils import PDB
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *
import config


def generate_topology():
	p = PDB()
	recipe_gen = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()
	input_list = list(config.param_dic['input']['molecules'].values())
	_ = p.sanitize(input_list)

	for i, mol in enumerate(config.param_dic['input']['molecules']):
		input_strc = config.param_dic['input']['molecules'][mol]
		recipe = recipe_gen.generate(recipe_file=config.param_dic['topology_recipe'],
		                             molecule_id=mol,
		                             protonation_dic={},
		                             prefix_folder='begin/',
		                             out_suffix='')
		jobs.delegate(job_num=i+1, job_id='generate', recipe_str=recipe, input_model=input_strc)
	cns.run(jobs)

	# retrieve OUTPUT somehow
	output = retrieve_output(jobs)
	pass


if __name__ == '__main__':
	init()
	generate_topology()
