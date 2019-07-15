# convert the analysis scripts to python
import os
from haddock.workflows.scoring.analysis.ana import Ana
from haddock.workflows.scoring.cns.engine import CNS
from haddock.workflows.scoring.cns.input import RecipeGenerator
from haddock.workflows.scoring.pdb.utils import PDB
from haddock.workflows.scoring.worker.distribution import JobCreator
import haddock.workflows.scoring.config as config
# from haddock.workflows.scoring.config import load_parameters


def run_scoring(ensamble_f):

	# param_dic = load_parameters()
	# input_f = param_dic['input']['ensemble']
	# md5_dictionary = extract_md5()
	pdb = PDB()
	recipe = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()

	# model_list = split_models()
	# sanitize(model_list)

	# Prepare PDBs
	pdb.prepare(ensamble_f)

	recipe.initialize()
	recipe.load()  # load the recipe defined in json

	jobs.delegate(recipe, model_list)

	# cns.run(jobs)


def run_analysis():

	ana = Ana()

	ana.retrieve_structures()

	ana.extract_energies()

	ana.calculate_haddock_score()

	# ana.cluster(cutoff=0.75, threshold=1)
	# ana.cluster(cutoff=0.75, threshold=2)
	ana.cluster(cutoff=0.75, threshold=4)

	# ana.cluster(cutoff=0.65, threshold=1)
	# ana.cluster(cutoff=0.65, threshold=2)
	ana.cluster(cutoff=0.65, threshold=4)

	ana.run_fastcontact()

	# ana.run_dfire()

	# ana.run_dockq()

	ana.output()


if __name__ == '__main__':
	input_f = config.param_dic['input']['ensemble']
	run_scoring(input_f)
	# run_analysis()