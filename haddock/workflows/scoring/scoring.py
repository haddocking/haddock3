# convert the analysis scripts to python
from haddock.workflows.scoring.analysis.ana import Ana
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.structure.utils import PDB
from haddock.modules.worker.distribution import JobCreator
import haddock.workflows.scoring.config as config


def run_scoring(ensamble_f):

	pdb = PDB()
	recipe_gen = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()

	# Prepare PDBs
	pdb.prepare(ensamble_f)

	# Generate the Recipe
	recipe = recipe_gen.generate()

	jobs.delegate(recipe, pdb.model_list)

	cns.run(jobs)


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