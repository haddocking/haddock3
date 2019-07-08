# convert the analysis scripts to python
from haddock.workflows.scoring.analysis.ana import Ana
from haddock.workflows.scoring.cns.engine import CNS
from haddock.workflows.scoring.cns.input import RecipeGenerator
from haddock.workflows.scoring.pdb.utils import extract_md5, split_models, sanitize
from haddock.workflows.scoring.worker.distribution import JobCreator
# from haddock.workflows.scoring.config import load_parameters


def run_scoring():

	# param_dic = load_parameters()
	# input_f = param_dic['input']['ensemble']
	# md5_dictionary = extract_md5()
	model_list = split_models()
	model_list = model_list
	sanitize(model_list)

	recipe = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()

	recipe.initialize()
	recipe.load()  # load the recipe defined in json

	jobs.delegate(recipe, model_list)

	cns.run(jobs)


def run_analysis():

	ana = Ana()

	ana.retrieve_structures()

	ana.extract_energies()

	ana.calculate_haddock_score()

	ana.run_fastcontact()

	ana.run_dfire()

	ana.run_dockq()

	ana.cluster(cutoff=0.75, threshold=1)
	ana.cluster(cutoff=0.75, threshold=2)
	ana.cluster(cutoff=0.75, threshold=4)

	ana.cluster(cutoff=0.65, threshold=1)
	ana.cluster(cutoff=0.65, threshold=2)
	ana.cluster(cutoff=0.65, threshold=4)


if __name__ == '__main__':
	# run_scoring()
	run_analysis()