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
	# print('-> Running fastcontact')
	# ana.run_fastcontact()
	# print('-> Running dfire')
	# ana.run_dfire()
	print('-> Creating contact matrix')
	ana.calculate_contact_matrix()
	print('-> Clustering')
	ana.cluster()

	# ana.calculate_rmsd()
	print(ana)
	pass

	# 4. Analyze
	#   4.1 Calculate HADDOCK Score
	#   4.2 Retrieve contacts
	#       4.2.1 Calculate FCC
	#       4.2.2 Cluster
	#   4.2 Run third-party software
	#       4.2 Profit (i-RMSD/l-RMSD/i-l-RMSD)
	#       4.3 Fastcontact
	#       4.4 Dfire
	#       4.5 Molprobity
	#       4.6 DockQ

	# 5. Output


if __name__ == '__main__':
	# run_scoring()
	run_analysis()