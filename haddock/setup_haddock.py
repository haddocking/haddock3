# setup the simulation
import sys
import toml
from datetime import datetime
from haddock.modules.functions import *
from haddock.modules.setup import Setup


def greeting():
	# start = datetime.now().replace(second=0, microsecond=0)
	start = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	python_version = sys.version
	print(f'''##############################################
#                                            #
#           Setup HADDOCK v3.0beta1          #
#                                            #
#             EXPERIMENTAL BUILD             #
#                                            #                               
##############################################

 Starting HADDOCK on {start}

 HADDOCK version: 3.0 beta 1
 Python {python_version}
''')


def adieu():
	end = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	salut = bye()
	print(f'''
 Your HADDOCK Run: {s.run_dir} has been correctly setup
 
 Finished at {end}
 
{salut}
''')


def pre_process(molecule_dic):

	p = PDB()

	treated_dic = p.treat_ensemble(molecule_dic)

	for mol in treated_dic:
		pdb_list = treated_dic[mol]
		p.sanitize(pdb_list)

	return treated_dic


if __name__ == '__main__':

	setup_dictionary = toml.load(sys.argv[1])

	greeting()

	s = Setup(setup_dictionary)
	s.prepare_folders()
	s.configure_recipes()

	_ = pre_process(setup_dictionary['molecules'])

	adieu()


