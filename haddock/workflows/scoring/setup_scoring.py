import sys
import toml
import argparse
from datetime import datetime
from haddock.modules.functions import *
from haddock.modules.setup import Setup


def greeting():
    now = datetime.now().replace(second=0, microsecond=0)
    python_version = sys.version
    return (f'''##############################################
#                                            #
#         Starting HADDOCK v3.0beta1         #
#                                            #
#             EXPERIMENTAL BUILD             #
#                                            #
#              Scoring Workflow              #
#                                            #
##############################################

Starting HADDOCK on {now}

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


def pre_process(raw_molecule_dic):

    p = PDB()

    # 1. Add or assign ChainID/SegID
    chainseg_dic = p.fix_chainseg(raw_molecule_dic)

    # 2. Check if it is an ensemble and split it
    ensemble_dic = p.treat_ensemble(chainseg_dic)

    # 3. Clean problematic parts
    clean_molecule_dic = p.sanitize(ensemble_dic)

    # 4. Match chains
    # matched_molecule_dic = p.match(clean_molecule_dic)

    return clean_molecule_dic


if __name__ == '__main__':

    print(greeting())

    parser = argparse.ArgumentParser(description='Setup your HADDOCK Scoring run')
    parser.add_argument("run_file", help="The run file containing the parameters of your scoring run (.toml)")
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_usage()
        exit()

    setup_dictionary = toml.load(sys.argv[1])

    greeting()

    s = Setup(setup_dictionary)
    s.prepare_folders()
    s.configure_recipes()

    processed_molecules = pre_process(setup_dictionary['molecules'])

    adieu()
