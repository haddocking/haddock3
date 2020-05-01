import os
import json
import toml
from haddock.utils.files import get_full_path
from haddock.modules.worker.recipe import RecipeComposer

etc_folder = get_full_path('haddock', 'etc')


class Setup:

    def __init__(self, setup_dic):
        self.setup_dic = setup_dic

        if not self.check_setup_dic(setup_dic):
            raise SetupError('+ Error: Please configure your setup file')

        self.protocol_path = get_full_path('haddock', 'protocols')

        with open(f'{etc_folder}/default.json', 'r') as fh:
            self.default_recipes = json.load(fh)
        fh.close()

        run_id = None
        try:
            run_id = self.setup_dic['identifier']['run']
        except KeyError:
            raise SetupError('+ ERROR: Run identifier not provided')

        for mol in self.setup_dic['molecules']:
            target_mol = self.setup_dic['molecules'][mol]
            if not os.path.isfile(target_mol):
                raise SetupError(f'+ ERROR: Molecule {target_mol} not found')

        self.run_dir = None
        if type(run_id) == int:
            self.run_dir = f"run{run_id}"
        elif type(run_id) == str:
            run_id = run_id.replace(' ', '-')
            self.run_dir = f"run-{run_id}"
        else:
            raise SetupError('+ ERROR: Run identifier can only be string or integer')

    @staticmethod
    def check_setup_dic(sdic):
        # Check if setup dic meets minimal requirements
        #  stage
        #  identifier
        #  molecules
        return True

    def prepare_folders(self):
        """ Create folder structure and copy/edit significant input files """

        if len(self.setup_dic['molecules']) >= 20:
            raise SetupError('+ ERROR: Too many molecules')

        # move input molecules to correct places
        if not os.path.isdir(self.run_dir):
            os.system(f'mkdir {self.run_dir}')
        else:
            print(f'+ WARNING: {self.run_dir} already present')
        # exit()

        data_dir = f'{self.run_dir}/data'
        if not os.path.isdir(data_dir):
            os.system(f'mkdir {data_dir}')

        with open(f'{self.run_dir}/data/run.toml', 'w') as setup_fh:
            _ = toml.dump(self.setup_dic, setup_fh)
        setup_fh.close()

        mol_dic = dict([(e, self.setup_dic['molecules'][e]) for e in self.setup_dic['molecules'] if 'mol' in e])
        for mol_id in mol_dic:
            molecule = self.setup_dic['molecules'][mol_id]

            if molecule == mol_id + '.pdb':
                raise SetupError("+ ERROR: Name mol_X.pdb not supported, please rename your molecule.")
            # Ensembles will be treated later
            if not os.path.isfile(molecule):
                raise SetupError(f'+ ERROR: {molecule} not found!')

            os.system(f'cp {molecule} {self.run_dir}/data/{mol_id}_1.pdb')

            self.setup_dic['molecules'][mol_id] = f'{self.run_dir}/data/{mol_id}_1.pdb'

        # ENHANCEMENT: Check if SegIDs defined in tbl files are on the structures/input
        try:
            ambig_fname = self.setup_dic['restraints']['ambig']
            os.system(f'cp {ambig_fname} {self.run_dir}/data/ambig.tbl')
        except KeyError:
            pass

        try:
            unambig_fname = self.setup_dic['restraints']['unambig']
            os.system(f'cp {unambig_fname} {self.run_dir}/data/unambig.tbl')
        except KeyError:
            pass

        try:
            hbond_fname = self.setup_dic['restraints']['hbond']
            os.system(f'cp {hbond_fname} {self.run_dir}/data/hbond.tbl')
        except KeyError:
            pass
        try:
            dihe_fname = self.setup_dic['restraints']['dihedrals']
            os.system(f'cp {dihe_fname} {self.run_dir}/data/dihe.tbl')
        except KeyError:
            pass

        # Q: What should this function return?
        return True

    def configure_recipes(self):
        """ Prepare recipes with parameters from setup dictionary """

        stage_dic = self.setup_dic['stage']
        for stage in stage_dic:
            recipe_id = stage_dic[stage]['recipe']

            if recipe_id == 'default':
                recipe_id = self.default_recipes[stage]
            else:
                # TODO: check if recipe is valid
                pass

            recipe_name = recipe_id.split('.')[0]
            recipe_params_file = f'{self.protocol_path}/{recipe_name}/{recipe_name}.json'

            if os.path.isfile(recipe_params_file):
                with open(recipe_params_file, 'r') as f:
                    default_parameter_dic = json.load(f)
                f.close()
            else:
                raise SetupError(f'+ ERROR: Default parameters not found for {recipe_id}')

            custom_parameter_dic = dict([(a, stage_dic[stage][a]) for a in stage_dic[stage] if a != 'recipe'])

            # TODO: Save the custom json file alongside the template
            #  find a way to combine default_parameter_dic with custom_parameter_dic and save it for reproductibility
            # parameter_dic = self.merge_parameters(default_parameter_dic, custom_parameter_dic)

            # Generate the recipe and place it in the appropriate place
            r = RecipeComposer(recipe_id, default_parameter_dic, custom_parameter_dic)
            complete_recipe = r.compose()

            if not os.path.isdir(f'{self.run_dir}'):
                os.mkdir(f'{self.run_dir}')

            if not os.path.isdir(f'{self.run_dir}/{stage}'):
                os.mkdir(f'{self.run_dir}/{stage}')
                if not os.path.isdir(f'{self.run_dir}/{stage}/template'):
                    os.mkdir(f'{self.run_dir}/{stage}/template')

            with open(f'{self.run_dir}/{stage}/template/{recipe_id}', 'w') as fh:
                fh.write(complete_recipe)
            fh.close()

        # Q: What should this function return?
        return True

class SetupError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None
        raise SystemExit('{0} '.format(self.message))
