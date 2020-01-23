import json
import os
import re
import random
from haddock.utils.files import get_full_path


class RecipeComposer:

    def __init__(self, recipe_file, default_params, custom_params):
        self.protocol_path = get_full_path('haddock', 'protocols')
        self.recipe_name = recipe_file.split('.')[0]
        self.recipe = f'{self.protocol_path}/{self.recipe_name}/{recipe_file}'

        # recipe_params = self.protocol_path + '/' + recipe_file.split('.')[0] + '.json'
        # if os.path.isfile(recipe_params):
        # 	with open(recipe_params, 'r') as f:
        # 		self.recipe_params = json.load(f)
        # 	f.close()
        # else:
        # 	print(f'+ ERROR: Default parameters not found for {recipe_file}')
        # 	exit()

        self.default_params = default_params
        self.custom_params = custom_params

    def compose(self):
        """ Load the CNS recipe identify which modules need to be loaded and compose the recipe body """
        module_regex = r'.*(\@|/)(.*\.cns)'

        _ = self.list_dependencies(self.recipe)

        with open(self.recipe) as f:
            target = f.readlines()
        f.close()

        new = []
        check = True
        module_check = [(i, e) for i, e in enumerate(target) if '.cns' in e if '!' not in e]
        while check:

            for idx, m in module_check:

                # module found, get name and load
                try:
                    m = re.findall(module_regex, target[idx])[0][-1]
                except IndexError:
                    print(target[idx])
                    exit()
                with open(f'{self.protocol_path}/{self.recipe_name}/{m}') as f:
                    module = f.readlines()
                f.close()

                new = target[:idx] + module + target[idx + 1:]

                break

            target = new
            module_check = [(i, e) for i, e in enumerate(target) if '.cns' in e if '!' not in e]
            if module_check:
                check = True
            else:
                check = False

        parameter_header = self.load_recipe_params()

        recipe = parameter_header + ''.join(new)

        return recipe

    @staticmethod
    def identify_modules(cns_file):
        """ Find all modules in this CNS file """
        module_list = []

        if not os.path.isfile(cns_file):
            print(f'+ ERROR: Module {cns_file} not found')
            exit()

        with open(cns_file) as f:
            for pos, l in enumerate(f.readlines()):
                if '!' not in l:
                    if '@' in l:
                        if '.cns' in l:
                            # module detected
                            module = l.split('@')[-1].split('/')[-1].split('\n')[0]
                            # module_name = f'{self.protocol_path}/{module}'
                            module_list.append(module)
                        # print(pos, module)
        return module_list

    def list_dependencies(self, target_f):
        """ List the CNS modular dependency of the target file """

        output_list = []

        target_f_name = target_f.split('/')[-1]

        # print(f'+ Creating dependency tree of {target_f_name}')
        # print(f'> {target_f_name}')

        # Lvl 1  Identify modules
        module_list = self.identify_modules(target_f)
        if module_list:
            tbw = ''
            for module in module_list:
                output_list.append(module)
                module_path = f'{self.protocol_path}/{self.recipe_name}/{module}'
                tbw += f'  |_ {module} \n'

                # Lvl 2 Dependencies
                # print(module)
                dependency_list = self.identify_modules(module_path)
                if dependency_list:
                    for dependency in dependency_list:
                        output_list.append(dependency)
                        dependency_path = f'{self.protocol_path}/{self.recipe_name}/{dependency}'
                        tbw += f'     |_ {dependency} \n'

                        # Lvl 3 Co-dependency
                        co_dependency_list = self.identify_modules(dependency_path)
                        if co_dependency_list:
                            for co_dependency in co_dependency_list:
                                output_list.append(co_dependency)
                                co_dependency_path = f'{self.protocol_path}/{self.recipe_name}/{co_dependency}'
                                tbw += f'        |_ {co_dependency} \n'

                                # Lvl 4 Co-co-dependency
                                co_co_dependency_list = self.identify_modules(co_dependency_path)
                                if co_co_dependency_list:
                                    for co_co_dependency in co_co_dependency_list:
                                        output_list.append(co_co_dependency)
                                        co_co_dependency_path = f'{self.protocol_path}/{self.recipe_name}/{co_co_dependency}'
                                        print(f'          |_ {co_co_dependency}')

                                        # Lvl 5 Co-co-co-dependency
                                        # You have gone too far...
                                        co_co_co_dependency_list = self.identify_modules(co_co_dependency_path)
                                        if co_co_co_dependency_list:
                                            for co_co_co_dependency in co_co_co_dependency_list:
                                                print(f'+ ERROR: Too many dependency levels {target_f_name} > '
                                                      f'{dependency} > {co_dependency} > {co_co_dependency} > {co_co_co_dependency}')
        # print(tbw)
        return output_list

    def load_recipe_params(self):

        param_header = '\n! Parameters\n'

        for param in self.default_params['params']:

            try:
                v = self.custom_params['params'][param]
            except KeyError:
                # custom definition not found, use default
                v = self.default_params['params'][param]

            if isinstance(v, bool):
                v = str(v).lower()
                param_header += f'eval (${param}={v})\n'

            elif isinstance(v, str):
                param_header += f'eval (${param}="{v}")\n'

            elif isinstance(v, int):
                param_header += f'eval (${param}={v})\n'

            elif isinstance(v, float):
                param_header += f'eval (${param}={v})\n'

            elif not v:
                # either 0 or empty string
                if isinstance(v, str):
                    v = '\"\"'
                    param_header += f'eval (${param}={v})\n'
                if isinstance(v, int):
                    v = 0.0
                    param_header += f'eval (${param}={v})\n'

        # if isinstance(v, int):
        # 	param_header += f'eval (${param}={v})\n'
        # elif isinstance(v, str):
        # 	param_header += f'eval (${param}="{v}")\n'

        if 'chain' in self.default_params:
            # load molecule specific things
            for mol in self.default_params['chain']:
                for param in self.default_params['chain'][mol]:
                    v = self.default_params['chain'][mol][param]
                    # this are LOGICAL, which means no quotes
                    param_header += f'eval (${param}_{mol}={v})\n'

        # evaluate($toppar.prot_segid_$nmol = $prot_segid_mol$nmol)
        # evaluate($toppar.fix_origin_$nmol =$fix_origin_mol$nmol)
        # evaluate($toppar.dna_$nmol =$dna_mol$nmol)
        # evaluate($toppar.cyclicpept_$nmol = $cyclicpept_mol$nmol)
        # evaluate($toppar.shape_$nmol = $shape_mol$nmol)
        # evaluate($toppar.cg_$nmol = $cg_mol$nmol)

        return param_header
