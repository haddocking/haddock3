"""Define common test variables."""
from pathlib import Path

from haddock.modules import modules_category


tests_path = Path(__file__).resolve().parents[0]
data_folder = Path(tests_path, 'data')
golden_data = Path(tests_path, 'golden_data')

# preprocessing files
broken_pdb = Path(data_folder, 'broken.pdb')
corrected_pdb = Path(data_folder, 'corrected.pdb')

residues_top = Path(data_folder, 'residues.top')

configs_data = Path(tests_path, 'configs')
emptycfg = Path(configs_data, 'empty.cfg')
haddock3_yaml_cfg_examples = Path(configs_data, 'yml_example.yml')
haddock3_yaml_converted = Path(configs_data, 'yaml2cfg_converted.cfg')
haddock3_yaml_converted_no_header = \
    Path(configs_data, 'yaml2cfg_converted_no_header.cfg')
clean_steps_folder = Path(tests_path, 'clean_output_data')

steptmp = Path(data_folder, "0_dummystep")

# defines which modules are already working
working_modules = [t for t in modules_category.items() if t[0] != 'topocg']
