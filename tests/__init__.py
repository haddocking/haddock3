"""Test module."""
from pathlib import Path

from haddock.modules import modules_category


tests_path = Path(__file__).resolve().parents[0]
golden_data = Path(tests_path, 'golden_data')
data_path = Path(tests_path, 'data')

configs_data = Path(tests_path, 'configs')
emptycfg = Path(configs_data, 'empty.cfg')
haddock3_yaml_cfg_examples = Path(configs_data, 'yml_example.yml')
haddock3_yaml_converted = Path(configs_data, 'yaml2cfg_converted.cfg')

# defines which modules are already working
working_modules = [t for t in modules_category.items() if t[0] != 'topocg']

adjusted_col_width_table = Path(data_path, 'adjusted_col_width_table.txt')
