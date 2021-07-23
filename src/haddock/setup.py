import toml
import shutil
import copy
from pathlib import Path
from haddock.error import SetupError
import logging

logger = logging.getLogger(__name__)


class Setup:
    """Sets up the folder structure for the simulation."""

    def __init__(self, workflow_f):
        self.params = toml.load(workflow_f)

    def validate(self):
        # ====
        # Major validation function
        # ===

        # check if steps are valid
        module_directory = Path(__file__).parent / 'modules'
        if 'order' not in self.params['input']:
            raise SetupError("Workflow does not specify the"
                             " order of execution")

        steps = copy.deepcopy(self.params['input']['order'])
        for step in steps:
            module_loc = module_directory / step / 'default.py'
            if not module_loc.exists():
                logger.warning(f'Module {step} is not valid, it will'
                               ' be EXCLUDED from execution')
                self.params['input']['order'].remove(step)

        try:
            params = self._create_folder_structure()
        except SetupError as se:
            logger.error(se)

        return params

    def _create_folder_structure(self):

        params = copy.deepcopy(self.params)
        input_params = params['input']
        project_dir = Path(input_params['project_dir'])
        input_params['project_dir'] = project_dir.absolute()

        data_dir = project_dir / 'data'
        begin_dir = project_dir / 'begin'

        if project_dir.exists():
            logger.warning('This code is not production ready')
            logger.warning(f'{project_dir} exists, it will be REMOVED')
            shutil.rmtree(project_dir)

        project_dir.mkdir()
        begin_dir.mkdir()
        data_dir.mkdir()

        for mol_identifier in input_params['molecules']:
            input_mol = Path(input_params['molecules'][mol_identifier])
            shutil.copy(input_mol, data_dir)

            begin_mol = (begin_dir / f'{mol_identifier}.pdb').absolute()
            shutil.copy(input_mol, begin_mol)

            input_params['molecules'][mol_identifier] = begin_mol

        return params
