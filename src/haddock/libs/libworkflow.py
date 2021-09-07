"""HADDOCK3 workflow logic"""
import logging
import importlib
import shutil
from pathlib import Path
from haddock.core.exceptions import HaddockError, StepError
# unused
# from haddock.core.defaults import MODULE_PATH_NAME, TOPOLOGY_PATH
from haddock.modules import modules_category
from haddock.libs.libutil import zero_fill

logger = logging.getLogger(__name__)


class WorkflowManager:
    """The WorkflowManager reads the workflow and executes them."""
    def __init__(self, workflow_params, start=0, **other_params):
        self.start = start
        # Create a workflow from a TOML file
        self.recipe = Workflow(workflow_params, **other_params)

    def run(self):
        """High level workflow composer"""
        for step in self.recipe.steps[self.start:]:
            step.execute()


class Workflow:
    """Represents a set of stages to be executed by HADDOCK"""
    def __init__(self, content, run_dir=None):
        # Create the list of steps contained in this workflow
        self.steps = []
        for num_stage, (stage_name, params) in enumerate(content.items()):
            logger.info(f"Reading instructions of [{stage_name}] step")

            try:
                _ = Step(
                    stage_name,
                    order=num_stage,
                    run_dir=run_dir,
                    **params,
                    )
                self.steps.append(_)

            except StepError as re:
                logger.error(f"Error found while parsing course {stage_name}")
                raise HaddockError from re


class Step:
    """Represents a Step of the Workflow."""

    def __init__(self, module_name, order=None, run_dir=None, **config_params):
        self.config = config_params
        self.module_name = module_name
        self.order = order

        self.working_path = Path(
            run_dir,
            zero_fill(self.order, digits=2) + "_" + self.module_name,
            )

    def execute(self):
        if self.working_path.exists():
            logger.warning(f"Found previous run ({self.working_path}),"
                           " removed")
            shutil.rmtree(self.working_path)
        self.working_path.resolve().mkdir(parents=True, exist_ok=False)

        # Import the module given by the mode or default
        module_name = ".".join([
            'haddock',
            'modules',
            modules_category[self.module_name],
            self.module_name
            ])
        module_lib = importlib.import_module(module_name)
        module = module_lib.HaddockModule(
            order=self.order,
            path=self.working_path)

        # Run module
        module.run(**self.config)
