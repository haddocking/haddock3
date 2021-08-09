"""HADDOCK3 workflow logic"""
import logging
import importlib
import shutil
from pathlib import Path
from haddock.error import HaddockError, StepError
# unused
# from haddock.defaults import MODULE_PATH_NAME, TOPOLOGY_PATH
from haddock.modules import modules_index
from haddock.libs.libutil import zero_fill

logger = logging.getLogger(__name__)


class WorkflowManager:
    """The WorkflowManager reads the workflow and executes them."""
    def __init__(self, workflow_params, start=0):
        self.start = start
        # Create a workflow from a TOML file
        self.recipe = Workflow(workflow_params)

    def run(self):
        """High level workflow composer"""
        for step in self.recipe.steps[self.start:]:
            step.execute()


class Workflow:
    """Represents a set of stages to be executed by HADDOCK"""
    def __init__(self, content):
        order = content['order']
        # Create the list of steps contained in this workflow
        self.steps = []
        for num_stage, stage in enumerate(order):
            try:
                logger.info(f"Reading instructions of [{stage}] step")
                is_substage = (len(stage.split('.')) == 2)
                if is_substage:
                    sub_stage, sub_id = stage.split('.')
                    self.steps.append(Step(sub_stage, num_stage))
                else:
                    self.steps.append(Step(content, stage))
            except StepError as re:
                logger.error(f"Error found while parsing course {stage}")
                raise HaddockError from re


class Step:
    """Represents a Step of the Workflow."""

    def __init__(self, stream, module_name):
        self.stream = stream
        self.module = module_name
        self.order = stream['order']
        self.raw_information = stream['stage'][module_name]
        self.module_index = self.order.index(self.module)
        self.working_path = Path(
            stream['run_dir'],
            zero_fill(self.module_index, digits=2) + "_" + self.module,
            )

        # ===================================================
        # self.mode = "default"
        # DISCUSSION: mode will change the behaviour of a module,
        #  shall we keep this or enforce that different
        #  behaviours = different modules?
        #
        #  for example:
        #   Module = rigid_body
        #    Mode = with-CM-restraints
        #    Mode = with-random-restraints
        #    Mode = with-ambig-restraints
        #  OR
        #   Module = rigid_body-CM
        #   Module = rigid_body-randair
        #   Module = rigid_body-ambig
        # =====
        # if "module" in self.raw_information and self.raw_information["module"]:
        #     self.mode = self.raw_information["module"]
        # else:
        #     self.mode = "default"
        # ===================================================

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
            modules_index[self.module],
            self.module
            ])
        module_lib = importlib.import_module(module_name)
        step_number = self.order.index(self.module)
        module = module_lib.HaddockModule(stream=self.stream,
                                          order=step_number,
                                          path=self.working_path)
        # Remove mode information as it is already used and won't be mapped
        self.raw_information.pop("mode", None)

        # Run module
        module.run(self.raw_information)
