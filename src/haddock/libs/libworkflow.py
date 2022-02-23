"""HADDOCK3 workflow logic."""
import importlib
import sys
from pathlib import Path
from time import time

from haddock import log
from haddock.core.exceptions import HaddockError, StepError
from haddock.gear.config_reader import get_module_name
from haddock.libs.libutil import (
    convert_seconds_to_min_sec,
    recursive_dict_update,
    zero_fill,
    )
from haddock.modules import (
    modules_category,
    non_mandatory_general_parameters_defaults,
    )


class WorkflowManager:
    """Read and execute workflows."""

    def __init__(self, workflow_params, start=0, **other_params):
        self.start = start
        # Create a workflow from a TOML file
        self.recipe = Workflow(workflow_params, **other_params)

    def run(self):
        """High level workflow composer."""
        for step in self.recipe.steps[self.start:]:
            step.execute()


class Workflow:
    """Represent a set of stages to be executed by HADDOCK."""

    def __init__(self, modules_parameters, **other_params):

        # filter out those parameters not belonging to the modules
        general_modules = {
            k: v
            for k, v in other_params.items()
            if k in non_mandatory_general_parameters_defaults
            }

        # Create the list of steps contained in this workflow
        self.steps = []
        _items = enumerate(modules_parameters.items())
        for num_stage, (stage_name, params) in _items:
            log.info(f"Reading instructions of [{stage_name}] step")

            # updates the module's specific parameter with global parameters
            # that are applicable to the modules. But keep priority to the local
            # level
            params_up = recursive_dict_update(general_modules, params)

            try:
                _ = Step(
                    get_module_name(stage_name),
                    order=num_stage,
                    **params_up,
                    )
                self.steps.append(_)

            except StepError as re:
                log.error(f"Error found while parsing course {stage_name}")
                raise HaddockError from re


class Step:
    """Represents a Step of the Workflow."""

    def __init__(
            self,
            module_name,
            order=None,
            **config_params,
            ):
        self.config = config_params
        self.module_name = module_name
        self.order = order

        self.working_path = \
            Path(zero_fill(self.order, digits=2) + "_" + self.module_name)

    def execute(self):
        """Execute simulation step."""
        self.working_path.resolve().mkdir(parents=False, exist_ok=False)

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
        start = time()
        try:
            module.run(**self.config)
        except KeyboardInterrupt:
            log.info("You have halted subprocess execution by hitting Ctrl+c")
            log.info("Exiting...")
            sys.exit(1)

        end = time()
        elapsed = convert_seconds_to_min_sec(end - start)
        module.log(f"took {elapsed}")
