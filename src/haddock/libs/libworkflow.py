"""HADDOCK3 workflow logic."""
import importlib
import sys
from pathlib import Path
from time import time

from haddock import log
from haddock.core.exceptions import HaddockError, StepError
from haddock.gear.clean_steps import clean_output
from haddock.gear.config_reader import get_module_name
from haddock.gear.zerofill import zero_fill
from haddock.libs.libtimer import convert_seconds_to_min_sec, log_time
from haddock.libs.libutil import recursive_dict_update
from haddock.modules import (
    modules_category,
    non_mandatory_general_parameters_defaults,
    )


class WorkflowManager:
    """Read and execute workflows."""

    def __init__(self, workflow_params, start=0, **other_params):
        self.start = start
        self.recipe = Workflow(workflow_params, start=0, **other_params)

    def run(self):
        """High level workflow composer."""
        for step in self.recipe.steps[self.start:]:
            step.execute()

    def clean(self):
        """Clean steps."""
        for step in self.recipe.steps:
            step.clean()


class Workflow:
    """Represent a set of stages to be executed by HADDOCK."""

    def __init__(self, modules_parameters, start=0, **other_params):

        # filter out those parameters not belonging to the modules
        general_modules = {
            k: v
            for k, v in other_params.items()
            if k in non_mandatory_general_parameters_defaults
            }

        # Create the list of steps contained in this workflow
        self.steps = []
        _items = enumerate(modules_parameters.items(), start=start)
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
        self.working_path = Path(zero_fill.fill(self.module_name, self.order))
        self.module = None

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
        self.module = module_lib.HaddockModule(
            order=self.order,
            path=self.working_path)

        # Run module
        start = time()
        try:
            self.module.update_params(**self.config)
            self.module.save_config(Path(self.working_path, "params.cfg"))
            self.module.run()
        except KeyboardInterrupt:
            log.info("You have halted subprocess execution by hitting Ctrl+c")
            log.info("Exiting...")
            sys.exit(1)

        end = time()
        elapsed = convert_seconds_to_min_sec(end - start)
        self.module.log(f"took {elapsed}")

    def clean(self):
        """Clean step output."""
        if self.module is None and self.config["clean"]:
            with log_time("cleaning output files took"):
                clean_output(self.working_path, self.config["ncores"])

        elif self.module is not None and self.module.params["clean"]:
            self.module.clean_output()
