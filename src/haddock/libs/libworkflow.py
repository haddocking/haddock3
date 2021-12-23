"""HADDOCK3 workflow logic."""
import importlib
import shutil
import sys
from pathlib import Path

from haddock import log
from haddock.core.exceptions import HaddockError, StepError
from haddock.gear.config_reader import get_module_name
from haddock.libs.libhpc import (
    HPCScheduler_CONCAT_DEFAULT,
    HPCWorker_QUEUE_DEFAULT,
    HPCWorker_QUEUE_LIMIT_DEFAULT,
    )
from haddock.libs.libutil import zero_fill
from haddock.modules import modules_category


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

    def __init__(
            self,
            content,
            ncores=None,
            #run_dir=None,
            cns_exec=None,
            config_path=None,
            mode='local',
            queue=HPCWorker_QUEUE_DEFAULT,
            concat=HPCScheduler_CONCAT_DEFAULT,
            queue_limit=HPCWorker_QUEUE_LIMIT_DEFAULT,
            relative_envvars=True,
            #self_contained=False,
            **others):
        # Create the list of steps contained in this workflow
        self.steps = []
        for num_stage, (stage_name, params) in enumerate(content.items()):
            log.info(f"Reading instructions of [{stage_name}] step")

            # uses gobal ncores parameter unless module-specific value
            # hasn't been used
            params.setdefault('ncores', ncores)
            params.setdefault('cns_exec', cns_exec)
            params.setdefault('config_path', config_path)
            params.setdefault('mode', mode)
            params.setdefault('queue', queue)
            params.setdefault('concat', concat)
            params.setdefault('queue_limit', queue_limit)
            params.setdefault('relative_envvars', relative_envvars)
            #params.setdefault('self_contained', self_contained)

            try:
                _ = Step(
                    get_module_name(stage_name),
                    order=num_stage,
                    #run_dir=run_dir or Path.cwd(),
                    **params,
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
            #run_dir=None,
            **config_params,
            ):
        self.config = config_params
        self.module_name = module_name
        self.order = order

        self.working_path = \
            Path(zero_fill(self.order, digits=2) + "_" + self.module_name)

    def execute(self):
        """Execute simulation step."""
        #if self.working_path.exists():
        #    log.warning(f"Found previous run ({self.working_path}), removed")
        #    shutil.rmtree(self.working_path)
        # this should run in the CWD
        # and the folder should have been removed
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
        try:
            module.run(**self.config)
        except KeyboardInterrupt:
            log.info("You have halted subprocess execution by hitting Ctrl+c")
            log.info("Exiting...")
            sys.exit()
