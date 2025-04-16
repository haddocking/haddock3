"""HADDOCK3 workflow logic."""
import importlib
import sys
from pathlib import Path
from time import time

from haddock import log
from haddock.clis.cli_analyse import main as cli_analyse
from haddock.clis.cli_traceback import main as cli_traceback
from haddock.core.exceptions import HaddockError, HaddockTermination, StepError
from haddock.core.typing import Any, ModuleParams, Optional
from haddock.gear.clean_steps import clean_output
from haddock.gear.config import get_module_name
from haddock.gear.zerofill import zero_fill
from haddock.libs.libtimer import convert_seconds_to_min_sec, log_time
from haddock.libs.libutil import recursive_dict_update
from haddock.modules import (
    modules_category,
    non_mandatory_general_parameters_defaults,
)


class WorkflowManager:
    """Read and execute workflows."""

    def __init__(
        self,
        workflow_params: ModuleParams,
        start: Optional[int] = 0,
        **other_params: Any,
    ) -> None:
        self.start = 0 if start is None else start
        self.recipe = Workflow(workflow_params, start=0, **other_params)
        # terminate is used to synchronize the `clean` option with the
        # `exit` module. If the `exit` module is removed in the future,
        # you can also remove and clean the `terminate` part here.
        self._terminated = None

    def run(self) -> None:
        """High level workflow composer."""
        for i, step in enumerate(self.recipe.steps[self.start :], start=self.start):
            try:
                step.execute()
            except HaddockTermination:
                self._terminated = i  # type: ignore
                break

    def clean(self, terminated: Optional[int] = None) -> None:
        """
        Clean steps.

        Parameters
        ----------
        terminated : int, None
            At which index of the workflow to stop the cleaning. If ``None``,
            uses the internal class configuration.
        """
        terminated = self._terminated if terminated is None else terminated
        for step in self.recipe.steps[:terminated]:
            step.clean()

    def postprocess(self, self_contained: bool = False) -> None:
        """Postprocess the workflow."""
        # is the workflow going to be cleaned?
        is_cleaned = self.recipe.steps[0].config['clean']
        # Is the workflow supposed to run offline
        offline = self.recipe.steps[0].config['offline']
        # running mode
        mode = self.recipe.steps[0].config['mode']
        # ncores
        ncores = self.recipe.steps[0].config['ncores']

        capri_steps: list[int] = []
        for step in self.recipe.steps:
            if step.module_name == "caprieval":
                capri_steps.append(step.order)  # type: ignore
        # call cli_analyse (no need for capri_dicts, it's all precalculated)
        cli_analyse(
            "./",
            capri_steps,
            top_clusters=10,
            format=None,
            scale=None,
            inter=False,
            is_cleaned=is_cleaned,
            offline=offline,
            mode=mode,
            ncores=ncores,
            self_contained=self_contained,
            )
        # call cli_traceback. If it fails, it's not a big deal
        try:
            cli_traceback("./", offline=offline)
        except Exception as e:
            log.warning(f"Error running traceback: {e}")


class Workflow:
    """Represent a set of stages to be executed by HADDOCK."""

    def __init__(
        self,
        modules_parameters: ModuleParams,
        start: Optional[int] = 0,
        **other_params: Any,
    ) -> None:
        if start is None:
            start = 0

        # filter out those parameters not belonging to the modules
        general_modules = {
            k: v
            for k, v in other_params.items()
            if k in non_mandatory_general_parameters_defaults
        }

        # Create the list of steps contained in this workflow
        self.steps: list[Step] = []
        _items = enumerate(modules_parameters.items(), start=start)
        for num_stage, (stage_name, params) in _items:
            stage_name = get_module_name(stage_name)
            log.info(f"Reading instructions step {num_stage}_{stage_name}")

            # updates the module's specific parameter with global parameters
            # that are applicable to the modules. But keep priority to the local
            # level
            params_up = recursive_dict_update(general_modules, params)

            try:
                _ = Step(
                    stage_name,
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
        module_name: str,
        order: Optional[int] = None,
        **config_params: Any,
    ) -> None:
        self.config = config_params
        self.module_name = module_name
        self.order = order
        self.working_path = Path(zero_fill.fill(self.module_name, self.order))  # type: ignore
        self.module = None

    def execute(self) -> None:
        """Execute simulation step."""
        self.working_path.resolve().mkdir(parents=False, exist_ok=False)

        # Import the module given by the mode or default
        module_name = ".".join(
            ["haddock", "modules", modules_category[self.module_name], self.module_name]
        )
        module_lib = importlib.import_module(module_name)
        self.module = module_lib.HaddockModule(order=self.order, path=self.working_path)

        # Run module
        start = time()
        try:
            self.module.update_params(**self.config)  # type: ignore
            self.module.save_config(Path(self.working_path, "params.cfg"))  # type: ignore
            self.module.run()  # type: ignore
        except KeyboardInterrupt:
            log.info("You have halted subprocess execution by hitting Ctrl+c")
            log.info("Exiting...")
            sys.exit(1)

        end = time()
        elapsed = convert_seconds_to_min_sec(end - start)
        self.module.log(f"took {elapsed}")  # type: ignore

    def clean(self) -> None:
        """Clean step output."""
        if self.module is None and self.config["clean"]:
            with log_time("cleaning output files took"):
                clean_output(self.working_path, self.config["ncores"])

        elif self.module is not None and self.module.params["clean"]:
            self.module.clean_output()
