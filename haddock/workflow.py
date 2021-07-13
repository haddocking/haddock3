"""HADDOCK3 workflow logic"""
import logging
import importlib
import shutil
from pathlib import Path
import toml
from haddock.error import HaddockError, RecipeError
from haddock.defaults import MODULE_PATH_NAME, TOPOLOGY_PATH


logger = logging.getLogger(__name__)


class WorkflowManager:
    """The WorkflowManager reads the recipes and executes them."""
    def __init__(self, recipe_path, start=0):
        path = Path(recipe_path)
        if not path.exists():
            raise HaddockError(f"Recipe {recipe_path} does not exist or cannot"
                               " be opened")
        self.start = start
        # Save current working path
        self.workspace = path.parent.absolute()
        # Create a recipe from a TOML file
        self.recipe = Recipe(path)

    def run(self):
        """High level workflow composer"""
        for step in self.recipe.steps[self.start:]:
            step.execute()


class Recipe:
    """Represents a set of stages to be executed by HADDOCK"""
    def __init__(self, recipe_path):
        content = toml.load(recipe_path)
        # Check if the order has been set
        if 'order' not in content:
            raise RecipeError("This recipe does not specify the order"
                              " of execution")

        logger.info(f"Recipe [{recipe_path}] contains {len(content['order'])}"
                    " steps")

        # Create the list of courses contained in this recipe
        self.steps = []
        for num_stage, stage in enumerate(content["order"]):
            try:
                logger.info(f"Reading instructions of [{stage}] step")
                is_substage = (len(stage.split('.')) == 2)
                if is_substage:
                    sub_stage, sub_id = stage.split('.')
                    self.steps.append(Step(sub_stage,
                                           num_stage,
                                           content["stage"][sub_stage][sub_id],
                                           recipe_path.parent))
                else:
                    self.steps.append(Step(stage,
                                           num_stage,
                                           content["stage"][stage],
                                           recipe_path.parent))
            except RecipeError as re:
                logger.error(f"Error found while parsing course {stage}")
                raise HaddockError from re


class Step:
    """Represents a step to be executed in the workflow."""

    def __init__(self, module_name, order, course_information, working_path):
        self.module = module_name
        self.order = order
        self.raw_information = course_information
        self.working_path = working_path
        if "mode" in self.raw_information and self.raw_information["mode"]:
            self.mode = self.raw_information["mode"]
        else:
            self.mode = "default"

    def execute(self):
        # Create course path structure
        if self.module == "topology":
            p = self.working_path / Path(TOPOLOGY_PATH)
        else:
            p = self.working_path / Path(f"{MODULE_PATH_NAME}{self.order}")
        if p.exists():
            logger.warning(f"Found previous run ({p}), removed")
            shutil.rmtree(p)
        p.absolute().mkdir(parents=True, exist_ok=False)

        # Import the module given by the mode or default
        module_name = f"haddock.modules.{self.module}.{self.mode}"
        module_lib = importlib.import_module(module_name)
        module = module_lib.HaddockModule(order=self.order, path=p.absolute())
        # Remove mode information as it is already used and won't be mapped
        self.raw_information.pop("mode", None)

        # Run module
        module.run(self.raw_information)
