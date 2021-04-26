import logging
import toml
import importlib
import shutil
from pathlib import Path
from haddock.error import HaddockError, RecipeError
from haddock.defaults import MODULE_PATH_NAME, TOPOLOGY_PATH


logger = logging.getLogger(__name__)


class Chef:
    """The Chef class is in charge of reading recipes and cooking"""
    def __init__(self, recipe_path, start=0):
        path = Path(recipe_path)
        if not path.exists():
            raise HaddockError(f"Recipe {recipe_path} does not exist or cannot be opened")
        self.start = start
        # Save current working path
        self.workspace = path.parent.absolute()
        # Create a recipe from a TOML file
        self.recipe = Recipe(path)

    def cook(self):
        """High level workflow composer"""
        for course in self.recipe.courses[self.start:]:
            course.cook()


class Recipe:
    """Represents a set of stages to be executed by HADDOCK"""
    def __init__(self, recipe_path):
        content = toml.load(recipe_path)
        # Check if the order has been set
        if 'order' not in content:
            raise RecipeError("This recipe does not specify the order of execution")

        # Basic check for defined stages by order
        for stage in content["order"]:
            if stage not in content["stage"]:
                raise RecipeError(f"{stage} is defined by the order variable but not defined in this recipe")

        logger.info(f"Recipe [{recipe_path}] contains {len(content['order'])} courses")

        # Create the list of courses contained in this recipe
        self.courses = []
        for num_stage, stage in enumerate(content["order"]):
            try:
                logger.info(f"Reading instructions of [{stage}] course")
                self.courses.append(Course(stage, num_stage, content["stage"][stage], recipe_path.parent))
            except RecipeError as re:
                logger.error(f"Error found while parsing course {stage}")
                raise HaddockError(str(re))


class Course:
    """A course represents a HADDOCK module"""
    def __init__(self, module_name, order, course_information, working_path):
        self.module = module_name
        self.order = order
        self.raw_information = course_information
        self.working_path = working_path

    def cook(self):
        try:
            self.flavour = self.raw_information["flavour"]
            if not self.flavour:
                self.flavour = "default"
        except KeyError:
            self.flavour = "default"

        # Create course path structure
        if self.module == "topology":
            p = self.working_path / Path(TOPOLOGY_PATH)
        else:
            p = self.working_path / Path(f"{MODULE_PATH_NAME}{self.order}")
        if p.exists():
            logger.warning(f"Found previous run ({p}), removed")
            shutil.rmtree(p)
        p.absolute().mkdir(parents=True, exist_ok=False)

        # Import the module given by the flavour or default
        module_name = f"haddock.modules.{self.module}.{self.flavour}"
        module_lib = importlib.import_module(module_name)
        module = module_lib.HaddockModule(order=self.order, path=p.absolute())
        # Remove flavour information as it is already used and won't be mapped
        self.raw_information.pop("flavour", None)

        # Run module
        module.run(self.raw_information)
