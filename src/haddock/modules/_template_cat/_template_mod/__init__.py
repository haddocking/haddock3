"""
Template module
===============

Use this template to implement your own HADDOCK3 module.

In this docstring, please write the documentation explaining what the module
does technically, scientifically, and its purpose.

Anything you want to write to describe the module, also do it here.
You should use restructureText syntax:

    https://docutils.sourceforge.io/docs/user/rst/quickstart.html
"""
# Import here what you need
from pathlib import Path
from haddock.core.typing import FilePath, Any
from haddock.core.defaults import MODULE_DEFAULT_YAML

# In case you need to import a Python library that is a run-time dependency,
# you should import it inside the `_run` method to avoid import errors for those
# users not using the new module, thus not having its dependencies installed.

# If your module does not use CNS, import the following
from haddock.modules import BaseHaddockModule

# If your module uses CNS import
from haddock.modules.base_cns_module import BaseCNSModule


# this is mandatory, don't erase nor edit these lines
RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


# this is the main class of the module. It should be named exactly as this.
class HaddockModule(BaseHaddockModule):
    # inherit from BaseCNSModule in case the module uses CNS
    """
    Here you can write any extra documentation you wish about the class.

    Use numpydoc syntax:

    https://numpydoc.readthedocs.io/en/latest/format.html
    """

    # this is mandatory
    name = RECIPE_PATH.name

    # this __init__ structure is mandatory, but you can extend it to the
    # module's needs as long as you keep the main structure. Surely,
    # *ignore and **everything can be edited to fit your needs.
    def __init__(
            self,
            order: int,
            path: Path,
            *ignore: Any,
            init_params: FilePath = DEFAULT_CONFIG,
            **everything: Any,
            ) -> None:

        # if your module uses CNS you might need to define where the main CNS
        # script is localted. See examples in `topoaa`, `emref`.
        #  else leave it out.
        # cns_script = Path(RECIPE_PATH, "cns", "main.cns")

        # use one of the following:
        # if the module does not use CNS:
        super().__init__(order, path, init_params)

        # if the module uses CNS:
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if the module is ready to use."""
        # here, you should write any code needed to confirm that all the
        # dependencies required by your module are installed.
        # this class method will be executed when HADDOCK3 starts.

        # if you module does not import any run-time dependency, just leave
        # this method blank
        return

    # here is where the module magic will happen
    def _run(self) -> None:
        # Import here any Python run-time dependencies that your module needs.

        # you can refer to other modules as examples to see how they perform.
        # Likely, the new module will need to load the output from the previous
        # module. For example:
        models_to_use = self.previous_io.retrieve_models()

        # once the module loads, you can access its parameters in
        # `self.params`. See more in the BaseHaddockModule class.
        #  `self.params` is a dictionary that contains both the general
        #  parameters such as the number of processors
        #  `self.params['ncores']` parameters as well as those
        #  defined in your `defaults.yml`

        # Add here all the magic that your module does. You can split this part
        # into many functions, classes, or add any extra files or subfolders to
        # the module's folder as you may need. You can even import other modules
        # if needed.

        # If you module creates models, the line below is almost mandatory.
        # see other modules, such as topoaa or mdref, emref to see how they
        # handle `output_models`.
        # the PDB references in `list_of_created_models` below must be instances
        # of the `libs.libontology.PDBFile` class.

        output_model_list = your_function(models_to_use, self.params)
        # output_model_list = [(pdb, psf, score), ...]

        list_of_created_models = []
        for element in output_model_list:
            pdb, psf, score = element
            # IMPORTANT: pass `file_name=Path.name`

            # alternatively to the strategy below, you can use the
            # `libs.libcns.prepare_expected_pdb` function if it better fits
            # your case. See `emref` and `mdref` for examples.
            pdb_object = PDBFile(
                Path(pdb).name,
                topology=TopologyFile(Path(psf).name, path="."),
                path=".")

            pdb_object.score = score
            list_of_created_models.append(pdb_object)

        # final section
        self.output_models = list_of_created_models
        self.export_io_models()
        # in case your module considers possible tolerance for generated models,
        # you can use:
        # self.export_io_models(faulty_tolerance=self.params["tolerance"])


# Finally, the haddock module's class inherit from BaseHaddockModule. It is
# important that you go through the BaseHaddockModule interface to understand
# it's functionalities. If you find any spots of missing documentation, let us
# know. Thanks! @ The Haddock Team
