Welcome!

To develop your own HADDOCK3 module follow the patterns used in other modules.

1. If your module belongs to a new category, create a folder for that category
   under the `modules` folder.
1. Else, create a folder for your module inside its relevant category.
   The name of that folder is the name of the module, i.e., the name used to call
   the module in the haddock3 configuration files and used throughout the code base.
1. Copy the `__init__.py` file here to the new module's folder and edit it accordingly
   to the instructions there.
1. Do the same for the `defaults.yaml` file.
1. You can then add any extra files needed inside your module's folder in order to
   develop your module fully.
1. If your module requires any extra libraries, describe how to install those libraries
   in the `docs/INSTALL.md` file. Unless approved by the Haddock Team, do not add
   those dependencies to the `requirements.*` files.
1. HADDOCK3 has already many features related to IO, subprocess run, sending jobs,
   etc. Please, look around the `libs` folder for pertinent functions, but, above all,
   feel welcomed to reach out to us with any doubts.
1. Please write also tests for your module. Our testing machinery already
   tests for the common patterns, for example, inspecting the `defaults.yaml` file.
   But you should write any additional tests to ensure that your module works properly.
   See other examples in the `tests/` folder.
1. Finally, add an example of how to use your module in the `examples/` folder.
   The example should have a short sampling scheme. Name the config file ending with
   `-test.cfg`.

Thanks, and we are happy to help you! :-)
@ The Haddock Team
