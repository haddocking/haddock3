Modules
=======

HADDOCK3 allows users to compose modular simulation workflows. Workflows are
composed in *steps*, and each step is a HADDOCK3 module. There are modules for
**sampling**, **refinement**, **analysis**, etc.

.. toctree::
   :maxdepth: 2

   topology/index
   sampling/index
   refinement/index
   scoring/index
   analysis/index
   extras/index

Parent code for modules
-----------------------

.. automodule:: haddock.modules
   :members:
   :show-inheritance:
   :inherited-members:

Parent code for CNS modules
---------------------------

.. automodule:: haddock.modules.base_cns_module
   :members:
   :show-inheritance:
   :inherited-members:

General Default parameters
--------------------------

General default parameters can be defined in the main section of the
configuration file, but can also be defined for each individual module (step).
In the later case, overriding the general definition.

.. include:: general_module_params.rst
