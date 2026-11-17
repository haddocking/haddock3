Modules
=======

HADDOCK3 allows users to compose modular simulation workflows. Workflows are
composed in *steps*, and each step is a HADDOCK3 module. There are modules for
**sampling**, **refinement**, **analysis**, etc.

.. toctree::
   :maxdepth: 2

   ../../src/modules/topology/index
   ../../src/modules/sampling/index
   ../../src/modules/refinement/index
   ../../src/modules/scoring/index
   ../../src/modules/analysis/index
   ../../src/modules/extras/index

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

.. include:: ../../src/modules/general_module_params.rst
