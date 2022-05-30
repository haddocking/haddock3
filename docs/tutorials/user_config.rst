Workflow configuration files
============================

Users can configure HADDOCK3 workflows through human-readable text files. The
configuration file is divided into two sections. The first contains the general
**mandatory** and **optional** parameters. The second section contains the
configuration parameters for each **step** in the workflow. The workflow steps
are defined by a title header with the module's name to execute followed by its
parameters.

.. code:: toml

    # first section: general parameters =======
    # mandatory parameters
    run_dir = "run1-test"

    molecules =  [
        "data/mol1.pdb",
        "data/mol2.pdb"
        ]

    # optional parameters
    ncores = 40

    # second section: step parameters =======
    # each workflow step executes a HADDOCK3 module
    # the step's header is the name of the module
    # the step's paramereters are defined below its header
    [topoaa]
    autohis = false
    [topoaa.mol1]
    nhisd = 0
    nhise = 1
    hise_1 = 75
    [topoaa.mol2]
    nhisd = 1
    hisd_1 = 76
    nhise = 1
    hise_1 = 15

    [rigidbody]
    tolerance = 20
    ambig_fname = "data/mol1.tbl"
    sampling = 20

    [caprieval]
    reference_fname = "data/mol1-mol2.pdb"

.. note::

    If you are a developer and wish to indagate further, the HADDOCK3
    configuration files follow the TOML syntax. However they are not TOML files.
    HADDOCK3 implements its own parser with additional features not covered by
    TOML but that were needed for this project. See :ref:`Workflow Configuration
    file reader`.


HADDOCK3 will use the default values for those parameters not specified by the
user.

To obtain the list of all parameters for each module, you can use the
``haddock-cfg`` command-line. For example, to list all parameters from the
*All-atom topology* module (``topoaa``):

.. code:: bash

    haddock-cfg -m topoaa

    # or to write the output to a file
    haddock-cfg -m topoaa > topoaa_defaults.cfg

    # obtain help on `haddoc-cfg` CLI
    haddock-cfg -h

The HADDOCK3 workflow will follow the order of the steps defined in the
configuration file, reading from top to bottom. For the example above, that
would be: ``topoaa``, ``rigidbody``, ``caprieval``. Therefore, to edit/add/remove
steps in a HADDOCK3 workflow, you simply need to add/remove/edit the text blocks
referring to each module.

Moreover, you can define the general **optional** parameters in the global scope of the
workflow or specifically for each module. Optional parameters defined in the
general scope of the configuration files affect all modules. While those
described inside a step will affect only that step. For example:

.. code:: toml

    run_dir = "run1-test"

    molecules =  [
        "data/mol1.pdb",
        "data/mol2.pdb"
        ]

    # all steps will use at most 40 CPU threads
    ncores = 40

    [topoaa]
    autohis = false

    [rigidbody]
    # while rigidbody is blocked to use only 10 CPU threads
    ncores = 10
    ambig_fname = "data/mol1.tbl"
    sampling = 20

    [caprieval]
    reference_fname = "data/mol1-mol2.pdb"

Inside the `examples
<https://github.com/haddocking/haddock3/tree/main/examples>`_ subfolders you will
a panoply of examples of workflow configuration files (``.cfg``).

Here is a list of all available :ref:`Modules`.

Finally, if you are a developer and wish to use HADDOCK3 as a library to read
and write configuration files please see the related Python modules:

* :ref:`Workflow configuration file reader`
* :ref:`Workflow configuration file writer`

For example, to read a workflow configuration file:

.. code:: python

    from pathlib import Path

    from haddock.gear.config_reader import read_config

    config_path = Path("path", "to", "config.cfg")
    config_dict = read_config(config_path)
