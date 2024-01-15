Workflow configuration files
============================

Users can configure HADDOCK3 workflows through human-readable text files. The
configuration file is divided into two sections. The first contains the general
**mandatory** and **optional** parameters. The second section contains the
configuration parameters for each **step** in the workflow. The workflow steps
are defined by a title header with the module's name to execute followed by its
parameters. The ``parameter = value`` pair definitions follow strictly the TOML
syntax. However, contrarily to TOML, in HADDOCK3 configuration files can have
headers with repeated names.

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

    If you are a developer and wish to know more details, see the
    :py:mod:`haddock.gear.config` module.


Each module does have default parameters defined and HADDOCK3 will use those
default values unless those are specified and changed by the user in the
configuration file.

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

    # each .job will produce 5 (or less) models
    mode = "batch"
    concat = 5

    [topoaa]
    autohis = false

    [rigidbody]
    ambig_fname = "data/mol1.tbl"
    sampling = 20

    [caprieval]
    reference_fname = "data/mol1-mol2.pdb"

    [flexref]
    # flexref jobs are bound to generate one model per job
    # therefore, we can specify the 'concat' parameter specifically for flexref
    concat = 1
    ambig_fname = "data/mol1.tbl"

    [caprieval]
    reference_fname = "data/mol1-mol2.pdb"

Inside the `examples
<https://github.com/haddocking/haddock3/tree/main/examples>`_ subfolders you will
a panoply of examples of workflow configuration files (``.cfg``).

Here is a list of all available :ref:`Modules`.

Finally, if you are a developer and wish to use HADDOCK3 as a library to read
and write configuration files please see the related Python modules:

* :ref:`Configuration file I/O`

For example, to read a workflow configuration file:

.. code:: python

    from pathlib import Path

    from haddock.gear.config import load

    config_path = Path("path", "to", "config.cfg")
    config_dict = load(config_path)

Compatibility with TOML
-----------------------

HADDOCK3 actually uses `TOML <https://pypi.org/project/toml/>`_ to read the
configuration files. However, some additional features are introduced to enhance
user experience. For example, the capability to repeat header names. Hence,
HADDOCK3 can also read pure TOML files as workflow configuration files. For
those cases, repeated modules should have a trailing integer in their
definition:

.. code:: toml

    [topoaa]
    # parameters here...

    [rigidbody]
    # parameters here...

    [caprieval]
    # parameters here...

    [flexref]
    # parameters here...

    ['caprieval.2']
    # parameters here...

The exact trailing integer is irrelevant as long as headers are not repeated.
HADDOCK3 will normalize the integers when reading the config file. For example,
the example below is the same as the above example. What matters is the
order in which the steps are presented in the configuration file.

.. code:: toml

    [topoaa]
    # parameters here...

    [rigidbody]
    # parameters here...

    ['caprieval.10']  # <- mind the 10 here
    # parameters here...

    [flexref]
    # parameters here...

    ['caprieval.2']
    # parameters here...
