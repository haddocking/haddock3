### Using haddock3 with module

This directory provides an example module for using haddock3 with the `module` command if installed on your system.
`module` is part of the lmod software (refer to documentation below for details). In short, it allows to enable and disable software at wish on a server.

This directory contains an example module file for haddock3. It assumes that haddock3 has been installed using miniconda3. The module will source the conda environment and activate haddock3. To use it in your enviroment, someone with admin rights will have to edit it and place it in the proper location.

To configure the module script for your system, edit the provided `haddock3/3.0.0.lua` file and change the following line to specifiy the installation location of miniconda3:

```
local conda_dir = "/trinity/home/software/miniconda3"

```

Then copy the entire haddock3 directory into the `modulefiles` directory of your lmod installation. E.g on our system this is: `/opt/ohpc/pub/modulefiles/` as lmod was installed as part of the open HPC software stack (see for example https://openhpc.github.io/cloudwg/) 



#### Documentation

Lmod Web Sites

*  Documentation:    http://lmod.readthedocs.org
*  Github:           https://github.com/TACC/Lmod
*  Sourceforge:      https://lmod.sf.net
*  TACC Homepage:    https://www.tacc.utexas.edu/research-development/tacc-projects/lmod
