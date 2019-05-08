![HADDOCK3](docs/media/HADDOCK3-logo.png)

The official repo of the new modular BioExcel2 version of HADDOCK.

**ATTENTION: This repository is under heavy development and may change abruptly.**

## Container setup for HADDOCK2.4
### Singularity
*[Install instructions](https://www.sylabs.io/guides/3.2/user-guide/installation.html)*

#### Build image from bootstrap file

Note: Change the username in line 22 of `haddock3/containers/Haddock.def` and `NUMJOB` according to your config.

```
git clone 
cd haddock3/containers
sudo singularity build --sandbox Haddock.img Haddock.def
```

#### Run a shell within a container
```
sudo singularity shell --writable Haddock.img
```

#### Test
```
singularity run Haddock.img
```

```
##############################################################################
#                                                                            #
# Starting HADDOCK2.4                                                        #
#                                                                            #
#         N-components version of HADDOCK (current maximum is 20)            #
#                                                                            #
#   Copyright 2003-2018 Alexandre Bonvin, Utrecht University.                #
#   Originally adapted from Aria 1.2 from Nilges and Linge, EMBL.            #
#   All rights reserved.                                                     #
#   This code is part of the HADDOCK software and governed by its            #
#   license. Please see the LICENSE file that should have been included      #
#   as part of this package.                                                 #
#                                                                            #
##############################################################################

Starting HADDOCK on: 2019-05-08 14:14:47

HADDOCK version:   2.4
Python version: 2.7.15 (default, Aug 22 2018, 13:28:29)
[GCC 6.3.0]
PYTHONPATH system variable contains:
['/software/haddock/Haddock', '/home/vagrant', '/software/haddock', '/usr/lib/python27.zip', '/usr/lib/python2.7', '/usr/lib/python2.7/plat-linux2', '/usr/lib/python2.7/lib-tk', '/usr/lib/python2.7/lib-old', '/usr/lib/python2.7/lib-dynload', '/usr/lib/python2.7/site-packages']
there is no run.cns OR run.param in your current directory.
=> HADDOCK stopped
##############################################################################
Finishing HADDOCK on:  2019-05-08 14:14:47
Au revoir.
Tot ziens.
Bye bye.
```

### Docker
```
git clone 
cd haddock3/containers
curl -O https://surfdrive.surf.nl/files/index.php/s/zi7DSzACww5B3jQ/download
mv download haddock2.4-img.tar.gz
gzip -d haddock2.4-img.tar.gz
 
```

#### Test container

```
cd docker
docker load -i haddock2.4-img.tar
docker run --tmpfs /tmp -it haddock2.4 /bin/bash
docker run -it --tmpfs /tmp haddock2.4 /bin/tcsh
```

Inside the container, type `haddock2.4`:

```
root@ffc004adbae3:/# haddock2.4 

##############################################################################
#                                                                            #
# Starting HADDOCK2.4                                                        #
#                                                                            #
#         N-components version of HADDOCK (current maximum is 20)            #
#                                                                            #
#   Copyright 2003-2018 Alexandre Bonvin, Utrecht University.                #
#   Originally adapted from Aria 1.2 from Nilges and Linge, EMBL.            #
#   All rights reserved.                                                     #
#   This code is part of the HADDOCK software and governed by its            #
#   license. Please see the LICENSE file that should have been included      #
#   as part of this package.                                                 #
#                                                                            #
##############################################################################
    
Starting HADDOCK on: 2019-04-26 10:39:44
 
HADDOCK version:   2.4
Python version: 2.7.15rc1 (default, Nov 12 2018, 14:31:15) 
[GCC 7.3.0]
PYTHONPATH system variable contains:
['/software/haddock/Haddock', '/', '/software/haddock', '/usr/lib/python2.7', '/usr/lib/python2.7/plat-x86_64-linux-gnu', '/usr/lib/python2.7/lib-tk', '/usr/lib/python2.7/lib-old', '/usr/lib/python2.7/lib-dynload', '/usr/local/lib/python2.7/dist-packages', '/usr/lib/python2.7/dist-packages']
there is no run.cns OR run.param in your current directory.
=> HADDOCK stopped
##############################################################################
Finishing HADDOCK on:  2019-04-26 10:39:44
Au revoir.
Tot ziens.
Bye bye.
```