# CNS Installation

## Download CNS

Downloading [CNS (Crystallography and NMR System)](http://cns-online.org/v1.3/) requires a license, which academic users can [request for free](http://cns-online.org/cns_request/). You will be sent a password for the download via e-mail.
Using the link provided in the same e-mail, download CNS ("Source installer (compilation required on your system)") to the directory in which you want to install it (for example `software`). Navigate to this directory, then uncompress and extract the archive:

```bash
cd ~/software/
gunzip cns_solve_1.3_all.tar.gz
tar -xvf cns_solve_1.3_all.tar
```

Go to the CNS directory:

```bash
cd cns_solve_1.3
```

If you are using csh/tcsh, edit `cns_solve_env` to point to the correct location (the current directory):

```
# CHANGE THE NEXT LINE TO POINT TO THE LOCATION OF THE CNSsolve DIRECTORY

            setenv CNS_SOLVE PATH/TO/cns_solve_1.3
```

If you are using sh/bash, edit `.cns_solve_env_sh`:

```
# CHANGE THE NEXT LINE TO POINT TO THE LOCATION OF THE CNSsolve DIRECTORY

        CNS_SOLVE=PATH/TO/cns_solve_1.3
```

## Add HADDOCK routines

HADDOCK requires some adjustments in CNS to function properly. Download the [set of routines to add](https://www.dropbox.com/s/wliubqovuusqdvr/cns.tgz?dl=0). Uncompress the archive and move all of its files into the CNS source directory (`PATH/TO/cns_solve_1.3/source/`).

## Recompile CNS

Make sure you are inside the main `cns_solve_1.3` directory.

If a suitable compiler is installed on your system, the following command should start the compilation of CNS:

```bash
make install
```

```bash
# For csh/tcsh:
source cns_solve_env
# For sh/bash:
source .cns_solve_env_sh
```
