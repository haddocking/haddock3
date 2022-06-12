# CNS Installation

## 1 Download CNS

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

## 2 Add HADDOCK routines

HADDOCK requires some adjustments in CNS to function properly. Download the [set of routines to add](https://www.dropbox.com/s/wliubqovuusqdvr/cns.tgz?dl=0). Uncompress the archive and move all of its files into the CNS source directory (`PATH/TO/cns_solve_1.3/source/`).

## 3 Compile CNS

#### Start compilation

Make sure you are inside the main `cns_solve_1.3` directory.

If a suitable compiler is installed on your system, the following command should start the compilation of CNS:

```bash
make install
```

If needed, you can specify the compiler you want to use by adding the `compiler` option, for example:

```bash
make install compiler=gfortran
```

You should see the progress of the compilation, followed by a confirmation that an executable has been created. For example (with gfortran):

```
created executable file cns_solve-2206031450.exe
```

If this worked, you can continue to the [next step](#4-Source-CNS).

#### Install a compiler

If no suitable compiler is found, you will get an error message:

```
Error: no suitable compiler (ifort, pgf95, ifort_i4, ifort_mp or pgf95_mp) found
```

In this case, you will have to install a compiler. For more information, check the [CNS Installation page](http://cns-online.org/v1.3/).
One Fortran compiler that is easy (and free) to install is `gfortran`. It is part of the GNU Compiler Collection (GCC), which can be installed with one command if you have homebrew:

```bash
brew install gcc
```

Alternatively, gfortran can be downloaded as a [separate binary](https://gcc.gnu.org/wiki/GFortranBinaries).

#### Makefile

Check whether a corresponding `Makefile` is provided for your system in `PATH/TO/cns_solve_1.3/instlib/machine/supported/YOUR-SYSTEM/`. If no such file (eg. `Makefile.header.2.gfortran`) exists, create one manually.

Then try [`make install`](#Start-compilation) again, as described above.

## 4 Source CNS

Lastly, CNS needs to be loaded in the current shell before it can be used.

For csh/tcsh:

```bash
source cns_solve_env
```

For sh/bash:

```bash
source .cns_solve_env_sh
```

To avoid having to repeat this step each time you open a terminal window to run CNS/HADDOCK, you may want to add the relevant `source` line to your `~/.cshrc` / `~/.tcshrc` / `~/.bashrc` / `~/.bash_profile`.

## 5 Check installation

To check that CNS has been installed, start the program (from any working directory) with the following command:

```bash
cns
```

If the installation was successful, CNS will open in your terminal with a header including the following information:

```
============================================================
|                                                          |
|            Crystallography & NMR System (CNS)            |
|                         CNSsolve                         |
|                                                          |
============================================================
 Version: 1.3 at patch level U
 Status: Special UU release with Rg, paramagnetic
         and Z-restraints (A. Bonvin, UU 2013)
============================================================
```

The command to exit CNS is `STOP`.

You now have a suitable version of CNS to [install HADDOCK3](INSTALL.md).
