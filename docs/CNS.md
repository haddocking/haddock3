# CNS Installation

The computational engine for most of HADDOCK3's modules is the [Crystallography and NMR System (CNS)](http://cns-online.org/v1.3/).
Once you have [cloned the HADDOCK3 repository](INSTALL.md), please follow these instructions to compile a CNS executable that can be used to run HADDOCK.

## 1 Download CNS

Downloading CNS requires a license, which academic users can [request for free](http://cns-online.org) (click on the `Download` menu and fill the form). You will be sent a password for the download via e-mail.
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

            setenv CNS_SOLVE /PATH/TO/cns_solve_1.3
```

If you are using sh/bash, edit `.cns_solve_env_sh`:

```
# CHANGE THE NEXT LINE TO POINT TO THE LOCATION OF THE CNSsolve DIRECTORY

        CNS_SOLVE=/PATH/TO/cns_solve_1.3
```

## 2 Check that your current OS is recognized by the CNS installation scripts

From the CNS installation directory, type the following command:

```bash
./bin/getarch
```

For an unknown system you might see as output:

```
~/software/cns_solve_1.31-UU> ./bin/getarch
unknown-arm64-Darwin
```

The current list of supported architectures can be found by looking at the content of the `instlib/machine/supported` directory.
In those directories, you will find Makefiles for various compilers, including some for gfortran.

For convenience, we are providing with HADDOCK3 updated CNS files for the supported systems in which the gfortran Makefiles have already been edited.
Simply copy these to your CNS installation directory with:

```bash
cp -r ~/software/haddock3/varia/cns1.3/instlib ./ 
```

Also copy an updated `getarch` script that will recognize Mac M1/2/3/ processors:

```bash
cp ~/software/haddock3/varia/cns1.3/bin/getarch ./bin 
```

Make sure to specify the correct location of your haddock3 installation as the above command is only an example.


The recommended gfortran Makefile options are:

```
# fortran options
F77 = gfortran
F77STD = -fdefault-integer-8 -w -fallow-argument-mismatch
F77OPT = -O3 $(CNS_MALIGN_I86) -funroll-loops -ffast-math -march=native -mtune=native
F77FLAGS = $(F77STD) $(F77OPT) $(EXT_F77FLAGS) $(F77BUG)

# C options
CC = gcc
CPP = g++
CCFLAGS = -O -DINTEGER='long int' -DCNS_ARCH_TYPE_$(CNS_ARCH_TYPE) $(EXT_CCFLAGS)

# link options
LD = gfortran
LDFLAGS = -w $(EXT_LDFLAGS) -static-libgfortran
```

If your OS/hardware is unknown, this is not per se a problem, but it means the compiler options might not be properly defined.
For example, the default CNS package will not recognize the Mac M1/2/3/ processors (but after copying the CNS files distributed with haddock3 it should!). To solve that edit the `./bin/get_arch` script and add the following lines after line 179 (this should be added after a line containin `exit 0`):

```
    arm64:Darwin:*:*)
    echo mac-arm64-darwin
    exit 0 ;;
```

If you OS/hardware is not in the list of supported systems, you will need to define it in the `getarch` script as just explained and create a directory under `instlib/machine/supported` with the same name as the one you defined in `getarch`. Copy then and edit if needed the files from another supported system into the directory you created.


## 3 Add HADDOCK routines

HADDOCK requires some adjustments in CNS to function properly.
The set of routines to add is provided in the haddock3 repository under [`varia/cns1.3/`](/varia/cns1.3/README.md).
Copy all of these files into the CNS source directory (`/PATH/TO/cns_solve_1.3/source/`), for example:

```bash
cp ~/software/haddock3/varia/cns1.3/* source/
```

Make sure to specify the correct location of your haddock3 installation as the above command is only an example.


## 4 Compile CNS


Make sure you are inside the main `cns_solve_1.3` directory. 
Source first the `cns_solve_env` (for csh/tcsh) or `.cns_solve_env_sh` script (for bash).

For csh/tcsh:

```bash
source ./cns_solve_env
```

For bash:

```bash
source ./.cns_solve_env_sh
```


To compile CNS we recommend using the gfortran compiler. The following command should start the compilation of CNS:

```bash
make install compiler=gfortran
```

You should see the progress of the compilation, followed by a confirmation that an executable has been created. For example (with gfortran):

```
created executable file cns_solve-2206031450.exe
```

If this does not work, check the [troubleshooting section](#Troubleshooting) below.


## 5 Source CNS

Lastly, CNS needs to be loaded in the current shell before it can be used.

For csh/tcsh:

```bash
source ./cns_solve_env
```

For sh/bash:

```bash
source ./.cns_solve_env_sh
```

## 6 Check installation

To check that CNS has been installed, start the program with the following command:

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

You now have a suitable CNS executable to finish [installing HADDOCK3](INSTALL.md).

***


## 7 Installing the created CNS executable into haddock3

You can find the freshly compiled CNS executable a newly created directory in the CNS installation directory named after the architecture of your system. E.g. for a M1/2/3/ Mac Arm processor it will be `mac-arm64-darwin`. In the `source` directory in that directory you can find the compiled CNS executable.

```bash
ls mac-arm64-darwin/source/*exe 
> mac-arm64-darwin/source/cns_solve-2405221019.exe
```

You can then copy this executable into the `haddock3/bin` directory, renaming it simply to `cns`



# Troubleshooting

## Compiler

If no suitable compiler is found, you will get an error message:

```
Error: no suitable compiler (ifort, pgf95, ifort_i4, ifort_mp or pgf95_mp) found
```

In this case, you will have to install a compiler. For more information, check the [CNS Installation page](http://cns-online.org/v1.3/installation/frame.html).
One Fortran compiler that is easy (and free) to install is `gfortran`. It is part of the GNU Compiler Collection (GCC), which can be installed with one command if you have homebrew:

```bash
brew install gcc
```

Alternatively, gfortran can be downloaded as a [separate binary](https://gcc.gnu.org/wiki/GFortranBinaries).

## Makefile

If a suitable compiler is installed but no corresponding `Makefile.header` is found, `make install` will give an error as well. For example (with gfortran):

```
Error: Makefile template for compiler gfortran is not available
```

Makefile headers for each supported system-compiler combination are provided in `instlib/machine/supported/` (for unsupported systems, see `instlib/machine/unsupported/`). If you already tried `make install`, the appropriate directory for your system (eg. `mac-intel-darwin`) will have been copied to the main `cns_solve_1.3` directory.
If this directory does not contain the desired file (eg. `Makefile.header.2.gfortran`), you can create one manually.

For example, for `gfortran`, create a new file named `Makefile.header.6.gfortran` inside the directory containing Makefile headers for your system. Save the following lines in this file:

```
# fortran options
F77 = gfortran
F77STD = -fdefault-integer-8 -w -fallow-argument-mismatch
F77OPT = -O3 $(CNS_MALIGN_I86) -funroll-loops -ffast-math -march=native -mtune=native
F77FLAGS = $(F77STD) $(F77OPT) $(EXT_F77FLAGS) $(F77BUG)

# C options
CC = gcc
CPP = g++
CCFLAGS = -O -DINTEGER='long int' -DCNS_ARCH_TYPE_$(CNS_ARCH_TYPE) $(EXT_CCFLAGS)

# link options
LD = gfortran
LDFLAGS = -w $(EXT_LDFLAGS) -static-libgfortran
```

Note the `-fallow-argument-mismatch` flag. From gfortran version 10 onward, type mismatches by default give an error rather than a warning. This flag returns warnings instead, so that CNS compilation can successfully complete. If a `Makefile.header.x.gfortran` is provided for your system but compilation fails with (many) `Error: Type mismatch`, you may need to add this flag to the provided file.
