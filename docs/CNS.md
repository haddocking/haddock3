# CNS Installation

The computational engine for most of HADDOCK3's modules is the [Crystallography and NMR System (CNS)](http://cns-online.org/v1.3/).
Once you have [cloned the HADDOCK3 repository](INSTALL.md), please follow these instructions to compile a CNS executable that can be used to run HADDOCK.

## 1. Downloading CNS

Downloading CNS requires a license, which academic users can [request for free](http://cns-online.org).

To Download CNS, go to their [website](http://cns-online.org) and  click on the `Download` menu and fill the form. You will be sent a password for the download via e-mail - keep in mind only education e-mails are allowed for this.

You will recieve a login/password in the e-mail, once you have those go to the [download page](http://cns-online.org/download/v1.3/cns_solve_1.3_all.tar.gz) and use the provided login/password to download the CNS package.

> The CNS website is not currently using SSL certificates, this mean your download might be blocked or it might give you a warning, double check in your browser.

Once the archive `cns_solve_1.3_all.tar.gz` has been downloaded, move it to a `$HOME/software` directory (or any other directory of your choice) and uncompress it:

```bash
mkdir ~/software
mv ~/Downloads/cns_solve_1.3_all.tar.gz ~/software/
tar -zxvf cns_solve_1.3_all.tar.gz
```

The command above will have create a `cns_solve_1.3` directory in the `~/software` directory, containing the CNS source code.

## 2. Configuring CNS

Go to the CNS directory:

```bash
cd cns_solve_1.3
```

Inside this directory there will be a `.cns_solve_env.sh` file.

> Mind the `.` in front of of the filename, this means that this file is "hidden", you can only see it if you type `ls -a`

You need to edit this file to point to the current directory.

You will need to change line 18 of this file to point to the current directory.

You can do that either by opening this file in a code editor such as VScode or vim, or by running the command below:

```bash
sed -i 's@_CNSsolve_location_@'"$PWD"'@' .cns_solve_env_sh
```

Check if the substitution worked:

```bash
$ sed -n 18p .cns_solve_env_sh
```

## 3. Patch CNS with HADDOCK custom files

For HADDOCK to use CNS, some modifications were made to the CNS source code. We provide those modifications together with the haddock3 source files.

Copy the CNS files distributed with HADDOCK3 to the CNS directory:

```bash
cp -rv ~/software/haddock3/varia/cns1.3/bin ~/software/cns_solve_1.3/bin
cp -rv ~/software/haddock3/varia/cns1.3/instlib/* ~/software/cns_solve_1.3/
cp -v ~/software/haddock3/varia/cns1.3/*.f ~/software/cns_solve_1.3/source/
cp -v ~/software/haddock3/varia/cns1.3/*.inc ~/software/cns_solve_1.3/source/
```

> You cannot use CNS in HADDOCK without this patched files!

## 4. Compiling CNS

Before installing CNS, you need to check wether your operating system is supported by CNS.

The list of supported architecthures can be seen with the following command:

```bash
$ ls ~/software/cns_solve_1.3/instlib/machine/supported
arm-aarch64-linux  intel-x86_64bit-linux  linux  mac-arm-darwin  mac-intel-darwin
```

To check what is *your* architecthure, type the following command:

```bash
$ ~/software/cns_solve_1.3/bin/getarch
intel-x86_64bit-linux
```

> You might get something like `unknown-arm64-Darwin`, in that case installation is still possible, however it means the compiler options might not be properly defined.

Several `Makefile` headers are already provided - here we will focus on the ones that use `gfortran` as the compiler.

See them with:

```bash
ls ~/software/cns_solve_1.3/instlib/machine/supported/*/Makefile*gfortran*
```

The recommended gfortran Makefile options are defined in the `Makefile.header.7.gfortran` file. If you are using a different system, you can copy this file to the main CNS directory:

```makefile
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
<!--
For example, the default CNS package will not recognize the Mac M1/2/3/ processors (but after copying the CNS files distributed with haddock3 it should!). To solve that edit the `./bin/get_arch` script and add the following lines after line 179 (this should be added after a line containin `exit 0`):

```
    arm64:Darwin:*:*)
    echo mac-arm64-darwin
    exit 0 ;;
``` -->
<!--
If you OS/hardware is not in the list of supported systems, you will need to define it in the `getarch` script as just explained and create a directory under `instlib/machine/supported` with the same name as the one you defined in `getarch`. Copy then and edit if needed the files from another supported system into the directory you created.

## 3. Add HADDOCK routines

HADDOCK requires some adjustments in CNS to function properly.
The set of routines to add is provided in the haddock3 repository under [`varia/cns1.3/`](/varia/cns1.3/README.md).
Copy all of these files into the CNS source directory (`/PATH/TO/cns_solve_1.3/source/`), for example:

```bash
cp ~/software/haddock3/varia/cns1.3/* source/
```

Make sure to specify the correct location of your haddock3 installation as the above command is only an example. -->


Before you start with the compilation, you will need to source the `.cns_solve_env_sh` file, to do that run the command:


```bash
source ~/software/cns_solve_1.3/.cns_solve_env_sh
```

In case your system is not recognized refer to the troubleshooting section below.

To compile CNS we recommend using the gfortran compiler. The following command should start the compilation of CNS:

```bash
make install compiler=gfortran
```

You should see the progress of the compilation, followed by a confirmation that an executable has been created. For example (with gfortran):

```
created executable file cns_solve-2206031450.exe
```

If this does not work, check the [troubleshooting section](#Troubleshooting) below.

## 4. Source CNS

Lastly, CNS needs to be loaded in the current shell before it can be used.

For csh/tcsh:

```bash
source ./cns_solve_env
```

For sh/bash:

```bash
source ./.cns_solve_env_sh
```

## 5. Check installation

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

## 6. Installing the created CNS executable into haddock3

You can find the freshly compiled CNS executable a newly created directory in the CNS installation directory named after the architecture of your system. E.g. for a M1/2/3/ Mac Arm processor it will be `mac-arm64-darwin`. In the `source` directory in that directory you can find the compiled CNS executable.

```bash
ls mac-arm64-darwin/source/*exe
> mac-arm64-darwin/source/cns_solve-2405221019.exe
```

You can then copy this executable into the `haddock3/bin` directory, renaming it simply to `cns`

# Troubleshooting

## 1. System not recognized

If the installation fails because your OS/hardware are not recognized, most likely the `bin/getarch` script needs to be updated.

From withon your CNS installation directory type the following command to check if your system is recognized:

```bash
./bin/getarch
```

For an unknown system you might see as output:

```
~/software/cns_solve_1.31-UU> ./bin/getarch
unknown-arm64-Darwin
```

If your OS/hardware is unknown, this is not per se a problem, but it means the compiler options might not be properly defined.
For example, the default CNS package will not recognize the Mac M1/2/3/ processors (but after copying the CNS files distributed with haddock3 it should!).
To solve that edit the `./bin/get_arch` script and add lines such as your system be recognized (e.g. after line 179 - this should be added after a line containin `exit 0`). Here is an example to add support for the Mac M1/2/3 processors:

```
    arm64:Darwin:*:*)
    echo mac-arm64-darwin
    exit 0 ;;
```

You should then create a directory under `instlib/machine/supported` with the same name as the one you defined in `getarch`. Copy then and edit if needed the files from another supported system into the directory you created.

## 2. Compiler

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

## 3. Makefile

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
