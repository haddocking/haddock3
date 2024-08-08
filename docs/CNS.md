# CNS Installation

The computational engine for most of HADDOCK3's modules is the [Crystallography and NMR System (CNS)](http://cns-online.org/v1.3/).
Once you have [cloned the HADDOCK3 repository](INSTALL.md), please follow these instructions to compile a CNS executable that can be used to run HADDOCK.

Before starting, make sure you download the system requirements, they are `gcc`, `gfortran` and `csh`. Luckily they are easily available from the package manager, so if you are using Ubuntu (for example):

```bash
sudo apt-get install gcc gfortran csh
```

> Make sure `gfortran` is in your `$PATH`, check it with `which gfortran`.

## 1. Downloading CNS

Downloading CNS requires a license, which academic users can [request for free](http://cns-online.org).

To Download CNS, go to their [website](http://cns-online.org) and  click on the `Download` menu and fill the form. You will be sent a password for the download via e-mail - keep in mind only education e-mails are allowed for this.

You will receive a login/password in the e-mail, once you have those go to the [download page](http://cns-online.org/download/v1.3/cns_solve_1.3_all.tar.gz) and use the provided login/password to download the CNS package - the archive you need is: `cns_solve_1.3_all.tar.gz`.

> The CNS website is not currently using SSL certificates, this mean your download might be blocked or it might give you a warning, double check in your browser.

Once the archive `cns_solve_1.3_all.tar.gz` has been downloaded, move it to a `$HOME/software` directory (or any other directory of your choice) and decompress it:

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

Inside this directory there will be a `cns_solve_env` file.

You need to edit this file to point to the current directory.

You will need to change line 18 of this file to point to the current directory.

You can do that either by opening this file in a code editor such as VScode or vim, or by running the command below:

```bash
sed -i 's@_CNSsolve_location_@'"$PWD"'@' cns_solve_env
```

> **Note:** On a Mac, the default sed version will give an error with the above command. Install instead `gnu-sed` using `brew`, and then replace `sed` by `gsed` in the above command.

Check if the substitution worked:

```bash
sed -n 18p cns_solve_env
```

The expected output will be (with the path to your CNS installation directory):

```shell
        setenv CNS_SOLVE '/Absolute/path/where/cns_solve/is/cns_solve_1.3'
```

## 3. Patch CNS with HADDOCK custom files

For HADDOCK to use CNS, some modifications were made to the CNS source code. We provide those modifications together with the haddock3 source files.

> ❗ You cannot use CNS in HADDOCK without these patched files❗

Copy the CNS files distributed with HADDOCK3 to the CNS directory:

```bash
\cp -r ~/software/haddock3/varia/cns1.3/[bis]*  ~/software/cns_solve_1.3/
```

> **Note:** Make sure to use the correct path to your haddock3 location in the above command.

## 4. Compiling CNS

Before installing CNS, you need to check wether your operating system is supported by CNS.

The list of supported architectures can be seen with the following command:

```bash
ls ~/software/cns_solve_1.3/instlib/machine/supported
arm-aarch64-linux  intel-x86_64bit-linux  linux  mac-arm-darwin  mac-intel-darwin
```

To check what is *your* architecture, type the following command:

```bash
~/software/cns_solve_1.3/bin/getarch
```

> You might get something like `unknown-...`, in that case installation is still possible, however it means the compiler options might not be properly defined - check the Troubleshooting section

Several `Makefile` headers are already provided - here we will use on the one that uses `gfortran` as the compiler.

Now it is time to compile CNS, you can do that with the following command:

```bash
make install compiler=gfortran
```

You should see the progress of the compilation, followed by a confirmation that an executable has been created.

The executable will be located in a new architecture-specific directory.

For example under Linux (Intel or AMD processors) it will be under `~software/cns_solve_1.3/intel-x86_64bit-linux/source/cns_solve-XXXXXXX.exe`
And on a Mac with arm processors (M1/2/3) is will be under `~software/cns_solve_1.3/mac-arm64-darwin/source/cns_solve-XXXXXXX.exe`

> The `XXXXXXX` will be a number specific to your compilation.

Check if it is working by executing it:

```bash
~/software/cns_solve_1.3/intel-x86_64bit-linux/source/cns_solve-2405221906.exe
```

If you see the CNS prompt, the compilation was successful! 🎉

```text
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

The command to exit CNS is `stop`!

If this does not work, check the [troubleshooting section](#troubleshooting) below.

> **Note:** If you want to use CNS as standalone application, you should first source the `cns_solve_env` file (csh) or `.cns_solve_env_sh` file (bash). CNS will then be defined in you path and can be simply started using `cns` as command. If using bash, make sure to first edit the `cns_solve_env` variable in the `.cns_solve_env_sh` file to define the installation directory as described above.

## 5. Installing the created CNS executable into haddock3

You can find the freshly compiled CNS executable a newly created directory in the CNS installation directory named after the architecture of your system. E.g. for a M1/2/3/ Mac Arm processor it will be `mac-arm64-darwin`. In the `source` directory in that directory you can find the compiled CNS executable.

```bash
ls mac-arm64-darwin/source/*exe
> mac-arm64-darwin/source/cns_solve-2405221019.exe
```

You can then copy this executable into the `haddock3/bin` directory, renaming it simply to `cns`.

## Troubleshooting

### 1. System not recognized

If the installation fails because your OS/hardware are not recognized, most likely the `bin/getarch` script needs to be updated.

From within your CNS installation directory type the following command to check if your system is recognized:

```text
./bin/getarch
```

For an unknown system you might see as output:

```text
~/software/cns_solve_1.31-UU> ./bin/getarch
unknown-...
```

If your OS/hardware is unknown, this is not per se a problem, but it means the compiler options might not be properly defined.
For example, the default CNS package will not recognize the Mac M1/2/3/ processors (but after copying the CNS files distributed with haddock3 it should!).
To solve that edit the `./bin/get_arch` script and add lines such as your system be recognized (e.g. after line 179 - this should be added after a line containing `exit 0`). Here is an example to add support for the Mac M1/2/3 processors:

```shell
    arm64:Darwin:*:*)
    echo mac-arm64-darwin
    exit 0 ;;
```

You should then create a directory under `instlib/machine/supported` with the same name as the one you defined in `getarch`. Copy then and edit if needed the files from another supported system into the directory you created.

### 2. Compiler

If no suitable compiler is found, you will get an error message:

```text
Error: no suitable compiler (ifort, pgf95, ifort_i4, ifort_mp or pgf95_mp) found
```

In this case, you will have to install a compiler. For more information, check the [CNS Installation page](http://cns-online.org/v1.3/installation/frame.html).
One Fortran compiler that is easy (and free) to install is `gfortran`. It is part of the GNU Compiler Collection (GCC), which can be installed with one command if you have homebrew:

```bash
brew install gcc
```

Alternatively, gfortran can be downloaded as a [separate binary](https://gcc.gnu.org/wiki/GFortranBinaries).

### 3. Makefile

If a suitable compiler is installed but no corresponding `Makefile.header` is found, `make install` will give an error as well. For example (with gfortran):

```text
Error: Makefile template for compiler gfortran is not available
```

Makefile headers for each supported system-compiler combination are provided in `instlib/machine/supported/` (for unsupported systems, see `instlib/machine/unsupported/`). If you already tried `make install`, the appropriate directory for your system (eg. `mac-intel-darwin`) will have been copied to the main `cns_solve_1.3` directory.
If this directory does not contain the desired file (eg. `Makefile.header.2.gfortran`), you can create one manually.

For example, for `gfortran`, create a new file named `Makefile.header.6.gfortran` inside the directory containing Makefile headers for your system. Save the following lines in this file:

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

Note the `-fallow-argument-mismatch` flag. From gfortran version 10 onward, type mismatches by default give an error rather than a warning. This flag returns warnings instead, so that CNS compilation can successfully complete. If a `Makefile.header.x.gfortran` is provided for your system but compilation fails with (many) `Error: Type mismatch`, you may need to add this flag to the provided file.
