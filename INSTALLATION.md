# Installation instructions

## Requirements
 * Crystallography & NMR System (CNS)
    * Make a request on the [CNS website](http://cns-online.org/v1.3/), 
    * It must be compiled with the contents of `cns1.3/` (provided)
 * Python 3.7.x
 * gcc 4 (or higher)
 
## Optional: Third-party software
These are used in the Analysis step
 * [Dcomplex](http://compbio.iupui.edu/group/4/pages/downloads)
 * [Fastcontact](http://structure.pitt.edu/software/FastContact/)
 * [DockQ](https://github.com/bjornwallner/DockQ/)


 
Edit `haddock3/haddock/etc/haddock3.ini` accordingly

 # Step-by-step
 
 Download and uncompress the [latest stable release](https://github.com/haddocking/haddock3/releases).
 

 ```bash
$ cd haddock3
$ pip install -r requirements.txt

# Add the paths to the HADDOCK3 and CNS executable
$ vim bin/activate
$ source bin/activate

# add it to your enviroment
$ cat bin/activate >> ~/.bashrc
or
$ cat bin/activate >> ~/.zshrc
or
$ cat bin/activate.csh >> ~/.cshrc

$ cd haddock/src
$ make

# done (:
```

check the [examples documentation](examples/README.md) for instructions on how to execute HADDOCK 3 