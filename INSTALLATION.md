# Installation instructions

## Requirements
 * Crystallography & NMR System (CNS)
    * Make a request on the [CNS website](http://cns-online.org/v1.3/), download the software and copy the contents of `/cns1.3` to CNS's source folder before compilation.
 * Python 3.7.x
 * gcc 4 (or higher)
 
## Optional: Third-party software
 * [Dcomplex](http://compbio.iupui.edu/group/4/pages/downloads)
 * [Fastcontact](http://structure.pitt.edu/software/FastContact/)
 * [DockQ](https://github.com/bjornwallner/DockQ/)
 * [ClustalO](http://www.clustal.org/omega/)
 * [ProFit](http://www.bioinf.org.uk/software/profit/)
 
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

## Running the example

```bash
$ cd examples/protein-protein

# change the number of processors (nproc) according to your system
$ vim run.toml

$ haddock3.py run.toml
```

## Results
Execution time will vary depending on your system and number of processors used. 

Two files will be created at `analysis/`, `fcc_0.6-4.out` with cluster-based results and `ss.stats`, single-structure information for each model.

Example cluster results:
```
Cluster 2 [-120.57, -127.92] -> (228) 0 5 54 70 1 13 10 102 
...
```
Cluster numbers are given according to how populated they are (Cluster 1 > Cluster 2) and sorted according to their mean value. First value inside `[ ]` is the mean score over all cluster elements, second element is the mean over the top 4 cluster elements.
Cluster element inside `( )` is the closest structure to the center of the cluster, each subsequent element is ranked according to its HADDOCK-score. 

Example single-structure results:

```
model ranking haddock-score total bonds angles improper dihe vdw elec air cdih coup sani vean dani desolv bsa cluster-0.6-4_name cluster-0.6-4_internal_ranking cluster-0.6-4_overall_ranking
water_refinement/complex_itw_000000.pdb 1 -152.997 -537.556 0.000 0.000 0.000 0.000 -49.403 -498.758 10.605 0.000 0.000 nan 0.000 0.000 -4.904 1544.040 2 1 1 
water_refinement/complex_itw_000005.pdb 2 -150.768 -495.467 0.000 0.000 0.000 0.000 -51.989 -454.567 11.089 0.000 0.000 nan 0.000 0.000 -8.975 1575.740 2 2 1
...
``` 