# Running the example

The current example is a simplistic protein-protein docking using the true-interface as ambiguous restraints.

```bash
$ cd examples/protein-protein

# change the number of processors (nproc) according to your system
$ vim run.toml

$ haddock3.py run.toml
```

## Results
Execution time will vary depending on your system and number of processors used. 

The main analysis consists of the HADDOCK-score, its components and the FCC clustering of the complexes. Extra steps can be performed by third-party software, if installed and configured.

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
water_refinement/complex_itw_000000.pdb 1 -152.997 -537.556 0.000 0.000 0.000 0.000 -49.403 -498.758 10.605 0.000 0.000 0.000 0.000 0.000 -4.904 1544.040 2 1 1 
water_refinement/complex_itw_000005.pdb 2 -150.768 -495.467 0.000 0.000 0.000 0.000 -51.989 -454.567 11.089 0.000 0.000 0.000 0.000 0.000 -8.975 1575.740 2 2 1
...
``` 