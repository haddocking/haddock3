# Continuing and editing runs

The modularity of HADDOCK3 allows users to take successful steps
(modules) of a workflow as starting points for new runs. There are
several ways users can continue or expand previous runs. The index below
lists those possibilities:

1. Restart a previous run from a given step
1. Expanding a previous run from a given step
1. Editing and expanding a previous run from a given step
1. Copy a run to a new folder
1. Extend a run with a partial config
1. Use a module as a seed for a new run
1. Final considerations

# Restart a previous run from a given step

This option allows you to restart a run from a specific step. It is useful when
you want to repeat part of the run or to continue a broken run without needing
to rerun all successful steps.

Let's imagine you execute a run with six steps (modules).

```
haddock3 my-haddock3-config-file.cfg
```

If the run completes successfully, the run directory shows as follows:

```
0_topoaa/
1_rigidbody/
2_caprieval/
3_seletop/
4_flexref/
5_caprieval/
data
log
```

Now, imagine that later in your analysis, you notice that you used wrong parameters for the flexible refinement (`4_flexref`) stage, setting `randremoval=true` instead of `randremoval=false`. Therefore, you want to repeat the run only for `4_flexref` onwards without repeating the previous successful steps. For that, you should do the following:

```
haddock3 my-haddock3-config-file.cfg --restart 4
```

This operation will delete modules' folders `4_` onwards and repeat those steps
following the information in the configuration file. You must run the `haddock3`
command from within the same directory as you did for the original run.

# Expanding a previous run from a given step

The `--restart` functionality explained above allows some exciting variations.
We will see here how to expand a previous run with additional modules.

Run 1 completed successfully, but now we want to complement our run with a
`mdref` step and perform a final CAPRI score evaluation. For that, add the
parameters for the new modules to the end of the original configuration file.
You must keep the rest of the configuration file the same.

Navigate back to the same folder you executed the previous runs and:

```
haddock3 expanded-configuration-file.cfg --restart 6
```

In this case, the `--restart` option takes the value `6` because we want to
continue the run from step 6, the new `mdref` step.

# Editing and expanding a previous run from a given step

Following the same rationale from the previous example, we can edit and expand a
previous run. For that, edit the original configuration file by altering
(add/edit/remove) the information for the new modules you want to execute. For
example, delete steps `4_` onwards and perform a `mdref` only. You should delete
the parameters for those modules in the original configuration file and replace
them by the parameters for a `mdref` execution. In the newly edited
configuration file, the`mdref` is now the `4_` module, so you should restart
from four.

```
haddock3 edit-config-file.cfg --restart 4
```

This setup would result in:

```
0_topoaa/
1_rigidbody/
2_caprieval/
3_seletop/
4_mdref
data
log
```

As you can see, we have reused the first four modules, up to `3_seletop`, and
replaced the `4_flexref` and `5_caprieval` with a `4_mdref` module. You could
execute as many modules as you wish after `4_mdref`.

# Copy a run to a new folder

Importantly to know, you can execute all previous examples on copies of the
original run. There is no need to rewrite the initial run if you wish to keep
it. For that, copy the original run to a new folder:

```
cp -r run1 run2
```

In the configuration file, edit the `run_dir` parameter to `run2` so it
matches with the newly created folder. Then, proceed by editing the
configuration file as discussed in the previous examples. Finally,

```
haddock3 new-config-file-run2.cfg --restart N
```

where `N` is the step you wish to start from, as explained previously. Notably,
you do not need to copy the whole folder when copying the original run. Instead,
you can copy only the modules folders you wish to maintain plus the `data`
folder.

# Extend a run with a partial config

HADDOCK3 allows you to extend a run without previous knowledge of the run
configuration and using just a partial configuration file.
For example, you have several runs for which you **now** want to
perform a CAPRI evaluation at the last step. How to do this?

Create a new configuration file containing just the `caprieval` module, or
whatever and many modules you desire. For example:

```toml
[caprieval]
reference_fname = "reference-structure.pdb"
```

Now, run the configuration file by giving the `--restart` flag the path of the
previous run directory:

```
haddock3 caprieval.cfg --restart <previous-run-directory>
```

In this case, the configuration file does not need to contain any information
about the molecules used as input or the previous modules. HADDOCK3 will just
execute the new modules on top of the previous run.

# Use a module as a seed for a new run

Imagine you want to take a specific module of a successful run as the starting
point of a new run. Copy the `data`, the `topoaa`, and the desired module to a
new folder. You can also use the `topoaa` as the initial module. At this point,
you can start the new run using both strategies we have detailed so far:

1) edit the original configuration file keeping only the information for the
modules you used as seed and editing the `run_dir` parameter to the new folder,
and execute:

```
haddock3 new-config.cfg --restart 2
```

If your seed module is directly the `topoaa` module, use `--restart 1`.

2) Create a new configuration file containing only the parameters for the
modules you want to run and execute:

```
haddock3 new-config.cfg --restart <run_dir>
```

where the `run_dir` is the new directory containing the seeding modules.

# Final considerations

When extending or continuing a run from a new configuration file, all the
modules' folders and their content will be updated to the corresponding position
and leading zeros. For example, the new takes steps 0 and 4 from a run with more
than 10 modules (note the leading double digit numbers):

```
00_topoaa/
04_flexref/
data/
```

If you extend this run with a series of `caprieval`, `emref`, `caprieval`, and
`clustfcc`, the final output will be:

```
0_topoaa/
1_flexref/
2_caprieval/
3_emref/
4_caprieval/
5_clustfcc/
```

Where `0_topoaa` and `1_flexref` are the original folders renamed. The files
within these folders were also edited to ensure paths compatibility and the
consistency of paths and the run as a whole.
