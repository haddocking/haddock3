# Continuing and editing runs

The modularity of HADDOCK3 allows users to take successful steps
(modules) of a workflow as starting points for new runs. There are
several ways users can continue or expand previous runs. The index below
lists those possibilities:

1. Restart, editing, and extending a previous run from a given step
1. Copy a run to a new folder, and continue from there
    1. Use a single step as a starting point for a new run
1. Extend a run with new modules
1. Final considerations

# Restart, editing, and extending a previous run from a given step

This option allows you to restart a run from a specific step. It is
useful when you want to repeat part of the run or to continue a broken
run without needing to rerun all successful steps.

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

Now, imagine that later in your analysis, you notice that you used the
wrong parameters for the flexible refinement (`4_flexref`) stage,
setting `randremoval=true` instead of `randremoval=false`. Therefore,
you want to repeat the run only for `4_flexref` onwards without
repeating the previous successful steps. For that, edit the
configuration file according to your new preferences, and run:

```
haddock3 my-haddock3-config-file.cfg --restart 4
```

This operation will delete modules' folders `4_` onwards and repeat
those steps following the configurations in the (new) configuration
file. You must run the `haddock3` command from within the same directory
as you did for the original run.

As you can see, this option allows you to **restart** or **repeat** a
run from a specific step using the **original** configuration file or to
**modify** a run by updating the parameters for the steps you want to
rerun.

Likewise, we can use the `--restart` functionality to **extend** a run
with additional steps. For example, Imagine `run 1` completed
successfully. Still, you now want to complement it with an `emref` step,
perform a final CAPRI score evaluation, and finalize the run by
clustering with FCC (`clustfcc` module). For that, add the parameters
for the new modules to the end of the original configuration file. The
original steps' order in the configuration file must be kept the same
for HADDOCK3 to understand *the previous* with the *new* modules.

Navigate back to the same folder you executed the previous runs and:

```
haddock3 expanded-configuration-file.cfg --restart 6
```

In this case, the `--restart` option takes the value `6` because we want to
continue the run from step 6, the new `emref` step. This would result in:

```
0_topoaa/
1_rigidbody/
2_caprieval/
3_seletop/
4_flexref/
5_caprieval/
6_emref/
7_caprieval/
8_clustfcc/
data
log
```

As you can see, you can use `--restart` to edit and rerun modules, to
add new modules, **or do both at the same time**.

# Copy a run to a new folder, and continue from there

You can execute all previous examples on copies of the original run.
There is no need to rewrite the initial run if you wish to keep it. For
that, copy the original run to a new folder:

```
cp -r run1 run2
```

You can also copy the run to different folders because HADDOCK3 treats
internal paths as relative paths. Therefore, modules can communicate
with each other regardless of the absolute location of the run
directory. The same applies to the installation requirements. Provided
the installation process was the same, you can copy HADDOCK3 run over
different machines. 

In the configuration file, edit the `run_dir` parameter to `run2` to
match the newly created folder. `run1` and `run2` are arbitrary names.

**Important:** The `run_dir` path must point to the new run directoy
**relatively** from where you will run the configuration file.
Otherwise, use an absolute path.

Then, proceed by editing the configuration file as discussed in the
previous examples. Finally, run:

```
haddock3 new-config-file-run2.cfg --restart N
```

where `N` is the step you wish to start from, as explained previously.

Notably, you do not need to copy the whole folder when copying the
original run. Instead, you can copy only the steps' folders you wish to
maintain plus the `data` folder. For example, if the original run has
ten modules, but you wish to keep only the first three, you can copy
only those first three modules to the new run folder.

## Use a single step as a starting point for a new run

Following the example above, you can copy a single step and use it as a
starting point for a new run. For example, taking the result of an
extensive rigid body sampling and exploring different refinement
strategies separately.

To use a single successful step (module) as a starting point for a new
run, copy the `data`, the `topoaa`, and the desired step folder to a new
run directory folder. The selected step does **not** need to be
consecutive to the topology module. It can be any step.  You can also
use the `topoaa` as the initial module. For example,

```bash
mkdir run2
cp -r run/data run1/0_topoaa run1/4_flexref/ run2
```

At this point, you can start the new run by editing the original
configuration file. You must keep the information for the modules you
want to use as starting points for the new run and edit the `run_dir`
parameter to match the new run directory name. Finally, execute:

```
haddock3 new-config.cfg --restart 2
```

If your starting step was the `topoaa` module, use `--restart 1`.

# Extend a run with new modules

On top of the previous examples using the  `--restart N` option, where
`N` is an integer, and variations of the original configuration file,
HADDOCK3 allows you to extend a run from a configuration file defining
only the **new** modules.

For example, imagine you have several completed runs for which you
**now** want to perform a CAPRI evaluation at the last step, or apply a
new scoring or clustering module recently released or a different new
minimization protocol.

How can you execute new modules on top of successful runs without
referring to their original configuration file?

Create a new configuration file containing **only** the `caprieval`
module, or whatever new and many modules you desire. For example:

```toml
[caprieval]
reference_fname = "reference-structure.pdb"
```

Now, run the configuration file by giving the `--restart` flag the path
of the previous run directory, instead of an integer `N`:

```
haddock3 caprieval.cfg --restart <previous-run-directory>
```

In this case, the configuration file does not need to contain any
information about the molecules used as input or the previous modules.
Instead, HADDOCK3 will execute the new modules following the run's last
module.

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
