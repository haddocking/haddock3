# Editing and extending runs

The modularity of HADDOCK3 allows restarting halted runs from a certain step,
extending previous runs with additional steps, and using successful steps as
starting points for new runs. Let's see some examples.

## Restarting a run

Use the `--restart` flag to restart a HADDOCK3 run. For example, imagine you
start a run:

```
haddock3 my-run-config.cfg
```

Unfortunately, something went wrong in step 4 (the fourth module defined in the
configuration file).  You want now to continue the run from where it stopped.
Use the `--restart` followed by `N`, where `N` is an integer corresponding to
the step you which to restart from. **Remember**, haddock steps are 0-indexed
(0, 1, 2, 3, ...); so, `--restart 3` is the fourth step.

```
haddock3 my-run-config.cfg --restart 3
```

**Important:** As soon you run the above command, The `--restart` option will
**delete** step folders from `3` onward (inclusive), that is, `3_`, `4_`,
`...`, will be deleted.

## Restarting a run with modified parameters

You can profit from the `--restart` option described above to (re)run modules
with modified parameters. For example, if you noticed a parameter was set to
`true` when you wanted it to `false`, or if you wish to increase the sampling
for a given refinement module. Simply edit the configuration file modifying the
parameters for the modules you want to (re)run, in this case, from the fourth
module onward.

```
haddock3 my-edited-config.cfg --restart 3
```

**Important:** As soon you run the above command, The `--restart` option will
**delete** step folders from `3` onward (inclusive), that is, `3_`, `4_`,
`...`, will be deleted.

## Restarting a run with additional steps

You can profit from the `--restart` option to add steps to a run. For example,
a run has **five** steps and is completed successfully. You now realize you want
to perform an additional refinement (`mdref`) followed by CAPRI evaluation
(`caprieval`) and FCC clustering (`clustfcc`). Add the parameters for those
modules to the **original** configuration file, keeping all the original
parameters unchanged, and run haddock:

```
haddock3 config-with-additional-step.cfg --restart 5
```

We use `--restart 5` because the first new step is the sixth in the new
workflow: the five initial steps and three new steps. Remember `--restart` is
0-indexed.


## Extend a run

You can extend a successful run with additional steps. For that, prepare a
configurtation file containing only the new steps you wish to execute on top of
the previously successful run. In these cases, you **don't** need to define the
`run_dir` and `molecules` parameters in this new configuration file because
they will be ignored. You **can** define the other general parameters like
`ncores`, `mode`, etc. To extend a run with additional modules:

```
haddock3 new-steps.cfg --extend-run <run_dir>
```

## Starting new runs from successful steps

You can also start an independent run from a successful step of a previous run.
Consider the following successful run:

```
run1/
|--- 0_topoaa/
|--- 1_rigidbody/
|--- 2_caprieval/
|--- 3_seletop/
|--- 4_flexref/
|--- (etc...)
|--- data/
```

Now, you want to start a new independent run from the `4_flexref` step.

First, you need to copy the step to a new run folder using our `haddock3-copy` CLI.

```
haddock3-copy -r run1 -m 0 4 -o run2
```

The above command copies steps `0_topoaa` and `4_flexref` to a new `run2`
directory. It also updates the path references in the step folders. Resulting in:

```
run2/
|--- 0_topoaa/
|--- 1_flexref/
|--- data/
     |--- 0_topoaa/
     |--- 1_flexref/
```

Do **not** use the bash `cp` command to emulate this operation because there
several internal aspects treated by `haddock3-copy` that wouldn't be treated
by `cp`.

**Note:** If the new run uses CNS-dependent modules, you **also need** to copy
the folder corresponding to the initial topology creation (the `topoaa` module).

Second, create a configuration file for the new run containing **only** the
parameters of the new modules you wish to execute after the `flexref`. You
**don't** need to define the `run_dir` and `molecules` parameters in this new
configuration file because they will be ignored. You **can** define the other
general parameters like `ncores`, `mode`, etc. For example:

```toml
ncores = 40

[emref]
tolerance = 20
ambig_fname = "path/to/air.tbl"

[caprieval]
reference_fname = "path/to/reference.pdb"
```

Following the example, to start the new run:

```
haddock3 my-new-config.cfg --extend-run run2
```

## Additional considerations

You can now combine the functionalities explained above in many different
variations.

Below is some quick Q&A regarding these options:

1. Can I copy run directories to different folders on my computer?

Yes. HADDOCK3 uses relative paths inside the run directory.

2. Can I copy run directories to different computers?

Yes, provided you have installed HADDOCK3 following the [INSTALL](INSTALL.md)
instructions in the two systems. You can copy a run directory (or some of its
steps) to a different computer/system and (re)run it using the `--restart` and
`--extend-run` options.

3. Modules are 0-indexed. Is this related to Python being 0-indexed also?

No. It is because we defined the `topoaa` module as zero.

4. Where can I find additional help?

You can find additional help by running the command: `haddock3 -h` and reading
the parameters' explanations. Otherwise, ask us in the ["issues"
forum](https://github.com/haddocking/haddock3/issues).
