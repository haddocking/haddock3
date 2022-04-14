# Editing and extending runs

The modularity of HADDOCK3 allows restarting halted runs from a certain step,
extending previous runs with additional steps, and using successful steps as
starting points for new runs. Let's see some examples.

## Restarting a run

Use the'- restart' flag to restart a HADDOCK3 run. For example, imagine you
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

## Restarting a run with modified parameters

You can profit from the `--restart` option described above to (re)run modules
with modified parameters. For example, if you noticed a parameter was set to
`true` when you wanted it to `false`, or if you wish to increase the sampling
for a given refinement module. Simply edit the configuration file modifying the
parameters for the modules you want to (re)run, in this case, from the fourth
module onwards.

```
haddock3 my-run-config.cfg --restart 3
```

## Restarting a run with additional steps

You can profit from the `--restart` option to add steps to a run.  For example,
a run has **five** steps and is completed successfully. You now realize you want
to perform an additional refinement (`mdref`) followed by CAPRI evaluation
(`caprieval`) and FCC clustering (`clustfcc`). Add the parameters for those
modules to the configuration file and run haddock:

```
haddock3 my-edited-config.cfg --restart 5
```

We use `--restart 5` because the first new step is the sixth in the new
workflow: the five initial steps and three new steps. Remember `--restart` is
0-indexed.

## Starting new runs from previous modules

You can also start a new independent run from a successful step of a previous
run. First, create a new run directory and copy the desired step folder from the
successful run, for example, `4_flexref`. **Note:** If the new run uses
CNS-dependent modules, you **also** need to copy the folder corresponding to the
initial topology creation (the `topoaa` module).

Create a configuration file for the new run. **Note:** In this case, you
**don't** need to define the `run_dir` and `molecules` parameters because these
will be taken from the initial step folder. You should not define the parameters
for the module you took as the starting point.

To start the new run:

```
haddock3 my-new-config.cfg --restart-from-run <path-to-run-dir>
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
``--restart-from-run`` options.

3. Modules are 0-indexed. Is this related to Python being 0-indexed also?

No. It is because we defined the `topoaa` module as zero.

4. Where can I find additional help?

You can find additional help by running the command: `haddock3 -h` and reading
the parameters' explanations. Otherwise, ask us in the ["issues"
forum](https://github.com/haddocking/haddock3/issues).
