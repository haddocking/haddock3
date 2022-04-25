# Usage

Before using HADDOCK3, install HADDOCK3 following the [installation
instructions](INSTALL.md).

## Run HADDOCK3

This is a minimal description of how to run HADDOCK3. We are currently
creating detailed documentation and tutorials. The minimal HADDOCK3
execution command is:

```bash
haddock3 <CONFIG WORKFLOW>
```

Where the `<CONFIG WORKFLOW>` is a haddock3 configuration file
describing the simulation workflow.

You can inspect all `haddock3` options with:

```bash
haddock3 -h
```

Inside the HADDOCK3 main folder, there is an `examples/` folder. There
you will find several subfolder with examples for specific docking
scenarios. Please experiment running those configuration files ending in
`-test.cfg`. For example:

```bash
cd examples/docking-protein-protein
haddock3 docking-protein-protein-test.cfg
```

Feel free to open the `cfg` files with your preferred test editor, these
are plain text files.

To get the list of all possible parameters for each module:

```bash
haddock3-cfg -h
haddock3-cfg -m <MODULE NAME>
haddock3-cfg -m rigidbody
```

We are actively working towards expanding our documentation pages.
Thanks for using HADDOCK3! For any question please [open an issue
here](https://github.com/haddocking/haddock3/issues).
