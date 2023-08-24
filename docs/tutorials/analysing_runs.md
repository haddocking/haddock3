# Analysing runs

HADDOCK3 allows to analyse different steps of the workflow, even after it has been completed.

## haddock3-analyse

The `haddock3-analyse` command is the main tool for the analysis of one or more workflow steps.

```
haddock3-analyse -r my-run-folder -m 2 5 6
```

Here `my-run-folder` is the run directory and 2, 5, and 6 are the steps that you want to analyse.

The command will inspect the folder, looking for the existing models. If the selected module is 
a `caprieval` module, `haddock3-analyse` simply loads the `capri_ss.tsv` and `capri_clt.tsv` files
produced by the `caprieval` module. Otherwise, `haddock3-analyse` runs a `caprieval` analysis of the models.
You can provide some [caprieval-specific parameters](https://github.com/haddocking/haddock3/blob/main/src/haddock/modules/analysis/caprieval/defaults.yaml)
using the following syntax:

```
haddock3-analyse -r my-run-folder -m 2 5 6 -p reference_fname my_ref.pdb receptor_chain F
```

Here the `-p` key tells the code that you are about to insert caprieval parameters, whose name should match the parameter name.

Another parameter that can be specified is `top_cluster`, which defines how many of the first N clusters will be considered in the analysis.
This value is set to 10 by default.

```
haddock3-analyse -r my-run-folder -m 2 5 6 --top_cluster 12
```

This number is meaningless when dealing with models with no cluster information, that is, models that have never been clustered before.

By default `haddock3-analyse` produces [plotly](https://plotly.com/python/) plots in the html `format`, but the user can select 
one of the formats available [here](https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_image.html), 
while also adjusting the resolution with the `scale` parameter:

```
haddock3-analyse -r my-run-folder -m 2 5 6 --format pdf --scale 2.0
```

### The analysis folder

After running `haddock3-analyse` you can check the content of the `analysis` directory in your run folder.
If everything went succesfully, one of the above commands could have produced an analysis folder structured as

```
my-run-folder/
|--- analysis/
     |--- 2_caprieval_analysis
     |--- 5_seletopclusts_analysis
     |--- 6_flexref_analysis
```

Each subfolder contains all the analysis plots related to that specific step of the workflow.

By default `haddock3-analyse` produces a set of scatter plots that compare each HADDOCK energy term 
(i.e., the HADDOCK score and its components) to the different metrics used to evaluate the quality of a model,
such as the interface-RMSD, Fnat, DOCKQ, and so on. An example is available [here](../figs/irmsd_score.html).

For each of the energy component and the metrics mentioned above `haddock3-analyse` produces also a box plot, in which each cluster 
is considered separately. An example is available [here]().

### The report

Scatter plots, box plots, CAPRI statistics and an interactive visualization of the models is available in the `report.html` file, present
in each analysis subfolder. In order to visualize the models it is necessary to start a local server at the end of the `haddock3-analyse` run,
following the indications provided in the log file:

```
[2023-08-24 10:09:09,552 cli_analyse INFO] View the results in analysis/12_caprieval_analysis/report.html
[2023-08-24 10:09:09,552 cli_analyse INFO] To view structures or download the structure files, in a terminal run the command
`python -m http.server --directory /haddock3/examples/docking-antibody-antigen/run1-CDR-acc-cltsel-test`.
By default, http server runs on `http://0.0.0.0:8000/`. Open the link
http://0.0.0.0:8000/analysis/12_caprieval_analysis/report.html in a web browser.
```

Launch this command to open the report:
```
python -m http.server --directory path-to-my-run
```

In the browser you can navigate to each analysis subfolder and open the `report.html` file. If you are not interested in
visualizing the models, you can simply open the `report.html` file in a standard browser.

## haddock3-traceback

HADDOCK3 is highly customisable and modular, as the user can introduce several refinement, clustering, and scoring steps in a workflow.
Quantifying the impact of the different modules is important while developing a novel docking protocol. The `haddock3-traceback` command
is developed to assist the user in this task, as it allows to "connect" all the models generated in a HADDOCK3 workflow:

```
haddock3-traceback my-run-folder
```

`haddock3-traceback` creates a traceback subfolder within the `my-run-folder` directory, containing a `traceback.tsv` table:

```
00_topo1     00_topo2        01_rigidbody            01_rigidbody_rank       04_seletopclusts        04_seletopclusts_rank   06_flexref      06_flexref_rank 
4G6K.psf     4I1B.psf        rigidbody_10.pdb        3                       cluster_1_model_1.pdb   1                       flexref_1.pdb   2       
4G6K.psf     4I1B.psf        rigidbody_11.pdb        10                      cluster_1_model_2.pdb   3                       flexref_3.pdb   1       
4G6K.psf     4I1B.psf        rigidbody_18.pdb        4                       cluster_2_model_1.pdb   2                       flexref_2.pdb   4      
4G6K.psf     4I1B.psf        rigidbody_20.pdb        15                      cluster_2_model_2.pdb   4                       flexref_4.pdb   3       
```

In this table each row represents a model that has been produced by the workflow. The (typically) two used topologies are reported first,
and then each module has its own column, containing the name and rank of the model at that stage. As an example, in the first row of the
table above `rigidbody_10.pdb` is ranked 3rd at the `rigidbody` stage. Then, it becomes `cluster_1_model_1.pdb` (ranked 1st) after 
the `seletopclusts` module. This model is then refined in `flexref_1.pdb`, which turns out to be the 2nd best model at the end of the workflow.

The table can be easily parsed and used to evaluate the impact of different refinement steps on the different models.

## The postprocess option

You may want to run the `haddock3-analyse` and `haddock3-traceback` commands by default at the end of the workflow.
The `postprocess` option of a standard HADDOCK3 configuration (.cfg) file is devoted to this task. At first, it forces HADDOCK3 
to execute `haddock3-analyse` on all the `caprieval` folders found in the workflow, therefore loading data present in the CAPRI tables.
Second, it executes the `haddock3-traceback` command.

To activate this, just set `postprocess` to `true` at the beginning of your configuration file:

```
 ====================================================================
# This is a HADDOCK3 configuration file

# directory in which the docking will be done
run_dir = "my-run-folder"

# postprocess the run
postprocess = true

...
```

You can find additional help by running the command: `haddock3-analyse -h` and `haddock3-traceback -h` and reading
the parameters' explanations. Otherwise, ask us in the ["issues" forum](https://github.com/haddocking/haddock3/issues).
