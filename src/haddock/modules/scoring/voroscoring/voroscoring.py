"""Voro scoring class.

This class holds all the machinery to perform scoring of input pdb models using
ftdmp voro-mqa-all tool.
For more information, please check: https://github.com/kliment-olechnovic/ftdmp

It is a third party module, and requires the appropriate set up and intallation
for it to run without issue.
"""

import os
import subprocess
import glob
import time

from random import randint

from haddock import log
from haddock.core.typing import Any, Generator, Path, Union
from haddock.libs.libio import working_directory
from haddock.libs.libontology import NaN, PDBFile


# Defines the SLURM job template
# Notes: Please feel free to modify the #SBATCH entries to fit your needs/setup
SLURM_HEADER_GPU = """#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
"""

SLURM_HEADER_CPU = """#SBATCH -J hd3-voroscoring-cpu
#SBATCH --partition haddock
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
"""

# Job template
VOROMQA_CFG_TEMPLATE = """#!/bin/bash
{HEADER}
#SBATCH -J {JOBNAME}

# Where to do the work
WORKDIR="{WORKDIR}"

# Name of the outputfile (.ssv for space separated values)
OUTPUT_FNAME="voro_scores.ssv"

# Define Constants
CONDA_INSTALL_DIR="{CONDA_INSTALL_DIR}"
CONDA_ENV_NAME="{CONDA_ENV_NAME}"
FTDMP_INSTALL_DIR="{FTDMP_INSTALL_DIR}"
VOROMQA_SCRIPT="ftdmp-qa-all"

# Define workflow variables
OUTPUT_FPATH="$WORKDIR/$OUTPUT_FNAME"
PDB_LIST_PATH="{PDB_LIST_PATH}"
OUT_MSG="Output file is here: $OUTPUT_FPATH"

# 1. Setup enviroments
# Load the gnu13 module...
# NOTE: specific to haddock-team users...
# module load gnu13
# Activate conda env
source "$CONDA_INSTALL_DIR/bin/activate"
conda activate $CONDA_ENV_NAME
echo "conda env: $CONDA_PREFIX"

# 2. Setup run directory
# Create working directory
mkdir -p $WORKDIR

# 3. Run voro-mqa (model quality assessment)
# Go to ftdmp install directory
cd $FTDMP_INSTALL_DIR
echo "Directory: $PWD"
# run voro-mqa
echo "./$VOROMQA_SCRIPT --conda-path $CONDA_INSTALL_DIR --conda-env $CONDA_ENV_NAME --workdir '$WORKDIR' --rank-names 'protein_protein_voromqa_and_global_and_gnn_no_sr' < $PDB_LIST_PATH > $OUTPUT_FPATH"
./$VOROMQA_SCRIPT --conda-path $CONDA_INSTALL_DIR --conda-env $CONDA_ENV_NAME --workdir $WORKDIR --rank-names 'protein_protein_voromqa_and_global_and_gnn_no_sr' --output-redundancy-threshold 1.0 < $PDB_LIST_PATH > $OUTPUT_FPATH
# Let the magic happen..

# 4. Analyze results
# Print final ouput file
echo $OUT_MSG
"""  # noqa : E501


class VoroMQA():
    """The Haddock3 implementation of voro-mqa-all as a python class."""

    def __init__(
            self,
            models: list[PDBFile],
            workdir: Union[str, Path],
            params: dict[str, Any],
            output: Union[str, Path] = "voroscoring_voro.tsv",
            ):
        """Init of the VoroMQA class.

        Parameters
        ----------
        models : list[PDBFile]
            List of input PDB files to be scored.
        workdir : Union[str, Path]
            Where to do the process.
        params : dict[str, Any]
            Config file parameters
        output : Path, optional
            Name of the generated file, by default Path("voroscoring_voro.tsv")
        """
        self.models = models
        self.workdir = workdir
        self.params = params
        self.output = Path(output)

    def run(self):
        """Process class logic."""
        # Obtain absolute paths
        self.workdir = Path(self.workdir).resolve()
        all_pdbs = [
            str(Path(mdl.path, mdl.file_name).resolve())
            for mdl in self.models
            ]
        # Loop over batches
        for bi, batch in enumerate(self.batched(all_pdbs, size=300)):
            # Run slurm
            self.run_voro_batch(batch, batch_index=bi + 1)
        # Recombine all batches output files
        scores_fpath = self.recombine_batches()
        log.info(f"Generated output file: {scores_fpath}")

    def run_voro_batch(
            self,
            pdb_filepaths: list[str],
            batch_index: int = 1,
            ) -> None:
        """Preset and launch predictions on subset of pdb files.

        Parameters
        ----------
        pdb_filepaths : list[str]
            List of absolute path to the PDBs to score
        batch_index : int, optional
            Index of the batch, by default 1
        """
        # Create workdir
        batch_workdir = Path(self.workdir, f"batch_{batch_index}")
        batch_workdir.mkdir(parents=True)

        # Create list of pdb files
        pdb_files_list_path = Path(batch_workdir, "pdbs.list")
        pdb_files_list_path.write_text(os.linesep.join(pdb_filepaths))

        # Format config file
        batch_cfg = VOROMQA_CFG_TEMPLATE.format(
            HEADER=SLURM_HEADER_CPU,
            CONDA_INSTALL_DIR=self.params["conda_install_dir"],
            CONDA_ENV_NAME=self.params["conda_env_name"],
            FTDMP_INSTALL_DIR=self.params["ftdmp_install_dir"],
            JOBNAME=f"hd3_voro_b{batch_index}",
            WORKDIR=batch_workdir,
            PDB_LIST_PATH=pdb_files_list_path,
            )
        
        # Write it
        batch_cfg_fpath = Path(batch_workdir, "vorobatchcfg.job")
        batch_cfg_fpath.write_text(batch_cfg)

        # Launch script
        self._launch_computation(batch_workdir, batch_cfg_fpath)
        #initdir = os.getcwd()
        #os.chdir(batch_workdir)
        #log.info(f"sbatch {batch_cfg_fpath}")
        #subprocess.run(f"sbatch {batch_cfg_fpath}", shell=True)
        #os.chdir(initdir)

    def _launch_computation(self, batch_workdir: str, batch_cfg_fpath: str) -> None:
        """Execute a given script from working directory.

        Parameters
        ----------
        batch_workdir : str
            Path to working directory
        batch_cfg_fpath : str
            Script to execute
        """
        exec_tool = "sbatch" if self.params["mode"] == "batch" else "bash"
        cmd_ = f"{exec_tool} {batch_cfg_fpath}"
        with working_directory(batch_workdir):
            log.info(cmd_)
            subprocess.run(cmd_, shell=True)

    def recombine_batches(self) -> str:
        """Recombine batches output file in a single one.

        Returns
        -------
        finale_output_fpath : str
            Filepath of the recombined scores
        """
        # Wait for all results to be obtained
        batches_result_paths = self.wait_for_termination()
        # Loop over them
        all_predictions: list[dict[str, str]] = []
        combined_header: list[str] = []
        for batch_results in batches_result_paths:
            # Read voro results
            with open(batch_results, 'r') as filin:
                header = filin.readline().strip().split(' ')
                for head in header:
                    if head not in combined_header:
                        combined_header.append(head)
                for line in filin:
                    s_ = line.strip().split(' ')
                    all_predictions.append({
                        head: s_[header.index(head)]
                        for head in header
                        })

        # Sort all batches entries
        sorted_entries = sorted(
            all_predictions,
            key=lambda k: float(k[self.params["metric"]]),
            reverse="_energy" not in self.params["metric"],
            )

        # Write final output file
        finale_output_fpath = f"{self.workdir}/{self.output}"
        with open(finale_output_fpath, "w") as filout:
            file_header = '\t'.join(combined_header)
            filout.write(file_header + os.linesep)
            for entry in sorted_entries:
                ordered_data = [
                    entry[h] if h in entry.keys() else '-'
                    for h in combined_header
                    ]
                line = '\t'.join(ordered_data)
                filout.write(line + os.linesep)
        return finale_output_fpath
    
    def wait_for_termination(self, wait_time: float = 60) -> list[Path]:
        """Wait until all results are accessible.

        Parameters
        ----------
        wait_time : int, optional
            Time in second between every termination checks, by default 60

        Returns
        -------
        output_files : list[Path]
            List of voro scores results for every batches.
        """
        batches_dirpath = glob.glob(f"{self.workdir}/batch_*/")
        log.info(
            f"Waiting for {len(batches_dirpath)} "
            "voro-mqa prediction batch(es) to finish..."
            )
        while True:
            try:
                output_files: list[Path] = []
                for batch_dir in batches_dirpath:
                    expected_outputfile = Path(batch_dir, "voro_scores.ssv")
                    assert expected_outputfile.exists()
                    assert expected_outputfile.stat().st_size != 0
                    output_files.append(expected_outputfile)
            except AssertionError:
                log.info(f"Waiting {wait_time} sec more...")
                time.sleep(wait_time)
            else:
                log.info(
                    "VoroMQA results are accessible: "
                    f"{len(output_files)} batch(es)"
                    )
                return output_files

    @staticmethod
    def batched(
            entries: list[str],
            size: int = 300,
            ) -> Generator[list[str], None, None]:
        """Generate batches of defined size.

        Parameters
        ----------
        entries : list[str]
            List of pdb files.
        size : int, optional
            Maximum size in every batch, by default 300

        Yields
        ------
        batch : Generator[list[str], None, None]
            List of pdb files <= size.
        """
        batch = []
        for pdb in entries:
            batch.append(pdb)
            if len(batch) == size:
                yield batch
                batch = []
        if batch:
            yield batch


def update_models_with_scores(
        voro_scoring_fname: Union[str, Path],
        models: list[PDBFile],
        metric: str = "jury_score",
        ) -> list[PDBFile]:
    """Update PDBfiles with computed scores.

    Parameters
    ----------
    output_fname : Union[str, Path]
        Path to the file where to access scoring data.
    models : list[PDBFile]
        List of PDBFiles to be updated.
    metric : str, optional
        Name of the metric to be retrieved, by default "jury_score"

    Returns
    -------
    models : list[PDBFile]
        The updated list of PDBfiles now holding the score and rank attributes.
    """
    scores_mapper: dict[str, float] = {}
    ranking_mapper: dict[str, int] = {}
    rank: int = 0
    # Read output file
    with open(voro_scoring_fname, 'r') as filin:
        for i, line in enumerate(filin):
            s_ = line.strip().split('\t')
            # Extract header
            if i == 0:
                header = s_
                continue
            # Extract data
            modelpath = str(s_[header.index("ID")])
            score = float(s_[header.index(metric)])
            # Only extract model filename
            model_filename = modelpath.split('/')[-1]
            # Reverse score if not an energy
            if "_energy" not in metric:
                score = -score
            # Hold score
            scores_mapper[model_filename] = score
            rank += 1
            ranking_mapper[model_filename] = rank

    # Compute rankings
    #ranking_mapper = {
    #    model_filename: rank
    #    for rank, model_filename in enumerate(sorted(scores_mapper), start=1)
    #    }

    # Loop over input models
    for model in models:
        # Add score and rank as attribute
        if model.file_name in scores_mapper.keys():
            model.score = scores_mapper[model.file_name]
            model.rank = ranking_mapper[model.file_name]
        # In some cases computation may fail
        else:
            # Go for (garlic cheese) naans
            model.score = NaN
            model.rank = NaN
        model.ori_name = model.file_name
    return models
