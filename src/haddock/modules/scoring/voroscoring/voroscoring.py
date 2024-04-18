import os
import subprocess
import glob
import time

from random import randint

from haddock import log
from haddock.core.typing import Any, Generator, Path, Union
from haddock.libs.libontology import PDBFile


VOROMQA_CFG_TEMPLATE = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:{GPUID}
#SBATCH --mem-per-gpu=1GB

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
# Load the gnu12 module...
# NOTE: specific to tintin users...
module load gnu12
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
"""


class VoroMQA():

    def __init__(
            self,
            models: list[Union[str, Path, PDBFile]],
            workdir: Union[str, Path],
            params: dict[str, Any],
            output: Path = Path("voroscoring.tsv"),
            ):
        self.models = models
        self.workdir = workdir
        self.params = params
        self.output = output

    def run(self):
        # Obtain absolute paths
        self.workdir = Path(self.workdir).resolve()
        all_pdbs = [
            Path(mdl.path, mdl.file_name).resolve()
            for mdl in self.models
            ]
        # Loop over batches
        for bi, batch in enumerate(self.batched(all_pdbs, size=300)):
            # Run slurm
            self.run_voro_batch(
                batch,
                batch_index=bi,
                gpuid=bi % self.params['nb_gpus'],
                )
        # Recombine all batches output files
        scores_fpath = self.recombine_batches(self.workdir)
        log.info(f"Generated output file: {scores_fpath}")

    def run_voro_batch(
            self,
            pdb_filepaths: list[Union[str, Path]],
            batch_index: int = 1,
            gpuid: int = -1,
            ) -> None:
        # Create workdir
        batch_workdir = Path(self.workdir, f"batch_{batch_index}")
        batch_workdir.mkdir(parents=True)

        # Create list of pdb files
        pdb_files_list_path = Path(batch_workdir, "pdbs.list")
        pdb_files_list_path.write_text(os.linesep.join(pdb_filepaths))

        # Get GPU id
        if gpuid < 0:
            gpuid = randint(0, self.params["nb_gpus"] - 1)

        # Format config file
        batch_cfg = VOROMQA_CFG_TEMPLATE.format(
            CONDA_INSTALL_DIR=self.params["conda_install_dir"],
            CONDA_ENV_NAME=self.params["conda_env_name"],
            FTDMP_INSTALL_DIR=self.params["ftdmp_install_dir"],
            GPUID=gpuid,
            WORKDIR=batch_workdir,
            PDB_LIST_PATH=pdb_files_list_path,
            )
        
        # Write it
        batch_cfg_fpath = Path(batch_workdir, "vorobatchcfg.job")
        batch_cfg_fpath.write_text(batch_cfg)

        # Launch slurm
        initdir = os.getcwd()
        os.chdir(batch_workdir)
        log.info(f"sbatch {batch_cfg_fpath}")
        subprocess.run(f"sbatch {batch_cfg_fpath}", shell=True)
        os.chdir(initdir)
    
    def recombine_batches(self) -> str:
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
                filout.write(line+os.linesep)
        return finale_output_fpath
    
    def wait_for_termination(
            self,
            wait_time: int = 60,
            ) -> list[Union[str, Path]]:
        batches_dirpath = glob.glob(f"{self.workdir}/batch_*/")
        while True:
            try:
                output_files: list[Union[str, Path]] = []
                for batch_dir in batches_dirpath:
                    expected_outputfile = Path(batch_dir, "voro_scores.ssv")
                    assert expected_outputfile.exists()
                    assert expected_outputfile.stat().st_size != 0
                    output_files.append(expected_outputfile)
            except AssertionError as _e:
                log.info(f"Waiting {wait_time} sec...")
                time.sleep(wait_time)
            else:
                return output_files

    @staticmethod
    def batched(entries: str, size: int = 300) -> Generator[list, None, None]:
        batch = []
        for pdb in entries:
            batch.append(pdb)
            if len(batch) == size:
                yield batch
                batch = []
        yield batch

def update_models_with_scores(
        output_fname: Union[str, Path],
        models: list[PDBFile],
        metric: str = "jury_score",
        ) -> tuple[list[PDBFile], dict[str, dict[str, float]]]:
    scores_mapper: dict[str, float] = {}
    # Read output file
    with open(output_fname, 'r') as filin:
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
            # Hold score
            scores_mapper[model_filename] = score

    # Compute rankings
    ranking_mapper = {
        model_filename: rank
        for rank, model_filename in enumerate(
            sorted(
                scores_mapper,
                reverse="_energy" not in metric,
                ),
            start=1,
            )
        }

    data_mapper = {
        model_filename: {
            "score": scores_mapper[model_filename],
            "rank": ranking_mapper[model_filename],
            }
        }

    # Loop over input models
    for model in models:
        # only modify the model score
        model.score = data_mapper[model.file_name]["score"]
        model.rank = data_mapper[model.file_name]["rank"]

    return models, data_mapper

def write_models_scores(
        models_scores: dict[str, dict[str, Union[float, int]]],
        filepath: Union[str, Path],
        ) -> None:
    header = ("structure", "original_name", "md5", "score", "rank", )
    with open(filepath, 'w') as filout:
        filout.write('\t'.join(header) + os.linesep)
        # sort models by keys
        sorted_models = sorted(
            models_scores,
            key=lambda k: models_scores[k]['rank'],
            )
        for modelname in sorted_models:
            scores = models_scores[modelname]
            newline_dt = f"{modelname}\t{modelname}\t-\t{scores['score']}\t{scores['rank']}"  # noqa : E501
            filout.write(newline_dt + os.linesep)