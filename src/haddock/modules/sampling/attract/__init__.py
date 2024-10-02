"""attract docking module
==================================

This module performs a fragment-based single-stranded (ss) RNA-protein docking 
using ATTRACT docking engine. This docking approach was developed to tackle the 
flexibility of ssRNA. Its core idea is to split ssRNA chain into overlapping 
trinucleotides (fragments), and dock them onto the rigid receptor separately, 
assembling the fragments back into the whole chain models afterwards.

#todo 
add short description of the protocol, including CG, NAlib, sampling, 
scoring (two ways) and assembly + possible restraints. 
."""
import os
#import sys
import tempfile
import subprocess
import shutil
import shlex
from pathlib import Path
from haddock.modules.sampling.attract.attractmodule import (
    rename_and_coarse_grain,
    process_rna_file,
)

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath, Any 
#from haddock.libs import libpdb
#from haddock.libs.libio import working_directory
from haddock.libs.libontology import Format, PDBFile
#from haddock.libs.libutil import check_subprocess
from haddock.modules import BaseHaddockModule

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)

class HaddockModule(BaseHaddockModule):
    """ATTRACT module."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 initial_params: FilePath = DEFAULT_CONFIG) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm that ATTRACT and its environment variables are properly set."""
        try:
            attract_dir = os.environ['ATTRACTDIR']
            attract_tools = os.environ['ATTRACTTOOLS']
            nalib = os.environ['LIBRARY']
            randsearch = os.environ['RANDSEARCH']
        except KeyError as e:
            raise EnvironmentError(f"Required environment variable not found: {e}")

        attract_exec = Path(attract_dir, 'attract')
        result = subprocess.run([str(attract_exec)], capture_output=True, text=True)
        
        if "Too few arguments" not in result.stderr:
            raise RuntimeError('ATTRACT is not installed properly')

        # pass paths to _run for further use
        cls.attract_dir = attract_dir
        cls.attract_tools = attract_tools
        cls.nalib = nalib

    def _run(self) -> None:
        """Execute module.
        Currently:
        1. Converts protein and RNA in ATTRACT coarse-grain 
        2. Splits RNA into overlapping fragments
        3. Creates required by ATTRACT files: motif.list, boundfrag.list, nalib
        4. Passes initial all-atom protein and RNA to next module
        """

        # Get the models generated in previous step
        models: list[PDBFile] = [
            p for p in self.previous_io.output
            if p.file_type == Format.PDB ]

        # Check that exactly two models are provided
        # practically attract needs protein structure and RNA *sequence* 
        # but at this stage it's more practical to ask for RNA structure     
        if len(models) !=2 : 
            _msg = "ATTRACT requires exactly two molecules"
            self.finish_with_error(_msg)
        
        # Copy each model to the working directory
        for model in models:
            src_model = Path(model.path, model.file_name)
            dest_model = Path(os.getcwd(), model.file_name)
            shutil.copyfile(src_model, dest_model)
       
        # Ensure we have exactly protein and RNA molecules
        model_1 = models[0]
        model_2 = models[1]
    
        attracttools = self.attract_tools
        attrac_reduce_path = Path(attracttools, 'reduce.py')

        _, label_1 = rename_and_coarse_grain(model_1.file_name, attrac_reduce_path)
        _, label_2 = rename_and_coarse_grain(model_2.file_name, attrac_reduce_path)
        
        if {label_1, label_2} != {'protein', 'rna'}:
            _msg = "ATTRACT requires protein and RNA molecules as input"
            self.finish_with_error(_msg)

        # Add required by ATTRACT files: 
        log.info("Preparing docking directory")
        
        nalib = self.nalib
        cmd = f"ln -s {nalib} nalib"
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        err = p.stderr.decode('utf-8')

        process_rna_file('rna-aar.pdb')

        tmp_dir = tempfile.mkdtemp(dir=os.getcwd())
        # will be used during the docking
        shutil.rmtree(tmp_dir)

        list_of_created_models = []     
        created_models = ['protein-aa.pdb','rna-aa.pdb']
        for model in created_models:
            pdb_object = PDBFile(Path(model).name,  path=".")
            list_of_created_models.append(pdb_object)
       
        self.output_models = list_of_created_models
        self.export_io_models()

            