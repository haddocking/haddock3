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
import subprocess
import shutil
import shlex
from pathlib import Path

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
        """Confirm that ATTRACT is installed."""
        attract_dir = os.environ['ATTRACTDIR']
        attract_exec = Path(attract_dir, 'attract')
        cmd = f'{attract_exec}'
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        err = p.stderr.decode('utf-8')
        if "Too few arguments" in err and "usage:" in err:
            # ATTRACT is installed correctly 
            pass
        else:
            raise Exception('ATTRACT is not installed properly')

    def determine_molecule_type(self, pdb_file, lines_to_check=10):
        """Check sevelal last lines of pdb file to label it 
        as rna, dna or protein."""
        dna_residues = {'DA', 'DG', 'DC', 'DT'}  
        rna_residues = {'A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC'}  
        with open(pdb_file, 'r') as f:
            lines = f.readlines()[-lines_to_check:]
        for line in reversed(lines):
            if line.startswith("ATOM"):
                residue_name = line[17:20].strip()
                if residue_name in dna_residues:
                    return 'dna'
                elif residue_name in rna_residues:
                    return 'rna'
        return 'protein'

    def rename_and_coarse_grain(self, file_name, new_name, reduce_path, is_rna=False):
        """Rename models to 'protein-aa.pdb' and 'rna-aa.pdb' and convert 
        them to ATTRACT coarse-grained representation."""
        os.rename(file_name, new_name)
        cmd = ['python', reduce_path, '--rna', new_name] if is_rna else ['python', reduce_path, new_name]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def get_frag_lists(self, file):
        """Get list of fragments and list of motifs from rna-aar.pdb."""
        tmp = []
        seq = ''
        # extract sequence from structure
        for line in file:
            if line.startswith("ATOM"):
                residue_name = line[17:20].strip()[-1]
                residue_id = int(line[22:26].strip())
                if [residue_name, residue_id] not in tmp:
                    tmp.append([residue_name, residue_id])
                    seq = seq + residue_name
        # get info for motif.list and boundfrag.list
        prev_frag = seq[:3]
        count = 1 
        boundfrag_list = [[count, prev_frag]]
        motif_list = [prev_frag]
        for x in seq[3:]:
            count += 1 
            next_frag = prev_frag[1:] + x
            prev_frag = next_frag
            boundfrag_list.append([count, next_frag])
            if prev_frag not in motif_list: motif_list.append(next_frag)
        return motif_list, boundfrag_list

    def get_frag_pdbs(self, file):
        """Extract fragments from full rna"""
        first_res_id = int(file[0][22:26].strip())
        last_res_id = int(file[-1][22:26].strip())
        frag_total = last_res_id - first_res_id - 1 
        fragmens=[]
        for curr_residue_id in range(first_res_id, first_res_id+frag_total):
            curr_frag = []
            for line in file:
                if line.startswith("ATOM"):
                    residue_id = int(line[22:26].strip())
                    if residue_id in [curr_residue_id, curr_residue_id+1, curr_residue_id+2]:
                        curr_frag.append(line)
                    elif residue_id > curr_residue_id+2:
                        fragmens.append([curr_frag])
                        break
        fragmens.append([curr_frag])
        return fragmens

    def _run(self) -> None:
        """Execute module.
            
            # todo 
            1. check library + 
            2. process input +
            3. make folder +
            4. run docking 
            5. run lrmsd (optional)
            6. HIPPO (optional)
            7. Assembly ?
            8. Scoring ?
            9. SeleTopChains ? """

        # Check if $LIBRARY exists 
        try:
            nalib = os.environ['LIBRARY']
        except:
            _msg = "Environment variable $LIBRARY not found."
            self.finish_with_error(_msg)
        log.info("NAlib found")

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
        label_1 = self.determine_molecule_type(model_1.file_name) 
        label_2 = self.determine_molecule_type(model_2.file_name) 
        labels = {label_1, label_2}
        
        if labels != {'protein', 'rna'}:
            _msg = "ATTRACT requires protein and RNA molecules as input"
            self.finish_with_error(_msg)

        # Convert each molecule to corse-grained representation
        # using $ATTRACTTOOLD/reduce.py 
        log.info("Converting to coarse-grain representation")
        try:
            attracttools = os.environ['ATTRACTTOOLS']
            attrac_reduce_path = Path(attracttools, 'reduce.py')

            if label_1 == 'protein':
                self.rename_and_coarse_grain(model_1.file_name, 'protien-aa.pdb', attrac_reduce_path, is_rna=False)
                self.rename_and_coarse_grain(model_2.file_name, 'rna-aa.pdb', attrac_reduce_path, is_rna=True)
            else:
                self.rename_and_coarse_grain(model_1.file_name, 'rna-aa.pdb', attrac_reduce_path, is_rna=True)
                self.rename_and_coarse_grain(model_2.file_name, 'protien-aa.pdb', attrac_reduce_path, is_rna=False)
        except:
            _msg = "Convertation to coarse-grain representation failes"
            self.finish_with_error(_msg)
        
        # Add required by ATTRACT files:
        log.info("Preparing docking directory")
        # 1.link fragment library          
        cmd = f"ln -s {nalib} nalib"
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        err = p.stderr.decode('utf-8')
        # 2. get motif.list and boundfrag.list
        f = open('rna-aar.pdb','r')
        file = f.readlines()
        f.close()
        motif, boundfrag = self.get_frag_lists(file)
        with open('boundfrag.list', 'w') as f:
            for item in boundfrag:
                f.write(f"{item[0]} {item[1]}\n")
        with open('motif.list', 'w') as f:
            for item in motif:
                f.write(f"{item}\n")
        # 3. create fragXr.pdb (optional)
        fragments = self.get_frag_pdbs(file)
        for i in range(1, len(fragments)+1):
            with open(f"frag{i}r.list", 'w') as f:
                for fragment in fragments[i-1]:
                    for line in fragment:
                        f.write(line)
        
        # Run docking 
        #log.info("Running ATTRACT")   
        
        
        # ???
        list_of_created_models = []     
        created_models = ['protien-aa.pdb','rna-aa.pdb']
        for model in created_models:
            pdb_object = PDBFile(Path(model).name,  path=".")
            list_of_created_models.append(pdb_object)

            self.output_models = list_of_created_models
            self.export_io_models()