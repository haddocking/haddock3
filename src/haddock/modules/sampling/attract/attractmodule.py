"""Set if functions related to [attract]"""

import os
import sys 
import subprocess
import shutil
import shlex
from pathlib import Path


def rename_and_coarse_grain(pdb_file: str, reduce_path: str, lines_to_check: int = 10):
    """
    Check several last lines of pdb file to label it as RNA, DNA or protein; 
    Then reduce it to ATTRACT coarse grained representation.

    Args:
    pdb_file (str): Path to the PDB file.
    lines_to_check (int): Number of lines to check from the end of the file.
    reduce_path (str); Path to $ATTRACTTOOLD/reduce.py 

    Returns:
    subprocess.CompletedProcess: Result of the subprocess run. 
    (ATTRACT CG atom mapping)
    molecule_type (str): Label of the input file
    """
    dna_residues = {'DA', 'DG', 'DC', 'DT'}
    rna_residues = {'A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC'}
    
    old_path = Path(pdb_file)
    
    with old_path.open('r') as f:
        lines = f.readlines()[-lines_to_check:]
    
    for line in lines:
        if line.startswith("ATOM"):
            residue_name = line[17:20].strip()
            if residue_name in dna_residues:
                molecule_type = 'dna'
            elif residue_name in rna_residues:
                molecule_type = 'rna'
            else:
                molecule_type = 'protein'
    
        new_path = old_path.with_name(f'{molecule_type}-aa.pdb')
        old_path.replace(new_path)
        
        cmd = [sys.executable, reduce_path]
        if molecule_type in ('dna', 'rna'):
            cmd.append(f'--{molecule_type}')
        cmd.append(str(new_path))
        # this works in haddock3 env 
        return subprocess.run(cmd, capture_output=True, text=True), molecule_type

def process_rna_file(pdb_file: str = 'rna-aar.pdb'):
    """
    Process PDB file to generate RNA fragment and motif lists, and extract individual fragments.
    Save generated files. 
    
    Args:
    pdb_file (str): Path to reduced RNA PDB file.    

    Note: 
    This script works on all-atom pdb as well
    """
    enumer_nucl = []
    sequence = ''
    fragments = []
    residue_buffer = []  
    current_residue = [] 
    last_residue_id = None

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                residue_name = line[17:20].strip()[-1]    
                residue_id = int(line[22:26].strip())

                # Get sequence
                if [residue_name, residue_id] not in enumer_nucl:
                    enumer_nucl.append([residue_name, residue_id])
                    sequence += residue_name
              
                # Get fragments
                if last_residue_id is None or residue_id == last_residue_id:
                    current_residue.append(line)
                else:
                    if current_residue:
                        residue_buffer.append(current_residue)
                        if len(residue_buffer) == 3:
                            fragments.append(residue_buffer[:])
                            residue_buffer.pop(0) 
                    current_residue = [line]
                last_residue_id = residue_id 

        # Add the last residue if any
        if current_residue:
            residue_buffer.append(current_residue)
            if len(residue_buffer) == 3:
                fragments.append(residue_buffer[:])
                
    # Generate motif and boundfrag lists
    prev_frag = sequence[:3]
    count = 1
    boundfrag_list = [[count, prev_frag]]
    motif_list = [prev_frag]

    for x in sequence[3:]:
        count += 1
        next_frag = prev_frag[1:] + x
        prev_frag = next_frag
        boundfrag_list.append([count, next_frag])
        if prev_frag not in motif_list: 
            motif_list.append(next_frag)

    # Save boundfrag_list
    with open('boundfrag.list', 'w') as f:
        for item in boundfrag_list:
            f.write(f"{item[0]} {item[1]}\n")
    # Save motif_list
    with open('motif.list', 'w') as f:
        for item in motif_list:
            f.write(f"{item}\n")
    # Save each fragments
    for i, fragment in enumerate(fragments, start=1):
        with open(f"frag{i}r.pdb", 'w') as frag_file:
            for nucleotide in fragment:
                for atomline in nucleotide:
                    if not atomline.endswith('\n'):
                        atomline += '\n'
                    frag_file.write(atomline)
    return