from pathlib import Path

import pytest
import string
import numpy as np

from haddock.modules.analysis.caprieval.capri import load_contacts
from haddock.libs.libontology import PDBFile


def calc_fnat_with_caprieval(model: Path, native: Path) -> float:
    model_pdb = PDBFile(model)
    native_pdb = PDBFile(native)

    model_contacts = load_contacts(model_pdb)
    native_contacts = load_contacts(native_pdb)

    intersection = native_contacts & model_contacts

    fnat = len(intersection) / float(len(model_contacts))

    return fnat


@pytest.fixture
def calc_fnat():
    return calc_fnat_with_caprieval


def generate_synthetic_pdb(
    n_atoms: int,
    output_path: Path,
    n_chains: int = 2,
    seed: int = 42,
) -> None:
    """Generate a synthetic PDB file for testing.

    Creates a multi-chain protein structure with CA atoms positioned randomly
    in space. Chains are spatially separated to enable interchain contacts.

    Parameters
    ----------
    n_atoms : int
        Total number of CA atoms to generate across all chains
    output_path : Path
        Path where PDB file will be written
    n_chains : int, optional
        Number of chains to generate (default: 2)
    seed : int, optional
        Random seed for reproducibility (default: 42)
    """
    np.random.seed(seed)

    # Standard amino acids
    amino_acids = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
        "LEU",
        "MET",
        "ASN",
        "PRO",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "VAL",
        "TRP",
        "TYR",
    ]

    # Chain identifiers
    chain_ids = string.ascii_uppercase[:n_chains]

    # Distribute atoms across chains
    atoms_per_chain = [n_atoms // n_chains] * n_chains
    for i in range(n_atoms % n_chains):
        atoms_per_chain[i] += 1

    with open(output_path, "w") as f:
        # Write header
        f.write("REMARK   Synthetic PDB generated for integration testing\n")
        f.write(f"REMARK   Total atoms: {n_atoms}, Chains: {n_chains}\n")

        atom_num = 1

        for chain_idx, chain_id in enumerate(chain_ids):
            n_residues = atoms_per_chain[chain_idx]

            if n_residues == 0:
                continue

            # Offset each chain spatially (25 Å apart)
            chain_offset_x = chain_idx * 25.0

            last_resname = "ALA"
            last_res_num = 1

            for res_num in range(1, n_residues + 1):
                # Random coordinates within a reasonable range
                x = np.random.uniform(0, 50) + chain_offset_x
                y = np.random.uniform(0, 50)
                z = np.random.uniform(0, 50)

                # Random amino acid
                resname = np.random.choice(amino_acids)

                # Write
                f.write(
                    f"ATOM  {atom_num:5d}  CA  {resname:3s} {chain_id}{res_num:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"  1.00 20.00           C  \n"
                )
                atom_num += 1

                # Track last residue for TER record
                last_resname = resname
                last_res_num = res_num

            # Write TER record after each chain
            f.write(
                f"TER   {atom_num:5d}      {last_resname:3s} {chain_id}{last_res_num:4d}\n"
            )
            atom_num += 1

        # END record
        f.write("END\n")
