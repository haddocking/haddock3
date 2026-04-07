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
    """Generate a synthetic PDB file for testing."""
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
    per_chain = [n_atoms // n_chains] * n_chains
    for i in range(n_atoms % n_chains):
        per_chain[i] += 1

    with open(output_path, "w") as f:
        # Write header
        f.write("REMARK   Synthetic PDB generated for integration testing\n")
        f.write(f"REMARK   Total atoms: {n_atoms}, Chains: {n_chains}\n")

        atom_num = 1

        for chain_idx, chain_id in enumerate(chain_ids):
            n_residues = per_chain[chain_idx]

            if n_residues == 0:
                continue

            last_resname = "ALA"
            last_res_num = 1

            for res_num in range(1, n_residues + 1):
                # Random coordinates within a reasonable range
                x = np.random.uniform(0, 50)
                y = np.random.uniform(0, 50)
                z = np.random.uniform(0, 50)

                # Random amino acid
                resname = np.random.choice(amino_acids)

                for _ in range(1, 6):
                    # Add small offset for each atom in the residue
                    x_atom = x + np.random.uniform(-0.5, 0.5)
                    y_atom = y + np.random.uniform(-0.5, 0.5)
                    z_atom = z + np.random.uniform(-0.5, 0.5)

                    f.write(
                        f"ATOM  {atom_num:5d}  X   {resname:3s} {chain_id}{res_num:4d}    "
                        f"{x_atom:8.3f}{y_atom:8.3f}{z_atom:8.3f}"
                        f"  1.00 20.00           C  \n"
                    )
                    atom_num += 1

                # Track last residue for TER record
                last_resname = resname
                last_res_num = res_num

            # Write TER record after each chain
            f.write(
                f"TER   {atom_num:5d}      {last_resname:3s} {chain_id}{last_res_num:4d}    \n"
            )
            atom_num += 1

        # END record
        f.write("END\n")
