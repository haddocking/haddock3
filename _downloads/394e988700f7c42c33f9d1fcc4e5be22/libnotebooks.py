"""
Helper functions for HADDOCK3 notebooks
"""

import gzip
import os
from io import StringIO

from Bio.PDB import PDBIO, PDBParser, Superimposer

# NOTE: `py3Dmol` is an optional dependency, it might not be installed
#  this is why we have the conditional import here
try:
    import py3Dmol
except ImportError:
    py3Dmol = None


def load_pdb_file(file_path):
    """."""
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return None

    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as f:
            return f.read()
    else:
        with open(file_path, "r") as f:
            return f.read()


def pdb_string_to_structure(pdb_string, structure_id):
    """."""
    parser = PDBParser(QUIET=True)
    pdb_io = StringIO(pdb_string)
    structure = parser.get_structure(structure_id, pdb_io)
    return structure


def structure_to_pdb_string(structure):
    """."""
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    output = StringIO()
    pdb_io.save(output)
    return output.getvalue()


# Full alignement - user-given chains used
# FIXME: This function has too many arguments, needs to be refactored
def align_full(
    pdb_path1,
    pdb_path2,
    chains=["A", "B"],
    width=800,
    height=600,
    model1_colors={"A": "red", "B": "orange", "C": "pink"},
    model2_colors={"A": "blue", "B": "green", "C": "lime"},
    atom_types=["P", "C1", "CA"],
    show_labels=False,
    show_per_chain_rmsd=True,
):
    """."""
    """
    Align two protein structures using all specified chains and visualize with py3Dmol.

    Parameters:
    -----------
    pdb_path1 : str
        Path to the first (reference) PDB file
    pdb_path2 : str
        Path to the second PDB file to align to the first
    chains : list, default ['A', 'B']
        List of chain IDs to include in alignment
    width : int, default 800
        Viewer width in pixels
    height : int, default 600
        Viewer height in pixels
    model1_colors : dict, default {'A': 'red', 'B': 'orange'}
        Colors for chains in model 1
    model2_colors : dict, default {'A': 'blue', 'B': 'green'}
        Colors for chains in model 2
    atom_types : list, default ['CA', 'P', 'C1']
        Atom types to use for alignment (['CA'] or ['CA', 'CB', 'N', 'C'])
    show_labels : bool, default False
        Whether to show descriptive labels
    show_per_chain_rmsd : bool, default True
        Whether to calculate and display per-chain RMSD values

    Returns:
    --------
    py3Dmol.view object

    Example:
    --------
    align_full_molecule('model1.pdb.gz', 'model2.pdb.gz')
    """

    def get_atoms_from_chains(structure, chain_ids, atom_types):
        atoms = []
        chain_info = {}

        for model in structure:
            for chain_id in chain_ids:
                if chain_id in model:
                    chain = model[chain_id]
                    chain_atoms = []
                    for residue in chain:
                        for atom_type in atom_types:
                            if atom_type in residue:
                                atoms.append(residue[atom_type])
                                chain_atoms.append(residue[atom_type])
                    chain_info[chain_id] = len(chain_atoms)

        return atoms, chain_info

    # Create viewer
    view = py3Dmol.view(width=width, height=height)

    # Load PDB files
    model_1_data = load_pdb_file(pdb_path1)
    model_2_data = load_pdb_file(pdb_path2)

    if not (model_1_data and model_2_data):
        print("Failed to load one or both PDB files")
        return view, None, {}

    overall_rmsd = None
    per_chain_rmsd = {}

    try:
        # Parse structures
        struct1 = pdb_string_to_structure(model_1_data, "model1")
        struct2 = pdb_string_to_structure(model_2_data, "model2")

        # Get atoms from all specified chains
        atoms_1, chain_info_1 = get_atoms_from_chains(struct1, chains, atom_types)
        atoms_2, chain_info_2 = get_atoms_from_chains(struct2, chains, atom_types)

        print(f"Atoms for alignment - Model 1: {chain_info_1}, Total: {len(atoms_1)}")
        print(f"Atoms for alignment - Model 2: {chain_info_2}, Total: {len(atoms_2)}")
        print(
            "Model 1: "
            + "; ".join(
                [
                    f"chain {chain} in {color}"
                    for chain, color in model1_colors.items()
                    if chain in chains
                ]
            )
        )
        print(
            "Model 2: "
            + "; ".join(
                [
                    f"chain {chain} in {color}"
                    for chain, color in model2_colors.items()
                    if chain in chains
                ]
            )
        )

        if len(atoms_1) > 0 and len(atoms_2) > 0:
            # Align using all atoms from specified chains
            min_atoms = min(len(atoms_1), len(atoms_2))
            ref_atoms = atoms_1[:min_atoms]
            alt_atoms = atoms_2[:min_atoms]

            print(f"Using {min_atoms} atom pairs for alignment")

            # Perform superimposition
            sup = Superimposer()
            sup.set_atoms(ref_atoms, alt_atoms)
            overall_rmsd = sup.rms

            print(f"Whole molecule alignment RMSD: {overall_rmsd:.3f} Å")

            # Apply transformation to all atoms in structure 2
            sup.apply(struct2.get_atoms())

            # Calculate per-chain RMSD if requested
            if show_per_chain_rmsd:
                for chain_id in chains:
                    try:
                        chain_atoms_1, _ = get_atoms_from_chains(
                            struct1, [chain_id], atom_types
                        )
                        chain_atoms_2, _ = get_atoms_from_chains(
                            struct2, [chain_id], atom_types
                        )
                        if len(chain_atoms_1) > 0 and len(chain_atoms_2) > 0:
                            min_chain = min(len(chain_atoms_1), len(chain_atoms_2))
                            sup_chain = Superimposer()
                            sup_chain.set_atoms(
                                chain_atoms_1[:min_chain], chain_atoms_2[:min_chain]
                            )
                            per_chain_rmsd[chain_id] = sup_chain.rms
                            print(f"Chain {chain_id} RMSD: {sup_chain.rms:.3f} Å")
                    except Exception as e:
                        print(f"Could not calculate RMSD for chain {chain_id}: {e}")

            # Convert back to PDB strings
            aligned_pdb_1 = structure_to_pdb_string(struct1)
            aligned_pdb_2 = structure_to_pdb_string(struct2)

            # Add models to viewer
            view.addModel(aligned_pdb_1, "pdb")
            view.addModel(aligned_pdb_2, "pdb")

        else:
            print(
                "Could not find sufficient atoms for alignment, adding original models"
            )
            view.addModel(model_1_data, "pdb")
            view.addModel(model_2_data, "pdb")

    except Exception as e:
        print(f"Alignment failed: {e}")
        view.addModel(model_1_data, "pdb")
        view.addModel(model_2_data, "pdb")

    # Apply styling
    view.setStyle({"model": 0}, {"cartoon": {}})
    view.setStyle({"model": 1}, {"cartoon": {}})

    for chain, color in model1_colors.items():
        if chain in chains:
            view.addStyle(
                {"model": 0, "chain": chain},
                {"cartoon": {"color": color, "opacity": 0.9}},
            )

    for chain, color in model2_colors.items():
        if chain in chains:
            view.addStyle(
                {"model": 1, "chain": chain},
                {"cartoon": {"color": color, "opacity": 0.6}},
            )

    # Add labels if requested
    if show_labels:
        view.addLabel(
            "Model 1 (Reference)",
            {
                "position": {"x": -20, "y": 20, "z": 0},
                "backgroundColor": "darkred",
                "fontColor": "white",
            },
        )
        view.addLabel(
            "Model 2 (Aligned)",
            {
                "position": {"x": 20, "y": 20, "z": 0},
                "backgroundColor": "darkgreen",
                "fontColor": "white",
            },
        )
        view.addLabel(
            "Full Molecule Alignment",
            {
                "position": {"x": 0, "y": -20, "z": 0},
                "backgroundColor": "navy",
                "fontColor": "white",
            },
        )
        if overall_rmsd:
            view.addLabel(
                f"Overall RMSD: {overall_rmsd:.3f} Å",
                {
                    "position": {"x": 0, "y": -30, "z": 0},
                    "backgroundColor": "purple",
                    "fontColor": "white",
                },
            )

            # Add per-chain RMSD labels
            y_offset = -40
            for chain_id, rmsd_val in per_chain_rmsd.items():
                view.addLabel(
                    f"Chain {chain_id}: {rmsd_val:.3f} Å",
                    {
                        "position": {"x": 0, "y": y_offset, "z": 0},
                        "backgroundColor": "gray",
                        "fontColor": "white",
                        "fontSize": 10,
                    },
                )
                y_offset -= 8

    view.zoomTo()

    return view
