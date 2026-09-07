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


# Full alignement of an ensemble onto a reference - user-given chains used
# FIXME: This function has too many arguments, needs to be refactored
def align_full_ens(
    ref_path,
    ensemble_path,
    chains=["A", "B"],
    width=800,
    height=600,
    ref_colors={"A": "red", "B": "orange", "C": "pink"},
    ensemble_color="lightgrey",
    atom_types=["P", "C1", "CA"],
    show_labels=False,
    show_per_model_rmsd=True,
    animate=False,
    interval=1000,
    show_model_number=False,
):
    """."""
    """
    Align an ensemble of models onto a reference structure and visualize
    with py3Dmol.

    Each model of the multi-model ensemble PDB file is superimposed
    independently onto the reference structure using the specified chains
    and atom types.

    When ``animate`` is True, the aligned ensemble is added as animation
    frames (one conformer shown at a time, cycling) while the reference
    structure stays permanently displayed.

    Parameters:
    -----------
    ref_path : str
        Path to the reference PDB file (single model)
    ensemble_path : str
        Path to the multi-model ensemble PDB file whose models are aligned
        to the reference
    chains : list, default ['A', 'B']
        List of chain IDs to include in alignment
    width : int, default 800
        Viewer width in pixels
    height : int, default 600
        Viewer height in pixels
    ref_colors : dict, default {'A': 'red', 'B': 'orange', 'C': 'pink'}
        Colors for the chains of the reference structure
    ensemble_color : str, default 'lightgrey'
        Color used for every aligned ensemble model
    atom_types : list, default ['P', 'C1', 'CA']
        Atom types to use for alignment (['CA'] or ['CA', 'CB', 'N', 'C'])
    show_labels : bool, default False
        Whether to show descriptive labels
    show_per_model_rmsd : bool, default True
        Whether to calculate and display the RMSD of each aligned model
    animate : bool, default False
        If True, add the aligned ensemble as animation frames (one
        conformer displayed at a time, cycling) with the reference kept
        static, instead of showing all conformers superimposed at once
    interval : int, default 1000
        Delay in milliseconds between animation frames (only used when
        ``animate`` is True)
    show_model_number : bool, default False
        If True (and ``animate`` is True), display the model number of the
        conformer currently shown as a label that updates with each frame

    Returns:
    --------
    py3Dmol.view object

    Example:
    --------
    align_full_ens('reference.pdb', 'ensemble.pdb.gz')
    """

    def get_atoms_from_chains(entity, chain_ids, atom_types):
        """Collect atoms from the given chains of a model (or structure)."""
        atoms = []
        chain_info = {}

        # A Structure iterates over its Models; a Model iterates over its
        # Chains. Normalise so we always iterate over models.
        if hasattr(entity, "get_models"):
            models = entity.get_models()
        else:
            models = [entity]

        for model in models:
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
    ref_data = load_pdb_file(ref_path)
    ensemble_data = load_pdb_file(ensemble_path)

    if not (ref_data and ensemble_data):
        print("Failed to load the reference and/or the ensemble PDB file")
        return view

    per_model_rmsd = {}
    n_ensemble_models = 0
    aligned_models = []

    try:
        # Parse structures
        ref_struct = pdb_string_to_structure(ref_data, "reference")
        ens_struct = pdb_string_to_structure(ensemble_data, "ensemble")

        # Reference atoms are taken from its first model
        ref_model = next(ref_struct.get_models())
        ref_atoms, ref_chain_info = get_atoms_from_chains(ref_model, chains, atom_types)
        print(
            f"Atoms for alignment - Reference: {ref_chain_info}, "
            f"Total: {len(ref_atoms)}"
        )
        print(
            "Reference: "
            + "; ".join(
                [
                    f"chain {chain} in {color}"
                    for chain, color in ref_colors.items()
                    if chain in chains
                ]
            )
        )

        if len(ref_atoms) == 0:
            print("Could not find any reference atoms for alignment")
            view.addModel(ref_data, "pdb")
            view.addModel(ensemble_data, "pdb")
        else:
            # Reference is the first model added to the viewer
            view.addModel(structure_to_pdb_string(ref_struct), "pdb")

            # Align each ensemble model independently onto the reference
            for model in ens_struct.get_models():
                model_id = model.id
                ens_atoms, _ = get_atoms_from_chains(model, chains, atom_types)

                if len(ens_atoms) == 0:
                    print(f"Model {model_id}: no atoms for alignment, skipping")
                    continue

                min_atoms = min(len(ref_atoms), len(ens_atoms))
                sup = Superimposer()
                sup.set_atoms(ref_atoms[:min_atoms], ens_atoms[:min_atoms])
                per_model_rmsd[model_id] = sup.rms

                # Apply the transformation to every atom of this model only
                sup.apply(model.get_atoms())
                aligned_models.append(model)

                if show_per_model_rmsd:
                    print(
                        f"Model {model_id}: RMSD {sup.rms:.3f} Å "
                        f"({min_atoms} atom pairs)"
                    )

                if not animate:
                    # Show all conformers at once: add each aligned model as
                    # its own viewer model (models 1..N).
                    view.addModel(_model_to_pdb_string(model), "pdb")

                n_ensemble_models += 1

            if animate and aligned_models:
                # Concatenate the aligned models into a single multi-MODEL
                # PDB string and add them as animation frames of one model.
                # 3Dmol.js cycles through frames while the reference (model 0)
                # stays static.
                frames_pdb = "\n".join(_model_to_pdb_string(m) for m in aligned_models)
                view.addModelsAsFrames(frames_pdb, "pdb")

                if show_model_number:
                    # Attach one label per frame. 3Dmol.js only shows a label
                    # whose `frame` matches the frame currently displayed, so
                    # the number updates as the animation cycles. Frames are
                    # 0-indexed; display the model's own serial number.
                    for frame_idx, m in enumerate(aligned_models):
                        model_no = getattr(m, "serial_num", None) or frame_idx + 1
                        view.addLabel(
                            f"Model {model_no}",
                            {
                                "frame": frame_idx,
                                "useScreen": True,
                                "position": {"x": 10, "y": 10, "z": 0},
                                "backgroundColor": "black",
                                "backgroundOpacity": 0.7,
                                "fontColor": "white",
                                "fontSize": 14,
                            },
                        )

            print(f"Aligned {n_ensemble_models} ensemble model(s) onto reference")

    except Exception as e:
        print(f"Alignment failed: {e}")
        view.addModel(ref_data, "pdb")
        view.addModel(ensemble_data, "pdb")

    # Apply styling: model 0 is the reference, models 1..N are the ensemble
    view.setStyle({"model": 0}, {"cartoon": {}})
    for chain, color in ref_colors.items():
        if chain in chains:
            view.addStyle(
                {"model": 0, "chain": chain},
                {"cartoon": {"color": color, "opacity": 0.9}},
            )

    # Colour the ensemble with a single colour. When animating, all frames
    # live in one model (index 1); otherwise every conformer is its own
    # model (indices 1..N).
    if animate:
        ensemble_model_indices = [1] if n_ensemble_models else []
    else:
        ensemble_model_indices = range(1, n_ensemble_models + 1)

    for model_idx in ensemble_model_indices:
        view.setStyle(
            {"model": model_idx},
            {"cartoon": {"color": ensemble_color, "opacity": 0.6}},
        )

    # Add labels if requested
    if show_labels:
        view.addLabel(
            "Reference",
            {
                "position": {"x": -20, "y": 20, "z": 0},
                "backgroundColor": "darkred",
                "fontColor": "white",
            },
        )
        view.addLabel(
            f"Ensemble ({n_ensemble_models} models)",
            {
                "position": {"x": 20, "y": 20, "z": 0},
                "backgroundColor": "gray",
                "fontColor": "white",
            },
        )
        if per_model_rmsd:
            mean_rmsd = sum(per_model_rmsd.values()) / len(per_model_rmsd)
            view.addLabel(
                f"Mean RMSD: {mean_rmsd:.3f} Å",
                {
                    "position": {"x": 0, "y": -20, "z": 0},
                    "backgroundColor": "purple",
                    "fontColor": "white",
                },
            )

    view.zoomTo()

    # Start cycling through the ensemble frames (reference stays static)
    if animate and n_ensemble_models:
        view.animate({"loop": "forward", "reps": 0, "interval": interval})

    return view


def _model_to_pdb_string(model):
    """Serialise a single Bio.PDB Model to a PDB-formatted string."""
    pdb_io = PDBIO()
    pdb_io.set_structure(model)
    output = StringIO()
    pdb_io.save(output)
    return output.getvalue()
