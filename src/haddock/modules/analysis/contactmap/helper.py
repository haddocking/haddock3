import os

import psutil


def get_available_memory() -> float:
    """
    Get the available system memory in GB.

    Returns
    -------
    float
        Available memory in gigabytes (GB).
    """
    return psutil.virtual_memory().available / (1024**3)


def get_necessary_memory(models: list) -> float:
    """
    Estimate the memory required to compute contact maps.

    Calculates memory based on the largest model's estimated atom count.
    The memory estimate assumes a distance matrix of size (atoms x atoms) * 8 bytes.

    Parameters
    ----------
    models : list
        List of model objects with file_name attribute.

    Returns
    -------
    float
        Estimated memory requirement in gigabytes (GB).
    """
    atoms = 0
    for model in models:
        try:
            # use `getsize` because its faster than reading the whole PDB
            file_size = os.path.getsize(model.file_name)
            # estimate: ~10 bytes per atom in PDB file
            estimated_atoms = file_size // 10
            atoms = max(atoms, estimated_atoms)
        except Exception:
            # skip this model
            continue

    if atoms == 0:
        # Fall back to a safe default if we can't estimate
        atoms = 10000

    matrix_size_bytes = (atoms * atoms) * 8
    matrix_size_gb = matrix_size_bytes / (1024**3)

    return matrix_size_gb
