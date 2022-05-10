import numpy as np

from haddock.libs.libontology import PDBFile


def calc_rmsd(V, W):
    """Calculate the RMSD from two vectors."""
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)


def kabsch(P, Q):
    """Find the rotation matrix using Kabsch algorithm."""
    # Covariance matrix
    P = np.array(P)
    Q = np.array(Q)
    C = np.dot(np.transpose(P), Q)
    # use SVD
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def centroid(X):
    """Get the centroid."""
    X = np.array(X)
    return X.mean(axis=0)


def load_coords(pdb_f, atoms, filter_resdic=None, numbering_dic=None):
    """Load coordinates from PDB."""
    coord_dic = {}
    chain_dic = {}
    idx = 0
    if isinstance(pdb_f, PDBFile):
        pdb_f = pdb_f.rel_path
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords = np.asarray([x, y, z])
                if numbering_dic:
                    try:
                        resnum = numbering_dic[chain][resnum]
                    except KeyError:
                        # this residue is not matched, and so it should
                        #  not be considered
                        # self.log(
                        #     f"WARNING: {chain}.{resnum}.{atom_name}"
                        #     " was not matched!"
                        #     )
                        continue
                # identifier = f"{chain}.{resnum}.{atom_name}"
                identifier = (chain, resnum, atom_name)
                if atom_name not in atoms[resname]:
                    continue
                if chain not in chain_dic:
                    chain_dic[chain] = []
                if filter_resdic:
                    # Only retrieve coordinates from the filter_resdic
                    if (
                            chain in filter_resdic
                            and resnum in filter_resdic[chain]
                            ):
                        coord_dic[identifier] = coords
                        chain_dic[chain].append(idx)
                        idx += 1
                else:
                    # retrieve everything
                    coord_dic[identifier] = coords
                    chain_dic[chain].append(idx)
                    idx += 1
    chain_ranges = {}
    for chain in chain_dic:
        min_idx = min(chain_dic[chain])
        max_idx = max(chain_dic[chain])
        chain_ranges[chain] = (min_idx, max_idx)
    return coord_dic, chain_ranges
