import numpy as np
from Bio.PDB import PDBParser


def generate_pose_initial_params(pose_file, param_count):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pose", pose_file)

    coords = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:

                    if atom.element == "H":
                        continue

                    coords.append(atom.coord)

    coords = np.array(coords)

    if len(coords) < 2:
        return np.zeros(param_count)

    # compute pairwise distances
    dists = []

    for i in range(len(coords)):
        for j in range(i+1,len(coords)):

            dist = np.linalg.norm(coords[i]-coords[j])

            dists.append(dist)

    dists = np.array(dists)

    # convert distances → amplitudes
    params = np.exp(-dists)

    if len(params) == 0:
        params = np.zeros(param_count)

    # resize to match ansatz parameter count
    params = np.resize(params,param_count)

    return params