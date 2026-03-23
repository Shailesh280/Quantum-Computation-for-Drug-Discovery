from Bio.PDB import PDBParser
import numpy as np
import os

VALID_ELEMENTS = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br"}

def build_ligand_atom_string_from_pose(
    pose_file,
    pocket_file=None,
    max_atoms=12
):

    parser = PDBParser(QUIET=True)

    ligand_structure = parser.get_structure("ligand", pose_file)

    ligand_atoms = []

    # -------------------------------------------------
    # Extract all valid ligand atoms
    # -------------------------------------------------

    for atom in ligand_structure.get_atoms():

        element = atom.element.strip()

        if not element:
            name = atom.get_name().strip()

            if name[:2].capitalize() in VALID_ELEMENTS:
                element = name[:2].capitalize()
            else:
                element = name[0].upper()

        if element not in VALID_ELEMENTS:
            continue

        ligand_atoms.append((element, atom.coord))

    if not ligand_atoms:
        raise ValueError("No ligand atoms extracted")

    # -------------------------------------------------
    # If pocket exists → rank ligand atoms by proximity
    # -------------------------------------------------

    if pocket_file and os.path.exists(pocket_file):

        pocket_structure = parser.get_structure("pocket", pocket_file)
        pocket_coords = [atom.coord for atom in pocket_structure.get_atoms()]

        ranked_atoms = []

        for element, coord in ligand_atoms:

            min_dist = min(
                np.linalg.norm(coord - p) for p in pocket_coords
            )

            ranked_atoms.append((element, coord, min_dist))

        # Sort by distance to pocket
        ranked_atoms.sort(key=lambda x: x[2])

        ligand_atoms = [(e, c) for e, c, _ in ranked_atoms]

    # -------------------------------------------------
    # Limit atom count
    # -------------------------------------------------

    ligand_atoms = ligand_atoms[:max_atoms]

    # Debug
    print(f"[DEBUG] Ligand atoms selected: {len(ligand_atoms)}")

    atom_strings = [
        f"{element} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}"
        for element, coord in ligand_atoms
    ]

    return "; ".join(atom_strings)