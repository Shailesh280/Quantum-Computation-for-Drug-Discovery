from Bio.PDB import PDBParser


# Allowed elements for quantum simulation
VALID_ELEMENTS = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br"}


# ==========================================================
# Extract Ligand Coordinates From Docking Pose
# ==========================================================

def build_ligand_atom_string_from_pose(pose_file, max_atoms=10):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pose", pose_file)

    atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:

                    element = atom.element.strip()

                    if not element:
                        name = atom.get_name().strip()

                        if name[:2].capitalize() in VALID_ELEMENTS:
                            element = name[:2].capitalize()
                        else:
                            element = name[0].upper()

                    if element not in VALID_ELEMENTS:
                        continue

                    coord = atom.coord

                    atoms.append(
                        f"{element} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}"
                    )

                    if len(atoms) >= max_atoms:
                        return "; ".join(atoms)

    if not atoms:
        raise ValueError("No valid atoms extracted from docking pose")

    return "; ".join(atoms)