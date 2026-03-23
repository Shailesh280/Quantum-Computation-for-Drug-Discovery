from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np


class PocketSelect(Select):

    def __init__(self, pocket_atoms):
        self.pocket_coords = {tuple(atom.coord) for atom in pocket_atoms}

    def accept_atom(self, atom):
        return tuple(atom.coord) in self.pocket_coords

def extract_pocket(protein_pdb, ligand_pose, output_file, cutoff=8.0):
    parser = PDBParser(QUIET=True)

    protein = parser.get_structure("protein", protein_pdb)
    ligand = parser.get_structure("ligand", ligand_pose)

    ligand_atoms = [atom.coord for atom in ligand.get_atoms()]

    pocket_atoms = []

    for atom in protein.get_atoms():
        coord = atom.coord

        for l in ligand_atoms:
            dist = np.linalg.norm(coord - l)

            if dist < cutoff:
                pocket_atoms.append(atom)
                break

    io = PDBIO()
    io.set_structure(protein)
    io.save(output_file, PocketSelect(pocket_atoms))

    return output_file