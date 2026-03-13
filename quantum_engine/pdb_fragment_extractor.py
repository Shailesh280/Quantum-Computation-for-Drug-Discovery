import os
import numpy as np
from Bio.PDB import PDBList, PDBParser


VALENCE = {
    "C": 4,
    "N": 3,
    "O": 2,
    "S": 2
}

ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16
}

MAX_HEAVY_ATOMS = 6


# -----------------------------------
# PDB Download
# -----------------------------------

def download_pdb(pdb_id, save_dir="pdb_files"):

    os.makedirs(save_dir, exist_ok=True)

    pdbl = PDBList()

    file_path = pdbl.retrieve_pdb_file(
        pdb_id,
        pdir=save_dir,
        file_format="pdb"
    )

    return file_path


# -----------------------------------
# Bond Estimation
# -----------------------------------

def estimate_bonds(atoms):

    bonds = {i: 0 for i in range(len(atoms))}

    for i, (elem1, coord1) in enumerate(atoms):
        for j, (elem2, coord2) in enumerate(atoms):

            if i >= j:
                continue

            distance = np.linalg.norm(coord1 - coord2)

            if distance < 1.8:
                bonds[i] += 1
                bonds[j] += 1

    return bonds


# -----------------------------------
# Hydrogen Capping
# -----------------------------------

def cap_with_hydrogens(atoms):

    capped_atoms = atoms.copy()

    bonds = estimate_bonds(atoms)

    directions = [
        np.array([1,0,0]),
        np.array([-1,0,0]),
        np.array([0,1,0]),
        np.array([0,-1,0]),
        np.array([0,0,1]),
        np.array([0,0,-1]),
    ]

    for i, (element, coord) in enumerate(atoms):

        if element not in VALENCE:
            continue

        expected = VALENCE[element]
        current = bonds[i]

        missing = expected - current

        for m in range(missing):

            direction = directions[m % len(directions)]
            direction = direction / np.linalg.norm(direction)

            h_coord = coord + direction * 1.1

            capped_atoms.append(("H", h_coord))

    return capped_atoms


# -----------------------------------
# Geometry Stabilization
# -----------------------------------

def stabilize_geometry(atoms, min_dist=0.8):

    atoms = atoms.copy()

    for _ in range(10):

        moved = False

        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):

                elem_i, coord_i = atoms[i]
                elem_j, coord_j = atoms[j]

                d_vec = coord_j - coord_i
                dist = np.linalg.norm(d_vec)

                if dist < min_dist:

                    direction = d_vec / (dist + 1e-6)

                    shift = (min_dist - dist) * 0.5

                    coord_i = coord_i - direction * shift
                    coord_j = coord_j + direction * shift

                    atoms[i] = (elem_i, coord_i)
                    atoms[j] = (elem_j, coord_j)

                    moved = True

        if not moved:
            break

    return atoms


# -----------------------------------
# Charge & Spin
# -----------------------------------

def compute_charge_and_spin(atoms):

    total_electrons = sum(ATOMIC_NUMBERS[a[0]] for a in atoms)

    charge = 0

    if total_electrons % 2 != 0:

        for i in range(len(atoms)-1, -1, -1):

            if atoms[i][0] == "H":
                atoms.pop(i)
                total_electrons -= 1
                break

    spin = 0

    return charge, spin


# -----------------------------------
# Residue Fragment Extraction
# -----------------------------------

def extract_residue_fragment(
    pdb_id,
    chain_id,
    residue_name,
    residue_number,
    atom_names
):

    pdb_file = download_pdb(pdb_id)

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(pdb_id, pdb_file)

    heavy_atoms = []

    for model in structure:
        for chain in model:

            if chain.id != chain_id:
                continue

            for residue in chain:

                if (
                    residue.get_resname() == residue_name
                    and residue.get_id()[1] == residue_number
                ):

                    for atom in residue:

                        if atom.get_name() in atom_names:

                            element = atom.element

                            if not element:
                                element = atom.get_name()[0].upper()

                            if element == "H":
                                continue

                            heavy_atoms.append(
                                (element, atom.coord)
                            )

                            if len(heavy_atoms) >= MAX_HEAVY_ATOMS:
                                break

    if not heavy_atoms:
        raise ValueError("Residue not found.")

    atoms = cap_with_hydrogens(heavy_atoms)

    atoms = stabilize_geometry(atoms)

    charge, spin = compute_charge_and_spin(atoms)

    atom_string = "; ".join(
        f"{elem} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}"
        for elem, coord in atoms
    )

    return atom_string, charge, spin