import os
import subprocess
from Bio.PDB import PDBList, PDBParser, PDBIO, Select


# --------------------------------------------
# Download PDB
# --------------------------------------------

def download_pdb(pdb_id, folder="data/pdb"):

    os.makedirs(folder, exist_ok=True)

    pdb_path = os.path.join(folder, f"{pdb_id}.pdb")

    if os.path.exists(pdb_path):
        print(f"PDB {pdb_id} already exists")
        return pdb_path

    print(f"Downloading PDB {pdb_id}...")

    pdbl = PDBList()

    downloaded = pdbl.retrieve_pdb_file(
        pdb_id,
        pdir=folder,
        file_format="pdb"
    )

    os.rename(downloaded, pdb_path)

    return pdb_path


# --------------------------------------------
# Remove waters + ligands
# --------------------------------------------

class ProteinOnlySelect(Select):

    def accept_residue(self, residue):

        # Remove water molecules
        if residue.get_resname() == "HOH":
            return False

        # Keep only standard protein residues
        if residue.id[0] != " ":
            return False

        return True


def clean_protein(pdb_path):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    clean_path = pdb_path.replace(".pdb", "_clean.pdb")

    io = PDBIO()
    io.set_structure(structure)
    io.save(clean_path, ProteinOnlySelect())

    return clean_path


# --------------------------------------------
# Convert to PDBQT with hydrogens
# --------------------------------------------

import subprocess
import os


def convert_to_pdbqt(pdb_path):

    pdbqt_folder = "data/pdbqt"
    os.makedirs(pdbqt_folder, exist_ok=True)

    name = os.path.basename(pdb_path).replace(".pdb", "")
    pdbqt_path = os.path.join(pdbqt_folder, f"{name}.pdbqt")

    if os.path.exists(pdbqt_path):
        return pdbqt_path

    # Convert with Open Babel
    subprocess.run([
        "obabel",
        pdb_path,
        "-O",
        pdbqt_path,
        "-xr",   # remove torsion tree → rigid receptor
        "-h"     # add hydrogens
    ], check=True)

    return pdbqt_path


# --------------------------------------------
# Full pipeline
# --------------------------------------------

def prepare_protein_for_docking(pdb_id):

    pdb_path = download_pdb(pdb_id)

    clean_path = clean_protein(pdb_path)

    pdbqt_path = convert_to_pdbqt(clean_path)

    return pdbqt_path