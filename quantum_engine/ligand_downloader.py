import os
import requests
import subprocess


def download_ligand(drug_name):

    ligand_folder = "data/ligands"
    os.makedirs(ligand_folder, exist_ok=True)

    sdf_path = os.path.join(ligand_folder, f"{drug_name}.sdf")
    pdbqt_path = os.path.join(ligand_folder, f"{drug_name}.pdbqt")

    if os.path.exists(pdbqt_path):
        return pdbqt_path

    print(f"Downloading ligand {drug_name} from PubChem...")

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/SDF"

    r = requests.get(url)

    if r.status_code != 200:
        raise RuntimeError(f"Failed to download ligand {drug_name}")

    with open(sdf_path, "wb") as f:
        f.write(r.content)

    # convert to pdbqt
    subprocess.run([
        "obabel",
        sdf_path,
        "-O",
        pdbqt_path,
        "-h"
    ], check=True)

    return pdbqt_path