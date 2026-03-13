from vina_engine import run_vina

result = run_vina(
    receptor="/mnt/d/Major project/Source code/pdb_files/protein.pdbqt",
    ligand="/mnt/d/Major project/Source code/pdb_files/ligand.pdbqt"
)

print("\nBest binding energy:", result["affinity"], "kcal/mol")
print("Structure file:", result["structure_file"])