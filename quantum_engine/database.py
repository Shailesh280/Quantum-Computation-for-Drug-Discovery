DISEASE_TARGETS = {
    "HIV_Protease": {
        "pdb_id": "1HVR",
        "chain_id": "A",
        "residue_name": "ASP",
        "residue_number": 25,
        "atom_names": ["CG", "OD1", "OD2"]
    }
}

# Hydrogen-capped drug fragments (closed shell)
DRUGS = {
    "Amide_Inhibitor": {
        "atom_string": (
            "C 0 0 0; "
            "O 0 0 1.2; "
            "N 1.1 0 0; "
            "H 0 1 0; "
            "H 1 1 0"
        )
    }
}
