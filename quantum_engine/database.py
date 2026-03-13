# =========================================================
# Disease Targets
# =========================================================

DISEASE_TARGETS = {

    "HIV_Protease": {
        "disease_name": "HIV / AIDS",
        "pdb_id": "1HVR",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "ASP",
                "residue_number": 25,
                "atom_names": ["CG", "OD1", "OD2"]
            },
            {
                "chain_id": "A",
                "residue_name": "ILE",
                "residue_number": 50,
                "atom_names": ["CG1", "CD1"]
            }
        ],

        "drugs": [
            "Amide_Inhibitor",
            "Hydroxyl_Inhibitor"
        ]
    },


    "SARS_CoV2_Main_Protease": {
        "disease_name": "COVID-19",
        "pdb_id": "6LU7",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 41,
                "atom_names": ["ND1", "NE2"]
            },
            {
                "chain_id": "A",
                "residue_name": "CYS",
                "residue_number": 145,
                "atom_names": ["SG"]
            }
        ],

        "drugs": [
            "Amide_Inhibitor",
            "Carbonyl_Fragment"
        ]
    },


    "COX2_Inflammation": {
        "disease_name": "Inflammation / Pain",
        "pdb_id": "5IKR",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "ARG",
                "residue_number": 120,
                "atom_names": ["CZ", "NH1", "NH2"]
            },
            {
                "chain_id": "A",
                "residue_name": "TYR",
                "residue_number": 355,
                "atom_names": ["OH"]
            }
        ],

        "drugs": [
            "Carboxyl_Fragment",
            "Aromatic_Fragment"
        ]
    },


    "Acetylcholinesterase": {
        "disease_name": "Alzheimer's Disease",
        "pdb_id": "4EY7",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "SER",
                "residue_number": 203,
                "atom_names": ["OG"]
            },
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 447,
                "atom_names": ["ND1", "NE2"]
            }
        ],

        "drugs": [
            "Amide_Inhibitor",
            "Ester_Fragment"
        ]
    },


    "Carbonic_Anhydrase": {
        "disease_name": "Glaucoma / Altitude Sickness",
        "pdb_id": "1CA2",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 94,
                "atom_names": ["ND1", "NE2"]
            },
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 96,
                "atom_names": ["ND1", "NE2"]
            }
        ],

        "drugs": [
            "Sulfonamide_Fragment",
            "Hydroxyl_Inhibitor"
        ]
    }

}


# =========================================================
# Drug Fragments (Quantum Friendly)
# =========================================================

DRUGS = {

    "Amide_Inhibitor": {
        "atom_string": (
            "C 0 0 0; "
            "O 0 0 1.2; "
            "N 1.1 0 0; "
            "H 0 1 0; "
            "H 1 1 0"
        )
    },

    "Hydroxyl_Inhibitor": {
        "atom_string": (
            "O 0 0 0; "
            "H 0 0 0.96"
        )
    },

    "Carbonyl_Fragment": {
        "atom_string": (
            "C 0 0 0; "
            "O 0 0 1.2"
        )
    },

    "Carboxyl_Fragment": {
        "atom_string": (
            "C 0 0 0; "
            "O 0 1.2 0; "
            "O 0 -1.2 0; "
            "H 1 0 0"
        )
    },

    "Aromatic_Fragment": {
        "atom_string": (
            "C 0 0 0; "
            "C 1.4 0 0; "
            "C 2.1 1.2 0; "
            "H 0 1 0"
        )
    },

    "Ester_Fragment": {
        "atom_string": (
            "C 0 0 0; "
            "O 0 0 1.2; "
            "O 1.2 0 0; "
            "H 1.2 1 0"
        )
    },

    "Sulfonamide_Fragment": {
        "atom_string": (
            "S 0 0 0; "
            "O 0 1.4 0; "
            "O 0 -1.4 0; "
            "N 1.2 0 0; "
            "H 1.8 0.5 0"
        )
    }

}