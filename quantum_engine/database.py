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
                "residue_number": 25
            },
            {
                "chain_id": "A",
                "residue_name": "ILE",
                "residue_number": 50
            }
        ],

        # Real HIV protease inhibitors
        "drugs": [
            "ritonavir",
            "indinavir",
            "saquinavir"
        ]
    },


    "SARS_CoV2_Main_Protease": {

        "disease_name": "COVID-19",
        "pdb_id": "6LU7",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 41
            },
            {
                "chain_id": "A",
                "residue_name": "CYS",
                "residue_number": 145
            }
        ],

        "drugs": [
            "nirmatrelvir",
            "boceprevir",
            "lopinavir"
        ]
    },


    "COX2_Inflammation": {

        "disease_name": "Inflammation / Pain",
        "pdb_id": "5IKR",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "ARG",
                "residue_number": 120
            },
            {
                "chain_id": "A",
                "residue_name": "TYR",
                "residue_number": 355
            }
        ],

        "drugs": [
            "celecoxib",
            "rofecoxib",
            "valdecoxib"
        ]
    },


    "Acetylcholinesterase": {

        "disease_name": "Alzheimer's Disease",
        "pdb_id": "4EY7",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "SER",
                "residue_number": 203
            },
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 447
            }
        ],

        "drugs": [
            "donepezil",
            "rivastigmine",
            "galantamine"
        ]
    },


    "Carbonic_Anhydrase": {

        "disease_name": "Glaucoma / Altitude Sickness",
        "pdb_id": "1CA2",

        "active_sites": [
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 94
            },
            {
                "chain_id": "A",
                "residue_name": "HIS",
                "residue_number": 96
            }
        ],

        "drugs": [
            "acetazolamide",
            "dorzolamide",
            "methazolamide"
        ]
    }

}