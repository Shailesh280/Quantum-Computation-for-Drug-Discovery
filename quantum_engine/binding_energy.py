import os
from Bio.PDB import PDBParser
from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL
from .pose_fragment_builder import build_ligand_atom_string_from_pose
from .database import DISEASE_TARGETS
from .vina_engine import run_vina
from .pose_parameter_initializer import generate_pose_initial_params
from .pdb_utils import prepare_protein_for_docking
import subprocess
from .ligand_downloader import download_ligand
from .pocket_extractor import extract_pocket
from qiskit.utils import algorithm_globals
import numpy as np
import shutil

# ==========================================================
# Distance-Based Atom Selection
# ==========================================================

def select_closest_atoms(protein_atoms, ligand_coords, max_atoms=10):

    ranked_atoms = []

    for atom in protein_atoms:
        coord = atom.coord

        min_dist = min(
            np.linalg.norm(coord - l) for l in ligand_coords
        )

        ranked_atoms.append((atom, min_dist))

    ranked_atoms.sort(key=lambda x: x[1])

    selected_atoms = [atom for atom, _ in ranked_atoms[:max_atoms]]

    return selected_atoms


# ==========================================================
# Closed Shell Enforcement
# ==========================================================

def enforce_closed_shell(atom_string, basis, charge):

    problem = build_problem(atom_string, basis=basis, charge=charge)

    num_alpha, num_beta = problem.num_particles
    total = num_alpha + num_beta

    if total % 2 != 0:
        charge += 1
        problem = build_problem(atom_string, basis=basis, charge=charge)

    return problem, charge


# ==========================================================
# VQE Worker
# ==========================================================

def run_vqe_for_pose(data):

    complex_atom_string, pose_file = data

    print("Running VQE for pose:", pose_file, flush=True)

    problem, _ = enforce_closed_shell(complex_atom_string, "sto3g", 0)

    initial_params = generate_pose_initial_params(
        pose_file,
        param_count=20
    )

    result = compute_energies(
        problem,
        label="POSE",
        pose_initial=initial_params
    )

    if result.get("solver_error"):
        print("VQE solver error for pose:", pose_file, flush=True)
        return None

    print("VQE completed for pose:", pose_file, flush=True)

    return result


# ==========================================================
# Main Pipeline
# ==========================================================

def compute_binding_from_selection(
    disease_key,
    drug_name,
    separation_distance=3.0,
    basis="sto3g"
):

    try:

        # =====================================================
        # Disease Info
        # =====================================================

        disease_info = DISEASE_TARGETS[disease_key]
        pdb_id = disease_info["pdb_id"]

        print("Selected disease:", disease_key, flush=True)
        print("Selected drug:", drug_name, flush=True)

        # =====================================================
        # Prepare receptor + ligand
        # =====================================================

        receptor_file = prepare_protein_for_docking(pdb_id)
        ligand_file = download_ligand(drug_name)

        # =====================================================
        # Run Docking
        # =====================================================

        print("Starting Vina docking...", flush=True)

        docking_result = run_vina(
            receptor=receptor_file,
            ligand=ligand_file,
            active_sites=disease_info["active_sites"]
        )

        vina_best_affinity = docking_result["best_affinity"]
        vina_top_poses = docking_result["top_poses"]

        print("Docking completed. Best affinity:", vina_best_affinity, flush=True)

        # =====================================================
        # Build Quantum Complexes (NO REDUNDANCY)
        # =====================================================

        pose_complexes = []

        for pose in vina_top_poses:

            pose_file = pose["structure_file"]
            pocket_file = f"data/docking_results/pocket_{os.path.basename(pose_file)}"

            # Extract pocket
            extract_pocket(receptor_file, pose_file, pocket_file)

            parser = PDBParser(QUIET=True)

            # Protein atoms
            structure = parser.get_structure("pocket", pocket_file)
            all_protein_atoms = list(structure.get_atoms())

            # Ligand atoms
            ligand_structure = parser.get_structure("ligand", pose_file)
            ligand_coords = [atom.coord for atom in ligand_structure.get_atoms()]

            # Select closest atoms
            selected_atoms = select_closest_atoms(
                all_protein_atoms,
                ligand_coords,
                max_atoms=10
            )

            protein_atom_strings = []
            for atom in selected_atoms:
                element = atom.element.strip()
                if not element:
                    element = atom.get_name()[0].upper()

                x, y, z = atom.coord
                protein_atom_strings.append(
                    f"{element} {x:.6f} {y:.6f} {z:.6f}"
                )

            protein_atom_string = "; ".join(protein_atom_strings)

            # Ligand fragment
            ligand_atom_string = build_ligand_atom_string_from_pose(
                pose_file,
                pocket_file,
                max_atoms=8
            )

            # Complex
            complex_atom_string = protein_atom_string + "; " + ligand_atom_string

            # Store everything (IMPORTANT)
            pose_complexes.append({
                "pose_file": pose_file,
                "protein_atom_string": protein_atom_string,
                "ligand_atom_string": ligand_atom_string,
                "complex_atom_string": complex_atom_string
            })

        # =====================================================
        # Run VQE for Each Pose
        # =====================================================

        results = []

        for pose_data in pose_complexes:

            result = run_vqe_for_pose(
                (pose_data["complex_atom_string"], pose_data["pose_file"])
            )

            if result and "vqe_energy_hartree" in result:
                results.append({
                    "pose_data": pose_data,
                    "result": result
                })

        if not results:
            return {"error": "All VQE runs failed"}

        print("All VQE runs completed", flush=True)

        # =====================================================
        # Select Best Pose (Quantum-based)
        # =====================================================

        best_entry = min(
            results,
            key=lambda r: r["result"]["vqe_energy_hartree"]
        )

        best_data = best_entry["pose_data"]

        # Reuse stored fragments (NO REBUILDING)
        protein_atom_string = best_data["protein_atom_string"]
        ligand_atom_string = best_data["ligand_atom_string"]
        complex_atom_string = best_data["complex_atom_string"]

        best_pose_file = os.path.basename(best_data["pose_file"])

        # =====================================================
        # Final VQE Calculations
        # =====================================================

        protein_problem, _ = enforce_closed_shell(protein_atom_string, basis, 0)
        protein_result = compute_energies(protein_problem, label="PROTEIN")
        E_protein = protein_result["vqe_energy_hartree"]

        ligand_problem, _ = enforce_closed_shell(ligand_atom_string, basis, 0)
        ligand_result = compute_energies(ligand_problem, label="LIGAND")
        E_ligand = ligand_result["vqe_energy_hartree"]

        complex_problem, _ = enforce_closed_shell(complex_atom_string, basis, 0)
        complex_result = compute_energies(complex_problem, label="COMPLEX")
        E_complex = complex_result["vqe_energy_hartree"]

        # =====================================================
        # Binding Energy
        # =====================================================

        binding_energy = E_complex - (E_protein + E_ligand)
        binding_energy_kcal = binding_energy * HARTREE_TO_KCAL

        verdict = (
            "Favorable Interaction"
            if float(vina_best_affinity) < 0
            else "Weak Interaction"
        )
        pocket_filename = f"pocket_{best_pose_file}"

        return {
            "disease": disease_key,
            "drug": drug_name,
            "best_quantum_pose": best_pose_file,
            "pocket_file": pocket_filename,
            "protein_energy": float(E_protein),
            "ligand_energy": float(E_ligand),
            "complex_energy": float(E_complex),
            "binding_energy": float(binding_energy),
            "binding_energy_kcal": float(binding_energy_kcal),
            "verdict": verdict,
            "circuit_metrics": {
                "qubit_count": complex_result["num_qubits"],
                "active_orbitals": complex_result["num_active_orbitals"],
                "active_electrons": complex_result["num_active_electrons"],
                "circuit_depth": complex_result["circuit_structure"]["depth"],
                "gate_count": complex_result["circuit_structure"]["gate_count"]
            },
            "convergence_data": complex_result["convergence_history"]
        }

    except Exception as e:
        print("Binding computation failed:", str(e), flush=True)
        return {
            "error": f"Binding computation failed: {str(e)}"
        }