import os

from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL
from .pdb_fragment_extractor import extract_residue_fragment
from .pose_fragment_builder import build_ligand_atom_string_from_pose
from .database import DISEASE_TARGETS
from .vina_engine import run_vina
from .pose_parameter_initializer import generate_pose_initial_params


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
# Main Binding Function
# ==========================================================

def compute_binding_from_selection(
    disease_key,
    active_site_index,
    drug_name,
    basis="sto3g",
    separation_distance=None
):

    try:

        # =====================================================
        # Disease + Active Site
        # =====================================================

        disease_info = DISEASE_TARGETS[disease_key]
        active_site = disease_info["active_sites"][active_site_index]

        pdb_id = disease_info["pdb_id"]

        # =====================================================
        # Extract Target Fragment
        # =====================================================

        target_atom_string, target_charge, _ = extract_residue_fragment(

            pdb_id=pdb_id,
            chain_id=active_site["chain_id"],
            residue_name=active_site["residue_name"],
            residue_number=active_site["residue_number"],
            atom_names=active_site["atom_names"],
        )

        target_problem, target_charge = enforce_closed_shell(
            target_atom_string,
            basis,
            target_charge
        )

        symbols = target_problem.molecule.symbols
        coords = target_problem.molecule.coords

        target_atom_string = "; ".join(
            f"{sym} {x} {y} {z}"
            for sym, (x, y, z) in zip(symbols, coords)
        )

        print("Target fragment built", flush=True)

        # =====================================================
        # Run Vina Docking
        # =====================================================

        print("Starting Vina docking...", flush=True)

        docking_result = run_vina(
            receptor="pdb_files/protein.pdbqt",
            ligand="pdb_files/ligand.pdbqt"
        )

        vina_best_affinity = docking_result["best_affinity"]
        vina_top_poses = docking_result["top_poses"]

        print("Docking completed. Best affinity:", vina_best_affinity, flush=True)

        # =====================================================
        # Build Quantum Complexes
        # =====================================================

        pose_complexes = []

        for pose in vina_top_poses:

            pose_file = pose["structure_file"]

            ligand_atom_string = build_ligand_atom_string_from_pose(pose_file)

            complex_atom_string = target_atom_string + "; " + ligand_atom_string

            pose_complexes.append((complex_atom_string, pose_file))

        print("Quantum complexes built:", len(pose_complexes), flush=True)

        # =====================================================
        # Sequential VQE Execution
        # =====================================================

        results = []

        for complex_atom_string, pose_file in pose_complexes:

            result = run_vqe_for_pose((complex_atom_string, pose_file))

            if result:
                results.append({
                    "pose_file": pose_file,
                    "result": result
                })

        if not results:
            return {"error": "All VQE runs failed"}

        print("All VQE runs completed", flush=True)

        # =====================================================
        # Select Best Pose (Lowest Quantum Energy)
        # =====================================================

        best_entry = min(
            results,
            key=lambda r: r["result"]["vqe_energy_hartree"]
        )

        best_result = best_entry["result"]
        best_pose_file = best_entry["pose_file"]

        complex_energy = best_result["vqe_energy_hartree"]
        complex_energy_kcal = complex_energy * HARTREE_TO_KCAL

        convergence_history = best_result["convergence_history"]

        # =====================================================
        # Interaction Verdict
        # =====================================================

        verdict = (
            "Favorable Interaction"
            if vina_best_affinity < -2
            else "Weak Interaction"
        )

        # =====================================================
        # Clean Pose Filename
        # =====================================================

        best_pose_file = os.path.basename(best_pose_file)

        # =====================================================
        # Response (Quantum-focused)
        # =====================================================

        return {

            "disease": disease_key,
            "drug": drug_name,

            # Selected pose
            "best_quantum_pose": best_pose_file,

            # Quantum energy
            "complex_energy": float(complex_energy),
            "complex_energy_kcal": float(complex_energy_kcal),

            # Interaction verdict
            "verdict": verdict,

            # Circuit metrics
            "circuit_metrics": {
                "active_orbitals": 4,
                "active_electrons": 4,
                "qubit_count": best_result["num_qubits"],
                "circuit_depth": best_result["circuit_structure"]["depth"],
                "gate_count": best_result["circuit_structure"]["gate_count"],
            },

            # VQE convergence
            "convergence_data": convergence_history
        }

    except Exception as e:

        print("Binding computation failed:", str(e), flush=True)

        return {
            "error": f"Binding computation failed: {str(e)}"
        }