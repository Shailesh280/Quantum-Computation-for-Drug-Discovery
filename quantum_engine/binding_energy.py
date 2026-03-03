from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL
from .pdb_fragment_extractor import extract_residue_fragment
from .database import DISEASE_TARGETS, DRUGS


MAX_DRUG_HEAVY_ATOMS = 4


# ==========================================================
# Utility: Reduce Drug Fragment
# ==========================================================
def reduce_drug_fragment(atom_string):
    heavy_atoms = []

    for atom in atom_string.split(";"):
        parts = atom.strip().split()
        element = parts[0]

        if element == "H":
            continue

        heavy_atoms.append(atom.strip())

        if len(heavy_atoms) >= MAX_DRUG_HEAVY_ATOMS:
            break

    return "; ".join(heavy_atoms)


# ==========================================================
# Utility: Ensure Even Electron Count
# ==========================================================
def enforce_closed_shell(atom_string, basis, charge):
    """
    Adjust charge to guarantee even total electrons.
    """

    problem = build_problem(
        atom_string,
        basis=basis,
        charge=charge,
    )

    num_alpha, num_beta = problem.num_particles
    total_electrons = num_alpha + num_beta

    if total_electrons % 2 != 0:
        adjusted_charge = charge + 1

        problem = build_problem(
            atom_string,
            basis=basis,
            charge=adjusted_charge,
        )

        return problem, adjusted_charge

    return problem, charge


# ==========================================================
# Main Binding Function
# ==========================================================
def compute_binding_from_selection(
    disease_name,
    drug_name,
    separation_distance=3.0,
    basis="sto3g",
):

    try:

        # ------------------------------------------------
        # Target Fragment
        # ------------------------------------------------
        target_info = DISEASE_TARGETS[disease_name]

        target_atom_string, target_charge, _ = extract_residue_fragment(
            pdb_id=target_info["pdb_id"],
            chain_id=target_info["chain_id"],
            residue_name=target_info["residue_name"],
            residue_number=target_info["residue_number"],
            atom_names=target_info["atom_names"],
        )

        # ------------------------------------------------
        # Drug Fragment
        # ------------------------------------------------
        full_drug_atom_string = DRUGS[drug_name]["atom_string"]
        drug_atom_string = reduce_drug_fragment(full_drug_atom_string)
        drug_charge = 0

        # ------------------------------------------------
        # Shift Drug
        # ------------------------------------------------
        shifted_drug_atoms = []
        for atom in drug_atom_string.split(";"):
            parts = atom.strip().split()
            element = parts[0]
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3]) + separation_distance
            shifted_drug_atoms.append(f"{element} {x} {y} {z}")

        shifted_drug_string = "; ".join(shifted_drug_atoms)

        # ------------------------------------------------
        # Complex Construction
        # ------------------------------------------------
        complex_atom_string = target_atom_string + "; " + shifted_drug_string
        complex_charge = target_charge + drug_charge

        # ------------------------------------------------
        # Build Problems
        # ------------------------------------------------
        drug_problem, drug_charge = enforce_closed_shell(
            drug_atom_string,
            basis,
            drug_charge,
        )

        target_problem, target_charge = enforce_closed_shell(
            target_atom_string,
            basis,
            target_charge,
        )

        complex_problem, complex_charge = enforce_closed_shell(
            complex_atom_string,
            basis,
            complex_charge,
        )

        # ------------------------------------------------
        # Compute Energies
        # ------------------------------------------------
        drug_result = compute_energies(drug_problem)
        target_result = compute_energies(target_problem)
        complex_result = compute_energies(complex_problem)

        # -------------------------
        # Solver Error Check
        # -------------------------
        for label, result in [
            ("drug", drug_result),
            ("target", target_result),
            ("complex", complex_result),
        ]:
            if result.get("solver_error"):
                return {
                    "error": f"{label.capitalize()} calculation failed",
                    "details": result["solver_error"],
                }

        # ------------------------------------------------
        # Chemical Accuracy Gate
        # ------------------------------------------------
        accuracy_flags = [
            drug_result["chemical_accuracy_met"],
            target_result["chemical_accuracy_met"],
            complex_result["chemical_accuracy_met"],
        ]

        if not all(flag is True for flag in accuracy_flags):
            return {
                "error": "Chemical accuracy not achieved for all systems.",
                "accuracy_details": {
                    "drug_error_kcal": drug_result["energy_error_kcal_mol"],
                    "target_error_kcal": target_result["energy_error_kcal_mol"],
                    "complex_error_kcal": complex_result["energy_error_kcal_mol"],
                }
            }

        # ------------------------------------------------
        # Extract Energies (Hartree)
        # ------------------------------------------------
        drug_energy = drug_result["vqe_energy_hartree"]
        target_energy = target_result["vqe_energy_hartree"]
        complex_energy = complex_result["vqe_energy_hartree"]

        # ------------------------------------------------
        # Binding Energy
        # ------------------------------------------------
        binding = complex_energy - (drug_energy + target_energy)
        binding_kcal = binding * HARTREE_TO_KCAL

        if binding_kcal < 0:
            verdict = "Favorable Interaction"
        else:
            verdict = "Unfavorable Interaction"

        return {
            "disease": disease_name,
            "drug": drug_name,
            "binding_hartree": float(binding),
            "binding_kcal_mol": float(binding_kcal),
            "verdict": verdict,
            "drug_qubits": drug_result["num_qubits"],
            "target_qubits": target_result["num_qubits"],
            "complex_qubits": complex_result["num_qubits"],
            "drug_ansatz": drug_result["ansatz_used"],
            "target_ansatz": target_result["ansatz_used"],
            "complex_ansatz": complex_result["ansatz_used"],
            "drug_error_kcal": drug_result["energy_error_kcal_mol"],
            "target_error_kcal": target_result["energy_error_kcal_mol"],
            "complex_error_kcal": complex_result["energy_error_kcal_mol"],
        }

    except Exception as e:
        return {"error": f"Binding computation failed: {str(e)}"}