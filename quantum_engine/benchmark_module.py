import numpy as np

from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL


# ==========================================================
# Benchmark Molecules (Full-Space Validation)
# ==========================================================

BENCHMARK_MOLECULES = {
    "H2": {
        "atom_string": "H 0 0 0; H 0 0 0.74",
        "charge": 0,
    },
    "HeH+": {
        "atom_string": "He 0 0 0; H 0 0 0.77",
        "charge": 1,
    },
    "H3+": {
        "atom_string": (
            "H 0.0 0.0 0.0; "
            "H 0.9 0.0 0.0; "
            "H 0.45 0.779 0.0"
        ),
        "charge": 1,
    },
    "H4": {
        "atom_string": (
            "H 0 0 0; "
            "H 0 0 1.0; "
            "H 0 0 2.0; "
            "H 0 0 3.0"
        ),
        "charge": 0,
    },
}

def to_python_type(obj):

    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {k: to_python_type(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_python_type(v) for v in obj]

    return obj
# ==========================================================
# Benchmark Runner
# ==========================================================

def run_benchmark(selected_molecule="H2"):

    print("\n===================================================")
    print("Quantum Engine Benchmark Validation Suite")
    print("Full Configuration Interaction (FCI) Reference")
    print("===================================================\n")

    summary = {}

    if selected_molecule not in BENCHMARK_MOLECULES:
        return {
            selected_molecule: {
                "status": "FAILURE",
                "error": "Invalid molecule selection"
            }
        }

    molecules_to_run = {
        selected_molecule: BENCHMARK_MOLECULES[selected_molecule]
    }


    for name, data in molecules_to_run.items():

        print(f"--- Molecule: {name} ---\n")

        # --------------------------------------------------
        # Build Full-Space Problem
        # --------------------------------------------------

        try:
            problem = build_problem(
                atom_string=data["atom_string"],
                basis="sto3g",
                charge=data["charge"],
                spin=None,
            )
        except Exception as e:
            summary[name] = {
                "status": "FAILURE",
                "error": f"Problem build failed: {str(e)}"
            }
            continue

        # --------------------------------------------------
        # Run Benchmark Mode (Full Space)
        # --------------------------------------------------

        result = compute_energies(problem, mode="benchmark")

        if result.get("solver_error"):
            summary[name] = {
                "status": "FAILURE",
                "error": result["solver_error"]
            }
            continue

        exact_energy = result["exact_energy_hartree"]
        vqe_energy = result["vqe_energy_hartree"]

        if exact_energy is None or vqe_energy is None:
            summary[name] = {
                "status": "FAILURE",
                "error": "Exact or VQE energy unavailable"
            }
            continue

        # --------------------------------------------------
        # Energy Analysis
        # --------------------------------------------------

        delta_hartree = abs(vqe_energy - exact_energy)
        delta_kcal = delta_hartree * HARTREE_TO_KCAL

        # Correlation recovery relative to exact correlation scale
        correlation_recovery = max(
            0.0,
            100.0 * (1.0 - delta_hartree / max(abs(exact_energy), 1e-8))
        )

        # --------------------------------------------------
        # Validation Logic
        # --------------------------------------------------

        # Variational principle check (VQE should not go below exact)
        variational_valid = vqe_energy >= exact_energy - 1e-6

        # Chemical accuracy threshold (optional scientific validation)
        chemical_valid = delta_kcal < 1.0

        if variational_valid:
            validation_status = "SUCCESS"
        else:
            validation_status = "FAILURE"

        # --------------------------------------------------
        # Console Output (Professional)
        # --------------------------------------------------

        print(f"Exact Energy (FCI)  : {exact_energy:.8f} Ha")
        print(f"VQE Energy          : {vqe_energy:.8f} Ha")
        print(f"Energy Deviation    : {delta_hartree:.8f} Ha")
        print(f"Deviation (kcal/mol): {delta_kcal:.6f}")
        print(f"Correlation Recovery: {correlation_recovery:.4f}%")
        print(f"Qubits Used         : {result['num_qubits']}")
        print(f"Ansatz Used         : {result['ansatz_used']}")
        print(f"Validation Status   : {validation_status}\n")

        # --------------------------------------------------
        # Structured JSON for Frontend
        # --------------------------------------------------

        summary[name] = {
            "exact_energy_hartree": exact_energy,
            "vqe_energy_hartree": vqe_energy,
            "delta_hartree": delta_hartree,
            "delta_kcal_mol": delta_kcal,
            "correlation_recovery_percent": correlation_recovery,
            "variational_principle_satisfied": variational_valid,
            "chemical_accuracy_met": chemical_valid,
            "num_qubits": result["num_qubits"],
            "ansatz_used": result["ansatz_used"],
            "multi_start_runs": result.get("multi_start_runs"),
            "adaptive_maxiter": result.get("adaptive_maxiter"),
            "convergence_history": result.get("convergence_history"),
            "circuit_structure": result.get("circuit_structure"),
            "status": validation_status,
        }

    print("===================================================\n")

    return to_python_type(summary)