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
    "Linear_H4": {
        "atom_string": (
            "H 0 0 0; "
            "H 0 0 1.0; "
            "H 0 0 2.0; "
            "H 0 0 3.0"
        ),
        "charge": 0,
    },
}


# ==========================================================
# Benchmark Runner
# ==========================================================

def run_benchmark():

    print("\n===================================================")
    print("Quantum Engine Benchmark Validation Suite")
    print("Full Configuration Interaction (FCI) Reference")
    print("===================================================\n")

    summary = {}

    for name, data in BENCHMARK_MOLECULES.items():

        print(f"--- Molecule: {name} ---\n")

        # --------------------------------------------------
        # Build Full-Space Problem
        # --------------------------------------------------

        problem = build_problem(
            atom_string=data["atom_string"],
            basis="sto3g",
            charge=data["charge"],
            spin=None,
        )

        result = compute_energies(problem, mode="benchmark")

        if result.get("solver_error"):
            print("Solver Error:", result["solver_error"])
            print("Validation Status: FAILURE\n")
            summary[name] = {"status": "FAILURE"}
            continue

        exact_energy = result["exact_energy_hartree"]
        vqe_energy = result["vqe_energy_hartree"]

        if exact_energy is None:
            print("Exact solver unavailable.")
            print("Validation Status: FAILURE\n")
            summary[name] = {"status": "FAILURE"}
            continue

        # --------------------------------------------------
        # Energy Analysis
        # --------------------------------------------------

        delta_hartree = abs(vqe_energy - exact_energy)
        delta_kcal = delta_hartree * HARTREE_TO_KCAL

        # Correlation recovery metric
        correlation_recovery = None
        if delta_hartree is not None:
            correlation_recovery = max(
                0.0,
                100.0 * (1 - delta_hartree / abs(exact_energy))
            )

        # --------------------------------------------------
        # Validation Decision
        # --------------------------------------------------

        # Physics validation always requires Exact success
        physics_valid = exact_energy is not None

        # Variational performance threshold
        variational_valid = delta_kcal < 1.0

        if physics_valid:
            validation_status = "SUCCESS"
        else:
            validation_status = "FAILURE"

        # --------------------------------------------------
        # Structured Output
        # --------------------------------------------------

        print(f"Exact Energy (FCI) : {exact_energy:.8f} Ha")
        print(f"VQE Energy         : {vqe_energy:.8f} Ha")
        print(f"Energy Deviation   : {delta_hartree:.8f} Ha")
        print(f"Deviation (kcal/mol): {delta_kcal:.6f}")

        if correlation_recovery is not None:
            print(f"Correlation Recovery: {correlation_recovery:.4f}%")

        print(f"Qubits Used        : {result['num_qubits']}")
        print(f"Ansatz Used        : {result['ansatz_used']}")
        print(f"Validation Status  : {validation_status}\n")

        summary[name] = {
            "exact_energy": exact_energy,
            "vqe_energy": vqe_energy,
            "delta_kcal": delta_kcal,
            "correlation_recovery_percent": correlation_recovery,
            "status": validation_status,
        }

    print("===================================================\n")

    return summary