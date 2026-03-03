import numpy as np

from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import SLSQP
from qiskit.primitives import Estimator

from qiskit_nature import settings
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.mappers import QubitConverter
from qiskit_nature.second_q.transformers import (
    FreezeCoreTransformer,
    ActiveSpaceTransformer,
)

from qiskit.circuit.library import EfficientSU2


# ==================================================
# Global Scientific Settings (API Safe)
# ==================================================

settings.use_pauli_sum_op = False

HARTREE_TO_KCAL = 627.509
CHEMICAL_ACCURACY_KCAL = 1.0
CHEMICAL_ACCURACY_HARTREE = 0.001593

SAFE_QUBIT_LIMIT = 14
TARGET_ACTIVE_ORBITALS = 4


# ==================================================
# Active Electron Selection (Safe)
# ==================================================

def choose_active_electrons(total_electrons, max_active=4):

    active = min(max_active, total_electrons)
    active -= active % 2

    if active < 2 and total_electrons >= 2:
        active = 2

    if (total_electrons - active) % 2 != 0:
        active -= 2

    return max(0, active)


# ==================================================
# Main Solver
# ==================================================
def compute_energies(problem, mode="production"):

    print("\n==============================", flush=True)
    print(f"Starting compute_energies() | Mode: {mode}", flush=True)

    # ==================================================
    # MODE CONTROL
    # ==================================================

    if mode == "production":

        print("Applying FreezeCore...", flush=True)

        try:
            problem = FreezeCoreTransformer().transform(problem)
        except Exception:
            print("FreezeCore skipped.", flush=True)

        try:
            num_alpha, num_beta = problem.num_particles
            total_electrons = num_alpha + num_beta
            total_orbitals = problem.num_spatial_orbitals

            print(f"Electrons: {total_electrons}", flush=True)
            print(f"Orbitals: {total_orbitals}", flush=True)

            active_orbitals = min(TARGET_ACTIVE_ORBITALS, total_orbitals)
            active_electrons = choose_active_electrons(total_electrons, 4)

            print(f"Active electrons: {active_electrons}", flush=True)
            print(f"Active orbitals: {active_orbitals}", flush=True)

            problem = ActiveSpaceTransformer(
                num_electrons=active_electrons,
                num_spatial_orbitals=active_orbitals,
            ).transform(problem)

        except Exception as e:
            return {"solver_error": f"Active space reduction failed: {str(e)}"}

        multi_starts = 1

    elif mode == "benchmark":

        print("Full space mode — no FreezeCore, no ActiveSpace", flush=True)

        try:
            num_alpha, num_beta = problem.num_particles
            print(f"Electrons: {num_alpha + num_beta}", flush=True)
            print(f"Orbitals: {problem.num_spatial_orbitals}", flush=True)
        except Exception as e:
            return {"solver_error": f"Metadata extraction failed: {str(e)}"}

        multi_starts = 3  # stronger convergence for benchmark

    else:
        return {"solver_error": "Invalid mode. Use 'production' or 'benchmark'."}

    # ==================================================
    # Build Hamiltonian
    # ==================================================

    try:
        print("Building Hamiltonian...", flush=True)

        mapper = JordanWignerMapper()
        converter = QubitConverter(mapper)

        second_q_ops = problem.second_q_ops()
        main_op = second_q_ops[0]

        qubit_op = converter.convert(
            main_op,
            num_particles=problem.num_particles,
        )

        num_qubits = qubit_op.num_qubits
        print(f"Qubits required: {num_qubits}", flush=True)

        if num_qubits > SAFE_QUBIT_LIMIT:
            return {"solver_error": f"System requires {num_qubits} qubits."}

    except Exception as e:
        return {"solver_error": f"Hamiltonian construction failed: {str(e)}"}

    nuclear_repulsion = problem.nuclear_repulsion_energy
    print(f"Nuclear repulsion: {nuclear_repulsion}", flush=True)

    # ==================================================
    # Exact Solver (for validation)
    # ==================================================

    exact_energy = None

    if num_qubits <= 12:
        try:
            print("Running Exact solver...", flush=True)
            exact_solver = NumPyMinimumEigensolver()
            exact_result = exact_solver.compute_minimum_eigenvalue(qubit_op)
            eigenvalue = getattr(exact_result, "eigenvalue", exact_result)
            exact_energy = float(np.real(eigenvalue)) + nuclear_repulsion
            print(f"Exact energy: {exact_energy}", flush=True)
        except Exception:
            exact_energy = None

    # ==================================================
    # Ansatz (Unified)
    # ==================================================

    print("Constructing Ansatz...", flush=True)

    hf = HartreeFock(
        problem.num_spatial_orbitals,
        problem.num_particles,
        converter,
    )

    ansatz = UCCSD(
        num_spatial_orbitals=problem.num_spatial_orbitals,
        num_particles=problem.num_particles,
        qubit_converter=converter,
        initial_state=hf,
    )

    ansatz_name = "UCCSD"

    print(f"Ansatz used: {ansatz_name}", flush=True)

    estimator = Estimator(options={"shots": None})

    # ==================================================
    # Adaptive VQE
    # ==================================================

    best_energy = np.inf

    try:
        print("Starting VQE...", flush=True)

        # Adaptive iteration scaling
        param_count = ansatz.num_parameters
        adaptive_maxiter = max(300, 100 * param_count)

        print(f"Adaptive maxiter: {adaptive_maxiter}", flush=True)
        print(f"Multi-start runs: {multi_starts}", flush=True)

        for run in range(multi_starts):

            print(f"VQE Run {run+1}/{multi_starts}", flush=True)

            optimizer = SLSQP(maxiter=adaptive_maxiter)

            # small random initialization for stability
            initial_point = np.random.uniform(
                -0.05, 0.05, param_count
            )

            vqe = VQE(
                estimator=estimator,
                ansatz=ansatz,
                optimizer=optimizer,
                initial_point=initial_point,
            )

            result = vqe.compute_minimum_eigenvalue(qubit_op)
            eigenvalue = getattr(result, "eigenvalue", result)
            energy = float(np.real(eigenvalue)) + nuclear_repulsion

            print(f"Run energy: {energy}", flush=True)

            if energy < best_energy:
                best_energy = energy

        vqe_energy = best_energy
        print(f"Best VQE energy: {vqe_energy}", flush=True)

    except Exception as e:
        return {"solver_error": f"VQE failed: {str(e)}"}

    # ==================================================
    # Accuracy Check
    # ==================================================

    energy_error_hartree = None
    energy_error_kcal = None
    chemical_accuracy_met = None

    if exact_energy is not None:
        energy_error_hartree = abs(vqe_energy - exact_energy)
        energy_error_kcal = energy_error_hartree * HARTREE_TO_KCAL
        chemical_accuracy_met = energy_error_kcal <= CHEMICAL_ACCURACY_KCAL

        print(f"Energy error (kcal/mol): {energy_error_kcal}", flush=True)
        print(f"Chemical accuracy met: {chemical_accuracy_met}", flush=True)

    print("Finished compute_energies()", flush=True)
    print("==============================\n", flush=True)

    return {
        "mode": mode,
        "ansatz_used": ansatz_name,
        "num_qubits": num_qubits,
        "exact_energy_hartree": exact_energy,
        "vqe_energy_hartree": vqe_energy,
        "energy_error_hartree": energy_error_hartree,
        "energy_error_kcal_mol": energy_error_kcal,
        "chemical_accuracy_met": chemical_accuracy_met,
    }