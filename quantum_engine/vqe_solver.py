import numpy as np
import time

from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import L_BFGS_B
from qiskit.primitives import Estimator

from qiskit_nature import settings
from qiskit_nature.second_q.mappers import JordanWignerMapper, QubitConverter
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.transformers import (
    FreezeCoreTransformer,
    ActiveSpaceTransformer,
)
import numpy as np

settings.use_pauli_sum_op = False


HARTREE_TO_KCAL = 627.509
CHEMICAL_ACCURACY_KCAL = 1.0

SAFE_QUBIT_LIMIT = 14
TARGET_ACTIVE_ORBITALS = 4


def debug(msg):
    print(f"[DEBUG] {msg}", flush=True)


def choose_active_electrons(total_electrons):

    active = min(6, total_electrons)

    if active % 2 != 0:
        active -= 1

    if active < 2:
        active = 2

    return active


def compute_energies(problem, label="SYSTEM", benchmark=False, pose_initial=None):

    debug("==============================")
    debug(f"{label} Starting compute_energies")

    total_start = time.time()

    # -----------------------------
    # Freeze core
    # -----------------------------
    if not benchmark:
        try:
            problem = FreezeCoreTransformer().transform(problem)
            debug("Freeze core applied")
        except Exception as e:
            debug(f"Freeze core skipped: {e}")
    else:
        debug("Benchmark mode → skipping freeze core")

    num_alpha, num_beta = problem.num_particles
    total_electrons = num_alpha + num_beta
    total_orbitals = problem.num_spatial_orbitals

    # -----------------------------
    # Active space selection
    # -----------------------------

    if benchmark:

        active_orbitals = total_orbitals
        active_electrons = total_electrons

    else:

        active_orbitals = min(TARGET_ACTIVE_ORBITALS, total_orbitals)
        active_electrons = choose_active_electrons(total_electrons)

    debug(
        f"Active space → electrons={active_electrons} orbitals={active_orbitals}"
    )

    if not benchmark:
        problem = ActiveSpaceTransformer(
            num_electrons=active_electrons,
            num_spatial_orbitals=active_orbitals,
        ).transform(problem)
    else:
        debug("Active space transformation skipped")

    # -----------------------------
    # Hamiltonian construction
    # -----------------------------

    mapper = JordanWignerMapper()
    converter = QubitConverter(mapper)

    fermionic_op = problem.second_q_ops()[0]

    qubit_op = converter.convert(
        fermionic_op,
        num_particles=problem.num_particles
    )

    num_qubits = qubit_op.num_qubits

    debug(f"Hamiltonian built | qubits = {num_qubits}")

    if num_qubits > SAFE_QUBIT_LIMIT:

        return {"solver_error": f"{num_qubits} qubits required"}

    # -----------------------------
    # Exact solver (reference)
    # -----------------------------

    exact_energy = None

    if num_qubits <= 12:

        try:

            debug("Running exact solver")

            exact_solver = NumPyMinimumEigensolver()

            exact_result = exact_solver.compute_minimum_eigenvalue(
                qubit_op
            )
            if benchmark:
                exact_energy = exact_result.eigenvalue.real + problem.nuclear_repulsion_energy
            else:
                exact_energy = exact_result.eigenvalue.real
            debug(f"Exact energy = {exact_energy}")

        except Exception as e:

            debug(f"Exact solver failed: {e}")

    # -----------------------------
    # Hartree–Fock reference state
    # -----------------------------

    hf = HartreeFock(
        problem.num_spatial_orbitals,
        problem.num_particles,
        converter
    )

    # -----------------------------
    # UCCSD Ansatz
    # -----------------------------

    ansatz = UCCSD(
        num_spatial_orbitals=problem.num_spatial_orbitals,
        num_particles=problem.num_particles,
        qubit_converter=converter,
        initial_state=hf
    )

    estimator = Estimator()

    param_count = ansatz.num_parameters

    debug(f"Ansatz parameters = {param_count}")

    # -----------------------------
    # Optimizer configuration
    # -----------------------------

    if benchmark:

        optimizer = L_BFGS_B(maxiter=1000, ftol=1e-9)

    else:

        optimizer = L_BFGS_B(maxiter=400, ftol=1e-6)

    # -----------------------------
    # Parameter initialization
    # -----------------------------

    initial_point = np.zeros(ansatz.num_parameters)

    # -----------------------------
    # Convergence tracking
    # -----------------------------

    convergence = []

    def callback(eval_count, params, energy, std):

        energy_val = energy.real

        convergence.append(
            {
                "iteration": int(eval_count),
                "energy": energy_val
            }
        )

    # -----------------------------
    # Run VQE
    # -----------------------------

    debug("Starting VQE")

    start = time.time()

    vqe = VQE(
        estimator=estimator,
        ansatz=ansatz,
        optimizer=optimizer,
        initial_point=initial_point,
        callback=callback
    )

    result = vqe.compute_minimum_eigenvalue(qubit_op)
    if benchmark:
        vqe_energy = result.eigenvalue.real + problem.nuclear_repulsion_energy
    else:
        vqe_energy = result.eigenvalue.real 
    debug(f"VQE finished in {time.time() - start:.2f}s")

    # -----------------------------
    # Accuracy metrics
    # -----------------------------

    energy_error_hartree = None
    energy_error_kcal = None
    chemical_accuracy_met = None

    if exact_energy is not None:

        energy_error_hartree = abs(vqe_energy - exact_energy)

        energy_error_kcal = energy_error_hartree * HARTREE_TO_KCAL

        chemical_accuracy_met = (
            energy_error_kcal < CHEMICAL_ACCURACY_KCAL
        )

    # -----------------------------
    # Circuit metadata
    # -----------------------------

    circuit = ansatz.decompose().decompose()

    circuit_metadata = {

        "qubits": int(num_qubits),
        "depth": int(circuit.depth()),
        "gate_count": int(circuit.size()),
        "parameters": int(ansatz.num_parameters),
    }

    debug(f"{label} energy = {vqe_energy}")
    debug(f"Total runtime {time.time() - total_start:.2f}s")

    # -----------------------------
    # Return results
    # -----------------------------

    return {

        "ansatz_used": "UCCSD",

        "num_qubits": num_qubits,

        "exact_energy_hartree": exact_energy,

        "vqe_energy_hartree": vqe_energy,

        "energy_error_hartree": energy_error_hartree,

        "energy_error_kcal_mol": energy_error_kcal,

        "chemical_accuracy_met": chemical_accuracy_met,

        "circuit_structure": circuit_metadata,

        "convergence_history": convergence

    }