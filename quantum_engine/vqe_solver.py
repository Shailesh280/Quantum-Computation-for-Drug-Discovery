import numpy as np

from qiskit_algorithms import VQE, NumPyMinimumEigensolver
from qiskit_algorithms.optimizers import SLSQP
from qiskit.primitives import Estimator

from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.mappers import QubitConverter
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

from qiskit.circuit.library import EfficientSU2


# -----------------------------
# Constants
# -----------------------------
HARTREE_TO_KCAL = 627.5
CHEMICAL_ACCURACY_HARTREE = 0.0016

EXACT_THRESHOLD = 10
UCCSD_THRESHOLD = 12


def compute_energies(problem):

    # -----------------------------
    # Active Space Reduction (Small Systems)
    # -----------------------------
    if problem.num_spatial_orbitals > 4:
        try:
            print("Applying active space reduction...")
            transformer = ActiveSpaceTransformer(
                num_electrons=2,
                num_spatial_orbitals=2
            )
            problem = transformer.transform(problem)
        except Exception as e:
            print("Active space reduction skipped:", e)

    # -----------------------------
    # Qubit Mapping
    # -----------------------------
    mapper = JordanWignerMapper()
    converter = QubitConverter(mapper)

    second_q_ops = problem.second_q_ops()
    main_op = second_q_ops[0]
    qubit_op = converter.convert(main_op, num_particles=problem.num_particles)

    num_qubits = qubit_op.num_qubits
    num_particles = problem.num_particles
    num_spatial_orbitals = problem.num_spatial_orbitals

    print(f"\nDetected {num_qubits} qubits after reduction.")

    # -----------------------------
    # Exact Solver (Safe Only)
    # -----------------------------
    exact_energy = None

    if num_qubits <= EXACT_THRESHOLD:
        try:
            print("Running exact solver...")
            exact_solver = NumPyMinimumEigensolver()
            exact_result = exact_solver.compute_minimum_eigenvalue(qubit_op)
            exact_energy = float(np.real(exact_result.eigenvalue))
            print("Exact energy:", exact_energy)
        except Exception as e:
            print("Exact solver failed:", str(e))
            exact_energy = None
    else:
        print("Skipping exact solver (qubit count too high).")

    # -----------------------------
    # Hartree-Fock Initial State
    # -----------------------------
    hf = HartreeFock(
        num_spatial_orbitals,
        num_particles,
        converter
    )

    # -----------------------------
    # Ansatz Selection
    # -----------------------------
    if num_qubits <= UCCSD_THRESHOLD:
        print("Using UCCSD ansatz.")
        ansatz = UCCSD(
            num_spatial_orbitals=num_spatial_orbitals,
            num_particles=num_particles,
            qubit_converter=converter,
            initial_state=hf,
        )
        ansatz_name = "UCCSD"
    else:
        print("Using EfficientSU2 ansatz.")
        ansatz = EfficientSU2(
            num_qubits=num_qubits,
            reps=3,
            entanglement="full"
        )
        ansatz_name = "EfficientSU2"

    print(f"Ansatz parameters: {ansatz.num_parameters}")
    print(f"Circuit depth: {ansatz.decompose().depth()}")

    # -----------------------------
    # VQE
    # -----------------------------
    convergence_history = []

    def callback(eval_count, parameters, mean, std):
        convergence_history.append(float(mean))

    optimizer = SLSQP(maxiter=200)
    estimator = Estimator()

    vqe = VQE(
        estimator,
        ansatz,
        optimizer,
        callback=callback
    )

    print("Starting VQE optimization...")
    vqe_result = vqe.compute_minimum_eigenvalue(qubit_op)
    vqe_energy = float(np.real(vqe_result.eigenvalue))
    print("VQE energy:", vqe_energy)

    # -----------------------------
    # Metrics
    # -----------------------------
    error = None
    error_kcal = None
    within_chemical_accuracy = None

    if exact_energy is not None:
        error = abs(vqe_energy - exact_energy)
        error_kcal = error * HARTREE_TO_KCAL
        within_chemical_accuracy = error < CHEMICAL_ACCURACY_HARTREE

    return {
        "ansatz_used": ansatz_name,
        "num_qubits": num_qubits,
        "num_parameters": ansatz.num_parameters,
        "circuit_depth": ansatz.decompose().depth(),

        "exact_energy": exact_energy,
        "vqe_energy": vqe_energy,
        "error": error,
        "error_kcal_mol": error_kcal,
        "within_chemical_accuracy": within_chemical_accuracy,

        "convergence_history": convergence_history
    }
