import numpy as np
from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL


def run_distance_sweep(
    atom_template,
    min_distance,
    max_distance,
    steps,
    basis="sto3g",
    charge=0,
    spin=0,
    active_electrons=None,
    active_orbitals=None,
    freeze_core=True,
):
    """
    Performs a potential energy surface (PES) scan
    over an interaction distance.

    atom_template must contain {distance} placeholder.
    """

    distances = np.linspace(min_distance, max_distance, steps)

    energies = []
    qubit_counts = []

    for d in distances:
        atom_string = atom_template.format(distance=d)

        problem = build_problem(
            atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            active_electrons=active_electrons,
            active_orbitals=active_orbitals,
            freeze_core=freeze_core,
        )

        result = compute_energies(problem)

        energies.append(result["vqe_energy"])
        qubit_counts.append(result["num_qubits"])

    energies = np.array(energies)

    # -------------------------
    # Equilibrium point
    # -------------------------
    min_index = np.argmin(energies)
    equilibrium_distance = distances[min_index]
    min_energy = energies[min_index]

    # -------------------------
    # Dissociation Energy
    # (energy difference between far distance and minimum)
    # -------------------------
    dissociation_energy = energies[-1] - min_energy
    dissociation_kcal = dissociation_energy * HARTREE_TO_KCAL

    # -------------------------
    # Curvature (stability indicator)
    # -------------------------
    curvature = None
    if 0 < min_index < len(energies) - 1:
        curvature = (
            energies[min_index - 1]
            - 2 * energies[min_index]
            + energies[min_index + 1]
        )

    return {
        "distances": distances.tolist(),
        "energies": energies.tolist(),

        "equilibrium_distance": float(equilibrium_distance),
        "min_energy": float(min_energy),

        "dissociation_energy_hartree": float(dissociation_energy),
        "dissociation_energy_kcal_mol": float(dissociation_kcal),

        "curvature_indicator": float(curvature) if curvature is not None else None,

        "qubit_counts": qubit_counts,
    }
