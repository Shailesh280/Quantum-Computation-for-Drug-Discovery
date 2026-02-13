import numpy as np
from .molecule_builder import build_problem
from .vqe_solver import compute_energies


def run_distance_sweep(
    atom_template,
    min_distance,
    max_distance,
    steps
):

    distances = np.linspace(min_distance, max_distance, steps)
    energies = []

    for d in distances:
        atom_string = atom_template.format(distance=d)
        problem = build_problem(atom_string)
        energy = compute_energies(problem)["vqe_energy"]
        energies.append(energy)

    min_energy = min(energies)
    equilibrium_distance = distances[energies.index(min_energy)]

    return {
        "distances": distances.tolist(),
        "energies": energies,
        "equilibrium_distance": float(equilibrium_distance),
        "min_energy": float(min_energy)
    }
