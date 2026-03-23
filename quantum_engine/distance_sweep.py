import numpy as np

from .molecule_builder import build_problem, auto_spin
from .vqe_solver import compute_energies, HARTREE_TO_KCAL
from .pdb_fragment_extractor import extract_residue_fragment
from .database import DISEASE_TARGETS


def enforce_closed_shell(atom_string, basis, charge, freeze_core,
                         active_electrons, active_orbitals):

    problem = build_problem(
        atom_string,
        basis=basis,
        charge=charge,
        spin=None,
    )

    num_alpha, num_beta = problem.num_particles
    total_electrons = num_alpha + num_beta

    if total_electrons % 2 != 0:
        charge += 1
        problem = build_problem(
            atom_string,
            basis=basis,
            charge=charge,
            spin=None,
        )

    return problem


def run_distance_sweep_from_selection(
    disease_name,
    drug_name,
    min_distance,
    max_distance,
    steps,
    basis="sto3g",
    active_electrons=None,
    active_orbitals=None,
    freeze_core=True,
):

    target_info = DISEASE_TARGETS[disease_name]

    target_atom_string, target_charge, _ = extract_residue_fragment(
        pdb_id=target_info["pdb_id"],
        chain_id=target_info["chain_id"],
        residue_name=target_info["residue_name"],
        residue_number=target_info["residue_number"],
        atom_names=target_info["atom_names"],
    )

    drug_atom_string = DRUGS[drug_name]["atom_string"]
    drug_charge = 0

    distances = np.linspace(min_distance, max_distance, steps)
    energies = []
    qubit_counts = []

    for d in distances:

        print(f"\n--- Distance: {d} Å ---", flush=True)

        shifted_drug_atoms = []
        for atom in drug_atom_string.split(";"):
            parts = atom.strip().split()
            element = parts[0]
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3]) + d
            shifted_drug_atoms.append(f"{element} {x} {y} {z}")

        shifted_drug_string = "; ".join(shifted_drug_atoms)

        complex_atom_string = target_atom_string + "; " + shifted_drug_string
        complex_charge = target_charge + drug_charge

        # Enforce closed shell
        problem = enforce_closed_shell(
            complex_atom_string,
            basis,
            complex_charge,
            freeze_core,
            active_electrons,
            active_orbitals,
        )

        result = compute_energies(problem, mode="production")

        if result.get("solver_error"):
            print("Solver error during sweep:", result["solver_error"], flush=True)
            energies.append(np.nan)
            qubit_counts.append(None)
            continue

        # Prefer exact energy if available
        energy = result.get("exact_energy_hartree")
        if energy is None:
            energy = result["vqe_energy_hartree"]

        energies.append(energy)
        qubit_counts.append(result["num_qubits"])

    energies = np.array(energies)

    min_index = np.nanargmin(energies)
    equilibrium_distance = distances[min_index]
    min_energy = energies[min_index]

    dissociation_energy = energies[-1] - min_energy
    dissociation_kcal = dissociation_energy * HARTREE_TO_KCAL

    curvature = None
    if 0 < min_index < len(energies) - 1:
        curvature = (
            energies[min_index - 1]
            - 2 * energies[min_index]
            + energies[min_index + 1]
        )

    return {
        "disease": disease_name,
        "drug": drug_name,
        "distances": distances.tolist(),
        "energies_hartree": energies.tolist(),
        "equilibrium_distance": float(equilibrium_distance),
        "min_energy_hartree": float(min_energy),
        "dissociation_energy_hartree": float(dissociation_energy),
        "dissociation_energy_kcal_mol": float(dissociation_kcal),
        "curvature_indicator": float(curvature) if curvature is not None else None,
        "qubit_counts": qubit_counts,
    }