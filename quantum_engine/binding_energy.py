from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL


def compute_binding_energy(drug_atom, active_atom, complex_atom):

    drug_problem = build_problem(drug_atom)
    active_problem = build_problem(active_atom)
    complex_problem = build_problem(complex_atom)

    drug_energy = compute_energies(drug_problem)["vqe_energy"]
    active_energy = compute_energies(active_problem)["vqe_energy"]
    complex_energy = compute_energies(complex_problem)["vqe_energy"]

    binding = complex_energy - (drug_energy + active_energy)

    return {
        "drug_energy": drug_energy,
        "active_energy": active_energy,
        "complex_energy": complex_energy,
        "binding_hartree": float(binding),
        "binding_kcal_mol": float(binding * HARTREE_TO_KCAL)
    }
