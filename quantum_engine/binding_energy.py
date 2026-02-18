from .molecule_builder import build_problem
from .vqe_solver import compute_energies, HARTREE_TO_KCAL


def compute_binding_energy(
    drug_atom,
    active_atom,
    complex_atom,
    basis="sto3g",
    charge=0,
    spin=0,
    active_electrons=None,
    active_orbitals=None,
    freeze_core=True,
):
    """
    Computes binding energy:

        E_binding = E_complex − (E_drug + E_active)

    All fragments use identical basis, active space,
    and freeze-core settings to ensure consistency.
    """

    # -------------------------
    # Build problems (consistent settings)
    # -------------------------
    drug_problem = build_problem(
        drug_atom,
        basis=basis,
        charge=charge,
        spin=spin,
        active_electrons=active_electrons,
        active_orbitals=active_orbitals,
        freeze_core=freeze_core,
    )

    active_problem = build_problem(
        active_atom,
        basis=basis,
        charge=charge,
        spin=spin,
        active_electrons=active_electrons,
        active_orbitals=active_orbitals,
        freeze_core=freeze_core,
    )

    complex_problem = build_problem(
        complex_atom,
        basis=basis,
        charge=charge,
        spin=spin,
        active_electrons=active_electrons,
        active_orbitals=active_orbitals,
        freeze_core=freeze_core,
    )

    # -------------------------
    # Compute energies
    # -------------------------
    drug_result = compute_energies(drug_problem)
    active_result = compute_energies(active_problem)
    complex_result = compute_energies(complex_problem)

    drug_energy = drug_result["vqe_energy"]
    active_energy = active_result["vqe_energy"]
    complex_energy = complex_result["vqe_energy"]

    # -------------------------
    # Binding Energy
    # -------------------------
    binding = complex_energy - (drug_energy + active_energy)
    binding_kcal = binding * HARTREE_TO_KCAL

    # -------------------------
    # Interpretation
    # -------------------------
    interaction_type = None
    if binding < 0:
        interaction_type = "Favorable (Exothermic Binding)"
    else:
        interaction_type = "Unfavorable (Endothermic)"

    # Strong vs Weak interaction estimate
    strength = None
    if binding_kcal < -20:
        strength = "Strong interaction"
    elif binding_kcal < -5:
        strength = "Moderate interaction"
    elif binding_kcal < 0:
        strength = "Weak interaction"
    else:
        strength = "No stable binding"

    return {
        "drug_energy": drug_energy,
        "active_energy": active_energy,
        "complex_energy": complex_energy,

        "binding_hartree": float(binding),
        "binding_kcal_mol": float(binding_kcal),

        "interaction_type": interaction_type,
        "interaction_strength": strength,

        "drug_qubits": drug_result["num_qubits"],
        "active_qubits": active_result["num_qubits"],
        "complex_qubits": complex_result["num_qubits"],

        "drug_ansatz": drug_result["ansatz_used"],
        "active_ansatz": active_result["ansatz_used"],
        "complex_ansatz": complex_result["ansatz_used"],
    }
