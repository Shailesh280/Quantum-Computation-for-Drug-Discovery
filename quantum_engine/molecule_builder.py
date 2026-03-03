from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit


# ---------------------------------------------------
# Atomic Numbers
# ---------------------------------------------------

ATOMIC_NUMBERS = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "Br": 35,
}


# ---------------------------------------------------
# Electron Counting Utility
# ---------------------------------------------------

def count_total_electrons(atom_string, charge=0):

    total_electrons = 0

    for atom in atom_string.split(";"):
        element = atom.strip().split()[0]

        if element not in ATOMIC_NUMBERS:
            raise ValueError(f"Unknown element: {element}")

        total_electrons += ATOMIC_NUMBERS[element]

    total_electrons -= charge
    return total_electrons


# ---------------------------------------------------
# Safe Spin Determination
# ---------------------------------------------------

def auto_spin(atom_string, charge=0):
    """
    Closed-shell if even electrons,
    doublet if odd.
    """

    total_electrons = count_total_electrons(atom_string, charge)

    if total_electrons % 2 == 0:
        return 0  # Singlet

    return 1  # Doublet


# ---------------------------------------------------
# Closed-Shell Enforcement
# ---------------------------------------------------

def enforce_closed_shell_charge(atom_string, charge):

    total_electrons = count_total_electrons(atom_string, charge)

    if total_electrons % 2 != 0:
        # Prefer removing one electron
        charge += 1

    return charge


# ---------------------------------------------------
# High-Accuracy Problem Builder
# ---------------------------------------------------

def build_problem(
    atom_string,
    basis="sto3g",
    charge=0,
    spin=None,
    enforce_closed_shell=True,
):
    """
    Builds full ElectronicStructureProblem.

    Design Principles:
    - No FreezeCore here
    - No ActiveSpace here
    - Full Hamiltonian preserved
    - Charge corrected if needed
    - Spin physically consistent
    """

    # ----------------------------------
    # Enforce Even Electron Count
    # ----------------------------------
    if enforce_closed_shell:
        charge = enforce_closed_shell_charge(atom_string, charge)

    # ----------------------------------
    # Determine Spin
    # ----------------------------------
    if spin is None:
        spin = auto_spin(atom_string, charge)

    # ----------------------------------
    # PySCF Driver Execution
    # ----------------------------------
    try:
        driver = PySCFDriver(
            atom=atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            unit=DistanceUnit.ANGSTROM,
        )

        problem = driver.run()

    except Exception as e:
        raise RuntimeError(f"PySCF driver failed: {str(e)}")

    # ----------------------------------
    # Final Electron Validation
    # ----------------------------------
    num_alpha, num_beta = problem.num_particles
    total_electrons = num_alpha + num_beta

    if total_electrons % 2 != 0:
        raise ValueError(
            "Electron count after PySCF is still odd. "
            "Check fragment construction."
        )

    return problem