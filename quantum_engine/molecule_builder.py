from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit

import time


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

    total_electrons = count_total_electrons(atom_string, charge)

    if total_electrons % 2 == 0:
        return 0   # singlet

    return 1       # doublet


# ---------------------------------------------------
# Closed-Shell Enforcement
# ---------------------------------------------------

def enforce_closed_shell_charge(atom_string, charge):

    total_electrons = count_total_electrons(atom_string, charge)

    if total_electrons % 2 != 0:
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
    # System Diagnostics
    # ----------------------------------

    atoms = [a.strip() for a in atom_string.split(";") if a.strip()]

    print("\n==============================", flush=True)
    print("Quantum System Construction", flush=True)
    print("==============================", flush=True)

    print("Basis:", basis, flush=True)
    print("Charge:", charge, flush=True)
    print("Spin:", spin, flush=True)

    print("Atom count:", len(atoms), flush=True)

    total_electrons = count_total_electrons(atom_string, charge)

    print("Total electrons:", total_electrons, flush=True)

    if len(atoms) > 40:
        print("⚠ WARNING: Large system detected!", flush=True)

    print("Starting PySCF driver...", flush=True)

    # ----------------------------------
    # PySCF Driver Execution
    # ----------------------------------

    try:

        t0 = time.time()

        driver = PySCFDriver(
            atom=atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            unit=DistanceUnit.ANGSTROM,
        )

        # Enable PySCF SCF iteration logs
        driver._pyscf_options = {"verbose": 4}

        problem = driver.run()

        t1 = time.time()

        print(f"PySCF finished in {t1 - t0:.3f} seconds", flush=True)

    except Exception as e:

        raise RuntimeError(f"PySCF driver failed: {str(e)}")

    # ----------------------------------
    # Final Electron Validation
    # ----------------------------------

    num_alpha, num_beta = problem.num_particles
    total_electrons = num_alpha + num_beta

    print("Final electrons:", total_electrons, flush=True)

    if total_electrons % 2 != 0:
        raise ValueError(
            "Electron count after PySCF is still odd. "
            "Check fragment construction."
        )

    print("Quantum chemistry problem built successfully.\n", flush=True)

    return problem