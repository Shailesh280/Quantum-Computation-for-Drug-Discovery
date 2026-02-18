from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer, FreezeCoreTransformer
from qiskit_nature.units import DistanceUnit


def build_problem(
    atom_string,
    basis="sto3g",
    charge=0,
    spin=0,
    active_electrons=None,
    active_orbitals=None,
    freeze_core=True
):
    """
    Build ElectronicStructureProblem with optional
    freeze-core and active-space reduction.
    """

    driver = PySCFDriver(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin,
        unit=DistanceUnit.ANGSTROM,
    )

    problem = driver.run()

    # -------------------------
    # Freeze Core
    # -------------------------
    if freeze_core:
        transformer = FreezeCoreTransformer()
        problem = transformer.transform(problem)

    # -------------------------
    # Active Space Reduction (Optional)
    # -------------------------
    if active_electrons is not None and active_orbitals is not None:
        transformer = ActiveSpaceTransformer(
            num_electrons=active_electrons,
            num_spatial_orbitals=active_orbitals,
        )
        problem = transformer.transform(problem)

    return problem
