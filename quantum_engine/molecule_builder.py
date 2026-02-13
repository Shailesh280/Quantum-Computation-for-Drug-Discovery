from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.problems import ElectronicStructureProblem
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.units import DistanceUnit


def build_problem(atom_string, basis="sto3g", charge=0, spin=0):
    driver = PySCFDriver(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin,
        unit=DistanceUnit.ANGSTROM,
    )
    problem = driver.run()
    return problem


def get_qubit_operator(problem):
    second_q_ops = problem.second_q_ops()
    hamiltonian = second_q_ops[0]

    mapper = JordanWignerMapper()
    qubit_op = mapper.map(hamiltonian)

    return qubit_op
