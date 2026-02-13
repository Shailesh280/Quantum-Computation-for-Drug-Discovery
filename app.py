from flask import Flask, request, jsonify
from flask_cors import CORS

from quantum_engine.molecule_builder import build_problem
from quantum_engine.vqe_solver import compute_energies
from quantum_engine.binding_energy import compute_binding_energy
from quantum_engine.distance_sweep import run_distance_sweep

app = Flask(__name__)
CORS(app)


@app.route("/compute", methods=["POST"])
def compute():
    data = request.json
    atom_string = data["atom_string"]

    problem = build_problem(atom_string)
    result = compute_energies(problem)

    return jsonify(result)


@app.route("/binding", methods=["POST"])
def binding():
    data = request.json

    result = compute_binding_energy(
        data["drug"],
        data["active_site"],
        data["complex"]
    )

    return jsonify(result)


@app.route("/distance_sweep", methods=["POST"])
def distance_sweep():
    data = request.json

    result = run_distance_sweep(
        data["atom_template"],
        data["min_distance"],
        data["max_distance"],
        data["steps"]
    )

    return jsonify(result)


if __name__ == "__main__":
    app.run(debug=True)
