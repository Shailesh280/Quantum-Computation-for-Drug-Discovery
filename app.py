from flask import Flask, request, jsonify
from flask_cors import CORS

from quantum_engine.molecule_builder import build_problem
from quantum_engine.vqe_solver import compute_energies
from quantum_engine.binding_energy import compute_binding_from_selection
from quantum_engine.distance_sweep import run_distance_sweep_from_selection

app = Flask(__name__)
CORS(app)

from quantum_engine.benchmark_module import run_benchmark


@app.route("/benchmark", methods=["POST"])
def benchmark():

    data = request.json
    molecule = data.get("molecule")

    if not molecule:
        molecule = "H2"

    result = run_benchmark(selected_molecule=molecule)

    return jsonify(result)


@app.route("/binding", methods=["POST"])
def binding():
    data = request.get_json()

    if not data:
        return jsonify({"error": "No JSON body received"}), 400

    disease_name = data.get("disease_name")
    drug_name = data.get("drug_name")
    separation_distance = data.get("separation_distance", 3.0)
    basis = data.get("basis", "sto3g")

    if not disease_name or not drug_name:
        return jsonify({
            "error": "disease_name and drug_name are required",
            "received": data
        }), 400

    result = compute_binding_from_selection(
        disease_name=disease_name,
        drug_name=drug_name,
        separation_distance=separation_distance,
        basis=basis,
    )

    return jsonify(result)


@app.route("/distance_sweep", methods=["POST"])
def distance_sweep():
    data = request.get_json()

    result = run_distance_sweep_from_selection(
        disease_name=data.get("disease_name"),
        drug_name=data.get("drug_name"),
        min_distance=data.get("min_distance"),
        max_distance=data.get("max_distance"),
        steps=data.get("steps"),
    )

    return jsonify(result)



if __name__ == "__main__":
    app.run(debug=False, use_reloader=False)
