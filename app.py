from flask import Flask, request, jsonify
from flask_cors import CORS
import os
from flask import send_from_directory

from quantum_engine.binding_energy import compute_binding_from_selection
from quantum_engine.distance_sweep import run_distance_sweep_from_selection
from quantum_engine.database import DISEASE_TARGETS
from quantum_engine.benchmark_module import run_benchmark

app = Flask(__name__)
CORS(app)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
POSE_DIR = os.path.join(BASE_DIR, "docking_results")

@app.route("/poses/<filename>")
def get_pose(filename):
    return send_from_directory(POSE_DIR, filename)
# =========================================================
# Benchmark Endpoint
# =========================================================

@app.route("/benchmark", methods=["POST"])
def benchmark():

    data = request.json
    molecule = data.get("molecule")

    if not molecule:
        molecule = "H2"

    result = run_benchmark(selected_molecule=molecule)

    return jsonify(result)


# =========================================================
# Get Available Diseases
# =========================================================

@app.route("/diseases", methods=["GET"])
def get_diseases():

    diseases = []

    for key, info in DISEASE_TARGETS.items():

        diseases.append({
            "key": key,
            "name": info["disease_name"]
        })

    return jsonify(diseases)


# =========================================================
# Get Active Sites + Drugs for Disease
# =========================================================

@app.route("/disease_details/<disease_key>", methods=["GET"])
def get_disease_details(disease_key):

    if disease_key not in DISEASE_TARGETS:
        return jsonify({"error": "Invalid disease"}), 400

    disease = DISEASE_TARGETS[disease_key]

    active_sites = []

    for i, site in enumerate(disease["active_sites"]):

        active_sites.append({
            "index": i,
            "label": f'{site["residue_name"]} {site["residue_number"]}'
        })

    return jsonify({
        "disease": disease["disease_name"],
        "active_sites": active_sites,
        "drugs": disease["drugs"]
    })


# =========================================================
# Binding Simulation
# =========================================================

@app.route("/binding", methods=["POST"])
def binding():

    data = request.get_json()

    disease_key = data.get("disease_key")
    active_site_index = data.get("active_site_index")
    drug_name = data.get("drug_name")

    separation_distance = data.get("separation_distance", 3.0)
    basis = data.get("basis", "sto3g")

    result = compute_binding_from_selection(
        disease_key=disease_key,
        active_site_index=active_site_index,
        drug_name=drug_name,
        separation_distance=separation_distance,
        basis=basis,
    )

    print("FINAL RESULT:", result, flush=True)

    import json
    return app.response_class(
        response=json.dumps(result, default=float),
        status=200,
        mimetype="application/json"
    )


# =========================================================
# Distance Sweep
# =========================================================

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


# =========================================================
# Run Server
# =========================================================

if __name__ == "__main__":
    app.run(debug=False, use_reloader=False)