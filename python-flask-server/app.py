
# imports
from flask import Flask, request, jsonify
import dcc.compute_utils as cutils 
import dcc.dcc_utils as dutils
import time 

# constants
DEBUG=False

# variables
app = Flask(__name__)
app.secret_key = "test_app_gpt"
map_gene_set_families = {}

logger = dutils.get_logger(__name__)
# p_value_cutoff = 0.3
P_VALUE_CUTOFF = 0.3
# p_value_cutoff = 0.05
MAX_NUMBER_GENE_SETS_FOR_COMPUTATION=100

def process_phenotypes(phenotype_list):
    """
    gets the gene scores from the compute module
    """
    # # For demo: just echo back uppercase names
    # return [p.upper() for p in phenotype_list]
    # log
    logger.info("got input phenotype list: {}".format(phenotype_list))

    # get the gene scores
    map_gene_scores = cutils.process_web_patient_phenotypes(list_input_phenotypes=phenotype_list, debug=True)

    # return
    return map_gene_scores


@app.route('/gene_scores', methods=['GET', 'POST'])
def handle_phenotypes():
    """
    Handles GET and POST requests with a 'phenotypes' input parameter.
    Converts it to a Python list, calls `process_phenotypes`, and returns result.
    """
    logs = []
    result = {"input": None, "phenotype_list": [], "output": None, "logs": logs}

    try:
        # --- 1. Retrieve input depending on method ---
        if request.method == 'POST':
            data = request.get_json(silent=True) or {}
            phenotypes_param = data.get("phenotypes", [])
        else:  # GET
            phenotypes_param = request.args.get("phenotypes", "")

        # --- 2. Convert input to a Python list ---
        phenotype_list = []
        if isinstance(phenotypes_param, list):
            # Already a list
            phenotype_list = [str(p).strip() for p in phenotypes_param if str(p).strip()]
        elif isinstance(phenotypes_param, str):
            # Comma-separated string like "diabetes,gout,asthma"
            phenotype_list = [p.strip() for p in phenotypes_param.split(",") if p.strip()]
        else:
            logs.append("Unsupported input format for 'phenotypes'.")

        result["input"] = phenotypes_param
        result["phenotype_list"] = phenotype_list

        # --- 3. Call downstream function ---
        if phenotype_list:
            output = process_phenotypes(phenotype_list)
            result["output"] = output
        else:
            logs.append("No valid phenotypes provided.")

    except Exception as e:
        logs.append(f"Error: {e}")

    # --- single return point ---
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8082)





    # # split the genes into list
    # if phenotypes:
    #     # build the phenotype list
    #     list_temp = phenotypes.split(",")
    #     list_select = []

    #     for value in list_temp:
    #         gene = value.strip()
    #         print("got phenotype: -{}-".format(gene))
    #         list_select.append(gene)

    #     # get the db connection
    #     conn = phenotype_utils.get_connection()

    #     # get the result diseases
    #     # map_disease = phenotype_utils.get_disease_score_for_phenotype_list(conn=conn, list_curies=list_select, log=False)
    #     list_disease = phenotype_utils.get_disease_score_sorted_list_for_phenotype_list(conn=conn, list_curies=list_select, log=False)
    #     print("got disease list size of: {}".format(len(list_disease)))

    #     # add to map
    #     # map_result['results'] = map_disease
    #     map_result['results'] = list_disease
