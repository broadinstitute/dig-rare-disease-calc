# imports
from flask import Flask, send_from_directory, render_template, request, flash

import dcc.startup_utils as sutils
import dcc.file_utils as futils
import dcc.matrix_utils as mutils 
import dcc.compute_utils as cutils 
import dcc.dcc_utils as dutils 
import dcc.sql_utils as sql_utils
import dcc.gui_utils as gutils

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

# in memory compute variables
map_conf = sutils.load_conf()

# load the data
db_file = map_conf.get('root_dir') +  map_conf.get('db_file')
logger.info("loading database file: {}".format(db_file))
sql_connection = sql_utils.db_sqlite_get_connection(db_path=db_file)

# map_gene_index, list_system_genes = futils.load_gene_file_into_map(file_path=map_conf.get('root_dir') + map_conf.get('gene_file'))
# load the genes - common data across gene sets
map_gene_index, list_system_genes, map_gene_ontology = sql_utils.db_load_gene_table_into_map(conn=sql_connection)
matrix_gene_sets, map_gene_set_index = mutils.load_geneset_matrix(map_gene_index=map_gene_index, 
                                                                  list_gene_set_files=map_conf.get('gene_set_files'), path_gene_set_files=map_conf.get('root_dir'), log=False)
(mean_shifts, scale_factors) = cutils._calc_X_shift_scale(X=matrix_gene_sets)

# load the phewas data
gene_phewas_bfs_in = map_conf.get('root_dir') +  "all_trait_gcat.trait_gene_combined.gt1.txt"
phenos, gene_pheno_Y, gene_pheno_combined_prior_Ys = mutils.read_gene_phewas_bfs(list_system_genes, map_gene_index, gene_phewas_bfs_in)

logger.info("================ Bayes App is UP! ===========================")

# test 
map_gene_set_families = sutils.load_gene_set_family_map(map_conf=map_conf, map_gene_index=map_gene_index, log=False)

logger.info("================ Test App is UP! ===========================")

@app.route("/heartbeat", methods=["GET"])
def heartbeat():
    # get the sizes of each gene set list
    map_gene_set_size = {}
    for key, value in map_gene_set_families.items():
        map_gene_set_size[key] = format(value.get_num_gene_sets(), ",")
       
    map_result = {'message': 'yes, I am up ;>', 'gene_sets': list(map_gene_set_families.keys()), 'gene_set_sizes': map_gene_set_size}

    return map_result


@app.route('/test_graph')
def index_graph():
    return send_from_directory('./static', 'factor_graph.html')  # Adjust 'index.html' to your main file

@app.route('/test_vis')
def index_vis():
    return send_from_directory('./static', 'vis.html')  # Adjust 'index.html' to your main file


@app.route("/query", methods=["POST"])
def post_genes():
    # initialize 
    map_result = {}
    list_input_genes = []
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)
    p_value_cutoff = P_VALUE_CUTOFF
    max_number_gene_sets = MAX_NUMBER_GENE_SETS_FOR_COMPUTATION
    gene_set_family_key = dutils.KEY_DEFAULT_GENE_SET_FAMILY

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get(dutils.KEY_REST_GENES)

    logger.info("got request: {} with gene inputs: {}".format(request.method, list_input_genes))
    logger.info("got gene inputs of size: {}".format(len(list_input_genes)))

    # get the p_value
    p_value_cutoff = process_numeric_value(json_request=data, name=dutils.KEY_REST_P_VALUE, cutoff_default=P_VALUE_CUTOFF)
    logger.info("using p_value: {}".format(p_value_cutoff))

    # get the max gene sets
    max_number_gene_sets = process_numeric_value(json_request=data, name=dutils.KEY_REST_MAX_NUMBER_GENE_SETS, cutoff_default=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION, is_float=False)
    logger.info("using max number gene sets: {}".format(max_number_gene_sets))

    # adding input to indicate the enrichment analysis type
    enrichment_analysis = process_string_value(json_request=data, name=dutils.KEY_REST_ENRICHMENT_ANALYSIS, default=dutils.DEFAULT_ENRICHMENT_ANALYSIS)
    str_message = "using enrichment analysis: {}".format(enrichment_analysis)
    logger.info(str_message)

    # adding input to indicate the factorization weight
    factorization_weight = process_string_value(json_request=data, name=dutils.KEY_REST_FACTORIZATION_WEIGHT, default=dutils.DEFAULT_FACTORIZATION_WEIGHT)
    str_message = "using factorization weight: {}".format(factorization_weight)
    logger.info(str_message)

    # get the gene set family name
    gene_set_family_key = process_string_value(json_request=data, name=dutils.KEY_REST_GENE_SET, default=dutils.KEY_DEFAULT_GENE_SET_FAMILY)
    exclude_controls = process_boolean_value(json_request=data, name=dutils.KEY_REST_EXCLUDE_CONTROLS, default=False)
    if not exclude_controls:
        gene_set_family_with_controls_key = gene_set_family_key + dutils.KEY_NEGATIVE_CONTROLS
        if gene_set_family_key in map_gene_set_families:
            gene_set_family_key = gene_set_family_with_controls_key
        else:
            str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_with_controls_key)
            logger.warning(str_message)
    logger.info("using input gene set family key: {}".format(gene_set_family_key))

    # get the gene set family object
    gene_set_family_object: sutils.GeneSetFamily = map_gene_set_families.get(gene_set_family_key)

    # make sure gene set family available
    if not gene_set_family_object:
        str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_key)
        logger.error(str_message)
        map_result = {"logs": [str_message]}

    else:
        # translate the genes into what the system can handle
        list_input_translated = sql_utils.db_get_gene_names_from_list(conn=sql_conn_query, list_input=list_input_genes)
        logger.info("got translated gene inputs of size: {}".format(len(list_input_translated)))

        # add the genes to the result
        # map_result['input_genes'] = list_input_genes
        if DEBUG:
            map_result['conf'] = map_conf

        # time
        start = time.time()

        # # compute
        # list_factor, list_factor_genes, list_factor_gene_sets, 
        #     gene_factor, gene_set_factor, map_gene_novelty, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=matrix_gene_sets, 
        #                                                                                             p_value=p_value_cutoff,
        #                                                                                             max_num_gene_sets=max_number_gene_sets,
        #                                                                                             list_gene=list_input_translated, 
        #                                                                                             list_system_genes=list_system_genes, 
        #                                                                                             map_gene_index=map_gene_index, 
        #                                                                                             map_gene_set_index=map_gene_set_index,
        #                                                                                             mean_shifts=mean_shifts, 
        #                                                                                             scale_factors=scale_factors,
        #                                                                                             log=True)

        # # time
        # end = time.time()

        # # format the data
        # # map_factors = cutils.group_factor_results(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes)
        # # map_result['data'] = map_factors
        # map_result = gutils.gui_build_results_map(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes, 
        #                                         map_gene_ontology=map_gene_ontology, list_input_gene_names=list_input_translated, map_gene_index=map_gene_index,
        #                                         matrix_gene_sets=matrix_gene_sets, map_gene_novelty=map_gene_novelty)

        list_factor, list_factor_genes, list_factor_gene_sets, \
            gene_factor, gene_set_factor, map_gene_factor_data, list_gene_set_p_values, logs_process = cutils.calculate_factors(
                                                                                                    matrix_gene_sets_gene_original=gene_set_family_object.matrix_gene_sets, 
                                                                                                    enrichment_analysis=enrichment_analysis,
                                                                                                    factorization_weight=factorization_weight,
                                                                                                    p_value=p_value_cutoff,
                                                                                                    max_num_gene_sets=max_number_gene_sets,
                                                                                                    list_gene=list_input_translated, 
                                                                                                    list_system_genes=list_system_genes, 
                                                                                                    map_gene_index=map_gene_index, 
                                                                                                    map_gene_set_index=gene_set_family_object.map_gene_set_index,
                                                                                                    mean_shifts=gene_set_family_object.mean_shifts, 
                                                                                                    scale_factors=gene_set_family_object.scale_factors,
                                                                                                    log=True)

        # time
        end = time.time()

        # format the data
        # map_factors = cutils.group_factor_results(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes)
        # map_result['data'] = map_factors
        map_result = gutils.gui_build_results_map(list_factor=list_factor, 
                                                    list_factor_gene_sets=list_factor_gene_sets, 
                                                    list_factor_genes=list_factor_genes, 
                                                    map_gene_ontology=map_gene_ontology, 
                                                    list_input_gene_names=list_input_translated, 
                                                    map_gene_index=map_gene_index,
                                                    matrix_gene_sets=gene_set_family_object.matrix_gene_sets, 
                                                    map_gene_factor_data=map_gene_factor_data)


        # add time
        logs_process.append("using gene set option: {}".format(gene_set_family_key))
        str_message = "total elapsed time is: {}s".format(end-start)
        logs_process.append(str_message)
        logs_process.append("code version is: {}".format(dutils.get_code_version()))
        map_result['logs'] = logs_process
        for row in logs_process:
            logger.info(row)

    # return
    return map_result


@app.route("/pigean", methods=["POST"])
def post_pigean_genes():
    # initialize 
    map_result = {}
    list_input_genes = []
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)
    p_value_cutoff = P_VALUE_CUTOFF
    list_logs = []

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got request: {} with gene inputs: {}".format(request.method, list_input_genes))
    str_message = "got gene inputs of size: {}".format(len(list_input_genes))
    logger.info(str_message)
    list_logs.append(str_message)

    # get the p_value
    p_value_cutoff = process_numeric_value(json_request=data, name='p_value', cutoff_default=P_VALUE_CUTOFF)
    str_message = "got p_value filter: {}".format(p_value_cutoff)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the max gene sets
    max_number_gene_sets = process_numeric_value(json_request=data, name='max_number_gene_sets', cutoff_default=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION, is_float=False)
    str_message = "got using max number of genes: {}".format(max_number_gene_sets)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate whether to generate factor label names
    is_generate_factor_labels = process_boolean_value(json_request=data, name=dutils.KEY_REST_GENERATE_FACTOR_LABELS, default=False)
    str_message = "using input whether generate factor labels (using LLM): {}".format(is_generate_factor_labels)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate whether to generate factor label names
    is_add_gene_scores = process_boolean_value(json_request=data, name=dutils.KEY_REST_ADD_GENE_SCORES, default=False)
    str_message = "using input whether to also calculate gene scores: {}".format(is_add_gene_scores)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate the enrichment analysis type
    enrichment_analysis = process_string_value(json_request=data, name=dutils.KEY_REST_ENRICHMENT_ANALYSIS, default=dutils.DEFAULT_ENRICHMENT_ANALYSIS)
    str_message = "using enrichment analysis: {}".format(enrichment_analysis)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate the factorization weight
    factorization_weight = process_string_value(json_request=data, name=dutils.KEY_REST_FACTORIZATION_WEIGHT, default=dutils.DEFAULT_FACTORIZATION_WEIGHT)
    str_message = "using factorization weight: {}".format(factorization_weight)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family name
    gene_set_family_key = process_string_value(json_request=data, name=dutils.KEY_REST_GENE_SET, default=dutils.KEY_DEFAULT_GENE_SET_FAMILY)
    exclude_controls = process_boolean_value(json_request=data, name=dutils.KEY_REST_EXCLUDE_CONTROLS, default=False)
    if not exclude_controls:
        gene_set_family_with_controls_key = gene_set_family_key + dutils.KEY_NEGATIVE_CONTROLS
        if gene_set_family_key in map_gene_set_families:
            gene_set_family_key = gene_set_family_with_controls_key
        else:
            str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_with_controls_key)
            logger.warning(str_message)
            list_logs.append(str_message)
    str_message = "using input gene set family key: {}".format(gene_set_family_key)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family object
    gene_set_family_object: sutils.GeneSetFamily = map_gene_set_families.get(gene_set_family_key)

    # make sure gene set family available
    if not gene_set_family_object:
        str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_key)
        logger.error(str_message)
        map_result = {"logs": [str_message]}

    else:
        # translate the genes into what the system can handle
        list_input_translated = sql_utils.db_get_gene_names_from_list(conn=sql_conn_query, list_input=list_input_genes)
        str_message = "got translated gene inputs of size: {}".format(len(list_input_translated))
        logger.info(str_message)
        list_logs.append(str_message)

        # add the genes to the result
        # map_result['input_genes'] = list_input_genes
        if DEBUG:
            map_result['conf'] = map_conf

        # time
        start = time.time()

        # compute
        # list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, \
        # gene_set_factor, map_gene_novelty, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=matrix_gene_sets, 
        #                                                                                                         p_value=p_value_cutoff,
        #                                                                                                         max_num_gene_sets=max_number_gene_sets,
        #                                                                                                         list_gene=list_input_translated, 
        #                                                                                                         list_system_genes=list_system_genes, 
        #                                                                                                         map_gene_index=map_gene_index, map_gene_set_index=map_gene_set_index,
        #                                                                                                         mean_shifts=mean_shifts, scale_factors=scale_factors,
        #                                                                                                         is_factor_labels_llm=is_generate_factor_labels,
        #                                                                                                         log=True)

        list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, \
        gene_set_factor, map_gene_novelty, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=gene_set_family_object.matrix_gene_sets, 
                                                                                                                enrichment_analysis=enrichment_analysis,
                                                                                                                factorization_weight=factorization_weight,
                                                                                                                p_value=p_value_cutoff,
                                                                                                                max_num_gene_sets=max_number_gene_sets,
                                                                                                                list_gene=list_input_translated, 
                                                                                                                list_system_genes=list_system_genes, 
                                                                                                                map_gene_index=map_gene_index, 
                                                                                                    map_gene_set_index=gene_set_family_object.map_gene_set_index,
                                                                                                    mean_shifts=gene_set_family_object.mean_shifts, 
                                                                                                    scale_factors=gene_set_family_object.scale_factors,
                                                                                                                is_factor_labels_llm=is_generate_factor_labels,
                                                                                                                log=True)
        # format the data
        # map_factors = cutils.group_factor_results(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes)
        # map_result['data'] = map_factors
        map_result = gutils.gui_build_pigean_app_results_map(list_input_genes=list_input_genes, list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, 
                                                            list_factor_genes=list_factor_genes, list_gene_set_p_values=list_gene_set_p_values)


        if is_add_gene_scores:
            # get the gene scores
            map_gene_scores, map_gene_set_scores, logs_gene_scores = cutils.calculate_gene_scores_map(matrix_gene_sets=gene_set_family_object.matrix_gene_sets, 
                                                                                                      list_input_genes=list_input_genes, 
                                                                                                      map_gene_index=map_gene_index, 
                                                                                                      map_gene_set_index=gene_set_family_object.map_gene_set_index,
                                                                                                      list_system_genes=list_system_genes,
                                                                                                      input_scale_factors=gene_set_family_object.scale_factors, 
                                                                                                      input_mean_shifts=gene_set_family_object.mean_shifts, 
                                                                                                      max_num_gene_sets=max_number_gene_sets, 
                                                                                                      log=True)


            # get the gene score gui elements
            map_result_gene_score = gutils.gui_build_gene_scores_app_results_map(map_gene_scores=map_gene_scores, map_gene_set_scores=map_gene_set_scores)

            # add the gene score data tologs and map result
            logs_process.append(logs_gene_scores)
            map_result.update(map_result_gene_score)

        # 20250129 - adding map for network graph data
        map_network_result = gutils.build_graph_node_edge_map(list_factor_input=list_factor, list_factor_gene_sets_input=list_factor_gene_sets, 
                                                            list_factor_genes_input=list_factor_genes)
        map_result.update(map_network_result)

        # time
        end = time.time()

        # add time
        str_message = "post /pigean total elapsed time is: {}s".format(end-start)
        logs_process.append(str_message)
        logs_process.append("code version is: {}".format(dutils.get_code_version()))

        # add the input to the logs
        logs_process = list_logs + logs_process
        map_result['logs'] = logs_process

        # print logs to app log
        for row in logs_process:
            logger.info(row)

    # return
    return map_result


@app.route("/translator_gene", methods=["POST"])
def post_translator_gene():
    # initialize 
    map_result = {}
    list_input_genes = []
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)
    p_value_cutoff = P_VALUE_CUTOFF
    list_logs = []

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got request: {} with gene inputs: {}".format(request.method, list_input_genes))
    str_message = "got gene inputs of size: {}".format(len(list_input_genes))
    logger.info(str_message)
    list_logs.append(str_message)

    # get the p_value
    p_value_cutoff = process_numeric_value(json_request=data, name='p_value', cutoff_default=P_VALUE_CUTOFF)
    str_message = "got p_value filter: {}".format(p_value_cutoff)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the max gene sets
    max_number_gene_sets = process_numeric_value(json_request=data, name='max_number_gene_sets', cutoff_default=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION, is_float=False)
    str_message = "got using max number of genes: {}".format(max_number_gene_sets)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate the enrichment analysis type
    enrichment_analysis = process_string_value(json_request=data, name=dutils.KEY_REST_ENRICHMENT_ANALYSIS, default=dutils.DEFAULT_ENRICHMENT_ANALYSIS)
    str_message = "using enrichment analysis: {}".format(enrichment_analysis)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate the factorization weight
    factorization_weight = process_string_value(json_request=data, name=dutils.KEY_REST_FACTORIZATION_WEIGHT, default=dutils.DEFAULT_FACTORIZATION_WEIGHT)
    str_message = "using factorization weight: {}".format(factorization_weight)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family name
    gene_set_family_key = process_string_value(json_request=data, name=dutils.KEY_REST_GENE_SET, default=dutils.KEY_DEFAULT_GENE_SET_FAMILY)
    exclude_controls = process_boolean_value(json_request=data, name=dutils.KEY_REST_EXCLUDE_CONTROLS, default=True)
    if not exclude_controls:
        gene_set_family_with_controls_key = gene_set_family_key + dutils.KEY_NEGATIVE_CONTROLS
        if gene_set_family_key in map_gene_set_families:
            gene_set_family_key = gene_set_family_with_controls_key
        else:
            str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_with_controls_key)
            logger.warning(str_message)
            list_logs.append(str_message)
    str_message = "using input gene set family key: {}".format(gene_set_family_key)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family object
    gene_set_family_object: sutils.GeneSetFamily = map_gene_set_families.get(gene_set_family_key)

    # make sure gene set family available
    if not gene_set_family_object:
        str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_key)
        logger.error(str_message)
        map_result = {"logs": [str_message]}

    else:
        # translate the genes into what the system can handle
        list_input_translated = sql_utils.db_get_gene_names_from_list(conn=sql_conn_query, list_input=list_input_genes)
        str_message = "got translated gene inputs of size: {}".format(len(list_input_translated))
        logger.info(str_message)
        list_logs.append(str_message)

        # add the genes to the result
        # map_result['input_genes'] = list_input_genes
        if DEBUG:
            map_result['conf'] = map_conf

        # time
        start = time.time()

        # compute
        list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, \
        gene_set_factor, map_gene_novelty, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=gene_set_family_object.matrix_gene_sets, 
                                                                                                                enrichment_analysis=enrichment_analysis,
                                                                                                                factorization_weight=factorization_weight,
                                                                                                                p_value=p_value_cutoff,
                                                                                                                max_num_gene_sets=max_number_gene_sets,
                                                                                                                list_gene=list_input_translated, 
                                                                                                                list_system_genes=list_system_genes, 
                                                                                                                map_gene_index=map_gene_index, 
                                                                                                    map_gene_set_index=gene_set_family_object.map_gene_set_index,
                                                                                                    mean_shifts=gene_set_family_object.mean_shifts, 
                                                                                                    scale_factors=gene_set_family_object.scale_factors,
                                                                                                                is_factor_labels_llm=True,
                                                                                                                log=True)
        # format the data
        # map_factors = cutils.group_factor_results(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes)
        # map_result['data'] = map_factors
        map_result = gutils.gui_build_translator_gene_results_map(list_input_genes=list_input_genes, list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, 
                                                            list_factor_genes=list_factor_genes, list_gene_set_p_values=list_gene_set_p_values)


        # time
        end = time.time()

        # add time
        str_message = "post /translator_gene total elapsed time is: {}s".format(end-start)
        logs_process.append(str_message)
        logs_process.append("code version is: {}".format(dutils.get_code_version()))

        # add the input to the logs
        logs_process = list_logs + logs_process
        map_result['logs'] = logs_process

        # print logs to app log
        for row in logs_process:
            logger.info(row)

    # return
    return map_result


@app.route("/gene_scores", methods=["POST"])
def post_gene_scores():
    # initialize 
    map_result = {}
    list_input_genes = []
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)
    p_value_cutoff = P_VALUE_CUTOFF
    list_logs = []

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got gene scores request: {} with gene inputs: {}".format(request.method, list_input_genes))
    str_message = "got gene inputs of size: {}".format(len(list_input_genes))
    logger.info(str_message)
    list_logs.append(str_message)

    # get the p_value
    p_value_cutoff = process_numeric_value(json_request=data, name='p_value', cutoff_default=P_VALUE_CUTOFF)
    str_message = "got using p_value: {}".format(p_value_cutoff)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the max gene sets
    max_number_gene_sets = process_numeric_value(json_request=data, name='max_number_gene_sets', cutoff_default=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION, is_float=False)
    str_message = "got using max number of gene sets: {}".format(max_number_gene_sets)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family name
    gene_set_family_key = process_string_value(json_request=data, name=dutils.KEY_REST_GENE_SET, default=dutils.KEY_DEFAULT_GENE_SET_FAMILY)
    exclude_controls = process_boolean_value(json_request=data, name=dutils.KEY_REST_EXCLUDE_CONTROLS, default=False)
    if not exclude_controls:
        gene_set_family_with_controls_key = gene_set_family_key + dutils.KEY_NEGATIVE_CONTROLS
        if gene_set_family_key in map_gene_set_families:
            gene_set_family_key = gene_set_family_with_controls_key
        else:
            str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_with_controls_key)
            logger.warning(str_message)
            list_logs.append(str_message)
    str_message = "using input gene set family key: {}".format(gene_set_family_key)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family object
    gene_set_family_object: sutils.GeneSetFamily = map_gene_set_families.get(gene_set_family_key)

    # make sure gene set family available
    if not gene_set_family_object:
        str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_key)
        logger.error(str_message)
        map_result = {"logs": [str_message]}

    else:
        # translate the genes into what the system can handle
        list_input_translated = sql_utils.db_get_gene_names_from_list(conn=sql_conn_query, list_input=list_input_genes)
        str_message = "got translated gene inputs of size: {}".format(len(list_input_translated))
        logger.info(str_message)
        list_logs.append(str_message)

        # add the genes to the result
        # map_result['input_genes'] = list_input_genes
        if DEBUG:
            map_result['conf'] = map_conf

        # time
        start = time.time()

        map_gene_scores, map_gene_set_scores, logs_process = cutils.calculate_gene_scores_map(matrix_gene_sets=gene_set_family_object.matrix_gene_sets, 
                                                                                              list_input_genes=list_input_genes, 
                                                                                              map_gene_index=map_gene_index, 
                                                                                              map_gene_set_index=gene_set_family_object.map_gene_set_index,
                                                                                              list_system_genes=list_system_genes,
                                                                                              input_scale_factors=gene_set_family_object.scale_factors, 
                                                                                              input_mean_shifts=gene_set_family_object.mean_shifts, log=True)

        # list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, \
        # gene_set_factor, map_gene_novelty, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=gene_set_family_object.matrix_gene_sets, 
        #                                                                                                         p_value=p_value_cutoff,
        #                                                                                                         max_num_gene_sets=max_number_gene_sets,
        #                                                                                                         list_gene=list_input_translated, 
        #                                                                                                         list_system_genes=list_system_genes, 
        #                                                                                                         map_gene_index=map_gene_index, 
        #                                                                                             map_gene_set_index=gene_set_family_object.map_gene_set_index,
        #                                                                                             mean_shifts=gene_set_family_object.mean_shifts, 
        #                                                                                             scale_factors=gene_set_family_object.scale_factors,
        #                                                                                                         is_factor_labels_llm=is_generate_factor_labels,
        #                                                                                                         log=True)
        # time
        end = time.time()

        # format the data
        # map_factors = cutils.group_factor_results(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes)
        # map_result['data'] = map_factors
        map_result = gutils.gui_build_gene_scores_app_results_map(map_gene_scores=map_gene_scores, map_gene_set_scores=map_gene_set_scores)


        # add time
        str_message = "post /gene_scores total elapsed time is: {}s".format(end-start)
        logs_process.append(str_message)
        logs_process.append("code version is: {}".format(dutils.get_code_version()))
        map_result['logs'] = logs_process

        # add the input to the logs
        logs_process = list_logs + logs_process

        # print logs to app log
        for row in logs_process:
            logger.info(row)

    # return
    return map_result


@app.route("/network_graph", methods=["POST"])
def post_network_graph():
    # initialize 
    map_result = {}
    list_input_genes = []
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)
    p_value_cutoff = P_VALUE_CUTOFF
    list_logs = []

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got gene scores request: {} with gene inputs: {}".format(request.method, list_input_genes))
    str_message = "got gene inputs of size: {}".format(len(list_input_genes))
    logger.info(str_message)
    list_logs.append(str_message)

    # get the p_value
    p_value_cutoff = process_numeric_value(json_request=data, name='p_value', cutoff_default=P_VALUE_CUTOFF)
    str_message = "got using p_value: {}".format(p_value_cutoff)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the max gene sets
    max_number_gene_sets = process_numeric_value(json_request=data, name='max_number_gene_sets', cutoff_default=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION, is_float=False)
    str_message = "got using max number of gene sets: {}".format(max_number_gene_sets)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family name
    gene_set_family_key = process_string_value(json_request=data, name=dutils.KEY_REST_GENE_SET, default=dutils.KEY_DEFAULT_GENE_SET_FAMILY)
    exclude_controls = process_boolean_value(json_request=data, name=dutils.KEY_REST_EXCLUDE_CONTROLS, default=False)
    if not exclude_controls:
        gene_set_family_with_controls_key = gene_set_family_key + dutils.KEY_NEGATIVE_CONTROLS
        if gene_set_family_key in map_gene_set_families:
            gene_set_family_key = gene_set_family_with_controls_key
        else:
            str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_with_controls_key)
            logger.warning(str_message)
            list_logs.append(str_message)
    str_message = "using input gene set family key: {}".format(gene_set_family_key)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate whether to generate factor label names
    is_generate_factor_labels = process_boolean_value(json_request=data, name=dutils.KEY_REST_GENERATE_FACTOR_LABELS, default=False)
    str_message = "using input whether generate factor labels (using LLM): {}".format(is_generate_factor_labels)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate the enrichment analysis type
    enrichment_analysis = process_string_value(json_request=data, name=dutils.KEY_REST_ENRICHMENT_ANALYSIS, default=dutils.DEFAULT_ENRICHMENT_ANALYSIS)
    str_message = "using enrichment analysis: {}".format(enrichment_analysis)
    logger.info(str_message)
    list_logs.append(str_message)

    # adding input to indicate the factorization weight
    factorization_weight = process_string_value(json_request=data, name=dutils.KEY_REST_FACTORIZATION_WEIGHT, default=dutils.DEFAULT_FACTORIZATION_WEIGHT)
    str_message = "using factorization weight: {}".format(factorization_weight)
    logger.info(str_message)
    list_logs.append(str_message)

    # get the gene set family object
    gene_set_family_object: sutils.GeneSetFamily = map_gene_set_families.get(gene_set_family_key)

    # make sure gene set family available
    if not gene_set_family_object:
        str_message = "got gene set family key which is not loaded: {}".format(gene_set_family_key)
        logger.error(str_message)
        map_result = {"logs": [str_message]}

    else:
        # translate the genes into what the system can handle
        list_input_translated = sql_utils.db_get_gene_names_from_list(conn=sql_conn_query, list_input=list_input_genes)
        str_message = "got translated gene inputs of size: {}".format(len(list_input_translated))
        logger.info(str_message)
        list_logs.append(str_message)

        # add the genes to the result
        # map_result['input_genes'] = list_input_genes
        if DEBUG:
            map_result['conf'] = map_conf

        # time
        start = time.time()

        list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, \
        gene_set_factor, map_gene_novelty, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=gene_set_family_object.matrix_gene_sets, 
                                                                                                                enrichment_analysis=enrichment_analysis,
                                                                                                                factorization_weight=factorization_weight,
                                                                                                                p_value=p_value_cutoff,
                                                                                                                max_num_gene_sets=max_number_gene_sets,
                                                                                                                list_gene=list_input_translated, 
                                                                                                                list_system_genes=list_system_genes, 
                                                                                                                map_gene_index=map_gene_index, 
                                                                                                    map_gene_set_index=gene_set_family_object.map_gene_set_index,
                                                                                                    mean_shifts=gene_set_family_object.mean_shifts, 
                                                                                                    scale_factors=gene_set_family_object.scale_factors,
                                                                                                                is_factor_labels_llm=is_generate_factor_labels,
                                                                                                                log=False)
        
        # easier to format the data for gene nmf call and thn extract the nodes/edges
        # get the /pigean data
        # TODO - abstract out the code that cobines the factor/gene/geneset data?
        # map_result_intermediate = gutils.gui_build_pigean_app_results_map(list_input_genes=list_input_genes, list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, 
        #                                                     list_factor_genes=list_factor_genes, list_gene_set_p_values=list_gene_set_p_values)

        # list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, \
        # gene_set_factor, map_gene_novelty, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=gene_set_family_object.matrix_gene_sets, 
        #                                                                                                         p_value=p_value_cutoff,
        #                                                                                                         max_num_gene_sets=max_number_gene_sets,
        #                                                                                                         list_gene=list_input_translated, 
        #                                                                                                         list_system_genes=list_system_genes, 
        #                                                                                                         map_gene_index=map_gene_index, 
        #                                                                                             map_gene_set_index=gene_set_family_object.map_gene_set_index,
        #                                                                                             mean_shifts=gene_set_family_object.mean_shifts, 
        #                                                                                             scale_factors=gene_set_family_object.scale_factors,
        #                                                                                                         is_factor_labels_llm=is_generate_factor_labels,
        #                                                                                                         log=True)
        # time
        end = time.time()

        # format the data
        # map_factors = cutils.group_factor_results(list_factor=list_factor, list_factor_gene_sets=list_factor_gene_sets, list_factor_genes=list_factor_genes)
        # map_result['data'] = map_factors
        map_result = gutils.build_graph_node_edge_map(list_factor_input=list_factor, list_factor_gene_sets_input=list_factor_gene_sets, 
                                                            list_factor_genes_input=list_factor_genes)

        # for debug, add pigean list
        map_result['list_factor'] = list_factor

        # add time
        str_message = "post /network_graph total elapsed time is: {}s".format(end-start)
        logs_process.append(str_message)
        logs_process.append("code version is: {}".format(dutils.get_code_version()))
        map_result['logs'] = logs_process

        # add the input to the logs
        logs_process = list_logs + logs_process

        # print logs to app log
        for row in logs_process:
            logger.info(row)

    # return
    return map_result


def log_and_add_to_list(logger, list_input, str_message):
    '''
    helper method to log a message and add to a log list
    '''
    # log message
    logger.info(str_message)

    # add to list
    if not list_input:
        list_input = []
    list_input.append(str_message)

    # return
    return list_input


@app.route("/phenotypes", methods=["POST"])
def post_phenotypes():
    # initialize 
    map_result = {}
    list_input_genes = []
    list_logs = []

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got phenotypes request: {} with gene inputs: {}".format(request.method, list_input_genes))
    str_message = "got gene inputs of size: {}".format(len(list_input_genes))
    logger.info(str_message)
    list_logs.append(str_message)

    # get max number of phenotypes
    max_number_phenotypes = process_numeric_value(json_request=data, name='max_number_phenotypes', cutoff_default=100, is_float=False)
    str_message = "got using max number of phenotypes: {}".format(max_number_phenotypes)
    logger.info(str_message)
    list_logs.append(str_message)

    # time
    start = time.time()

    # run phewas
    p_values, beta_tildes, ses = cutils.calculate_phewas(list_input_genes, list_system_genes, map_gene_index, phenos, gene_pheno_Y, gene_pheno_combined_prior_Ys)

    # build the results
    result_list = cutils.build_phewas_p_value_list(phenos, p_values, max_number_phenotypes)

    # time
    end = time.time()

    # log filter results
    str_message = "got phenotypes results of size: {} out of {} phenotypes".format(len(result_list), len(phenos))
    logger.info(str_message)
    list_logs.append(str_message)

    # add time
    str_message = "post /phenotypes total elapsed time is: {}s".format(end-start)
    list_logs.append(str_message)

    # add results to map
    map_result['phenotypes'] = result_list
    map_result['logs'] = list_logs

    # return
    return map_result


@app.route("/novelty_query", methods=["POST"])
def post_novelty_genes():
    '''
    will respond to the novelty score
    '''
    # initialize 
    map_result = {}
    list_input_genes = []
    p_value_cutoff = P_VALUE_CUTOFF

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got request: {} with gene inputs: {}".format(request.method, list_input_genes))
    logger.info("got gene inputs of size: {}".format(len(list_input_genes)))

    # get the p_value
    p_value_cutoff = process_numeric_value(json_request=data, name='p_value', cutoff_default=P_VALUE_CUTOFF)
    logger.info("got using p_value: {}".format(p_value_cutoff))

    # get the max gene sets
    max_number_gene_sets = process_numeric_value(json_request=data, name='max_number_gene_sets', cutoff_default=MAX_NUMBER_GENE_SETS_FOR_COMPUTATION)
    logger.info("got using p_value: {}".format(p_value_cutoff))

    # adding input to indicate the enrichment analysis type
    enrichment_analysis = process_string_value(json_request=data, name=dutils.KEY_REST_ENRICHMENT_ANALYSIS, default=dutils.DEFAULT_ENRICHMENT_ANALYSIS)
    logger.info("using enrichment analysis: {}".format(enrichment_analysis))

    # adding input to indicate the factorization weight
    factorization_weight = process_string_value(json_request=data, name=dutils.KEY_REST_FACTORIZATION_WEIGHT, default=dutils.DEFAULT_FACTORIZATION_WEIGHT)
    logger.info("using factorization weight: {}".format(factorization_weight))

    # get the calculated data
    map_gene_factor_data, list_input_translated = process_genes(list_input_genes=list_input_genes, p_value_cutoff=p_value_cutoff)

    # format the data
    map_result = gutils.gui_build_novelty_results_map(map_gene_ontology=map_gene_ontology, list_input_gene_names=list_input_translated, map_gene_index=map_gene_index,
                                              matrix_gene_sets=matrix_gene_sets, map_gene_factor_data=map_gene_factor_data)


    # return
    return map_result



@app.route("/curie_query", methods=["POST"])
def post_gene_curies():
    # initialize 
    list_input_genes = []
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)

    # get the input
    data = request.get_json()
    if data:
        list_input_genes = data.get('genes')

    logger.info("got request: {} with gene inputs: {}".format(request.method, list_input_genes))
    logger.info("got gene inputs of size: {}".format(len(list_input_genes)))

    # translate the genes into what the system can handle
    list_input_translated = sql_utils.db_get_gene_curies_from_list(conn=sql_conn_query, list_input=list_input_genes)

    # return
    return list_input_translated


def process_genes(list_input_genes, p_value_cutoff, enrichment_analysis=dutils.DEFAULT_ENRICHMENT_ANALYSIS, factorization_weight=dutils.DEFAULT_FACTORIZATION_WEIGHT, log=False):
    '''
    processes the input genes
    '''
    # initialize 
    sql_conn_query = sql_utils.db_sqlite_get_connection(db_path=db_file)

    # preprocess
    # translate the genes into what the system can handle
    list_input_translated = sql_utils.db_get_gene_names_from_list(conn=sql_conn_query, list_input=list_input_genes)
    logger.info("got translated gene inputs of size: {}".format(len(list_input_translated)))

    # do the calculations
    list_factor, list_factor_genes, list_factor_gene_sets, gene_factor, gene_set_factor, map_gene_factor_data, list_gene_set_p_values, logs_process = cutils.calculate_factors(matrix_gene_sets_gene_original=matrix_gene_sets, 
                                                                                                               enrichment_analysis=enrichment_analysis,
                                                                                                               factorization_weight=factorization_weight,
                                                                                                               p_value=p_value_cutoff,
                                                                                                               list_gene=list_input_translated, 
                                                                                                               list_system_genes=list_system_genes, 
                                                                                                               map_gene_index=map_gene_index, map_gene_set_index=map_gene_set_index,
                                                                                                               mean_shifts=mean_shifts, scale_factors=scale_factors,
                                                                                                               log=True)
    
    # log
    for row in logs_process:
        logger.info(row)

    # return
    return map_gene_factor_data, list_input_translated


def process_numeric_value(json_request, name, cutoff_default, is_float=True, log=False):
    '''
    will extract the p_value from the request; will return the default if none
    '''
    # initialize
    numeric_value = None

    # extract the value
    try:
        # Retrieve the float value from the request
        numeric_input = json_request.get(name)

        # Attempt to convert the input to a float
        if numeric_input is None:
            raise ValueError("No '{}' field provided in the request".format(name))

        if is_float:
            numeric_value = float(numeric_input)
        else:
            numeric_value = int(numeric_input)


    except ValueError as e:
        logger.error("got error for {}: {}".format(name, str(e)))
        numeric_value = cutoff_default

    except Exception as e:
        logger.error("got error for {}: {}".format(name, str(e)))
        numeric_value = cutoff_default

    # return
    return numeric_value


def process_string_value(json_request, name, default, log=False):
    '''
    will extract the string value given by the key
    '''
    # initialize
    value = None

    # extract the value
    # Retrieve the float value from the request
    value = json_request.get(name)

    # if null
    if not value:
        value = default

    # return
    return value

def process_boolean_value(json_request, name, default=False, log=False):
    '''
    Convert string values to boolean.
    '''
    # initialize
    value = None
    result = False

    # extract the value
    # Retrieve the float value from the request
    value = json_request.get(name)

    # if null
    if value:
        # if type boolean, return as is
        if isinstance(value, bool):
            result = value
        else:
           result = value.lower() in ['true', '1', 'yes', 'on']

    # return
    return result


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
