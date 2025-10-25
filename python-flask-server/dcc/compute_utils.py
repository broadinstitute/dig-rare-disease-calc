
# imports
import requests
import json
from typing import Dict, Any, List

# COMPUTE LOGIC
# - For each patient phenotype, hit this to get ranked list of phenotypes back, maybe take all top within 75% of top score?
#   - N Patient phenotypes -> M pigean phenotypes (with mapping)
# - For each pigean phenotype, query Pigean endpoint to get combined_D score (probability) and take all genes within 75% of top gene
# - Now we have N patient phenotypes -> M pigean phenotypes -> K genes
# - For each patient phenotype, take weighted mean of gene probabilities across all pigean phenotypes weighted 
#     by semantic match of pigean phenotype to patient phenotype, to get a list of gene probabilities per patient phenotypes (edited) 
#     (equivalent to normalizing semantic weights to sum to 1 per patient phenotype)
# - Finally, get two patient-level gene lists: one is max gene probability across all patient phenotypes, 
#   other is average probability (assume probability is 0.05 if gene is not on list for a patient phenotype) — allow for weighted average if user inputs weights per patient phenotype


# constants
URL_PHENOTYPE_ENBEDDINGS = "https://api.kpndataregistry.org/api/search/phenotypes"
KEY_DATA = 'data'
KEY_LOGS = 'logs'
KEY_ID = 'id'
KEY_DESCRIPTION = 'description'
KEY_SCORE = 'score'


# methods
def get_rest_phenotype_similarity(term: str, similarity_threshold: float=0.0, percent: float=75.0, log: bool=False):
    """
    Queries the KPN Data Registry phenotype API for a given search term.
    Returns a dictionary with id, description, group, and score.
    Ensures type safety and has a single return statement.
    """
    # initialize
    url = "https://api.kpndataregistry.org/api/search/phenotypes"
    params = {"q": term, "similarity_threshold": similarity_threshold}
    logs = []
    result = {KEY_DATA: [], KEY_LOGS: logs}  # initialize single return object

    try:
        # Perform GET request
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()  # Raise exception for HTTP errors

        # Parse JSON response
        json_data = response.json()

        # Extract and type-cast relevant fields
        entries = []
        for item in json_data.get(KEY_DATA, []):
            try:
                entry = {
                    KEY_ID: str(item.get(KEY_ID, "")),
                    KEY_DESCRIPTION: str(item.get(KEY_DESCRIPTION, "")),
                    KEY_SCORE: float(item.get(KEY_SCORE, 0.0))
                }

                entries.append(entry)

            except (ValueError, TypeError) as e:
                result.get(KEY_LOGS).append(f"Skipping malformed entry: {e}")

        # Filter within specified percent of top score
        if entries:
            top_score = max(e[KEY_SCORE] for e in entries)
            threshold = (float(percent) / 100.0) * top_score
            filtered = [e for e in entries if e[KEY_SCORE] >= threshold]
            result[KEY_DATA] = sorted(filtered, key=lambda x: x[KEY_SCORE], reverse=True)

            # log
            logs.append("using score cutoff of: {}".format(threshold))
            logs.append("filter: {} to cutoff of: {}".format(len(entries), len(result[KEY_DATA])))

        else:
            logs.append("No entries found in API response.")

    except requests.exceptions.RequestException as e:
        result.get(KEY_LOGS).append(f"Request error: {e}")

    except ValueError as e:
        result.get(KEY_LOGS).append(f"JSON decoding error: {e}")

    # single return point
    return result


def get_rest_genes_for_pigean_phenotype(term: str, sigma: int = 2, size: str = "Small", keep_frac: float = 0.75):
    """
    Query the PIGEAn gene-phenotype endpoint for a phenotype term and return genes whose
    combined_D probability is within `keep_frac` of the top gene's score.
    
    - Robust to varied field names: tries 'combined_D', 'combined_d', 'combined_proba', 'combined' (case-insensitive).
    - Single return statement.
    - Includes try/except and type casting.
    """
    url = "https://bioindex-dev.hugeamp.org/api/bio/query/pigean-gene-phenotype"
    params = {"q": f"{term},{sigma},{size}"}
    logs = []
    key_score = 'combined'
    result: Dict[str, Any] = {
        "q": [term, str(sigma), size],
        "score_key": None,
        "top_score": None,
        "threshold_score": None,
        "keep_fraction": float(keep_frac),
        KEY_LOGS: logs,
        "data": []  # list of {"gene": str, "score": float}
    }

    # Candidate score keys in priority order (case-insensitive handling)
    candidate_keys = ["combined_D", "combined_d", "combined_proba", "combined"]

    try:
        resp = requests.get(url, params=params, timeout=15)
        resp.raise_for_status()

        try:
            payload = resp.json()

        except ValueError as e:
            logs.append(f"JSON decoding error: {e}")

        else:
            # Get list of dicts
            items = payload.get("data", [])
            if not isinstance(items, list):
                items = []

            # # Determine actual score key present (case-insensitive)
            # score_key_found = None
            # if items:
            #     lower_keys = {k.lower(): k for k in items[0].keys()}
            #     for ck in candidate_keys:
            #         if ck.lower() in lower_keys:
            #             score_key_found = lower_keys[ck.lower()]
            #             break

            result["score_key"] = key_score

            # Collect (gene, score) with casting
            scores: List[Dict[str, Any]] = []
            for it in items:
                try:
                    gene = str(it.get("gene", "")).strip()
                    # Skip rows without a gene name
                    if not gene:
                        continue
                    raw_score = it.get(key_score, None)
                    # Cast to float, skip if invalid
                    if raw_score is None:
                        continue
                    score = float(raw_score)
                    scores.append({"gene": gene, "score": score})

                except (TypeError, ValueError):
                    logs.append(f"skipping item: {it}, decoding error: {e}")
                    # Skip malformed entries
                    continue

            # Compute top score and threshold, then filter
            if scores:
                top = max(x["score"] for x in scores)
                threshold = float(keep_frac) * top
                result["top_score"] = float(top)
                result["threshold_score"] = float(threshold)

                kept = [x for x in scores if x["score"] >= threshold]
                # Optional: sort descending by score for readability
                kept.sort(key=lambda x: x["score"], reverse=True)
                result["data"] = kept

    except requests.exceptions.RequestException as e:
        logs.append(f"Request error: {e}")

    # Single return point
    return result


def compute_weighted_gene_scores(
    phenotypes_json: Dict,
    pigean_results: Dict[str, Dict]
) -> Dict[str, List[Dict]]:
    """
    Combine per-phenotype PIGEAn gene scores into weighted mean probabilities,
    weighted by normalized semantic similarity of phenotype to patient phenotype.

    Args:
        phenotypes_json: dict containing {"data": [{"id": ..., "score": ...}, ...]}
        pigean_results: dict mapping phenotype_id -> {"data": [{"gene": ..., "score": ...}, ...]}

    Returns:
        dict with:
          {
            "data": [{"gene": "GENE1", "weighted_score": value}, ...],
            "logs": [...]
          }
    """
    logs = []
    result = {"data": [], "logs": logs}

    try:
        # --- Extract and normalize phenotype similarity scores ---
        phenotype_entries = phenotypes_json.get("data", [])
        if not phenotype_entries:
            logs.append("No phenotype entries provided.")
            return result

        # Normalize scores so they sum to 1 (for weighting)
        total_score = sum(float(p.get("score", 0.0)) for p in phenotype_entries)
        if total_score <= 0:
            logs.append("All phenotype scores are zero or invalid.")
            return result

        weights = {p["id"]: float(p["score"]) / total_score for p in phenotype_entries}
        logs.append(f"Normalized weights: {weights}")

        # --- Aggregate gene scores across phenotypes ---
        combined_scores = {}

        for phenotype_id, weight in weights.items():
            data_block = pigean_results.get(phenotype_id, {}).get("data", [])
            if not data_block:
                logs.append(f"No gene data found for phenotype '{phenotype_id}'.")
                continue

            for entry in data_block:
                try:
                    gene = str(entry.get("gene", "")).strip()
                    score = float(entry.get("score", 0.0))
                    if not gene:
                        continue

                    # Accumulate weighted score
                    combined_scores[gene] = combined_scores.get(gene, 0.0) + score * weight
                except (ValueError, TypeError):
                    logs.append(f"Skipping malformed gene entry in {phenotype_id}: {entry}")

        # --- Prepare sorted output ---
        result["data"] = [
            {"gene": g, "weighted_score": s}
            for g, s in sorted(combined_scores.items(), key=lambda x: x[1], reverse=True)
        ]

    except Exception as e:
        logs.append(f"Error computing weighted gene scores: {e}")

    return result


def process_patient_phenotypes(
    list_input_phenotypes: List[str], debug: bool=False) -> Dict[str, List[Dict]]:
    """
    Combine per-phenotype PIGEAn gene scores into weighted mean probabilities,
    weighted by normalized semantic similarity of phenotype to patient phenotype.

    Args:
        phenotypes_json: dict containing {"data": [{"id": ..., "score": ...}, ...]}
        pigean_results: dict mapping phenotype_id -> {"data": [{"gene": ..., "score": ...}, ...]}

    Returns:
        dict with:
          {
            "data": [{"gene": "GENE1", "weighted_score": value}, ...],
            "logs": [...]
          }
    """
    # initialize
    map_patient_pigean_phenotypes = {}
    map_pigean_phenotype_genes = {}
    map_patient_phenotype_gene_weights = {}

    # for each patient phenotype, get the similarity scores for pigean phenotypes
    for phenotype_patient in list_input_phenotypes:
        map_patient_pigean_phenotypes[phenotype_patient] = get_rest_phenotype_similarity(term=phenotype_patient)

    # debug log
    if debug:
        print("got patient phenotype to pigean map: \n{}".format(json.dumps(map_patient_pigean_phenotypes, indent=2)))

    # build a map of the distinct pigean phenotypes and their gene scores
    for row_patient in map_patient_pigean_phenotypes.values():
        for row_pigean in row_patient.get(KEY_DATA):
            pigean_phenotype_id = row_pigean.get(KEY_ID)
            if not map_pigean_phenotype_genes.get(pigean_phenotype_id):
                map_pigean_phenotype_genes[pigean_phenotype_id] = get_rest_genes_for_pigean_phenotype(term=pigean_phenotype_id)

    # for each patient phenotype, get the weighted gene scores
    for phenotype_patient in map_patient_pigean_phenotypes.keys():
        map_patient_phenotype_gene_weights[phenotype_patient] = compute_weighted_gene_scores(phenotypes_json=map_patient_pigean_phenotypes.get(phenotype_patient), 
                                                                                             pigean_results=map_pigean_phenotype_genes)
    # debug log
    if debug:
        print("got patient phenotype to gene map: \n{}".format(json.dumps(map_patient_phenotype_gene_weights, indent=2)))

    # TODO: return 2 lists
    return map_patient_phenotype_gene_weights


from typing import Dict, Any
import json

def combine_patient_gene_probabilities(phenotype_gene_map: Dict[str, Any], default_prob: float = 0.05) -> Dict[str, Any]:
    """
    Given a mapping of phenotype -> list of genes with weighted_score,
    returns two aggregated patient-level gene lists:
      1. max gene probability across all phenotypes
      2. average gene probability (assigning default_prob if gene missing)
    """

    logs = []
    result = {"max_genes": [], "avg_genes": [], "logs": logs}

    try:
        # --- Step 1: Collect phenotype-level data ---
        phenotype_genes = {}  # phenotype -> {gene: score}
        for phenotype, pdata in phenotype_gene_map.items():
            pgenes = {}
            for entry in pdata.get("data", []):
                try:
                    gene = str(entry.get("gene", "")).strip()
                    score = float(entry.get("weighted_score", 0.0))
                    if gene:
                        pgenes[gene] = score
                except (ValueError, TypeError):
                    logs.append(f"Skipping malformed entry in {phenotype}: {entry}")
            phenotype_genes[phenotype] = pgenes

        if not phenotype_genes:
            logs.append("No phenotype gene data found.")
            return result

        # --- Step 2: Determine all unique genes across phenotypes ---
        all_genes = set()
        for gmap in phenotype_genes.values():
            all_genes.update(gmap.keys())

        # --- Step 3: Compute per-gene aggregates ---
        max_scores = {}
        avg_scores = {}

        n_phenotypes = len(phenotype_genes)
        for gene in all_genes:
            # collect all phenotype scores (default if missing)
            scores = [phenotype_genes[p].get(gene, default_prob) for p in phenotype_genes]

            # compute metrics
            max_scores[gene] = max(scores)
            avg_scores[gene] = sum(scores) / n_phenotypes

        # --- Step 4: Sort results descending by score ---
        result["max_genes"] = [
            {"gene": g, "max_score": s}
            for g, s in sorted(max_scores.items(), key=lambda x: x[1], reverse=True)
        ]
        result["avg_genes"] = [
            {"gene": g, "avg_score": s}
            for g, s in sorted(avg_scores.items(), key=lambda x: x[1], reverse=True)
        ]

        logs.append(f"Processed {len(phenotype_genes)} phenotypes, {len(all_genes)} unique genes.")

    except Exception as e:
        logs.append(f"Error combining gene probabilities: {e}")

    # --- single return point ---
    return result


def process_web_patient_phenotypes(
    list_input_phenotypes: List[str], debug: bool=False) -> Dict[str, List[Dict]]:
    '''
    handles the list of phenotypes from the REST call
    '''
    # get the weighted gene list for each phenotype
    map_patient_result = process_patient_phenotypes(list_input_phenotypes=list_input_phenotypes, debug=debug)

    # get the two genes lists from the above result
    map_gene_scores = combine_patient_gene_probabilities(phenotype_gene_map=map_patient_result)

    # return
    return map_gene_scores


# main
if __name__ == "__main__":
    # list of phenotypes
    list_patient_phenotypes = ['diabetes', 'arthritis']

    # # get the weighted gene list for each phenotype
    # map_patient_result = process_patient_phenotypes(list_input_phenotypes=list_patient_phenotypes, debug=True)

    # # get the two genes lists from the above result
    # map_gene_scores = combine_patient_gene_probabilities(phenotype_gene_map=map_patient_result)

    map_gene_scores = process_web_patient_phenotypes(list_input_phenotypes=list_patient_phenotypes, debug=True)
    print("got final gene scores: \n{}".format(json.dumps(map_gene_scores, indent=2)))



    # # test the similarity embeddings REST
    # result_phenotypes = get_rest_phenotype_similarity(term='diabetes')

    # # print
    # print(json.dumps(result_phenotypes, indent=2))

    # # test the pigean call per phenotype
    # result_genes = {}
    # for row in result_phenotypes.get(KEY_DATA):
    #     result_genes[row.get('id')] = get_rest_genes_for_pigean_phenotype(term=row.get('id'))

    # # print
    # print(json.dumps(result_genes, indent=2))

    # # get the weighted list for a patient phentype
    # result_weighted_genes = compute_weighted_gene_scores(phenotypes_json=result_phenotypes, pigean_results=result_genes)

    # # print
    # print(json.dumps(result_weighted_genes, indent=2))



