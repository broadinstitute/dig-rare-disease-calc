
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
def get_rest_phenotype_similarity(term: str, similarity_threshold: float = 0.0):
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
        for item in json_data.get(KEY_DATA, []):
            try:
                entry = {
                    KEY_ID: str(item.get(KEY_ID, "")),
                    KEY_DESCRIPTION: str(item.get(KEY_DESCRIPTION, "")),
                    KEY_SCORE: float(item.get(KEY_SCORE, 0.0))
                }

                result["data"].append(entry)

            except (ValueError, TypeError) as e:
                result.get(KEY_LOGS).append(f"Skipping malformed entry: {e}")

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



# main
if __name__ == "__main__":
    # test the similarity embeddings REST
    result_phenotypes = get_rest_phenotype_similarity(term='diabetes')

    # print
    print(json.dumps(result_phenotypes, indent=2))

    # test the pigean call per phenotype
    result_genes = get_rest_genes_for_pigean_phenotype(term='IBD')

    # print
    print(json.dumps(result_genes, indent=2))


