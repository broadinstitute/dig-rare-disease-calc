
# imports
import requests
import json


# constants
URL_PHENOTYPE_ENBEDDINGS = "https://api.kpndataregistry.org/api/search/phenotypes"
KEY_DATA = 'data'
KEY_LOGS = 'logs'
KEY_ID = 'id'
KEY_DESCRIPTION = 'description'
KEY_SCORE = 'score'


# methods
def get_rest_phenotype_similarity(term: str = "diabetes", similarity_threshold: float = 0.0):
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


# main
if __name__ == "__main__":
    # test the similarity embeddings REST
    list_phenotypes = get_rest_phenotype_similarity(term='diabetes')

    # print
    print(json.dumps(list_phenotypes, indent=2))