import requests
import pandas as pd
import sys
import os
import json

def fetch_and_save_dataset_urls(collection_id, create_directory=False):
    domain_name = "cellxgene.cziscience.com"
    site_url = f"https://{domain_name}"
    api_url_base = f"https://api.{domain_name}"
    collection_path = f"/curation/v1/collections/{collection_id}"
    collection_url = f"{api_url_base}{collection_path}"
    
    res = requests.get(url=collection_url)
    res.raise_for_status()  # Raises a HTTPError if the response status code is 4XX or 5XX
    res_content = res.json()
    urls = [dataset['assets'][0]['url'] for dataset in res_content['datasets']]
    
    df = pd.DataFrame({'Dataset_URLs': urls})
    
    # option to create a directory for the collection_id
    if create_directory:
        directory_path = os.path.join(os.getcwd(), 'collections', collection_id)
        os.makedirs(directory_path, exist_ok=True)  # Create the directory if it doesn't exist
        csv_filename = os.path.join(directory_path, f'{collection_id}_dataset_urls.csv')
        json_filename = os.path.join(directory_path, f'{collection_id}_metadata.json')
    else:
        csv_filename = f'{collection_id}_dataset_urls.csv'
        json_filename = f'{collection_id}_metadata.json'
    
    df.to_csv(csv_filename, header=None, index=False)
    with open(json_filename, 'w') as json_file:
        json.dump(res_content, json_file, indent=4)
    
    print(f"Dataset URLs have been saved to {csv_filename}")
    print(f"Metadata has been saved to {json_filename}")

# example for testing
# collection_id = sys.argv[1]
# fetch_and_save_dataset_urls(collection_id, create_directory=True)
