#!/usr/bin/python

import requests
import pandas as pd
import sys

def fetch_and_save_dataset_urls(collection_id):
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
    csv_filename = f'{collection_id}_dataset_urls.csv'
    df.to_csv(csv_filename, header=None, index=False)
    
    print(f"Dataset URLs have been saved to {csv_filename}")

# Example usage:
collection_id = sys.argv[1]
fetch_and_save_dataset_urls(collection_id)
