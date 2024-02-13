#!/usr/bin/python

import requests
import sys

collection_id = sys.argv[1]
domain_name = "cellxgene.cziscience.com"
site_url = f"https://{domain_name}"
api_url_base = f"https://api.{domain_name}"

collection_path = f"/curation/v1/collections/{collection_id}"
collection_url = f"{api_url_base}{collection_path}"
res = requests.get(url=collection_url)
res.raise_for_status()
res_content = res.json()
print(res_content)
