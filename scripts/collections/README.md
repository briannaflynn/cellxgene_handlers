Example sub-directory

Each directory contains a file with {collection_id}_metadata.json, a JSON file containing all metadata such as publication information related to the collection and associated datasets. 

Files named {collection_id}_dataset_urls.csv contain a list of URLs to pull the associated datasets (in .h5ad format). 

The URL listings are used by ```iter_download.py``` and ```download_utility.py``` to pull the datasets from CELLXGENE.

If you're in a pinch, you can also use wget to pull data by the urls like so : ```wget your_url.h5ad```
