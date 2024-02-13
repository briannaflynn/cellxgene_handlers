# cellxgene_handlers
Scripts for processing high dimensional data from the CZI cell x gene dataset 

* [The main repo](https://github.com/chanzuckerberg/single-cell-curation/tree/main) 
* [Dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md)

Quickstart:

1) Use this command to get all available collections
```
curl -X 'GET' 'https://api.cellxgene.cziscience.com/curation/v1/collections?visibility=PUBLIC' -H 'accept: application/json' > collections.json
```

2) Extract collection IDs in collections.json with ```write_collection_ids.py```, and pass to  ```process_collections.py``` to get directory structure organized by collection ID with metadata (in JSON format) and URL information for subsequent datasets, with URL beginning in https://datasets.cellxgene.cziscience.com/ and ending in .h5ad.

3) Use ```iter_download.py``` to download the data. Optional: Use ```generate_parallel_job.py``` and ```split_for_launch.py``` for parallel download on TACC.

4) Use ```h5ad_summarizer.py``` to check out .h5ad formatted datasets.

   

   
