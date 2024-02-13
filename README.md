# cellxgene_handlers
Scripts for processing high dimensional data from the CZI cell x gene dataset 

* [The main repo](https://github.com/chanzuckerberg/single-cell-curation/tree/main) 
* [Dataset schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/4.0.0/schema.md)

Notes:

```curl -X 'GET' 'https://api.cellxgene.cziscience.com/curation/v1/collections?visibility=PUBLIC' -H 'accept: application/json' > collections.json```

1) Use this command to get all available collections

2) Use each collection ID in collections.json, and pass to  ```request_collection.py``` to get URL information with dataset URL beginning in https://datasets.cellxgene.cziscience.com/ and ending in .h5ad

3) Use wget to download data from each dataset URL, rename so that {collection-id}_{dataset_id}.h5ad to make things easier to organize later on

   
